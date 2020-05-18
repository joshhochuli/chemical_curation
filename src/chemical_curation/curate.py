from rdkit import Chem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import PandasTools

from rdkit import rdBase

import molvs.normalize
import molvs.fragment
import molvs.tautomer
import molvs.metal

rdBase.DisableLog('rdApp.*')

import os

import math #for rounding

import pandas

import logging
import pathlib

#list of atoms allowed for dragon descriptor calculatoin
dragon_allowed_atoms = set(["H","B","C","N","O","F","Al","Si","P","S","Cl","Cr",
    "Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Mo","Ag","Cd","In",
    "Sn","Sb","Te","I","Gd","Pt","Au","Hg","Ti","Pb","Bi"])

class Mol:
    
    @classmethod
    def from_rdkit_mol(cls, mol, precise_activities = None, imprecise_activities = None):
        inchi = process_mol(mol)
        return cls(inchi, precise_activities, imprecise_activities)


    @classmethod
    def from_mol_string(cls, mol_string, precise_activities = None, imprecise_activities = None):
        inchi = process_mol_string(mol_string)
        return cls(inchi, precise_activities, imprecise_activities)

    @classmethod
    def from_inchi(cls, inchi, precise_activities = None, imprecise_activities = None):
        inchi = process_inchi(inchi)
        return cls(inchi, precise_activities, imprecise_activities)

    @classmethod
    def from_smiles(cls, smiles, precise_activities = None, imprecise_activities = None):
        inchi = process_smiles(smiles)
        return cls(inchi, precise_activities, imprecise_activities)

    #activities should be a dict of name:value (e.g. {'kd':7})
    #assume nanomolar for concentrations
    def __init__(self, inchi, precise_activities, imprecise_activities):

        if precise_activities == None and imprecise_activities == None:
            raise Exception("no activity provided to Mol init())")

        self.precise_activities = precise_activities
        self.imprecise_activities = imprecise_activities
        self.inchi = inchi

    def has_activity(self, activity_name):
        if self.precise_activities and self.imprecise_activities:
            return activity_name in self.precise_activities or activity_name in self.imprecise_activities
        if self.precise_activities:
            return activity_name in self.precise_activities
        if self.imprecise_activities:
            return activity_name in self.imprecise_activities
        return False

    def has_precise_activity(self, activity_name):
        if self.precise_activities:
            return activity_name in self.precise_activities
        return False

    def get_precise_activity(self, activity_name):
        return self.precise_activities[activity_name]

    def has_imprecise_activity(self, activity_name):
        if self.imprecise_activities:
            return activity_name in self.imprecise_activities
        return False

    def get_imprecise_activity(self, activity_name):
        return self.imprecise_activities[activity_name]

    def __str__(self):
        return f"{self.inchi}\n Precise: {self.precise_activities}\n Imprecise: {self.imprecise_activities}\n"

    def __repr__(self):
        return self.__str__()

class ManualReviewException(Exception):
    pass

#needs original filename for manual review cases
#parses dataframe for mols with desired activities
def get_activities(df, original_filename, activity_fields, mol_field = "mol"):
    """
    Find the columns for the mol and each target; if more than one column for
    any of these, we have a problem and we quit.

    Filters the df to get only rows that have a valid activity value for at least
    one target. Finds the activities (separated into 'precise' and 'imprecise') for
    each molecule.

    Return a list of Mol objects and a list of molecules that require manual review.

    """

    # why is /targets/ renamed to /activity_fields/?

    # for each target, look through the columns of the passed dataframe and see
    # if any of them contain the target
    # make a dict with the targets as the keys and a list of all column names
    # containing that target as the items.
    # Could rewrite as:
    ## activity_hits = {target:[] for target in activity_fields}
    ## for target in activity_hits:
    ##     activity_hits[target] = [col for col in df.columns if target.lower() in col.lower()]
        
    activity_hits = {}
    for activity_field in activity_fields:
        activity_hits[activity_field] = []
        for col_name in df.columns:
            if activity_field.lower() in col_name.lower():
                # why is this check necessary, didn't we just put it in activity hits 3 lines ago
                if activity_field in activity_hits:
                    activity_hits[activity_field].append(col_name)

    # look through the columns again, this time looking to see if they contain the mol_field string
    # Should this require equality, not containment, and then we just tell them to pass the full name
    # of the mol column, with the default being "mol"?
    ## mol_hits = [col in df.columns if col.lower().contains(mol_field.lower())]
    mol_hits = []
    for col_name in df.columns:
        if mol_field.lower() in col_name.lower():
            mol_hits.append(col_name)

    # remove the targets that we couldn't find a column for.
    # If we find multiple possible columns for a target, give up.
    to_remove = []
    for activity_field, hits in activity_hits.items():
        if len(hits) == 0:
            logging.warning(f"Supplied activity field name '{activity_field}' not found in dataframe. Skipping.")
            to_remove.append(activity_field)
            continue
        if len(hits) > 1:
            logging.warning(f"Supplied activity field name '{activity_field}' is a substring of multiple fields.\n Hits: {activity_hits}\n Exiting.")
            # I don't think exit() is what we want here.
            exit()
        else:
            activity_hits[activity_field] = hits[0]
            print(f"Activity field found: {activity_hits[activity_field]}")

    for item in to_remove:
        # why not just remove it when you found it???
        activity_hits.pop(item)

    # If we find zero or many columns with the mol_field string, give up.
    # Otherwise, select the one column we found.
    if len(mol_hits) == 0:
        logging.warning(f"Supplied molecule field name '{mol_field}' not found in dataframe. Exiting.")
        exit()
    if len(mol_hits) > 1:
        logging.warning(f"Supplied molecule field name '{mol_field}' is a substring of multiple fields.\n Hits: {mol_hits}\n Exiting.")
        exit()
    else: 
        logging.info(f"Molecule field found: {mol_hits[0]}")
        mol_name = mol_hits[0]

    for_review = {}

    #get set of all rows with valid entries of any activity
    valid_cols = []

    # For each target, get a Series obj /u/ of length of df with True if the cell
    # contents at that position is not equal to ''
    # For each target, get a Series obj /v/ of length of df with True if the cell
    # contents at that position are not NaN
    # Take the logical AND of the two Series, which is a new Series /valid/
    # Append /valid/ to a list
    # Do we deal with converting datasets that may have been loaded with 'nan' or other common NaN pitfalls?
    for activity_name, df_activity_name in activity_hits.items():        
        #check for empty and null
        # True if not empty, False if empty string
        u = df[df_activity_name] != ''
        # True if not null, False if null
        v = pandas.notnull(df[df_activity_name])
        valid = (u & v)
        valid_cols.append(u & v)

    #find rows without any valid activities and report

    # Find the logical OR of the boolean Series for all the target columns
    # Get the rows from the df that are still False at the end
    # Report them
    final = valid_cols[0]
    for v in valid_cols[1:]:
        final = final | v
    not_valid = df[~final]

    if len(not_valid) > 0:
        for row in not_valid.iterrows():
            s = f"row {row[0]} in {original_filename} has no valid activities"
            if "input_files" not in for_review:
                for_review["input_files"] = s

            else:
                for_review["input_files"].append(s)

    # Get the rows from the df that had at least one True in the composite boolean Series
    trimmed_df = df[final]

    
    mols = []
    # /stats/ keeps track of the count of the number of mols with 'precise' and
    # 'imprecise' activities for each target
    # Initialize so that we don't have to do the checking below?
    ## stats = {target:{'precise':0, 'imprecise':0} for target in activity_hits}
    # or, why not reuse /activity_hits/?
    stats = {}

    # DON'T USE ITERROWS
    for i, row in trimmed_df.iterrows():
        # only used for error messages
        actual_row = i + 2 #index by one, account for header
        # /mol_name/ is the column that we found to contain the mols
        original_mol = row[mol_name]

        # Then for each mol we have a dict of precise and imprecise activities that we find
        precise = {}
        imprecise = {}
        
        # For each target, we again check that this row has a valid entry
        # Do we even need to do the valid checking above then? Won't an invalid
        # item just drop through this check
        for activity_name, df_activity_name in activity_hits.items():
            if pandas.notnull(row[df_activity_name]) and row[df_activity_name] != '':
                val = row[df_activity_name]
                # What kind of values are we expecting to see for 'imprecise' activities?
                # If it were a boolean represented as 1 or 0, that would still be considered
                # 'precise' by this conversion
                try:
                    val = float(val)
                    precise[activity_name] = val
                except:
                    imprecise[activity_name] = val

                if activity_name not in stats:
                    stats[activity_name] = {}
                    if precise:
                        stats[activity_name]['precise'] = 1
                        stats[activity_name]['imprecise'] = 0
                    else:
                        stats[activity_name]['imprecise'] = 1
                        stats[activity_name]['precise'] = 0
                else:
                    if precise:
                        stats[activity_name]['precise'] += 1
                    else:
                        stats[activity_name]['imprecise'] += 1

        # Is this conditional why mol_field wasn't a passed parameter?
        # I think there is probably a different way to deal with this (maybe look
        # at the type of the contents of the mol column?)

        # Create a Mol obj (our obj, not an rdkit Mol obj) from the mol_string
        # If creating the mol fails, add the info to the /for_review/ dict
        if mol_field == "mol":
            try:
                # rdkit.Chem.MolToMolBlock: creates an MDL molblock as a string
                # Does this deal with all the cases we may encounter?
                # Why do we do this conversion?
                mol = Mol.from_mol_string(precise_activities = precise,
                                          imprecise_activities = imprecise,
                                          mol_string = Chem.MolToMolBlock(original_mol))
                mols.append(mol)
            except Exception as e:
                s = f"{original_filename}, row {actual_row}, {str(e)}"
                # stop doing this; initialize the dict with the keys and an empty list as the item for each
                if activity_name in for_review:
                    for_review[activity_name].append(s)
                else:
                    for_review[activity_name] = [s]
                continue

        elif mol_field == "smiles":
            try:
                mol = Mol.from_smiles(smiles = original_mol, precise_activities = precise,
                        imprecise_activities = imprecise)
                mols.append(mol)
            except Exception as e:
                s = f"{original_filename}, {original_mol}, {str(e)}"
                if activity_name in for_review:
                    for_review[activity_name].append(s)
                else:
                    for_review[activity_name] = [s]
                continue
            
    # Return /mols/ (list), /stats/ (dict), /for_review/ (dict)
    return mols, stats, for_review

#if shared keys, form a list 
def extend_dict(a, b):
    """
    Absorb dict /b/ into dict /a/.

    Take in two dict objects /a/ and /b/. For each item in /b/, check if the
    item's key is in /a/. 

    If yes, then append the item's value to the value for that key in /a/
    (assumed to be a list). 

    If no, create the item in /a/.

    This operation is done in place.
    """

    for key, value in b.items():
        if key in a:
            # why do we not need to check the type of /value/ here but we do below?
            # couldn't we end up with a[key] = [x, y, [value]] ?
            # also why not just enforce that /value/ is a list at creation?
            a[key].append(value)
        else:
            if isinstance(value, list):
                a[key] = value
            else:
                a[key] = [value]

#performs structure standardization and mixture handling, outputs inchi
def process_mol_string(mol_string):

    mol = Chem.MolFromMolBlock(mol_string)
    processed_inchi = process_mol(mol)
    return processed_inchi

def process_inchi(inchi):
    mol = Chem.MolFromInchi(inchi)
    processed_inchi = process_mol(mol)
    return processed_inchi

def process_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ManualReviewException(f"SMILES string '{smiles}' does not parse successfully")
    processed_inchi = process_mol(mol)
    return processed_inchi

def check_valid_atoms(mol, allowed_list = dragon_allowed_atoms):
    for atom in mol.GetAtoms():
        s = atom.GetSymbol()
        if s not in allowed_list:
            raise ManualReviewException(f"'{s}' atom not in list of allowed atoms") 

def process_mol(mol):

    #removal of mixtures
    fragmenter_object = molvs.fragment.LargestFragmentChooser(prefer_organic = True)
    mol = fragmenter_object.choose(mol)
    if mol is None:
        logging.info("Mixture removal failed for molecule")

    #removal of inorganics
    if not molvs.fragment.is_organic(mol):
        raise ManualReviewException("Molecule is not organic")

    #removal of salts
    remover = SaltRemover.SaltRemover()
    mol = remover.StripMol(mol, dontRemoveEverything=True) #tartrate is listed as a salt? what do?
    if mol is None:
        raise ManualReviewException("Salt removal failed for molecule")

    #structure normalization
    normalizer = molvs.normalize.Normalizer(normalizations=molvs.normalize.NORMALIZATIONS,
            max_restarts = molvs.normalize.MAX_RESTARTS)
    mol = normalizer.normalize(mol)
    if mol is None:
        raise ManualReviewException("Normalization failed for molecule")

    #tautomer selection
    tautomerizer = molvs.tautomer.TautomerCanonicalizer(transforms=molvs.tautomer.TAUTOMER_TRANSFORMS, scores =
            molvs.tautomer.TAUTOMER_SCORES, max_tautomers=molvs.tautomer.MAX_TAUTOMERS)
    if mol is None:
        raise ManualReviewException("Tautomerization failed for molecule")

    #disconnect metals
    metal_remover = molvs.metal.MetalDisconnector()
    mol = metal_remover.disconnect(mol)
    if mol is None:
        raise ManualReviewException("Metal removal failed for molecule")

    #final check for only valid atoms
    check_valid_atoms(mol)

    inchi = Chem.MolToInchi(mol)

    return inchi

def imprecise_strings_equal(s1, s2):

    s1_sign = s1[0]
    s2_sign = s2[0]
    if s1_sign != s2_sign:
        return False

    s1_val = float(s1[1:])
    s2_val = float(s2[1:])

    if s1_val != s2_val:
        return False

    return True

def get_mols_from_files(filenames, targets, verbose = True):
    """
    Read each file into its own Pandas dataframe. File type is based on the file
    extension. Currently supported filetypes are .sdf, .smi, .csv, and .tsv.

    For each file, extract the mols, stats, and the molecules that require review.

    Bring cleaned mols from all files into one list, /all_mols/, and all mols
    requiring review into one dict, /all_for_review/.
    """

    all_for_review = {}
    all_mols = []

    for filename in filenames:
        logging.info(filename)

        # Determine the type of the filename by the extension
        file_ext = pathlib.Path(filename).suffix
        ## Mol_field should probably be a passable agument, defaulting to "mol"?
        mol_field = "mol"

        # Read file depending on file extension
        if file_ext == ".sdf":
            df = PandasTools.LoadSDF(filename, molColName = mol_field)
        elif file_ext in [".csv", ".tsv", ".smi"]:
            sep = ","
            if file_ext == ".tsv":
                sep = "\t"
            if file_ext == ".smi":
                mol_field = "smiles"
            df = pandas.read_csv(filename, sep = sep)
        else:
            # TODO Throw an error
            pass

        # Stats is never used?
        mols, stats, for_review = get_activities(df, original_filename = filename,
                                                 activity_fields = targets,
                                                 mol_field = mol_field)

        # Report the number of mols with activity for each target
        for target in targets:
            # We iterate over all the mols A LOT. Can that be reduced at all?
            # Also why did we make and return the /stats/ dict if we were just going to count
            # the stuff in /mols/ to get the same info???
            t = [x for x in mols if x.has_activity(target)]
            logging.info(f"{filename} {target} hits: {len(t)}")

        # Add the mols from the file to the list of all mols    
        all_mols.extend(mols)
        # Add the mols that require review from this file to the dict of all mols requiring review
        extend_dict(all_for_review, for_review)

    # Return the list of /all_mols/ that have at least one valid activity and the mols that need to be reviewed
    return all_mols, all_for_review

def get_unique_mols(mols, data_type, target):

    unique_mols = {}

    for mol in mols:
        if mol.inchi in unique_mols:
            if(data_type == "precise"):
                unique_mols[mol.inchi].append(mol.get_precise_activity(target))
            elif(data_type == "imprecise"):
                unique_mols[mol.inchi].append(mol.get_imprecise_activity(target))
        else:
            if(data_type == "precise"):
                unique_mols[mol.inchi] = [mol.get_precise_activity(target)]
            elif(data_type == "imprecise"):
                unique_mols[mol.inchi] = [mol.get_imprecise_activity(target)]

    return unique_mols

def deduplicate(inchis, activities):

    d = {}

    for i, inchi in enumerate(inchis):
        if inchi in d:
            d[inchi].append(activities[i])
        else:
            d[inchi] = [activities[i]]

    return d


#data_type should be either "precise" or "imprecise"
def deduplicate_mols(mols, data_type, target, review_threshold, verbose = True):

    unique_mols = get_unique_mols(mols, data_type, target)

    for_review = {}
    multiple_diff_count = 0
    multiple_same_count = 0
    review_count = 0
    total_removed = 0

    dedup = {}
    for key, value in unique_mols.items():

        if data_type == "precise":
            if(len(value) > 1):
                min_val = min(value)
                max_val = max(value)

                if max_val / min_val > review_threshold:
                    s = f'{key}, {target}, "values differ by more than a factor of {review_threshold}", {value}'

                    if target in for_review:
                        for_review[target].append(s)
                    else:
                        for_review[target] = [s]

                    review_count += 1
                    total_removed += len(value)

                else:
                    if all(x == value[0] for x in value): #if all same
                        multiple_same_count += 1
                    else:
                        multiple_diff_count += 1
                    avg_val = sum(value) / len(value)
                    dedup[key] = avg_val
                    total_removed += len(value) - 1

            else: #no clash
                dedup[key] = value[0]
        elif data_type == "imprecise":
            if(len(value) > 1): 

                #can't really compare values, just see if they are all the same
                first_val = value[0]
                different = False
                for i in range(2,len(value)):
                    if not imprecise_strings_equal(first_val, value[i]):
                        different = True
                        break
                if different:
                    s = f'{key}, {target}, "imprecise values do not match", {value}'

                    if target in for_review:
                        for_review[target].append(s)
                    else:
                        for_review[target] = [s]

                    review_count += 1
                    total_removed += len(value)
                else:
                    multiple_same_count += 1
                    dedup[key] = value[0]
                    total_removed += len(value) - 1
            else:
                dedup[key] = value[0]

    logging.debug(f"Length before deduplication: {len(mols)}")
    logging.debug(f"Length after deduplication: {len(dedup)}")
    logging.debug(f"Total removed: {total_removed}")
    
    if data_type == "precise":
        logging.debug(f"{multiple_same_count} molecules with exact duplicate activities found")
        logging.debug(f"{multiple_diff_count} molecules with different activities found and averaged")
        logging.debug(f"{review_count} molecules with value differing by more than a factor of {review_threshold} and flagged for review")
    elif data_type == "imprecise":
        logging.debug(f"{multiple_same_count} molecules with exact duplicate activities found")
        logging.debug(f"{review_count} molecules with different values flagged for review")

    return dedup, for_review


def main(filenames, output_dir, targets, review_threshold, verbosity):

    # define logger
    int_log_level = getattr(logging, verbosity.upper())
    if not isinstance(int_log_level, int):
        raise ValueError('Invalid log level: %s' % int_log_level)
    logging.basicConfig(level=int_log_level)

    logging.info(f'Beginning curation for {len(filenames)} files.')
    logging.info(f'Sending results to {output_dir}')

    output_ending = "curated"

    # aggregation
    all_mols, all_for_review = get_mols_from_files(filenames, targets)

    # deduplication, done by target
    for target in targets:
        logging.info(f"\nDeduplicating {target}")

        precise_mols = [x for x in all_mols if x.has_precise_activity(target)]
        imprecise_mols = [x for x in all_mols if x.has_imprecise_activity(target)]

        for data_type in ["imprecise", "precise"]:
            print('\n' + data_type)

            if(data_type == "precise"):
                mols = precise_mols
            elif(data_type == "imprecise"):
                mols = imprecise_mols

            dedup, for_review = deduplicate_mols(mols, data_type = data_type, target = target, review_threshold = review_threshold) 

            extend_dict(all_for_review, for_review)

            target_dir = f"{output_dir}/{target}/{data_type}"
            os.makedirs(target_dir, exist_ok = True)
            csv_filename = f"{target_dir}/{target}_{data_type}_{output_ending}.tsv"
            sdf_filename = f"{target_dir}/{target}_{data_type}_{output_ending}.sdf"
            csv = open(csv_filename, 'w')

            mols = []
            for key,value in dedup.items():
                csv.write(f"{key}\t{value}\n")
                mol = Chem.MolFromInchi(key)
                mol.SetProp("Activity",str(value))
                mols.append(mol)
            csv.close()

            sdfwriter = Chem.rdmolfiles.SDWriter(sdf_filename)
            for mol in mols:
                sdfwriter.write(mol)

        review_dir = f"{output_dir}/for_review"
        os.makedirs(review_dir, exist_ok = True)
        review_filename = f"{review_dir}/for_review.txt"
        f = open(review_filename, 'w')

        total_for_review = 0 
        for target, issues in all_for_review.items():
            f.write(target + '\n')
            print('\n' + target + '\n')
            for issue in issues:
                total_for_review += 1
                s = f"    {issue}\n"
                f.write(s)
                print(s, end = "")

    if total_for_review == 0:
        print("\nNo items needing manual review")
    else:
        print(f"\n{total_for_review} items needing manual review, stored at {review_filename}\n")


if __name__ == "__main__":
    main(filenames = ["/home/josh/git/cdk9_design/data/uncleaned/sdf/chembl_cdk9.sdf"],
         output_dir = "/home/josh/tmp/curation_test",
         targets = ["ic50", "ki", "kd"],
         review_threshold=10,
         verbosity='INFO')

    #filenames = ["/home/josh/git/chemical_curation/test/failures.smi",
    #"/home/josh/git/cdk9_design/data/uncleaned/sdf/chembl_cdk9.sdf"]
    # filenames = ["/home/josh/git/chemical_curation/test/failures.smi"]
    #filenames = ["/home/josh/git/cdk9_design/data/uncleaned/sdf/chembl_cdk9.sdf"]
