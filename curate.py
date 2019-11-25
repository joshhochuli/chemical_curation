from rdkit import Chem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import PandasTools

from rdkit import rdBase

from molvs import Standardizer

rdBase.DisableLog('rdApp.*')

import os

import math #for rounding

import pandas

class Mol:

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
        return cls(inchi, precise_activities, imprecise_activities, smiles = smiles)

    #activities should be a dict of name:value (e.g. {'kd':7})
    #assume nanomolar for concentrations
    def __init__(self, inchi, precise_activities, imprecise_activities):

        if not precise_activities and not imprecise_activities:
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
        return f"{self.input_smiles}\n Precise: {self.precise_activities}\n Imprecise: {self.imprecise_activities}\n"

    def __repr__(self):
        return self.__str__()

class ManualReviewException(Exception):
    pass

#needs original filename for manual review cases
#parses dataframe for mols with desired activities
def get_activities(df, original_filename, activity_fields, mol_field = "mol"):

    activity_hits = {}
    for activity_field in activity_fields:
        activity_hits[activity_field] = []
        for col_name in df.columns:
            if activity_field.lower() in col_name.lower():
                if activity_field in activity_hits:
                    activity_hits[activity_field].append(col_name)

    mol_hits = []
    for col_name in df.columns:
        if mol_field.lower() in col_name.lower():
            mol_hits.append(col_name)

    to_remove = []
    for activity_field, hits in activity_hits.items():
        if len(hits) == 0:
            print(f"Supplied activity field name '{activity_field}' not found in dataframe. Skipping.")
            to_remove.append(activity_field)
            continue
        if len(hits) > 1:
            print(f"Supplied activity field name '{activity_field}' is a substring of multiple fields.\n Hits: {activity_hits}\n Exiting.")
            exit()
        else:
            activity_hits[activity_field] = hits[0]
            print(f"Activity field found: {activity_hits[activity_field]}")

    for item in to_remove:
            activity_hits.pop(item)

    if len(mol_hits) == 0:
        print(f"Supplied molecule field name '{mol_field}' not found in dataframe. Exiting.")
        exit()
    if len(mol_hits) > 1:
        print(f"Supplied molecule field name '{mol_field}' is a substring of multiple fields.\n Hits: {mol_hits}\n Exiting.")
        exit()
    else: 
        print(f"Molecule field found: {mol_hits[0]}")
        mol_name = mol_hits[0]

    for_review = {}

    #get set of all rows with valid entries of any activity
    valid_cols = []

    for activity_name, df_activity_name in activity_hits.items():

        #check for empty and null
        u = df[df_activity_name] != ''
        v = pandas.notnull(df[df_activity_name])
        valid = (u & v)
        valid_cols.append(u & v)

    #find rows without any valid activities and report
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


    trimmed_df = df[final]

    mols = []
    stats = {}

    for i, row in trimmed_df.iterrows():
        actual_row = i + 2 #index by one, account for header
        original_mol = row["mol"]

        precise = {}
        imprecise = {}
        for activity_name, df_activity_name in activity_hits.items():
            if pandas.notnull(row[df_activity_name]) and row[df_activity_name] != '':
                val = row[df_activity_name]
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



        try:
            mol = Mol.from_mol_string(precise_activities = precise, imprecise_activities = imprecise, mol_string = Chem.MolToMolBlock(original_mol))
            mols.append(mol)
        except Exception as e:
            s = f"{original_filename}\trow {actual_row}\t{str(e)}"
            if activity_name in for_review:
                for_review[activity_name].append(s)
            else:
                for_review[activity_name] = [s]
            continue

    return mols, stats, for_review

#if shared keys, form a list 
def extend_dict(a, b):

    for key, value in b.items():
        if key in a:
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
    processed_inchi = process_mol(mol)
    return processed_inchi

def process_mol(mol):

    s = Standardizer()
    mol = s.standardize(mol)
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

    all_for_review = {}
    all_mols = []

    for filename in filenames:
        print(filename)

        if ".sdf" in filename:
            df = PandasTools.LoadSDF(filename, molColName = "mol")
            mols, stats, for_review = get_activities(df,
                                                     original_filename = filename,
                                                     activity_fields = targets)

        else:
            if ".csv" in filename:
                sep = ","
            if ".tsv" in filename:
                sep = "\t"

            df = pandas.read_csv(filename, sep = sep)

            mols, stats, for_review = get_activities(df, original_filename = filename, activity_fields = targets)

        for target in targets:
            t = [x for x in mols if x.has_activity(target)]
            print(f"{filename} {target} hits: {len(t)}")

        all_mols.extend(mols)
        extend_dict(all_for_review, for_review)

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

    if verbose:

        print(f"Length before deduplication: {len(mols)}")
        print(f"Length after deduplication: {len(dedup)}")
        print(f"Total removed: {total_removed}")
    
        if data_type == "precise":
            print(f"{multiple_same_count} molecules with exact duplicate activities found")
            print(f"{multiple_diff_count} molecules with different activities found and averaged")
            print(f"{review_count} molecules with value differing by more than a factor of {review_threshold} and flagged for review")
        elif data_type == "imprecise":
            print(f"{multiple_same_count} molecules with exact duplicate activities found")
            print(f"{review_count} molecules with different values flagged for review")

    return dedup, for_review


def main():

    targets = ["ic50", "ki", "kd"]
    output_ending = "curated"
    review_threshold = 10
    output_dir = "/home/josh/tmp/curation_test"


    filenames = ["/home/josh/git/cdk9_design/data/uncleaned/sdf/chembl_cdk9.sdf"]

    all_mols, all_for_review = get_mols_from_files(filenames, targets)

    for target in targets:
        print(f"\nDeduplicating {target}")

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
            for issue in issues:
                total_for_review += 1
                s = f"    {issue}\n"
                f.write(s + '\n')

        if total_for_review == 0:
            print("\nNo items needing manual review")
        else:
            print(f"\n{total_for_review} items needing manual review, stored at {review_filename}\n")


if __name__ == "__main__":
    main()
