from rdkit import Chem
from rdkit.Chem import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize

from rdkit import rdBase
rdBase.DisableLog('rdApp.*')

import pandas

class Mol:

    #activities should be a dict of name:value (e.g. {'kd':7})
    #assume nanomolar for concentrations
    def __init__(self, smiles, precise_activities, imprecise_activities):
        self.input_smiles = smiles
        self.precise_activities = precise_activities
        self.imprecise_activities = imprecise_activities

        try:
            self.inchi = process_smiles(smiles)
        except Exception as e:
            raise

    def has_activity(self, activity_name):
        return activity_name in self.precise_activities or activity_name in self.imprecise_activities

    def has_precise_activity(self, activity_name):
        return activity_name in self.precise_activities

    def get_precise_activity(self, activity_name):
        return self.precise_activities[activity_name]

    def has_imprecise_activity(self, activity_name):
        return activity_name in self.imprecise_activities

    def get_imprecise_activity(self, activity_name):
        return self.imprecise_activities[activity_name]


    def __str__(self):
        return f"{self.smiles} {self.activities}"

class ManualReviewException(Exception):
    pass

#todo: do without df.iterrows()
#needs original filename for manual review cases
def get_activities(df, original_filename, activity_fields, smiles_field = "SMILES"):

    activity_hits = {}
    for activity_field in activity_fields:
        activity_hits[activity_field] = []
        for col_name in df.columns:
            if activity_field.lower() in col_name.lower():
                if activity_field in activity_hits:
                    activity_hits[activity_field].append(col_name)

    smiles_hits = []
    for col_name in df.columns:
            if smiles_field.lower() in col_name.lower():
                smiles_hits.append(col_name)


    for activity_field, hits in activity_hits.items():
        if len(hits) == 0:
            print(f"Supplied activity field name '{activity_field}' not found in dataframe. Exiting.")
            exit()
        if len(hits) > 1:
            print(f"Supplied activity field name '{activity_field}' is a substring of multiple fields.\n Hits: {activity_hits}\n Exiting.")
            exit()
        else: 
            activity_hits[activity_field] = hits[0]
            print(f"Activity field found: {activity_hits[activity_field]}")

    if len(smiles_hits) == 0:
        print(f"Supplied SMILES field name '{activity_field}' not found in dataframe. Exiting.")
        exit()
    if len(smiles_hits) > 1:
        print(f"Supplied SMILES field name '{activity_field}' is a substring of multiple fields.\n Hits: {smiles_hits}\n Exiting.")
        exit()
    else: 
        print(f"SMILES field found: {smiles_hits[0]}")
        smiles_name = smiles_hits[0]


    #get set of all rows with valid entries of any activity
    valid_cols = []
    for activity_name, df_activity_name in activity_hits.items():

        v = pandas.notnull(df[df_activity_name])
        valid_cols.append(v)

    final = valid_cols[0]
    for v in valid_cols[1:]:
        final = final | v


    trimmed_df = df[final]

    mols = []
    for_review = []

    for i, row in df.iterrows():
        smiles = row[smiles_name]

        precise = {}
        imprecise = {}
        for activity_name, df_activity_name in activity_hits.items():
            if pandas.notnull(row[df_activity_name]):
                val = row[df_activity_name]
                try:
                    val = float(val)
                    precise[activity_name] = val
                except:
                    imprecise[activity_name] = val

        try:
            mol = Mol(smiles, precise, imprecise)
            mols.append(mol)
        except Exception as e:
            s = original_filename + ", " + str(e)
            for_review.append(s)
            continue

    return mols, for_review

def main():

    targets = ["ic50", "ki"]

    output_ending = "curated.csv"
    review_filename = "for_review.txt"

    all_for_review = []
    all_mols = []

    filenames = ["bindingdb_1286.tsv", "pdbbind_hits.csv"]
    for filename in filenames:
        print(filename)
        if ".csv" in filename:
            sep = ","
        if ".tsv" in filename:
            sep = "\t"

        df = pandas.read_csv(filename, sep = sep)

        mols, for_review = get_activities(df, original_filename = filename, activity_fields = targets)
        for target in targets:
            t = [x for x in mols if x.has_activity(target)]
            print(f"{filename} {target} hits: {len(t)}")
        all_mols.extend(mols)
        all_for_review.extend(for_review)


    dedup = {}
    for target in targets:
        dedup[target] = {}

        print(f"Deduplicating {target}")

        mols = [x for x in all_mols if x.has_precise_activity(target)]
        imp_mols = [x for x in all_mols if x.has_imprecise_activity(target)]
        print(f"Precise count for {target}: {len(mols)}")
        print(f"Imprecise count for {target}: {len(imp_mols)}")

        print(f"Length before deduplication: {len(mols)}")
        #deduplicate
        d = {}
        for mol in mols:
            if mol.inchi in d:
                d[mol.inchi].append(mol.get_precise_activity(target))
            else:
                d[mol.inchi] = [mol.get_precise_activity(target)]

        for key, value in d.items():
            if(len(value) > 1):
                min_val = min(value)
                max_val = max(value)

                if max_val / min_val > 10:
                    all_for_review.append(
                            f'{key}, "values differ by more than an order of magnitude", {value}')
                else:
                    avg_val = sum(value) / len(value)
                    dedup[target][key] = avg_val

            else: #no clash
                dedup[target][key] = value[0]

        print(f"Length after deduplication: {len(dedup[target])}")

    for target, d in dedup.items():
        final_filename = f"{target}_{output_ending}"
        f = open(final_filename, 'w')
        for key, value in d.items():
            f.write(f"{key}, {value}\n")
        f.close()

    f = open(review_filename, 'w')
    for line in all_for_review:
        f.write(line + '\n')



#performs structure standardization and mixture handling, outputs inchi
def process_smiles(smiles):

    mol = Chem.MolFromSmiles(smiles)

    salt_remover = SaltRemover.SaltRemover()
    #remove salts
    mol = salt_remover.StripMol(mol)

    #standardize
    start_inchi = Chem.MolToInchi(mol)
    rdMolStandardize.Cleanup(mol)
    end_inchi = Chem.MolToInchi(mol)
    if start_inchi != end_inchi:
        standardized += 1

    #check for mixtures
    n = Chem.rdmolops.GetMolFrags(mol)
    if len(n) != 1:
        raise ManualReviewException(f'{Chem.MolToSmiles(mol)}, "multiple fragments detected"')

    inchi = Chem.MolToInchi(mol)

    return inchi

if __name__ == "__main__":
    main()
