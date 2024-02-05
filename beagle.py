import pandas as pd
import numpy as np

###########################
#FOR THE IMPUTED VCF TABLE#
###########################

#Get column names
names = pd.read_csv('imputed.table', sep='\t', nrows=0).columns.tolist()
#Remove last name because it is "Unnamed" and a nonsense column
names = names[0:-1]

#Creates a dictionary of new column names, removing the extra annoying characters
new_names = {}

for ind,name in enumerate(names):
    #the CHROM label has more extra characters than the others
    if 'CHROM' in name:
        new_names[name] = name[5:]
    #index 0-7 are the field level labels
    elif 1 <= ind <= 7:
        new_names[name] = name[3:]
    #sample level labels starts at index 8
    elif ind > 7:
        split = name.split(']')
        split = split[1].split(':')
        new_names[name] = split[0]

#Read in the CSV as it is
imputed_df = pd.read_csv('imputed.table', sep='\t')
#Remove the last column because it is a nonsense column
imputed_df.drop(imputed_df.columns[-1], axis=1, inplace=True)
#Rename the columns with the dictionary created above
imputed_df.rename(columns = new_names, inplace=True)

#This will differ from the original below because this contains missing genotypes
replacement_dict = {'./.': -9,'0/0': 0, '0/1': 1, '1/0' : 1, '1/1': 2, '0|0': 0, '0|1': 1, '1|0' : 1, '1|1': 2}
#use the replacement dictinary to substitute genotypes for numerical values
imputed_df.iloc[:, 7:] = imputed_df.iloc[:, 7:].replace(replacement_dict)
            

 THE MASKED VCF TABLE#
##########################

#See notes above for what this code does
#Get column names
names = pd.read_csv('masked.table', sep='\t', nrows=0).columns.tolist()
#Remove last name because it is "Unnamed" and a nonsense column
names = names[0:-1]

new_names = {}

for ind,name in enumerate(names):
    #the CHROM label has more extra characters than the others
    if 'CHROM' in name:
        new_names[name] = name[5:]
    #index 0-7 are the field level labels
    elif 1 <= ind <= 7:
        new_names[name] = name[3:]
    #sample level labels starts at index 8
    elif ind > 7:
        split = name.split(']')
        split = split[1].split(':')
        new_names[name] = split[0]

#Read in the CSV as it is
masked_df = pd.read_csv('masked.table', sep='\t')
#Remove the last column because it is a nonsense column
masked_df.drop(masked_df.columns[-1], axis=1, inplace=True)
#Rename the columns with the dictionary created above
masked_df.rename(columns = new_names, inplace=True)
#This will differ from the original below because this contains missing genotypes
replacement_dict = {'./.' : -9, '0/0': 0, '0/1': 1, '1/0' : 1, '1/1': 2, '.|.' : -9, '0|0': 0, '0|1': 1, '1|0' : 1, '1|1': 2}
#use the replacement dictinary to substitute genotypes for numerical values
masked_df.iloc[:, 8:] = masked_df.iloc[:, 8:].replace(replacement_dict)

############################
#FOR THE ORIGINAL VCF TABLE#
############################

#See notes above for what this code does
#Get column names
names = pd.read_csv('original.table', sep='\t', nrows=0).columns.tolist()
#Remove last name because it is "Unnamed" and a nonsense column
names = names[0:-1]

new_names = {}

for ind,name in enumerate(names):
    #the CHROM label has more extra characters than the others
    if 'CHROM' in name:
        new_names[name] = name[5:]
    #index 0-7 are the field level labels
    elif 1 <= ind <= 7:
        new_names[name] = name[3:]
    #sample level labels starts at index 8
    elif ind > 7:
        split = name.split(']')
        split = split[1].split(':')
        new_names[name] = split[0]

#Read in the CSV as it is
original_df = pd.read_csv('original.table', sep='\t')
#Remove the last column because it is a nonsense column
original_df.drop(original_df.columns[-1], axis=1, inplace=True)
#Rename the columns with the dictionary created above
original_df.rename(columns=new_names, inplace=True)
#This will differ from the original below because this contains missing genotypes
replacement_dict = {'0/0': 0, '0/1': 1, '1/0' : 1, '1/1': 2, '0|0': 0, '0|1': 1, '1|0' : 1, '1|1': 2}
#Use the replacement dictionary to substitute genotypes for numerical values
original_df.iloc[:, 8:] = original_df.iloc[:, 8:].replace(replacement_dict)

# Extract unique values from imputed_df's 'ID' column
ids_to_keep = imputed_df['ID'].unique()

# Subset masked_df based on 'ID'
masked_df_subset = masked_df[masked_df['ID'].isin(ids_to_keep)].reset_index(drop=True)

# Subset original_df based on 'ID'
original_df_subset = original_df[original_df['ID'].isin(ids_to_keep)].reset_index(drop=True)

imputed_df = imputed_df.reset_index(drop=True)

# Initialize an empty list to store dictionaries for each iteration
rows_data = []

# Iterate through the rows of masked_df_subset
for index, row in masked_df_subset.iterrows():
    # Find columns where the value is -9 in masked_df_subset
    masked_columns = row[row == -9].index

    # Iterate through masked_columns
    for sample in masked_columns:
        # Extract values from masked_df_subset, imputed_df_subset, and original_df_subset
        masked_value = masked_df_subset.at[index, sample]
        imputed_value = imputed_df.iloc[index, imputed_df.columns.get_loc(sample)]
        original_value = original_df_subset.at[index, sample]

        # Append a dictionary for the current iteration
        rows_data.append({
            'ID': masked_df_subset.at[index, 'ID'],
            'Sample': sample,
            'MAF': masked_df_subset.at[index, 'AF'],
            'Masked': masked_value,
            'Imputed': imputed_value,
            'Original': original_value
        })

# Create a DataFrame from the list of dictionaries
result_df = pd.DataFrame(rows_data)

result_df.to_csv('B54_Ref54_Ne100_err05Maps_window200.tsv', sep='\t', index=False)
