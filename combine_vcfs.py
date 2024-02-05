import pandas as pd
import vcf

nugen_path = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr{chrom}/chr{chrom}.vcf' 

# Open the VCF file using pyvcf
vcf_reader = vcf.Reader(open(nugen_path, 'r'))

# Extract header information from VCF
vcf_header_lines = []
for line in vcf_reader._header_lines:
    if line.startswith('##'):
        vcf_header_lines.append(line)
    else:
        break

# Extract data from VCF and create a list of dictionaries
vcf_data = []
for record in vcf_reader:
    variant_dict = {
        'CHROM': record.CHROM,
        'POS': record.POS,
        'ID': record.ID,
        'REF': record.REF,
        'ALT': ",".join(map(str, record.ALT)),
        'QUAL': record.QUAL,
        'FILTER': '.',
        'INFO': '.',
        'FORMAT': 'GT'
    }
    
    # Extract genotype information
    for sample in record.samples:
        sample_name = sample.sample
        variant_dict[sample_name] = sample['GT']

    vcf_data.append(variant_dict)

# Create a Pandas DataFrame from the list of dictionaries
nugen_df = pd.DataFrame(vcf_data)

stitch_path = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Imputation/STITCH/results/FINAL/split/chr{chrom}/stitch.chr{chrom}.vcf'

# Open the VCF file using pyvcf
vcf_reader = vcf.Reader(open(stitch_path, 'r'))

# Extract data from VCF and create a list of dictionaries
vcf_data = []
for record in vcf_reader:
    variant_dict = {
        'CHROM': record.CHROM,
        'POS': record.POS,
        'ID': record.ID,
        'REF': record.REF,
        'ALT': ",".join(map(str, record.ALT)),
        'QUAL': record.QUAL,
        'FILTER': '.',
        'INFO': '.',
        'FORMAT': 'GT'
    }
    
    # Extract genotype information
    for sample in record.samples:
        sample_name = sample.sample
        variant_dict[sample_name] = sample['GT']

    vcf_data.append(variant_dict)

# Create a Pandas DataFrame from the list of dictionaries
stitch_df = pd.DataFrame(vcf_data)

common_columns = stitch_df.columns.intersection(nugen_df.columns)

stitch_df = stitch_df[common_columns]

# Reorder columns 4 to the end for nugen_df
nugen_df = nugen_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + list(nugen_df.columns[9:])]

# Reorder columns 4 to the end for stitch_df based on the order in nugen_df
stitch_df = stitch_df[['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + list(nugen_df.columns[9:])]

# Iterate through the rows and columns of nugen_df
for index, row in nugen_df.iterrows():
    for col in nugen_df.columns.difference(['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']):
        # Check if the genotype is './.'
        if row[col] == './.':
            # Find the corresponding spot in stitch_df
            stitch_value = stitch_df.loc[
                (stitch_df['ID'] == row['ID']) & (stitch_df[col] != './.'),
                col
            ].values
            # Replace the './.' genotype with the corresponding value from stitch_df if available
            if len(stitch_value) > 0:
                nugen_df.at[index, col] = stitch_value[0]

# Check for rows in stitch_df that are not in nugen_df and add them
stitch_df_ids_not_in_nugen = set(stitch_df['ID']) - set(nugen_df['ID'])
stitch_rows_to_add = stitch_df[stitch_df['ID'].isin(stitch_df_ids_not_in_nugen)].copy()
stitch_rows_to_add

nugen_df = pd.concat([nugen_df, stitch_rows_to_add], ignore_index=True)
nugen_df.fillna('.', inplace=True)
nugen_df.sort_values(by='POS', ascending=True, inplace=True)

output_vcf_path = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Imputation/STITCH/results/FINAL/split/chr{chrom}/stitch.chr{chrom}.vcf'

# Open the output VCF file for writing
with open(output_vcf_path, 'w') as output_vcf:
    # Write the header lines
    for line in vcf_reader._header_lines:
        output_vcf.write(str(line) + '\n')

    # Write the column headers line
    output_vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
    for sample in nugen_df.columns:
        if sample not in ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']:
            output_vcf.write(f"\t{sample}")
    output_vcf.write('\n')

    # Iterate through the rows of nugen_df and write formatted VCF lines
    for _, row in nugen_df.iterrows():
        # Format the VCF line
        vcf_line = f"{row['CHROM']}\t{row['POS']}\t{row['ID']}\t{row['REF']}\t{row['ALT']}\t{row['QUAL']}\t{row['FILTER']}\t{row['INFO']}\t{row['FORMAT']}"
        
        # Add genotype information
        for sample in nugen_df.columns:
            if sample not in ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']:
                vcf_line += f"\t{row[sample]}"

        # Write the formatted VCF line to the output file
        output_vcf.write(vcf_line + '\n')
