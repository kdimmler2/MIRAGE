itrs = ['set' + str(i) for i in range(1, config['iterations'] + 1)]
chrs = [str(i) for i in range(1,32)]
#chrs = ['1']
maf_bins = ['0_0.05', '0.05_0.1', '0.1_0.15', '0.15_0.2', '0.2_0.25', '0.25_0.3', '0.3_0.35', '0.35_0.4', '0.4_0.45', '0.45_0.5']
dr2_bins = ['0_0.15', '0.15_0.3', '0.3_0.45', '0.45_0.6', '0.6_0.75', '0.75_0.9', '0.9_1.0']


rule all:
    input:
        expand('0.1_masking_results/mask_assess/{itr}/reformat/imputed/split_tables', itr=itrs),
        expand('0.1_masking_results/mask_assess/{itr}/reformat/masked/split_tables', itr=itrs),
        expand('0.1_masking_results/mask_assess/{itr}/reformat/original/split_tables', itr=itrs),
        expand('0.1_masking_results/mask_assess/{itr}/reformat/reformatted_split_tables/reformatted_chunk_{i}.csv', itr=itrs, i=range(1, 101)),
        expand('0.1_masking_results/mask_assess/{itr}/dr2_table.txt', itr=itrs),
        expand('0.1_masking_results/mask_assess/{itr}/reformat/combined_table.txt', itr=itrs),
        expand('0.1_masking_results/mask_assess/{itr}/split_tables/', itr=itrs),
        expand('0.1_masking_results/mask_assess/{itr}/split_tables_maf/', itr=itrs),
        expand('0.1_masking_results/mask_assess/{itr}/split_tables_dr2/', itr=itrs),
        expand('0.1_masking_results/mask_assess/all_sets/chroms/chr{chr}/all_sets.csv', chr=chrs),
        expand('0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/all_sets.csv', maf_bin=maf_bins),
        expand('0.1_masking_results/mask_assess/all_sets/DR2/{dr2_bin}/all_sets.csv', dr2_bin=dr2_bins),
        expand('0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/disc_maf_genotype.csv', maf_bin=maf_bins),
        expand('0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/disc_maf_overall.csv', maf_bin=maf_bins),
        expand('0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/dr2_maf_genotype.csv', maf_bin=maf_bins),
        expand('0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/dr2_maf_overall.csv', maf_bin=maf_bins),
        expand('0.1_masking_results/mask_assess/all_sets/DR2/{dr2_bin}/final_metrics/disc_dr2_genotype.csv', dr2_bin=dr2_bins),
        expand('0.1_masking_results/mask_assess/all_sets/DR2/{dr2_bin}/final_metrics/disc_dr2_overall.csv', dr2_bin=dr2_bins),
        expand('0.1_masking_results/mask_assess/all_sets/chroms/chr{chr}/final_metrics/IQS.csv', chr=chrs),
        expand('0.1_masking_results/mask_assess/all_sets/chroms/chr{chr}/final_metrics/no_dupes.IQS.csv', chr=chrs),
        '0.1_masking_results/final_metrics/disc_maf_genotype.csv',
        '0.1_masking_results/final_metrics/disc_maf_overall.csv',
        '0.1_masking_results/final_metrics/dr2_maf_genotype.csv',
        '0.1_masking_results/final_metrics/dr2_maf_overall.csv',
        '0.1_masking_results/final_metrics/disc_dr2_genotype.csv',
        '0.1_masking_results/final_metrics/disc_dr2_overall.csv',
        '0.1_masking_results/final_metrics/IQS.csv',

rule split_imputed_table:
    input:
        imputed_full = '0.1_masking_results/mask_assess/{itr}/full_imputed.table'
    output:
        directory('0.1_masking_results/mask_assess/{itr}/reformat/imputed/split_tables') 
    params:
        num_chunks = 100
    resources:
        time   = 20,
        mem_mb = 6000
    run:
        # Open the input file
        with open(input.imputed_full, 'rt') as imp:
            # Read the header
            header = imp.readline().rstrip()

            # Read all remaining lines
            lines = imp.readlines()

            # Calculate the number of lines per chunk
            chunk_size = len(lines) // params.num_chunks
            remainder = len(lines) % params.num_chunks

            # Write the chunks
            for i in range(params.num_chunks):
                # Determine the lines for this chunk
                start_idx = i * chunk_size + min(i, remainder)
                end_idx = start_idx + chunk_size + (1 if i < remainder else 0)

                # Generate the output chunk file path
                chunk_file_path = '0.1_masking_results/mask_assess/{}/reformat/imputed/split_tables/chunk_{}.csv'.format(
                    wildcards.itr, i + 1
                )

                # Write the output chunk file
                os.makedirs(os.path.dirname(chunk_file_path), exist_ok=True)  # Ensure the directory exists
                with open(chunk_file_path, 'wt') as chunk_file:
                    chunk_file.write(header + '\n')
                    chunk_file.writelines(lines[start_idx:end_idx])

        print('Chunks written successfully for iteration: {}'.format(wildcards.itr))

rule split_masked_table:
    input:
        masked_full = '0.1_masking_results/mask_assess/{itr}/full_masked.table'
    output:
        directory('0.1_masking_results/mask_assess/{itr}/reformat/masked/split_tables') 
    params:
        num_chunks = 100
    resources:
        time   = 20,
        mem_mb = 6000
    run:
        # Open the input file
        with open(input.masked_full, 'rt') as msk:
            # Read the header
            header = msk.readline().rstrip()

            # Read all remaining lines
            lines = msk.readlines()

            # Calculate the number of lines per chunk
            chunk_size = len(lines) // params.num_chunks
            remainder = len(lines) % params.num_chunks

            # Write the chunks
            for i in range(params.num_chunks):
                # Determine the lines for this chunk
                start_idx = i * chunk_size + min(i, remainder)
                end_idx = start_idx + chunk_size + (1 if i < remainder else 0)

                # Generate the output chunk file path
                chunk_file_path = '0.1_masking_results/mask_assess/{}/reformat/masked/split_tables/chunk_{}.csv'.format(
                    wildcards.itr, i + 1
                )

                # Write the output chunk file
                os.makedirs(os.path.dirname(chunk_file_path), exist_ok=True)  # Ensure the directory exists
                with open(chunk_file_path, 'wt') as chunk_file:
                    chunk_file.write(header + '\n')
                    chunk_file.writelines(lines[start_idx:end_idx])

        print('Chunks written successfully for iteration: {}'.format(wildcards.itr))

rule split_original_table:
    input:
        original_full = '0.1_masking_results/mask_assess/{itr}/full_original.table'
    output:
        directory('0.1_masking_results/mask_assess/{itr}/reformat/original/split_tables')
    params:
        num_chunks = 100
    resources:
        time   = 20,
        mem_mb = 6000
    run:
        # Open the input file
        with open(input.original_full, 'rt') as org:
            # Read the header
            header = org.readline().rstrip()

            # Read all remaining lines
            lines = org.readlines()

            # Calculate the number of lines per chunk
            chunk_size = len(lines) // params.num_chunks
            remainder = len(lines) % params.num_chunks

            # Write the chunks
            for i in range(params.num_chunks):
                # Determine the lines for this chunk
                start_idx = i * chunk_size + min(i, remainder)
                end_idx = start_idx + chunk_size + (1 if i < remainder else 0)

                # Generate the output chunk file path
                chunk_file_path = '0.1_masking_results/mask_assess/{}/reformat/original/split_tables/chunk_{}.csv'.format(
                    wildcards.itr, i + 1
                )

                # Write the output chunk file
                os.makedirs(os.path.dirname(chunk_file_path), exist_ok=True)  # Ensure the directory exists
                with open(chunk_file_path, 'wt') as chunk_file:
                    chunk_file.write(header + '\n')
                    chunk_file.writelines(lines[start_idx:end_idx])

        print('Chunks written successfully for iteration: {}'.format(wildcards.itr))

rule reformat_chunk:
    input:
        masked = '0.1_masking_results/mask_assess/{itr}/reformat/masked/split_tables/chunk_{i}.csv', 
        original = '0.1_masking_results/mask_assess/{itr}/reformat/original/split_tables/chunk_{i}.csv',
        imputed = '0.1_masking_results/mask_assess/{itr}/reformat/imputed/split_tables/chunk_{i}.csv',
    output:
        reformatted_chunk = '0.1_masking_results/mask_assess/{itr}/reformat/reformatted_split_tables/reformatted_chunk_{i}.csv'
        #directory = directory('0.1_masking_results/mask_assess/{itr}/reformat/reformatted_split_tables'),
    resources:
        time = 30,
        mem_mb = 4000
    run:
        import pandas as pd
        import os

        # Mapping dictionary
        mapping_dict = {
            "0/0": 0,
            "0/1": 1,
            "1/0": 1,
            "1/1": 2,
            "0|0": 0,
            "0|1": 1,
            "1|0": 1,
            "1|1": 2,
            "./.": pd.NA
        }

        # Load data
        print(f"Loading data for chunk {wildcards.i}...")
        try:
            # Format chunk number to match {i}
            chunk_number = str(int(wildcards.i))

            # Construct file paths
           # masked_file = os.path.join(input.masked_dir, "chunk_" + chunk_number + ".csv")
           # original_file = os.path.join(input.original_dir, "chunk_" + chunk_number + ".csv")
           # imputed_file = os.path.join(input.imputed_dir, "chunk_" + chunk_number + ".csv")

            # Load the data into DataFrames
            masked = pd.read_csv(input.masked, sep="\t", index_col=False).replace("NA", pd.NA)
            original = pd.read_csv(input.original, sep="\t", index_col=False).replace("NA", pd.NA)
            imputed = pd.read_csv(input.imputed, sep="\t", index_col=False).replace("NA", pd.NA)

            print(f"Data loaded for chunk {wildcards.i}. Masked: {masked.shape}, Original: {original.shape}, Imputed: {imputed.shape}")
        except Exception as e:
            print(f"Error loading data for chunk {wildcards.i}: {e}")
            raise

        # Ensure MAF, AC, and QUAL are preserved
        print(f"Retaining MAF, AC, and QUAL for chunk {wildcards.i}...")
        metadata_columns = ["ID", "MAF", "AC", "QUAL"]  # Add MAF, AC, QUAL explicitly
        data_columns = [col for col in masked.columns if col not in metadata_columns]

        # Apply the mapping dictionary to all data columns
        print(f"Applying mapping dictionary to chunk {wildcards.i}...")
        for df in [masked, original]:
            df[data_columns] = df[data_columns].applymap(lambda x: mapping_dict.get(x, pd.NA))

        print(f"Mapping applied. Checking for NA values in masked DataFrame for chunk {wildcards.i}...")
        print(masked.isna().sum())

        # Process the data
        print(f"Processing chunk {wildcards.i}...")
        results = []
        for idx, row in masked.iterrows():
            variant_id = row["ID"]
            # Retain MAF, AC, QUAL values
            maf = row["MAF"]
            ac = row["AC"]
            qual = row["QUAL"]

            print(f"Processing row {idx}: Variant ID = {variant_id}, MAF = {maf}, AC = {ac}, QUAL = {qual}")
            for col in data_columns:
                if pd.isna(row[col]):
                    print(f"Found NA at row {idx}, column {col}")
                    try:
                        results.append({
                            "Variant": variant_id,
                            "Sample": col,
                            "MAF": maf,
                            "AC": ac,
                            "QUAL": qual,
                            "Original": original.loc[idx, col],
                            "Imputed": imputed.loc[idx, col]
                        })
                    except KeyError as e:
                        print(f"Error accessing Original or Imputed value at row {idx}, column {col}: {e}")

        print(f"Processing complete for chunk {wildcards.i}. Total results: {len(results)}")

        # Create DataFrame
        print(f"Creating DataFrame for chunk {wildcards.i}...")
        final_df = pd.DataFrame(results)

        # Step 1: Exclude metadata rows
        metadata_samples = ["CHROM", "POS", "REF", "ALT", "FILTER"]
        filtered_df = final_df[~final_df["Sample"].isin(metadata_samples)]

        # Step 2: Reorder columns to place MAF, AC, and QUAL after Sample
        column_order = ["Variant", "Sample", "MAF", "AC", "QUAL", "Original", "Imputed"]
        reordered_df = filtered_df[column_order]

        # Step 3: Save the results
        output_file = output.reformatted_chunk  # Replace with the desired output file path
        reordered_df.to_csv(output_file, index=False)
        print(f"Output saved to {output_file}")

rule combine_tables:
    input:
        reformatted_chunk = expand('0.1_masking_results/mask_assess/{itr}/reformat/reformatted_split_tables/reformatted_chunk_{i}.csv', itr=itrs, i=range(1, 101)), 
        input_directory = directory('0.1_masking_results/mask_assess/{itr}/reformat/reformatted_split_tables'),
    output:
        combined_table = '0.1_masking_results/mask_assess/{itr}/reformat/combined_table.txt',
    resources:
        time   = 20,
        mem_mb = 24000
    run:
        import pandas as pd
        import glob
        import os

        # Directory containing the reformatted chunk files
        reformatted_dir = input.input_directory  # Replace with your directory path
        output_file = output.combined_table  # The final combined output file

        # Get a list of all reformatted chunk files
        file_list = glob.glob(os.path.join(reformatted_dir, "reformatted_chunk_*.csv"))

        # Combine files
        print("Combining files...")
        combined_df = pd.concat(
            [pd.read_csv(file) for file in file_list], ignore_index=True
        )

        # Save the combined DataFrame
        combined_df.to_csv(output_file, index=False)
        print(f"Combined file saved to {output_file}. Total rows: {combined_df.shape[0]}")

rule get_DR2:
    input:
        imputed_vcf = '0.1_masking_results/FINAL/{itr}/imputed.vcf.gz',
    output:
        dr2_table = '0.1_masking_results/mask_assess/{itr}/dr2_table.txt', 
    resources:
        time   = 20,
        mem_mb = 2000
    shell:
        '''
            gatk VariantsToTable \
                -V {input.imputed_vcf} \
                -F ID \
                -F DR2 \
                -O {output.dr2_table}
        '''

rule add_DR2_combined:
    input:
        dr2_table = '0.1_masking_results/mask_assess/{itr}/dr2_table.txt',
        combined_table = '0.1_masking_results/mask_assess/{itr}/reformat/combined_table.txt'
    output:
        dr2_combined_table = '0.1_masking_results/mask_assess/{itr}/dr2_combined_table.txt'
    resources:
        time=60,
        mem_mb=60000
    run:
        import pandas as pd

        # Read the combined_table
        print(f"Reading combined_table: {input.combined_table}")
        try:
            combined_df = pd.read_csv(input.combined_table)
            print("Combined DataFrame (first few rows):")
            print(combined_df.head())
        except Exception as e:
            print(f"Error reading {input.combined_table}: {e}")
            raise

        # Read the dr2_table
        print(f"Reading dr2_table: {input.dr2_table}")
        try:
            dr2_df = pd.read_csv(input.dr2_table, sep='\t')
            print("DR2 DataFrame (first few rows):")
            print(dr2_df.head())
        except Exception as e:
            print(f"Error reading {input.dr2_table}: {e}")
            raise

        # Merge the DataFrames on Variant and ID
        print("Merging combined_table with dr2_table...")
        try:
            merged_df = combined_df.merge(dr2_df[["ID", "DR2"]], left_on="Variant", right_on="ID", how="left")
            merged_df = merged_df.drop(columns=["ID"])
            print("Merge successful. Merged DataFrame (first few rows):")
            print(merged_df.head())
        except Exception as e:
            print(f"Error during merge: {e}")
            raise

        # Write the merged DataFrame to the output file
        print(f"Writing merged DataFrame to: {output.dr2_combined_table}")
        try:
            merged_df.to_csv(output.dr2_combined_table, index=False)
            print(f"DR2 combined table written to: {output.dr2_combined_table}")
        except Exception as e:
            print(f"Error writing output file: {e}")
            raise

rule split_tables_IQS:
    input:
        dr2_combined_table = '0.1_masking_results/mask_assess/{itr}/dr2_combined_table.txt'
    output:
        split_table=directory('0.1_masking_results/mask_assess/{itr}/split_tables/') 
    resources:
        time=20,
        mem_mb=12000
    run:
        import os
        import pandas as pd

        # Read the input CSV
        print("Reading input CSV: {}".format(input.dr2_combined_table))
        try:
            df = pd.read_csv(input.dr2_combined_table)
            print("Input DataFrame (first few rows):")
            print(df.head())

            # Extract chromosome information from the Variant column
            df["Chromosome"] = df["Variant"].str.split("_").str[0]
            print("Extracted Chromosome column (unique values):")
            print(df["Chromosome"].unique())
        except Exception as e:
            print("Error reading input CSV: {}".format(e))
            raise

        # Split the DataFrame and save by chromosome
        for chr_num in range(1, 33):
            chr_name = "chr{}".format(chr_num)
            print("Processing {}...".format(chr_name))
            try:
                # Filter rows for the current chromosome
                chr_df = df[df["Chromosome"] == chr_name]
                print("Rows matching {}: {}".format(chr_name, len(chr_df)))

                # Skip if no rows match
                if chr_df.empty:
                    print("No rows found for {}. Skipping.".format(chr_name))
                    continue

                # Construct the output directory and file path
                output_dir = "{}/{}".format(output.split_table, chr_name)
                output_file = "{}/split_table.csv".format(output_dir)

                # Create the directory if it doesn't exist
                os.makedirs(output_dir, exist_ok=True)

                # Save the filtered DataFrame
                chr_df.to_csv(output_file, index=False)
                print("Written: {}".format(output_file))
            except Exception as e:
                print("Error processing {}: {}".format(chr_name, e))
                raise

rule split_tables_maf:
    input:
        dr2_combined_table='0.1_masking_results/mask_assess/{itr}/dr2_combined_table.txt'
    output:
        split_table=directory('0.1_masking_results/mask_assess/{itr}/split_tables_maf/')
    resources:
        time=20,
        mem_mb=12000
    run:
        import os
        import pandas as pd

        # Read the input CSV
        print("Reading input CSV: {}".format(input.dr2_combined_table))
        try:
            df = pd.read_csv(input.dr2_combined_table)
            print("Input DataFrame (first few rows):")
            print(df.head())
        except Exception as e:
            print("Error reading input CSV: {}".format(e))
            raise

        # Define MAF bins
        bin_edges = [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50]
        bin_labels = [f"[{bin_edges[i]}, {bin_edges[i+1]})" for i in range(len(bin_edges) - 1)]
        df["MAF_bin"] = pd.cut(df["MAF"], bins=bin_edges, labels=bin_labels, include_lowest=True)
        print("MAF bins created with labels:", bin_labels)

        # Split the DataFrame and save by MAF bin
        for maf_bin in bin_labels:
            print("Processing MAF bin {}...".format(maf_bin))
            try:
                # Filter rows for the current MAF bin
                bin_df = df[df["MAF_bin"] == maf_bin]
                print("Rows matching MAF bin {}: {}".format(maf_bin, len(bin_df)))

                # Skip if no rows match
                if bin_df.empty:
                    print("No rows found for MAF bin {}. Skipping.".format(maf_bin))
                    continue

                # Construct the output directory and file path
                output_dir = "{}/{}".format(output.split_table, maf_bin.replace(" ", "").replace("[", "").replace(")", "").replace(",", "_"))
                output_file = "{}/split_table.csv".format(output_dir)

                # Create the directory if it doesn't exist
                os.makedirs(output_dir, exist_ok=True)

                # Save the filtered DataFrame
                bin_df.to_csv(output_file, index=False)
                print("Written: {}".format(output_file))
            except Exception as e:
                print("Error processing MAF bin {}: {}".format(maf_bin, e))
                raise

rule split_tables_dr2:
    input:
        dr2_combined_table='0.1_masking_results/mask_assess/{itr}/dr2_combined_table.txt'
    output:
        split_table=directory('0.1_masking_results/mask_assess/{itr}/split_tables_dr2/')
    resources:
        time=20,
        mem_mb=12000
    run:
        import os
        import pandas as pd

        # Read the input CSV
        print("Reading input CSV: {}".format(input.dr2_combined_table))
        try:
            df = pd.read_csv(input.dr2_combined_table)
            print("Input DataFrame (first few rows):")
            print(df.head())
        except Exception as e:
            print("Error reading input CSV: {}".format(e))
            raise

        # Define MAF bins
        bin_edges = [0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.0]
        bin_labels = [f"[{bin_edges[i]}, {bin_edges[i+1]})" for i in range(len(bin_edges) - 1)]
        df["DR2_bin"] = pd.cut(df["DR2"], bins=bin_edges, labels=bin_labels, include_lowest=True)
        print("DR2 bins created with labels:", bin_labels)

        # Split the DataFrame and save by MAF bin
        for dr2_bin in bin_labels:
            print("Processing DR2 bin {}...".format(dr2_bin))
            try:
                # Filter rows for the current MAF bin
                bin_df = df[df["DR2_bin"] == dr2_bin]
                print("Rows matching DR2 bin {}: {}".format(dr2_bin, len(bin_df)))

                # Skip if no rows match
                if bin_df.empty:
                    print("No rows found for DR2 bin {}. Skipping.".format(dr2_bin))
                    continue

                # Construct the output directory and file path
                output_dir = "{}/{}".format(output.split_table, dr2_bin.replace(" ", "").replace("[", "").replace(")", "").replace(",", "_"))
                output_file = "{}/split_table.csv".format(output_dir)

                # Create the directory if it doesn't exist
                os.makedirs(output_dir, exist_ok=True)

                # Save the filtered DataFrame
                bin_df.to_csv(output_file, index=False)
                print("Written: {}".format(output_file))
            except Exception as e:
                print("Error processing DR2 bin {}: {}".format(dr2_bin, e))
                raise

rule gather_sets_chrom:
    input:
        split_tables_chrom=lambda wildcards: expand(
            '0.1_masking_results/mask_assess/{itr}/split_tables/chr{chr}/split_table.csv',
            itr=itrs,
            chr=[wildcards.chr]
        )
    output:
        all_sets='0.1_masking_results/mask_assess/all_sets/chroms/chr{chr}/all_sets.csv'
    resources:
        time=60,
        mem_mb=24000
    run:
        import pandas as pd

        # Initialize an empty list to store dataframes
        dfs = []

        # Loop over input files and read them into dataframes
        print("Reading split_tables files for chr{}:".format(wildcards.chr))
        for file in input.split_tables_chrom:
            print("Reading file: {}".format(file))
            try:
                df = pd.read_csv(file)
                dfs.append(df)
            except Exception as e:
                print("Error reading file {}: {}".format(file, e))

        # Concatenate all dataframes
        try:
            combined_df = pd.concat(dfs, ignore_index=True)
            print("Successfully concatenated split_tables for chr{}.".format(wildcards.chr))
        except Exception as e:
            print("Error during concatenation of split_tables for chr{}: {}".format(wildcards.chr, e))
            raise

        # Print the first few rows of combined_df for verification
        print("Final Combined DataFrame for chr{} (first few rows):".format(wildcards.chr))
        print(combined_df.head())

        # Write the combined dataframe to the output file
        try:
            combined_df.to_csv(output.all_sets, index=False)
            print("Combined table for chr{} written to: {}".format(wildcards.chr, output.all_sets))
        except Exception as e:
            print("Error writing combined table to file for chr{}: {}".format(wildcards.chr, e))
            raise

rule gather_sets_maf:
    input:
        split_tables_maf=lambda wildcards: expand(
            '0.1_masking_results/mask_assess/{itr}/split_tables_maf/{maf_bin}/split_table.csv',
            itr=itrs,
            maf_bin=[wildcards.maf_bin]
        )
    output:
        all_sets='0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/all_sets.csv'
    resources:
        time=60,
        mem_mb=24000
    run:
        import pandas as pd

        # Initialize an empty list to store dataframes
        dfs = []

        # Loop over input files and read them into dataframes
        print("Reading split_tables files for maf_bin: {}".format(wildcards.maf_bin))
        for file in input.split_tables_maf:
            print("Reading file: {}".format(file))
            try:
                df = pd.read_csv(file)
                dfs.append(df)
            except Exception as e:
                print("Error reading file {}: {}".format(file, e))

        # Concatenate all dataframes
        try:
            combined_df = pd.concat(dfs, ignore_index=True)
            print("Successfully concatenated split_tables for maf_bin: {}".format(wildcards.maf_bin))
        except Exception as e:
            print("Error during concatenation of split_tables for maf_bin {}: {}".format(wildcards.maf_bin, e))
            raise

        # Print the first few rows of combined_df for verification
        print("Final Combined DataFrame for maf_bin {} (first few rows):".format(wildcards.maf_bin))
        print(combined_df.head())

        # Write the combined dataframe to the output file
        try:
            combined_df.to_csv(output.all_sets, index=False)
            print("Combined table for maf_bin {} written to: {}".format(wildcards.maf_bin, output.all_sets))
        except Exception as e:
            print("Error writing combined table to file for maf_bin {}: {}".format(wildcards.maf_bin, e))
            raise

rule gather_sets_dr2:
    input:
        split_tables_dr2=lambda wildcards: expand(
            '0.1_masking_results/mask_assess/{itr}/split_tables_dr2/{dr2_bin}/split_table.csv',
            itr=itrs,
            dr2_bin=[wildcards.dr2_bin]
        )
    output:
        all_sets='0.1_masking_results/mask_assess/all_sets/DR2/{dr2_bin}/all_sets.csv'
    resources:
        time=60,
        mem_mb=24000
    run:
        import pandas as pd

        # Initialize an empty list to store dataframes
        dfs = []

        # Loop over input files and read them into dataframes
        print("Reading split_tables files for dr2_bin: {}".format(wildcards.dr2_bin))
        for file in input.split_tables_dr2:
            print("Reading file: {}".format(file))
            try:
                df = pd.read_csv(file)
                dfs.append(df)
            except Exception as e:
                print("Error reading file {}: {}".format(file, e))

        # Concatenate all dataframes
        try:
            combined_df = pd.concat(dfs, ignore_index=True)
            print("Successfully concatenated split_tables for dr2_bin: {}".format(wildcards.dr2_bin))
        except Exception as e:
            print("Error during concatenation of split_tables for dr2_bin {}: {}".format(wildcards.dr2_bin, e))
            raise

        # Print the first few rows of combined_df for verification
        print("Final Combined DataFrame for dr2_bin {} (first few rows):".format(wildcards.dr2_bin))
        print(combined_df.head())

        # Write the combined dataframe to the output file
        try:
            combined_df.to_csv(output.all_sets, index=False)
            print("Combined table for dr2_bin {} written to: {}".format(wildcards.dr2_bin, output.all_sets))
        except Exception as e:
            print("Error writing combined table to file for dr2_bin {}: {}".format(wildcards.dr2_bin, e))
            raise

rule disc_maf_metrics:
    input:
        all_sets='0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/all_sets.csv'
    output:
        disc_maf_genotype = '0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/disc_maf_genotype.csv',
        disc_maf_overall = '0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/disc_maf_overall.csv',
    resources:
        time   = 30,
        mem_mb = 24000
    run:
        import pandas as pd
        import numpy as np
        from scipy.stats import norm

        # Read input
        df = pd.read_csv(input.all_sets)

        # Add a column for the bin (based on the input file or context)
        df["Bin"] = wildcards.maf_bin  # Use the wildcard or manually specify the bin

        # Define conditions and corresponding values for Genotype
        conditions = [
            df["Original"] == 0,  # Original is 0
            df["Original"] == 1,  # Original is 1
            df["Original"] == 2   # Original is 2
        ]
        choices = ["Major", "Het", "Minor"]
        df["Genotype"] = np.select(conditions, choices, default="Unknown")

        # Calculate Discordance
        df["Discordance"] = (df["Original"] - df["Imputed"]).abs()

        # Function to calculate mean and confidence intervals
        def agg_with_ci(data, column):
            mean = data[column].mean()
            std = data[column].std()
            n = len(data[column])
            if n > 1:
                # 95% CI using normal approximation
                ci_range = 1.96 * (std / np.sqrt(n))
                lower_ci = mean - ci_range
                upper_ci = mean + ci_range
            else:
                # If n = 1, CIs cannot be calculated
                lower_ci, upper_ci = np.nan, np.nan
            return pd.Series({"mean": mean, "lower_ci": lower_ci, "upper_ci": upper_ci})

        # Calculate Discordance stats grouped by Genotype
        result_disc_genotype = (
            df.groupby("Genotype")
            .apply(lambda group: agg_with_ci(group, "Discordance"))
            .reset_index()
        )
        result_disc_genotype.rename(columns={"mean": "avg_Discordance"}, inplace=True)

        # Add the Bin column to the Genotype output
        result_disc_genotype["Bin"] = wildcards.maf_bin

        # Fill missing values with NaN explicitly
        result_disc_genotype["avg_Discordance"] = result_disc_genotype["avg_Discordance"].fillna(np.nan)
        result_disc_genotype["lower_ci"] = result_disc_genotype["lower_ci"].fillna(np.nan)
        result_disc_genotype["upper_ci"] = result_disc_genotype["upper_ci"].fillna(np.nan)

        # Save Genotype-specific results
        result_disc_genotype.to_csv(output.disc_maf_genotype, index=False, na_rep="NA")

        # Calculate Discordance stats without considering Genotype
        result_disc_overall = agg_with_ci(df, "Discordance").to_frame().T
        result_disc_overall.rename(columns={"mean": "avg_Discordance"}, inplace=True)

        # Add the Bin column to the Overall output
        result_disc_overall["Bin"] = wildcards.maf_bin

        # Fill missing values with NaN explicitly
        result_disc_overall["avg_Discordance"] = result_disc_overall["avg_Discordance"].fillna(np.nan)
        result_disc_overall["lower_ci"] = result_disc_overall["lower_ci"].fillna(np.nan)
        result_disc_overall["upper_ci"] = result_disc_overall["upper_ci"].fillna(np.nan)

        # Save Overall results
        result_disc_overall.to_csv(output.disc_maf_overall, index=False, na_rep="NA")

rule dr2_maf_metrics:
    input:
        all_sets='0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/all_sets.csv'
    output:
        dr2_maf_genotype = '0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/dr2_maf_genotype.csv',
        dr2_maf_overall = '0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/dr2_maf_overall.csv',
    resources:
        time   = 30,
        mem_mb = 24000
    run:
        import pandas as pd
        import numpy as np
        from scipy.stats import norm

        # Read input
        df = pd.read_csv(input.all_sets)

        # Add a column for the bin (based on the input file or context)
        df["Bin"] = wildcards.maf_bin  # Use the wildcard or manually specify the bin

        # Define conditions and corresponding values for Genotype
        conditions = [
            df["Original"] == 0,  # Original is 0
            df["Original"] == 1,  # Original is 1
            df["Original"] == 2   # Original is 2
        ]
        choices = ["Major", "Het", "Minor"]
        df["Genotype"] = np.select(conditions, choices, default="Unknown")

        # Function to calculate mean and confidence intervals
        def agg_with_ci(data, column):
            mean = data[column].mean()
            std = data[column].std()
            n = len(data[column])
            if n > 1:
                # 95% CI using normal approximation
                ci_range = 1.96 * (std / np.sqrt(n))
                lower_ci = mean - ci_range
                upper_ci = mean + ci_range
            else:
                # If n = 1, CIs cannot be calculated
                lower_ci, upper_ci = np.nan, np.nan
            return pd.Series({"mean": mean, "lower_ci": lower_ci, "upper_ci": upper_ci})

        # Calculate DR2 stats grouped by Genotype
        result_dr2_genotype = (
            df.groupby("Genotype")
            .apply(lambda group: agg_with_ci(group, "DR2"))
            .reset_index()
        )
        result_dr2_genotype.rename(columns={"mean": "avg_DR2"}, inplace=True)

        # Add the Bin column to the Genotype output
        result_dr2_genotype["Bin"] = wildcards.maf_bin

        # Fill missing values with NaN explicitly
        result_dr2_genotype["avg_DR2"] = result_dr2_genotype["avg_DR2"].fillna(np.nan)
        result_dr2_genotype["lower_ci"] = result_dr2_genotype["lower_ci"].fillna(np.nan)
        result_dr2_genotype["upper_ci"] = result_dr2_genotype["upper_ci"].fillna(np.nan)

        # Save Genotype-specific results
        result_dr2_genotype.to_csv(output.dr2_maf_genotype, index=False, na_rep="NA")

        # Calculate DR2 stats without considering Genotype
        result_dr2_overall = agg_with_ci(df, "DR2").to_frame().T
        result_dr2_overall.rename(columns={"mean": "avg_DR2"}, inplace=True)

        # Add the Bin column to the Overall output
        result_dr2_overall["Bin"] = wildcards.maf_bin

        # Fill missing values with NaN explicitly
        result_dr2_overall["avg_DR2"] = result_dr2_overall["avg_DR2"].fillna(np.nan)
        result_dr2_overall["lower_ci"] = result_dr2_overall["lower_ci"].fillna(np.nan)
        result_dr2_overall["upper_ci"] = result_dr2_overall["upper_ci"].fillna(np.nan)

        # Save Overall results
        result_dr2_overall.to_csv(output.dr2_maf_overall, index=False, na_rep="NA")

rule disc_dr2_metrics:
    input:
        all_sets='0.1_masking_results/mask_assess/all_sets/DR2/{dr2_bin}/all_sets.csv'
    output:
        disc_dr2_genotype = '0.1_masking_results/mask_assess/all_sets/DR2/{dr2_bin}/final_metrics/disc_dr2_genotype.csv',
        disc_dr2_overall = '0.1_masking_results/mask_assess/all_sets/DR2/{dr2_bin}/final_metrics/disc_dr2_overall.csv',
    resources:
        time   = 30,
        mem_mb = 48000
    run:
        import pandas as pd
        import numpy as np
        from scipy.stats import norm

        # Read input
        df = pd.read_csv(input.all_sets)

        # Add a column for the bin (based on the input file or context)
        df["Bin"] = wildcards.dr2_bin  # Use the wildcard or manually specify the bin

        # Define conditions and corresponding values for Genotype
        conditions = [
            df["Original"] == 0,  # Original is 0
            df["Original"] == 1,  # Original is 1
            df["Original"] == 2   # Original is 2
        ]
        choices = ["Major", "Het", "Minor"]
        df["Genotype"] = np.select(conditions, choices, default="Unknown")

        # Calculate Discordance
        df["Discordance"] = (df["Original"] - df["Imputed"]).abs()

        # Function to calculate mean and confidence intervals
        def agg_with_ci(data, column):
            mean = data[column].mean()
            std = data[column].std()
            n = len(data[column])
            if n > 1:
                # 95% CI using normal approximation
                ci_range = 1.96 * (std / np.sqrt(n))
                lower_ci = mean - ci_range
                upper_ci = mean + ci_range
            else:
                # If n = 1, CIs cannot be calculated
                lower_ci, upper_ci = np.nan, np.nan
            return pd.Series({"mean": mean, "lower_ci": lower_ci, "upper_ci": upper_ci})

        # Calculate Discordance stats grouped by Genotype
        result_disc_genotype = (
            df.groupby("Genotype")
            .apply(lambda group: agg_with_ci(group, "Discordance"))
            .reset_index()
        )
        result_disc_genotype.rename(columns={"mean": "avg_Discordance"}, inplace=True)

        # Add the Bin column to the Genotype output
        result_disc_genotype["Bin"] = wildcards.dr2_bin

        # Fill missing values with NaN explicitly
        result_disc_genotype["avg_Discordance"] = result_disc_genotype["avg_Discordance"].fillna(np.nan)
        result_disc_genotype["lower_ci"] = result_disc_genotype["lower_ci"].fillna(np.nan)
        result_disc_genotype["upper_ci"] = result_disc_genotype["upper_ci"].fillna(np.nan)

        # Save Genotype-specific results
        result_disc_genotype.to_csv(output.disc_dr2_genotype, index=False, na_rep="NA")

        # Calculate Discordance stats without considering Genotype
        result_disc_overall = agg_with_ci(df, "Discordance").to_frame().T
        result_disc_overall.rename(columns={"mean": "avg_Discordance"}, inplace=True)

        # Add the Bin column to the Overall output
        result_disc_overall["Bin"] = wildcards.dr2_bin

        # Fill missing values with NaN explicitly
        result_disc_overall["avg_Discordance"] = result_disc_overall["avg_Discordance"].fillna(np.nan)
        result_disc_overall["lower_ci"] = result_disc_overall["lower_ci"].fillna(np.nan)
        result_disc_overall["upper_ci"] = result_disc_overall["upper_ci"].fillna(np.nan)

        # Save Overall results
        result_disc_overall.to_csv(output.disc_dr2_overall, index=False, na_rep="NA")

rule IQS:
    input:
        all_sets='0.1_masking_results/mask_assess/all_sets/chroms/chr{chr}/all_sets.csv'
    output:
        IQS = '0.1_masking_results/mask_assess/all_sets/chroms/chr{chr}/final_metrics/IQS.csv',
    resources:
        time   = 720,
        mem_mb = 48000
    run:
        import pandas as pd
        import numpy as np
        from scipy.stats import norm

        # Read input
        df = pd.read_csv(input.all_sets)

        # Initialize results list
        results = []

        # Group by Variant
        for variant, group in df.groupby("Variant"):
            print(f"Processing Variant: {variant}")
            print(f"Group Data:\n{group[['Original', 'Imputed']]}\n")

            # Create a 3x3 n_ij matrix (AA, AB, BB)
            nij_matrix = np.zeros((3, 3))  # 3x3 for AA, AB, BB

            # Populate the n_ij matrix
            for true, dosage in zip(group["Original"], group["Imputed"]):
                n00 = max(0, 2 - dosage)
                n01 = max(0, 2 - 2 * abs(dosage - 1))
                n02 = max(0, dosage - 1)

                print(f"True Genotype: {true}, Dosage: {dosage}")
                print(f"Contributions: n00={n00}, n01={n01}, n02={n02}")

                if true in [0, 1, 2]:  # Ensure valid row indices
                    nij_matrix[true] += [n00, n01, n02]

                print(f"Updated n_ij matrix:\n{nij_matrix}\n")

            # Total sum of n_ij matrix
            n_total = np.sum(nij_matrix)
            print(f"Total n (sum of n_ij matrix): {n_total}\n")

            # Observed agreement (Po)
            po = np.trace(nij_matrix) / n_total if n_total > 0 else 0
            print(f"Observed agreement (P_o): {po}\n")

            # Marginals
            row_marginals = np.sum(nij_matrix, axis=1)  # Row totals
            col_marginals = np.sum(nij_matrix, axis=0)  # Column totals
            print(f"Row Marginals: {row_marginals}")
            print(f"Column Marginals: {col_marginals}\n")

            # Chance agreement (Pc)
            pc = np.sum((row_marginals * col_marginals) / (n_total ** 2)) if n_total > 0 else 0
            print(f"Chance agreement (P_c): {pc}\n")

            # Handle perfect concordance
            if po == 1 and pc == 1:
                print(f"Perfect concordance detected for Variant: {variant}")
                iqs = 1
            else:
                iqs = (po - pc) / (1 - pc) if pc < 1 else 0  # Avoid division by zero

            print(f"Imputation Quality Score (IQS): {iqs}\n")

            # Append results
            results.append({"Variant": variant, "IQS": iqs})

        # Create results DataFrame
        results_df = pd.DataFrame(results)

        # Merge MAF values into results_df
        results_df = results_df.merge(df[["Variant", "MAF"]].drop_duplicates(), on="Variant", how="left")

        # Save to CSV (optional)
        results_df.to_csv(output.IQS, index=False)

rule remove_dup_IQS:
    input:
        IQS = '0.1_masking_results/mask_assess/all_sets/chroms/chr{chr}/final_metrics/IQS.csv',
    output:
        IQS_no_dupes = '0.1_masking_results/mask_assess/all_sets/chroms/chr{chr}/final_metrics/no_dupes.IQS.csv',
    resources:
        time   = 60,
        mem_mb = 24000
    run:
        infile = open(input.IQS, 'rt')
        outfile = open(output.IQS_no_dupes, 'wt')

        header = infile.readline().rstrip()
        print(header, file=outfile)

        d = {}

        for line in infile:
            line = line.rstrip()
            split = line.split(',')
            if split[0] not in d:
                d[split[0]] = line

        for key,value in d.items():
            print(value, file=outfile)

rule combine_disc_maf:
    input:
        disc_maf_genotype = expand('0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/disc_maf_genotype.csv', maf_bin=maf_bins),
        disc_maf_overall = expand('0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/disc_maf_overall.csv', maf_bin=maf_bins),
    output:
        combined_table_genotype = '0.1_masking_results/final_metrics/disc_maf_genotype.csv',
        combined_table_overall = '0.1_masking_results/final_metrics/disc_maf_overall.csv',
    resources:
        time=60,
        mem_mb=24000
    run:
        import pandas as pd

        # Combine Genotype-Specific Files
        print("Combining input files for disc_maf_genotype:")
        dfs_genotype = []
        for file in input.disc_maf_genotype:
            print(f"Reading file: {file}")
            try:
                df = pd.read_csv(file)
                dfs_genotype.append(df)
            except Exception as e:
                print(f"Error reading file {file}: {e}")

        # Concatenate genotype-specific dataframes
        try:
            combined_genotype_df = pd.concat(dfs_genotype, ignore_index=True)
            print("Successfully concatenated genotype-specific input files.")
        except Exception as e:
            print(f"Error during concatenation of genotype-specific files: {e}")
            raise

        # Print the first few rows of the genotype-specific combined DataFrame
        print("Genotype-Specific Combined DataFrame (first few rows):")
        print(combined_genotype_df.head())

        # Write the genotype-specific combined dataframe to the output file
        try:
            combined_genotype_df.to_csv(output.combined_table_genotype, index=False)
            print(f"Genotype-specific combined table written to: {output.combined_table_genotype}")
        except Exception as e:
            print(f"Error writing genotype-specific combined table to file: {e}")
            raise

        # Combine Overall Files
        print("\nCombining input files for disc_maf_overall:")
        dfs_overall = []
        for file in input.disc_maf_overall:
            print(f"Reading file: {file}")
            try:
                df = pd.read_csv(file)
                dfs_overall.append(df)
            except Exception as e:
                print(f"Error reading file {file}: {e}")

        # Concatenate overall dataframes
        try:
            combined_overall_df = pd.concat(dfs_overall, ignore_index=True)
            print("Successfully concatenated overall input files.")
        except Exception as e:
            print(f"Error during concatenation of overall files: {e}")
            raise

        # Print the first few rows of the overall combined DataFrame
        print("Overall Combined DataFrame (first few rows):")
        print(combined_overall_df.head())

        # Write the overall combined dataframe to the output file
        try:
            combined_overall_df.to_csv(output.combined_table_overall, index=False)
            print(f"Overall combined table written to: {output.combined_table_overall}")
        except Exception as e:
            print(f"Error writing overall combined table to file: {e}")
            raise

rule combine_dr2_maf:
    input:
        dr2_maf_genotype = expand('0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/dr2_maf_genotype.csv', maf_bin=maf_bins),
        dr2_maf_overall = expand('0.1_masking_results/mask_assess/all_sets/MAF/{maf_bin}/final_metrics/dr2_maf_overall.csv', maf_bin=maf_bins),
    output:
        combined_table_genotype = '0.1_masking_results/final_metrics/dr2_maf_genotype.csv',
        combined_table_overall = '0.1_masking_results/final_metrics/dr2_maf_overall.csv',
    resources:
        time=60,
        mem_mb=24000
    run:
        import pandas as pd

        # Combine Genotype-Specific Files
        print("Combining input files for dr2_maf_genotype:")
        dfs_genotype = []
        for file in input.dr2_maf_genotype:
            print(f"Reading file: {file}")
            try:
                df = pd.read_csv(file)
                dfs_genotype.append(df)
            except Exception as e:
                print(f"Error reading file {file}: {e}")

        # Concatenate genotype-specific dataframes
        try:
            combined_genotype_df = pd.concat(dfs_genotype, ignore_index=True)
            print("Successfully concatenated genotype-specific input files.")
        except Exception as e:
            print(f"Error during concatenation of genotype-specific files: {e}")
            raise

        # Print the first few rows of the genotype-specific combined DataFrame
        print("Genotype-Specific Combined DataFrame (first few rows):")
        print(combined_genotype_df.head())

        # Write the genotype-specific combined dataframe to the output file
        try:
            combined_genotype_df.to_csv(output.combined_table_genotype, index=False)
            print(f"Genotype-specific combined table written to: {output.combined_table_genotype}")
        except Exception as e:
            print(f"Error writing genotype-specific combined table to file: {e}")
            raise

        # Combine Overall Files
        print("\nCombining input files for dr2_maf_overall:")
        dfs_overall = []
        for file in input.dr2_maf_overall:
            print(f"Reading file: {file}")
            try:
                df = pd.read_csv(file)
                dfs_overall.append(df)
            except Exception as e:
                print(f"Error reading file {file}: {e}")

        # Concatenate overall dataframes
        try:
            combined_overall_df = pd.concat(dfs_overall, ignore_index=True)
            print("Successfully concatenated overall input files.")
        except Exception as e:
            print(f"Error during concatenation of overall files: {e}")
            raise

        # Print the first few rows of the overall combined DataFrame
        print("Overall Combined DataFrame (first few rows):")
        print(combined_overall_df.head())

        # Write the overall combined dataframe to the output file
        try:
            combined_overall_df.to_csv(output.combined_table_overall, index=False)
            print(f"Overall combined table written to: {output.combined_table_overall}")
        except Exception as e:
            print(f"Error writing overall combined table to file: {e}")
            raise

rule combine_disc_dr2:
    input:
        disc_dr2_genotype = expand('0.1_masking_results/mask_assess/all_sets/DR2/{dr2_bin}/final_metrics/disc_dr2_genotype.csv', dr2_bin=dr2_bins),
        disc_dr2_overall = expand('0.1_masking_results/mask_assess/all_sets/DR2/{dr2_bin}/final_metrics/disc_dr2_overall.csv', dr2_bin=dr2_bins),
    output:
        combined_table_genotype = '0.1_masking_results/final_metrics/disc_dr2_genotype.csv',
        combined_table_overall = '0.1_masking_results/final_metrics/disc_dr2_overall.csv',
    resources:
        time=60,
        mem_mb=24000
    run:
        import pandas as pd

        # Combine Genotype-Specific Files
        print("Combining input files for disc_dr2_genotype:")
        dfs_genotype = []
        for file in input.disc_dr2_genotype:
            print(f"Reading file: {file}")
            try:
                df = pd.read_csv(file)
                dfs_genotype.append(df)
            except Exception as e:
                print(f"Error reading file {file}: {e}")

        # Concatenate genotype-specific dataframes
        try:
            combined_genotype_df = pd.concat(dfs_genotype, ignore_index=True)
            print("Successfully concatenated genotype-specific input files.")
        except Exception as e:
            print(f"Error during concatenation of genotype-specific files: {e}")
            raise

        # Print the first few rows of the genotype-specific combined DataFrame
        print("Genotype-Specific Combined DataFrame (first few rows):")
        print(combined_genotype_df.head())

        # Write the genotype-specific combined dataframe to the output file
        try:
            combined_genotype_df.to_csv(output.combined_table_genotype, index=False)
            print(f"Genotype-specific combined table written to: {output.combined_table_genotype}")
        except Exception as e:
            print(f"Error writing genotype-specific combined table to file: {e}")
            raise

        # Combine Overall Files
        print("\nCombining input files for disc_dr2_overall:")
        dfs_overall = []
        for file in input.disc_dr2_overall:
            print(f"Reading file: {file}")
            try:
                df = pd.read_csv(file)
                dfs_overall.append(df)
            except Exception as e:
                print(f"Error reading file {file}: {e}")

        # Concatenate overall dataframes
        try:
            combined_overall_df = pd.concat(dfs_overall, ignore_index=True)
            print("Successfully concatenated overall input files.")
        except Exception as e:
            print(f"Error during concatenation of overall files: {e}")
            raise

        # Print the first few rows of the overall combined DataFrame
        print("Overall Combined DataFrame (first few rows):")
        print(combined_overall_df.head())

        # Write the overall combined dataframe to the output file
        try:
            combined_overall_df.to_csv(output.combined_table_overall, index=False)
            print(f"Overall combined table written to: {output.combined_table_overall}")
        except Exception as e:
            print(f"Error writing overall combined table to file: {e}")
            raise

rule combine_IQS:
    input:
        IQS = expand('0.1_masking_results/mask_assess/all_sets/chroms/chr{chr}/final_metrics/no_dupes.IQS.csv', chr=chrs)
    output:
        combined_table = '0.1_masking_results/final_metrics/IQS.csv'
    resources:
        time=60,
        mem_mb=24000
    run:
        import pandas as pd

        # Initialize an empty list to store dataframes
        dfs = []

        # Loop over input files and read them into dataframes
        print("Combining input files for IQS:")
        for file in input.IQS:
            print(f"Reading file: {file}")
            try:
                df = pd.read_csv(file)
                dfs.append(df)
            except Exception as e:
                print(f"Error reading file {file}: {e}")

        # Concatenate all dataframes
        try:
            combined_df = pd.concat(dfs, ignore_index=True)
            print("Successfully concatenated all input files.")
        except Exception as e:
            print(f"Error during concatenation: {e}")
            raise

        # Print the first few rows of the combined DataFrame for verification
        print("Combined DataFrame (first few rows):")
        print(combined_df.head())

        # Write the combined dataframe to the output file, ensuring the header is written only once
        try:
            combined_df.to_csv(output.combined_table, index=False)
            print(f"Combined table written to: {output.combined_table}")
        except Exception as e:
            print(f"Error writing combined table to file: {e}")
            raise
