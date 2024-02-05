rule all:
	input:
#		expand('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr{chrom}/chr{chrom}.vcf', chrom=[str(i) for i in range(1, 32)] + ['X']),
#		expand('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr{chrom}/plink.log', chrom=[str(i) for i in range(1, 32)] + ['X']),
#		expand('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr{chrom}/plink.imiss', chrom=[str(i) for i in range(1, 32)] + ['X']),
#		expand('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr{chrom}/plink.lmiss', chrom=[str(i) for i in range(1, 32)] + ['X']),
#		expand('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Imputation/STITCH/results/FINAL/split/chr{chrom}/stitch.chr{chrom}.vcf', chrom=[str(i) for i in range(1, 32)] + ['X']),
#		expand('stitch/results/combine/chr{chrom}/chr{chrom}.vcf', chrom=[str(i) for i in range(1, 32)] + ['X']),
#		expand('stitch/results/combine/chr{chrom}/plink.log', chrom=[str(i) for i in range(1, 32)] + ['X']),
#		expand('stitch/results/combine/chr{chrom}/plink.lmiss', chrom=[str(i) for i in range(1, 32)] + ['X']),
#		expand('stitch/results/combine/chr{chrom}/plink.imiss', chrom=[str(i) for i in range(1, 32)] + ['X']),
#		'/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/nugen_rate.txt',
#		'stitch/results/combined_rate.txt',
#		'stitch/results/compare_rates.txt',
#		expand('stitch/results/combine/chr{chrom}/chr{chrom}.vcf.gz', chrom=[str(i) for i in range(1, 32)] + ['X']),
#		expand('stitch/results/combine/chr{chrom}/chr{chrom}.vcf.gz.tbi', chrom=[str(i) for i in range(1, 32)] + ['X']),
#		'stitch/results/vcf.list',
		'stitch/results/full.combined.vcf.gz',
		'stitch/results/full.combined.vcf.gz.tbi'


rule split_nugen:
	input:
		nugen = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/nugen.Jan2024.vcf.gz',
	output:
		split_nugen = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr{chrom}/chr{chrom}.vcf',
	threads: 4
	resources:
		time	= 60,
		mem_mb	= 24000,
		cpus	= 4
	shell:
		'''
			bcftools view \
			-r chr{wildcards.chrom} \
			-o {output.split_nugen} \
			{input.nugen}
		'''
rule plink_nugen:
	input:
		split_nugen = rules.split_nugen.output.split_nugen 
	output:
		logfile = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr{chrom}/plink.log',
		imiss = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr{chrom}/plink.imiss',
		lmiss = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr{chrom}/plink.lmiss',
	threads: 4
	resources:
		time    = 20,
		mem_mb  = 24000,
		cpus    = 4
	params:
		outdir = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr{chrom}/plink'
	shell:
		'''
			plink \
			--missing \
			--horse \
			--out {params.outdir} \
			--vcf {input.split_nugen}
		'''

rule split_stitch:
	input:
		stitch = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Imputation/STITCH/results/FINAL/for_beagle.vcf.gz',
	output:
		split_stitch = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Imputation/STITCH/results/FINAL/split/chr{chrom}/stitch.chr{chrom}.vcf',
	threads: 4
	resources:
		time    = 60,
		mem_mb  = 24000,
		cpus    = 4
	shell:
		'''
			bcftools view \
			-r chr{wildcards.chrom} \
			-o {output.split_stitch} \
			{input.stitch}
		'''
rule plink_stitch:
	input:
		split_stitch = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Imputation/STITCH/results/FINAL/split/chr{chrom}/stitch.chr{chrom}.vcf' 
	output:
		logfile = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Imputation/STITCH/results/FINAL/split/chr{chrom}/plink.log',
		imiss = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Imputation/STITCH/results/FINAL/split/chr{chrom}/plink.imiss', 
		lmiss = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Imputation/STITCH/results/FINAL/split/chr{chrom}/plink.lmiss',
	threads: 4
	resources:
		time    = 20,
		mem_mb  = 24000,
		cpus    = 4
	params:
		outdir = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Imputation/STITCH/results/FINAL/split/chr{chrom}/plink',
	shell:
		'''
			plink \
			--missing \
			--horse \
			--out {params.outdir} \
			--vcf {input.split_stitch}
		'''

rule combine:
	input:
		split_nugen = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr{chrom}/chr{chrom}.vcf', 
		split_stitch = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/Imputation/STITCH/results/FINAL/split/chr{chrom}/stitch.chr{chrom}.vcf',
	output:
		combined = 'stitch/results/combine/chr{chrom}/chr{chrom}.vcf' 
	threads: 4
	resources:
		time    = 60,
		mem_mb  = 24000,
		cpus    = 4
	run:
		import pandas as pd
		import vcf

		nugen_path = input.split_nugen 

		# Open the VCF file using pyvcf
		vcf_reader = vcf.Reader(open(nugen_path, 'r'))

		# Extract header information from VCF
		vcf_header_lines = [line for line in vcf_reader._header_lines if line.startswith('##')]

		# Extract data from VCF using list comprehension
		vcf_data = [
			{
				'CHROM': record.CHROM,
				'POS': record.POS,
				'ID': record.ID,
				'REF': record.REF,
				'ALT': ",".join(map(str, record.ALT)),
				'QUAL': record.QUAL,
				'FILTER': '.',
				'INFO': '.',
				'FORMAT': 'GT',
				**{sample.sample: sample['GT'] for sample in record.samples}
			}
			for record in vcf_reader
		]

		# Create a Pandas DataFrame directly from the list of dictionaries
		nugen_df = pd.DataFrame(vcf_data)

		count_all = (nugen_df == './.').sum().sum()
		count_all

		stitch_path = input.split_stitch 

		# Open the VCF file using pyvcf
		vcf_reader = vcf.Reader(open(stitch_path, 'r'))

		# Extract data from VCF using list comprehension
		vcf_data = [
			{
				'CHROM': record.CHROM,
				'POS': record.POS,
				'ID': record.ID,
				'REF': record.REF,
				'ALT': ",".join(map(str, record.ALT)),
				'QUAL': record.QUAL,
				'FILTER': '.',
				'INFO': '.',
				'FORMAT': 'GT',
				**{sample.sample: sample['GT'] for sample in record.samples}
			}
			for record in vcf_reader
		]

		# Create a Pandas DataFrame directly from the list of dictionaries
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

		count_all = (nugen_df == './.').sum().sum()
		count_all

		# Check for rows in stitch_df that are not in nugen_df and add them
		stitch_df_ids_not_in_nugen = set(stitch_df['ID']) - set(nugen_df['ID'])
		stitch_rows_to_add = stitch_df[stitch_df['ID'].isin(stitch_df_ids_not_in_nugen)].copy()
		stitch_rows_to_add

		nugen_df = pd.concat([nugen_df, stitch_rows_to_add], ignore_index=True)
		nugen_df.fillna('.', inplace=True)
		nugen_df.sort_values(by='POS', ascending=True, inplace=True)

		output_vcf_path = output.combined 

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

rule plink_combined:
	input:
		combined = rules.combine.output.combined 
	output:
		logfile = 'stitch/results/combine/chr{chrom}/plink.log', 
		imiss = 'stitch/results/combine/chr{chrom}/plink.imiss',
		lmiss = 'stitch/results/combine/chr{chrom}/plink.lmiss',
	threads: 4
	resources:
		time    = 20,
		mem_mb  = 24000,
		cpus    = 4
	params:
		outdir = 'stitch/results/combine/chr{chrom}/plink', 
	shell:
		'''
			plink \
			--missing \
			--horse \
			--out {params.outdir} \
			--vcf {input.combined}
		'''

rule compare_rate:
	input:
		combined_logs = sorted(expand(
		'stitch/results/combine/{chrom}/plink.log',
		chrom=['chr' + str(i) for i in range(1,31)] + ['chrX']
		)),
		nugen_logs = sorted(expand(
		'/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/{chrom}/plink.log',
		chrom=['chr' + str(i) for i in range(1,31)] + ['chrX']
		)),
	output:
		nugen_rate = '/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/nugen_rate.txt',
		combined_rate = 'stitch/results/combined_rate.txt', 
		compare_rates = 'stitch/results/compare_rates.txt',
	threads: 4
	resources:
		time    = 20,
		mem_mb  = 24000,
		cpus    = 4
	run:
		outfile = open(output.combined_rate, 'wt')

		for chrom in range(1,32):
			infile = open('stitch/results/combine/chr' + str(chrom) + '/plink.log', 'rt')
			for line in infile:
				line = line.rstrip()
				if 'genotyping rate' in line:
					split = line.split(' ')
					rate = split[4][2:4] + '.' + split[4][4]
					print('chr' + str(chrom) + '\t' + str(rate), file=outfile)
		infile = open('stitch/results/combine/chrX/plink.log', 'rt')
		for line in infile:
			line = line.rstrip()
			if 'genotyping rate' in line:
				split = line.split(' ')
				rate = split[4][2:4] + '.' + split[4][4]
				print('chrX' + '\t' + str(rate), file=outfile)

		outfile = open(output.nugen_rate, 'wt')

		for chrom in range(1,32):
			infile = open('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chr' + str(chrom) + '/plink.log', 'rt')
			for line in infile:
				line = line.rstrip()
				if 'genotyping rate' in line:
					split = line.split(' ')
					rate = split[4][2:4] + '.' + split[4][4]
					print('chr' + str(chrom) + '\t' + str(rate), file=outfile)
		infile = open('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/chrX/plink.log', 'rt')
		for line in infile:
			line = line.rstrip()
			if 'genotyping rate' in line:
				split = line.split(' ')
				rate = split[4][2:4] + '.' + split[4][4]
				print('chrX' + '\t' + str(rate), file=outfile)

		infile1 = open('/panfs/jay/groups/27/mccuem/dimml002/NuGEN/targets/split/nugen_rate.txt', 'rt')
		infile2 = open('stitch/results/combined_rate.txt', 'rt')
		outfile = open(output.compare_rates, 'wt')

		print('CHROM' + '\t' + 'nugen' + '\t' + 'combined', file=outfile)

		chroms = []
		nugen = []
		combined = []

		for line in infile1:
			line = line.rstrip()
			split = line.split('\t')
			chroms.append(split[0])
			nugen.append(split[1])

		for line in infile2:
			line = line.rstrip()
			split = line.split('\t')
			combined.append(split[1])

		for index,value in enumerate(nugen):
			print(chroms[index], nugen[index], combined[index],
					sep = '\t',
					file=outfile)

rule zip_index:
	input:
		combined = rules.combine.output.combined 
	output:
		zipped = 'stitch/results/combine/chr{chrom}/chr{chrom}.vcf.gz',
		tbi = 'stitch/results/combine/chr{chrom}/chr{chrom}.vcf.gz.tbi',
	threads: 4
	resources:
		time    = 20,
		mem_mb  = 24000,
		cpus    = 4
	shell:
		'''
			bgzip {input.combined}

			gatk IndexFeatureFile -I stitch/results/combine/chr{wildcards.chrom}/chr{wildcards.chrom}.vcf.gz

		'''

rule vcf_list:
	output:
		vcf_list = 'stitch/results/vcf.list',
	threads: 4
	resources:
		time    = 20,
		mem_mb  = 24000,
		cpus    = 4
	run:
		outfile = open(output.vcf_list, 'wt')

		for chrom in range(1,32):
			print('stitch/results/combine/chr' + str(chrom) + '/chr' + str(chrom) + '.vcf.gz', file=outfile)

		print('stitch/results/combine/chrX/chrX.vcf.gz', file=outfile)
	
rule MergeVcfs:
	input:
		zipped = expand('stitch/results/combine/chr{chrom}/chr{chrom}.vcf.gz', chrom=[str(i) for i in range(1, 32)] + ['X']),
		vcf_list = rules.vcf_list.output.vcf_list,
	output:
		full = 'stitch/results/full.combined.vcf.gz',
		full_tbi = 'stitch/results/full.combined.vcf.gz.tbi',
	threads: 4
	resources:
		time    = 20,
		mem_mb  = 24000,
		cpus    = 4
	shell:
		'''
			picard MergeVcfs \
				-I {input.vcf_list} \
				-O {output.full}

			gatk IndexFeatureFile -I {output.full}
		'''
