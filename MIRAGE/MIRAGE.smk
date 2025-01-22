import gzip
import random

#chroms = ['12', '20', '23', '27']
#chroms = ['27', '12', '23']
chroms = [str(i) for i in range(1,32)]
#chroms = ['1']

itrs = ['set' + str(i) for i in range(1, config['iterations'] + 1)]
#itrs = ['set1']

#bins_list = {"0.10_0.15": (0.1, 0.15),
#            "0.15_0.20": (0.15, 0.20)
#}

#bins_list = {
#    "0.01_0.02": (0.01, 0.02),
#    "0.02_0.05": (0.02, 0.05),
#    "0.05_0.10": (0.05, 0.10),
#    "0.10_0.15": (0.10, 0.15),
#    "0.15_0.20": (0.15, 0.20),
#    "0.20_0.25": (0.20, 0.25),
#    "0.25_0.30": (0.25, 0.30),
#    "0.30_0.35": (0.30, 0.35),
#    "0.35_0.40": (0.35, 0.40),
#    "0.40_0.45": (0.40, 0.45),
#    "0.45_0.50": (0.45, 0.50)
#}


bins_list = {
    "0_0.05": (0, 0.05),
    "0.05_0.15": (0.05, 0.15),
    "0.15_0.25": (0.15, 0.25),
    "0.25_0.35": (0.25, 0.35),
    "0.35_0.50": (0.35, 0.50)
}

bins = bins_list.keys()

rule all:
    input:
        expand('0.1_masking_results/split/original/chr{chrom}/chr{chrom}.vcf.gz', chrom=chroms),
        expand('0.1_masking_results/split/original/chr{chrom}/chr{chrom}.vcf.gz.tbi', chrom=chroms),
        expand('0.1_masking_results/split/original/prepped/chr{chrom}/chr{chrom}.no_missing.vcf.gz', chrom=chroms),
        expand('0.1_masking_results/split/original/prepped/chr{chrom}/chr{chrom}.no_missing.vcf.gz.tbi', chrom=chroms),
        expand('0.1_masking_results/split/original/prepped/chr{chrom}/chr{chrom}.no_missing.IDs.vcf.gz', chrom=chroms),
#        expand('0.1_masking_results/masking/maf_bins/chr{chrom}/{bin}/full/output.vcf.gz',chrom=chroms, bin=bins),
#        expand('0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/to_mask/output.vcf.gz', chrom=chroms, bin=bins, itr=itrs),
#        expand('0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/masked/output.vcf.gz', chrom=chroms, bin=bins, itr=itrs),
#        expand('0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/combine/isec/0000.vcf', chrom=chroms, bin=bins, itr=itrs),
#        expand('0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/vcf_list.txt', chrom=chroms, itr=itrs),
#        expand('0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/combine/combined.vcf.gz', chrom=chroms, itr=itrs),
        expand('0.1_masking_results/masking/{itr}/chr{chrom}/combine/combined.MAF.vcf.gz', chrom=chroms, itr=itrs),
#        expand('0.1_masking_results/imputed/{itr}/chr{chrom}/imputed.chr{chrom}.vcf.gz', chrom=chroms, itr=itrs),
        expand('0.1_masking_results/imputed/{itr}/chr{chrom}/imputed.chr{chrom}.MAF.vcf.gz', chrom=chroms, itr=itrs),
        expand('0.1_masking_results/FINAL/{itr}/imputed.vcf.gz', itr = itrs),
#        expand('0.1_masking_results/imputed/{itr}/imputed.list', itr=itrs),
#        expand('0.1_masking_results/FINAL/{itr}/imputed.vcf.gz', itr=itrs),
        expand('0.1_masking_results/mask_assess/{itr}/intersect/masked/chr{chrom}/combined.MAF.vcf.gz', chrom=chroms, itr=itrs),
        expand('0.1_masking_results/mask_assess/{itr}/chr{chrom}/imputed.table', itr=itrs, chrom=chroms),
        expand('0.1_masking_results/mask_assess/{itr}/full_imputed.table', itr=itrs),
        '0.1_masking_results/mask_assess/header.txt'
#        expand('0.1_masking_results/mask_assess/{itr}/chr{chrom}/mismatched_genos.txt', itr=itrs, chrom=chroms),

rule split_original:
    input:
        original_vcf = config['vcf_to_mask'],
        original_vcf_tbi = config['vcf_to_mask_tbi'],
    output:
        original_vcf_split = '0.1_masking_results/split/original/chr{chrom}/chr{chrom}.vcf.gz',
        original_vcf_split_tbi = '0.1_masking_results/split/original/chr{chrom}/chr{chrom}.vcf.gz.tbi',
    resources:
        time    = 30,
        mem_mb  = 4000,
    shell:
            '''
                bcftools view \
                    -r chr{wildcards.chrom} \
                    -Oz -o {output.original_vcf_split} \
                    {input.original_vcf}

                gatk IndexFeatureFile -I {output.original_vcf_split}
            '''

rule prep_original:
    input:
        original_vcf_split = '0.1_masking_results/split/original/chr{chrom}/chr{chrom}.vcf.gz',
        original_vcf_split_tbi = '0.1_masking_results/split/original/chr{chrom}/chr{chrom}.vcf.gz.tbi'
    output:
        original_vcf_prepped = '0.1_masking_results/split/original/prepped/chr{chrom}/chr{chrom}.no_missing.vcf.gz',
        original_vcf_prepped_tbi = '0.1_masking_results/split/original/prepped/chr{chrom}/chr{chrom}.no_missing.vcf.gz.tbi',
        #original_vcf_cleaned = '0.1_masking_results/split/original/prepped/chr{chrom}/chr{chrom}.cleaned.vcf.gz',
        #original_vcf_cleaned_tbi = '0.1_masking_results/split/original/prepped/chr{chrom}/chr{chrom}.cleaned.vcf.gz.tbi',
    resources:
        time = 60,
        mem_mb = 12000,
    params:
        directory = '0.1_masking_results/split/original/prepped'
    shell:
        r'''
        
            set -euo pipefail

           # echo "Starting bcftools view -f PASS"
           # bcftools view -f PASS {input.original_vcf_split} -Oz -o {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.temp.vcf.gz
           # echo "Finished bcftools view -f PASS"

           # echo "Starting bcftools norm"
           # bcftools norm -m -both {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.temp.vcf.gz -Oz -o {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.temp2.vcf.gz
           # echo "Finished bcftools norm"

            echo "Starting bcftools +fill-tags"
            bcftools +fill-tags {input.original_vcf_split} -- -t MAF | bcftools view -Oz -o {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.temp.vcf.gz
            echo "Finished bcftools +fill-tags"

           # echo "Starting bcftools annotate"
           # bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.temp3.vcf.gz -o {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.cleaned.vcf.gz
           # echo "Finished bcftools annotate"

           # echo "Starting GATK IndexFeatureFile for cleaned VCF"
           # gatk IndexFeatureFile -I {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.cleaned.vcf.gz
           # echo "Finished GATK IndexFeatureFile for cleaned VCF"

            echo "Starting bcftools view -g ^miss"
            bcftools filter -e 'AN < 3564' {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.temp.vcf.gz -Oz -o {output.original_vcf_prepped} 
            echo "Finished bcftools view -g ^miss"

            #echo "Starting bcftools view -i 'MAF>=0.01'"
            #bcftools view -i 'MAF>=0.01' {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.temp2.vcf.gz -Oz -o {output.original_vcf_prepped}
            #echo "Finished bcftools view -i 'MAF>=0.01'"

            echo "Starting gatk IndexFeatureFile for prepped VCF"
            gatk IndexFeatureFile -I {output.original_vcf_prepped}
            echo "Finished GATK IndexFeatureFile for prepped VCF"

            echo "Cleaning up temporary files"
            rm -f {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.temp.vcf.gz
            rm -f {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.temp2.vcf.gz
           # rm -f {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.temp3.vcf.gz
            #rm -f {params.directory}/chr{wildcards.chrom}/chr{wildcards.chrom}.temp4.vcf.gz
            echo "Cleanup complete"

        '''

rule annotate_original:
    input:
        original_vcf_prepped = '0.1_masking_results/split/original/prepped/chr{chrom}/chr{chrom}.no_missing.vcf.gz',
    output:
        original_annotated = '0.1_masking_results/split/original/prepped/chr{chrom}/chr{chrom}.no_missing.IDs.vcf.gz',
    resources:
        time = 60,
        mem_mb = 12000,
    params:
        directory = '0.1_masking_results/split/original/prepped/chr{chrom}'
    shell:
            '''
                bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o {output.original_annotated} {input.original_vcf_prepped} \

                gatk IndexFeatureFile -I {output.original_annotated} 

            '''


rule filter_maf_bins:
    input:
        original_vcf_prepped = '0.1_masking_results/split/original/prepped/chr{chrom}/chr{chrom}.no_missing.vcf.gz',
    output:
        maf_bin_vcf = '0.1_masking_results/masking/maf_bins/chr{chrom}/{bin}/full/output.vcf.gz',  
    params:
        maf_min = lambda wildcards: bins_list[wildcards.bin][0],
        maf_max = lambda wildcards: bins_list[wildcards.bin][1]
    shell:
        '''
        
            bcftools view -i 'MAF>={params.maf_min} & MAF<{params.maf_max}' {input.original_vcf_prepped} -Oz -o {output.maf_bin_vcf}

            tabix {output.maf_bin_vcf}
        
        '''


rule isolate_random_variants:
    input:
        maf_bin_vcf = '0.1_masking_results/masking/maf_bins/chr{chrom}/{bin}/full/output.vcf.gz',
    output:
        random_selection_vcf = '0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/to_mask/output.vcf.gz',
    params:
        directory = '0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/to_mask',
    resources:
        time    = 30,
        mem_mb  = 36000,
    shell:
        '''
        
            set -euo pipefail

            echo "Counting total variants"
            total_variants=$(bcftools view -H {input.maf_bin_vcf} | wc -l)
            echo "Total variants: $total_variants"

            echo "Calculating 30% of total variants"
            thirty_percent=$(echo "($total_variants * 0.3)/1" | bc)
            echo "30% of total variants: $thirty_percent"

            # Extract the header
            echo "Extracting header"
            bcftools view -h {input.maf_bin_vcf} > {params.directory}/temp_header.vcf
            echo "Header extracted"

            # Extract 30% of variants and sort them
            echo "Extracting and sorting 30% of variants"
            bcftools view -H {input.maf_bin_vcf} | shuf -n $thirty_percent | sort -k1,1 -k2,2n > {params.directory}/temp_variants.vcf
            echo "10% of variants extracted and sorted"

            # Combine header and sorted variants, and compress the result
            echo "Combining header and sorted variants"
            cat {params.directory}/temp_header.vcf {params.directory}/temp_variants.vcf | bgzip -c > {output.random_selection_vcf}
            echo "Combined and compressed"

            # Index the resulting VCF
            echo "Indexing the resulting VCF"
            bcftools index {output.random_selection_vcf}
            echo "Indexing complete"

            # Clean up temporary files
            echo "Cleaning up temporary files"
            rm {params.directory}/temp_header.vcf {params.directory}/temp_variants.vcf
            echo "Cleanup complete"
        
        '''

#rule get_clean_list:
#    input:
#        original_vcf_cleaned = expand('0.1_masking_results/split/original/prepped/chr{chrom}.cleaned.vcf.gz', chrom=chroms),
#    output:
#        cleaned_vcf_list = '0.1_masking_results/split/original/prepped/cleaned_vcf_list.txt',
#    resources:
#        time = 20,
#        mem_mb = 24000,
#    run:
#        outfile = open('0.1_masking_results/split/original/prepped/cleaned_vcf_list.txt', 'wt')
#
#        for num in range(1,32):
#            print('0.1_masking_results/split/original/prepped/chr' + str(num) + '.cleaned.vcf.gz', file=outfile) 
#
#rule concat_cleaned:
#    input:
#        cleaned_vcf_list = '0.1_masking_results/split/original/prepped/cleaned_vcf_list.txt',
#    output:
#        merged_cleaned_vcf = '0.1_masking_results/split/original/prepped/merged.cleaned.vcf.gz',
#        merged_cleaned_vcf_tbi = '0.1_masking_results/split/original/prepped/merged.cleaned.vcf.gz.tbi',
#    resources:
#        time = 60,
#        mem_mb = 24000,
#    shell:
#       '''
#            gatk MergeVcfs -I {input.cleaned_vcf_list} -O {output.merged_cleaned_vcf}
#
#            gatk IndexFeatureFile -I {output.merged_cleaned_vcf}
#        '''

rule mask:
    input:
        random_selection_vcf = '0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/to_mask/output.vcf.gz',
    output:
        masked_vcf = '0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/masked/output.vcf.gz',
        masked_tbi = '0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/masked/output.vcf.gz.tbi',
        #masked_vcf = '0.1_masking_results/masking/{itr}/chr{chrom}/masked.bin.{bin}.vcf.gz',
        #masked_tbi = '0.1_masking_results/masking/{itr}/chr{chrom}/masked.bin.{bin}.vcf.gz.tbi',
    resources:
        time    = 30,
        mem_mb  = 4000,
#    conda:
#        'mask.yaml' #this command for setGT only works with bcftools 1.17
    params:
        seed = lambda wildcards, input: random.randint(1, 100),
        mask_proportion = config['mask_proportion'],
        bcftools_directory = config['bcftools_directory']
    shell:
            '''
                {params.bcftools_directory} +setGT -o {output.masked_vcf} {input.random_selection_vcf} -- -t r:{params.mask_proportion} -s {params.seed} -n . 
                
                gatk IndexFeatureFile -I {output.masked_vcf}

            '''

rule isec:
    input:
        maf_bin_vcf = '0.1_masking_results/masking/maf_bins/chr{chrom}/{bin}/full/output.vcf.gz',
        masked_vcf = '0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/masked/output.vcf.gz',
    output:
        isec_vcfs = '0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/combine/isec/0000.vcf',
    resources:
        time = 60,
        mem_mb = 4000,
    params:
        directory = '0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/combine/isec/'
    shell:
       '''
            bcftools isec -p {params.directory} {input.maf_bin_vcf} {input.masked_vcf}
        '''

rule masked_list:
    input:
        isec_vcfs = expand('0.1_masking_results/masking/maf_bins/{itr}/chr{chrom}/{bin}/combine/isec/0000.vcf', chrom=chroms, bin=bins, itr=itrs),
    output:
        vcf_list = '0.1_masking_results/masking/{itr}/chr{chrom}/vcf_list.txt',
    resources:
        time = 10,
        mem_mb = 4000,
    run: 
        # Define the output file
        output_file = Path(output.vcf_list)

        # Create a list to hold the VCF paths
        vcf_paths = []

        # Loop through each bin and VCF file pattern
        for bin in bins:
            for vcf in ["0000.vcf", "0003.vcf"]:
                # Construct the VCF path using Pathlib
                vcf_path = "0.1_masking_results/masking/maf_bins/" + wildcards.itr + "/chr" + wildcards.chrom + "/" + bin + "/combine/isec/" + vcf
                vcf_paths.append(vcf_path)

        # Write the VCF paths to the output file
        with output_file.open("w") as f:
            for path in vcf_paths:
                f.write(path + "\n")

rule combine_masked:
    input:
        vcf_list = '0.1_masking_results/masking/{itr}/chr{chrom}/vcf_list.txt',
    output:
        combined_masked_vcf = '0.1_masking_results/masking/{itr}/chr{chrom}/combine/combined.vcf.gz',
    resources:
        time = 60,
        mem_mb = 12000,
    shell:
            '''
                java -jar /panfs/jay/groups/27/mccuem/dimml002/picard.jar MergeVcfs \
                    -I {input.vcf_list} \
                    -O {output.combined_masked_vcf}

                gatk IndexFeatureFile -I {output.combined_masked_vcf}
            '''

rule annotate_masked:
    input:
        combined_masked_vcf = '0.1_masking_results/masking/{itr}/chr{chrom}/combine/combined.vcf.gz',
    output:
        annotated_masked = '0.1_masking_results/masking/{itr}/chr{chrom}/combine/combined.MAF.vcf.gz',
    params:
        directory = '0.1_masking_results/masking/{itr}/chr{chrom}/combine'
    resources:
        time = 60,
        mem_mb = 12000,
    shell:
            '''
                bcftools +fill-tags {input.combined_masked_vcf} -Oz -o {params.directory}/temp.vcf.gz -- -t MAF,AC \

                bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o {output.annotated_masked} {params.directory}/temp.vcf.gz \

                gatk IndexFeatureFile -I {output.annotated_masked} \

                rm {params.directory}/temp.vcf.gz
            '''

rule split_phased:
    input:
        reference_vcf = config['reference_panel'],
        reference_vcf_tbi = config['reference_panel_tbi'] 
    output:
        reference_vcf_split = '0.1_masking_results/split/ref/chr{chrom}/chr{chrom}.vcf.gz',
        reference_vcf_split_tbi = '0.1_masking_results/split/ref/chr{chrom}/chr{chrom}.vcf.gz.tbi'
    resources:
        time    = 360,
        mem_mb  = 12000,
    shell:
            '''
                bcftools view \
                -r chr{wildcards.chrom} \
                -Oz -o {output.reference_vcf_split} \
                {input.reference_vcf}

                gatk IndexFeatureFile -I {output.reference_vcf_split}
            '''

rule beagle_impute:
    input:
        combined_masked_vcf = '0.1_masking_results/masking/{itr}/chr{chrom}/combine/combined.MAF.vcf.gz', 
        reference_vcf = '0.1_masking_results/split/ref/chr{chrom}/chr{chrom}.vcf.gz', 
        reference_vcf_tbi = '0.1_masking_results/split/ref/chr{chrom}/chr{chrom}.vcf.gz.tbi',
    output:
        imputed_vcf = '0.1_masking_results/imputed/{itr}/chr{chrom}/imputed.chr{chrom}.vcf.gz',
        imputed_tbi = '0.1_masking_results/imputed/{itr}/chr{chrom}/imputed.chr{chrom}.vcf.gz.tbi'
    params:
        prefix = lambda wildcards, output: output.imputed_vcf.rsplit('.',2)[0],
        breed = config['breed'],
        breed_ab = config['breed_ab'],
    threads: 8
    resources:
        time    = 360,
        mem_mb  = 24000,
    shell:
        '''
            java -jar beagle.27Jan18.7e1.jar \
            gtgl={input.combined_masked_vcf} \
            chrom=chr{wildcards.chrom} \
            map=recombination_maps/beagle_maps/BEAGLE_{params.breed}_ECA{wildcards.chrom}_map.txt \
            ref={input.reference_vcf} \
            ne=8460 \
            nthreads=8 \
            gprobs=true \
            impute=false \
            out={params.prefix}

            tabix {output.imputed_vcf} 
        '''

rule annotate_imputed:
    input:
        imputed_vcf = '0.1_masking_results/imputed/{itr}/chr{chrom}/imputed.chr{chrom}.vcf.gz',
    output:
        imputed_vcf_maf = '0.1_masking_results/imputed/{itr}/chr{chrom}/imputed.chr{chrom}.MAF.vcf.gz',
        imputed_tbi_maf = '0.1_masking_results/imputed/{itr}/chr{chrom}/imputed.chr{chrom}.MAF.vcf.gz.tbi',
    params:
        directory = '0.1_masking_results/imputed/{itr}/chr{chrom}'
    resources:
        time    = 20,
        mem_mb  = 4000,
    shell:
            '''
                bcftools +fill-tags {input.imputed_vcf} -Oz -o {params.directory}/temp.vcf.gz -- -t MAF,AC \

                bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o {output.imputed_vcf_maf} {params.directory}/temp.vcf.gz \

                gatk IndexFeatureFile -I {output.imputed_vcf_maf} \

                rm {params.directory}/temp.vcf.gz
            '''
rule imputed_list:
    input:
        imputed_vcfs = sorted(expand(
        '0.1_masking_results/imputed/{itr}/chr{chrom}/imputed.chr{chrom}.MAF.vcf.{ext}',
        chrom=[str(i) for i in range(1,31)], # NO CHROM M CORRECT?
        ext=['gz', 'gz.tbi'],
        itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]
        ))
    output:
        sorted_list = '0.1_masking_results/imputed/{itr}/imputed.list',
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    run:
        # drop indices from input
        for itr in range(1, config['iterations'] + 1):
            outfile = open('0.1_masking_results/imputed/set' + str(itr) + '/imputed.list', 'wt')
            for chrom in range(1,32):
                print('0.1_masking_results/imputed/set' + str(itr) + '/chr' + str(chrom) + '/imputed.chr' + str(chrom) + '.vcf.gz',
                    file = outfile)

rule combine_imputed:
    input:
        sorted_list = '0.1_masking_results/imputed/{itr}/imputed.list', 
        imputed_vcfs = sorted(expand(
        '0.1_masking_results/imputed/{itr}/chr{chrom}/imputed.chr{chrom}.vcf.{ext}', 
        chrom=[str(i) for i in range(1,31)], # NO CHROM M CORRECT?
        ext=['gz', 'gz.tbi'],
        itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]
        ))
    output:
        sorted_vcf = '0.1_masking_results/FINAL/{itr}/imputed.vcf.gz',
        sorted_tbi = '0.1_masking_results/FINAL/{itr}/imputed.vcf.gz.tbi' 
    threads: 4
    resources:
        time   = 120,
        mem_mb = 4000
    shell:
        '''
            bcftools concat \
            -Oz -o {output.sorted_vcf} \
            -f {input.sorted_list}

            gatk IndexFeatureFile -I {output.sorted_vcf}
        '''

rule get_intersect:
    input:
        imputed_vcf_maf = '0.1_masking_results/imputed/{itr}/chr{chrom}/imputed.chr{chrom}.MAF.vcf.gz',
        annotated_masked = '0.1_masking_results/masking/{itr}/chr{chrom}/combine/combined.MAF.vcf.gz',
        original_annotated = '0.1_masking_results/split/original/prepped/chr{chrom}/chr{chrom}.no_missing.IDs.vcf.gz',
    output:
        masked_intersect = '0.1_masking_results/mask_assess/{itr}/intersect/masked/chr{chrom}/combined.MAF.vcf.gz',
        original_intersect = '0.1_masking_results/mask_assess/{itr}/intersect/original/chr{chrom}/chr{chrom}.no_missing.IDs.vcf.gz',
    params:
        imp_d = '0.1_masking_results/imputed/{itr}/chr{chrom}',
    resources:
        time   = 60,
        mem_mb = 4000
    shell:
        '''
            bcftools query -f "%ID\n" {input.imputed_vcf_maf} > {params.imp_d}/imputed_IDs.txt \

            bcftools view --include ID=@{params.imp_d}/imputed_IDs.txt -Oz -o {output.masked_intersect} {input.annotated_masked} \

            tabix {output.masked_intersect} 
            
            bcftools view --include ID=@{params.imp_d}/imputed_IDs.txt -Oz -o {output.original_intersect} {input.original_annotated} \

            tabix {output.original_intersect}

        '''

rule to_table:
    input:
        masked_intersect = '0.1_masking_results/mask_assess/{itr}/intersect/masked/chr{chrom}/combined.MAF.vcf.gz',
        original_intersect = '0.1_masking_results/mask_assess/{itr}/intersect/original/chr{chrom}/chr{chrom}.no_missing.IDs.vcf.gz',
        imputed_vcf = '0.1_masking_results/imputed/{itr}/chr{chrom}/imputed.chr{chrom}.MAF.vcf.gz',
    output:
        imputed_table = '0.1_masking_results/mask_assess/{itr}/chr{chrom}/imputed.table',
        original_table = '0.1_masking_results/mask_assess/{itr}/chr{chrom}/original.table',
        masked_table = '0.1_masking_results/mask_assess/{itr}/chr{chrom}/masked.table',
    resources:
        time   = 30,
        mem_mb = 4000
    shell:
        '''
             bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%MAF\t%AC\t%QUAL\t%FILTER\t[%DS\t]\n' -H {input.imputed_vcf}  > {output.imputed_table}

             bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%MAF\t%AC\t%QUAL\t%FILTER\t[%GT\t]\n' -H {input.original_intersect}  > {output.original_table}

             bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%MAF\t%AC\t%QUAL\t%FILTER\t[%GT\t]\n' -H {input.masked_intersect}  > {output.masked_table}
    
        '''

rule get_header:
    input:
        imputed_table = '0.1_masking_results/mask_assess/set1/chr1/imputed.table' 
    output:
        header = '0.1_masking_results/mask_assess/header.txt'
    resources:
        time   = 30,
        mem_mb = 12000
    run:

        infile = open(input.imputed_table, 'rt')

        header = infile.readline().rstrip()

        split = header.split('\t')

        new_header = []

        for x in split:
            if x[0] == '#':
                new_header.append('CHROM')
            elif x[0] == '[' and 'DS' not in x:
                new_header.append(x.split(']')[1])
            else:
                new_header.append(x.split(']')[1].split(':')[0])

        # Output file path
        output_file = output.header

        with open(output_file, 'w') as f:
            f.write('\t'.join(new_header) + '\n')

rule combine_tables:
    input:
        imputed_tables = expand('0.1_masking_results/mask_assess/{itr}/chr{chrom}/imputed.table', itr=itrs, chrom=chroms),
        original_tables = expand('0.1_masking_results/mask_assess/{itr}/chr{chrom}/original.table', itr=itrs, chrom=chroms),
        masked_tables = expand('0.1_masking_results/mask_assess/{itr}/chr{chrom}/masked.table', itr=itrs, chrom=chroms),
        header = '0.1_masking_results/mask_assess/header.txt',
    output:
        imputed_full = '0.1_masking_results/mask_assess/{itr}/full_imputed.table',
        original_full = '0.1_masking_results/mask_assess/{itr}/full_original.table',
        masked_full = '0.1_masking_results/mask_assess/{itr}/full_masked.table',
    resources:
        time   = 30,
        mem_mb = 12000
    shell:
        '''
            # Combine imputed tables: Use the generated header as the first line
            cat {input.header} > {output.imputed_full}
            for file in {input.imputed_tables}; do
                tail -n +2 "$file" >> {output.imputed_full}
            done

            # Combine original tables: Use the generated header as the first line
            cat {input.header} > {output.original_full}
            for file in {input.original_tables}; do
                tail -n +2 "$file" >> {output.original_full}
            done

            # Combine masked tables: Use the generated header as the first line
            cat {input.header} > {output.masked_full}
            for file in {input.masked_tables}; do
                tail -n +2 "$file" >> {output.masked_full}
            done 
        '''
rule compare_genos:
    input:
        imputed_full = '0.1_masking_results/mask_assess/{itr}//full_imputed.table',
        original_full = '0.1_masking_results/mask_assess/{itr}/full_original.table',
        masked_full = '0.1_masking_results/mask_assess/{itr}/full_masked.table',
    output:
        mismatched_genos = '0.1_masking_results/mask_assess/{itr}/mismatched_genos.txt',
        remove_imputed = '0.1_masking_results/mask_assess/{itr}/remove_imputed.list',
        readme = '0.1_masking_results/mask_assess/{itr}/README.txt',
    resources:
        time   = 360,
        mem_mb = 24000
    run:
        infile1 = open(input.original_full, 'rt')
        infile2 = open(input.masked_full, 'rt')
        infile3 = open(input.imputed_full, 'rt')

        outfile1 = open(output.mismatched_genos, 'wt')
        outfile2 = open(output.remove_imputed, 'wt')
        outfile3 = open(output.readme, 'wt')

        #Dictionary to keep track of 
        org = {}

        total_genos = 0

        line = infile1.readline()
        samples_temp = line.split('\t')
        samples = []

        for sample in samples_temp:
            if '.' in sample:
                sample = sample.replace('.GT', '')
                samples.append(sample)

        sample_num = len(samples)
        ind = sample_num + 1

        total_variants = 0
        for line in infile1:
            line = line.rstrip()
            total_variants += 1
            split = line.split('\t')
            pos = split[0] + ':' + str(split[1])
            org[pos] = {}
            for i in range(3,ind):
                sample_ind = i - 3
                geno = split[i].replace('|', '/')
                al1,al2 = geno.split('/')
                alleles = []
                alleles.append(al1)
                alleles.append(al2)
                alleles.sort()
                geno = alleles[0] + '/' + alleles[1]
                org[pos][samples[sample_ind]] = geno
                total_genos += 1

        masked = {}

        line = infile2.readline()

        samples_temp = line.split('\t')
        samples = []

        for sample in samples_temp:
            if '.' in sample:
                sample = sample.replace('.GT', '')
                samples.append(sample)

        masked_genos = 0

        for line in infile2:
            line = line.rstrip()
            split = line.split('\t')
            pos = split[0] + ':' + str(split[1])
            if './.' in line:
                masked[pos] = {}
                for i in range(3,ind):
                    sample_ind = i - 3
                    if split[i] == './.':
                        split[i] = split[i].replace('|', '/')
                        masked[pos][samples[sample_ind]] = split[i]
                        masked_genos += 1

        imputed = {}

        line = infile3.readline()
        samples_temp = line.split('\t')
        samples = []

        for sample in samples_temp:
            if '.' in sample:
                sample = sample.replace('.GT', '')
                samples.append(sample)

        imputed_genos = 0

        for line in infile3:
            line = line.rstrip()
            split = line.split('\t')
            pos = split[0] + ':' + str(split[1])
            imputed[pos] = {}
            for i in range(3,ind):
                sample_ind = i - 3
                geno = split[i].replace('|', '/')
                al1,al2 = geno.split('/')
                alleles = []
                alleles.append(al1)
                alleles.append(al2)
                alleles.sort()
                geno = alleles[0] + '/' + alleles[1]
                imputed[pos][samples[sample_ind]] = geno
                imputed_genos += 1

        correct_snps = 0
        incorrect_snps = 0
        correct_indels = 0
        incorrect_indels = 0
        missing_variants = 0

        for variant,value in imputed.items():
            for sample,geno in value.items():
                if variant in org:
                    if org[variant][sample] == geno:
                        if len(geno) > 3:
                            correct_indels += 1
                        else:
                            correct_snps += 1
                    else:
                        if len(geno) > 3:
                            incorrect_indels += 1
                        else:
                            incorrect_snps += 1
                        print(variant, sample + ':' + org[variant][sample], sample + ':' + imputed[variant][sample],
                                sep='\t',
                                file=outfile1)
                else:
                    missing_variants += 1

            for variant,value in org.items():
                if variant in imputed:
                    chrom,pos = variant.split(':')
                    print(chrom + '\t' + pos, file=outfile2)


        missing_genos = total_genos - imputed_genos

        print('Total variants: ' + str(total_variants), file=outfile3)
        print('Total genotypes: ' + str(total_genos), file=outfile3)
        print('Masked genotypes: ' + str(masked_genos), file=outfile3)
        print('Imputed genotypes: ' + str(imputed_genos), file=outfile3)
        #print('Correctly imputed genotypes: ' + str(correct), file=outfile3)
        #print('Incorrectly imputed genotypes: ' + str(incorrect), file=outfile3)
        print('Missing variants: ' + str(missing_variants), file=outfile3)
        print('Missing genotype: ' + str(missing_genos) + '\n', file=outfile3)

        print('Percent Imputed: ' + str(imputed_genos/total_genos*100), file=outfile3)
        print('SNP Accuracy: ' + str(correct_snps/(correct_snps+incorrect_snps)*100), file=outfile3)
        print('Indel Accuracy: ' + str(correct_indels/(correct_indels+incorrect_indels)*100), file=outfile3)


rule remove_imputed:
    input:
        original_vcf = config['vcf_to_mask'],
        remove_imputed = '0.1_masking_results/mask_assess/{itr}/chr{chrom}/remove_imputed.list',
    output:
        no_imputed = '0.1_masking_results/mask_assess/{itr}/chr{chrom}/no_imputed.vcf.gz',
        no_imputed_tbi = '0.1_masking_results/mask_assess/{itr}/chr{chrom}/no_imputed.vcf.gz.tbi',
    threads: 4
    resources:
        time   = 30,
        mem_mb = 24000
    shell:
        '''
            bcftools view \
                -Oz \
                -o {output.no_imputed} \
                --min-ac 1 \
                -T {input.remove_imputed} \
                {input.original_vcf}

            gatk IndexFeatureFile -I {output.no_imputed}
 
        '''

rule combine_vcfs:
    input:
        no_imputed_chrom = '0.1_masking_results/mask_assess/{itr}/chr{chrom}/no_imputed.vcf.gz',
        imputed_vcf_chrom = '0.1_masking_results/FINAL/masked_{itr}/chr{chrom}/imputed.vcf.gz',
    output:
        full_vcf_chrom = '0.1_masking_results/mask_assess/{itr}/chr{chrom}/full.vcf.gz',
        full_tbi_chrom = '0.1_masking_results/mask_assess/{itr}/chr{chrom}/full.vcf.gz.tbi',
    threads: 4
    resources:
        time   = 30,
        mem_mb = 24000
    shell:
        '''
            bcftools concat \
                -Oz \
                -o {output.full_vcf_chrom} \
                --allow-overlaps \
                {input.imputed_vcf_chrom} {input.no_imputed_chrom}

            gatk IndexFeatureFile -I {output.full_vcf_chrom}
 
        '''

#
#rule plink_missing:
#    input:
#        full_vcf = '0.1_masking_results/mask_assess/{itr}/full.vcf.gz', 
#    output:
#        imiss = '0.1_masking_results/mask_assess/{itr}/plink.imiss',
#        lmiss = '0.1_masking_results/mask_assess/{itr}/plink.lmiss',
#        nosex = '0.1_masking_results/mask_assess/{itr}/plink.nosex',
#        logfile = '0.1_masking_results/mask_assess/{itr}/plink.log',
#    threads: 4
#    resources:
#        time   = 30,
#        mem_mb = 24000
#    conda:
#        'plink.yaml'
#    params:
#        directory = '0.1_masking_results/mask_assess/{itr}/plink'
#    shell:
#        '''
#            plink \
#                --missing \
#                --horse \
#                --vcf {input.full_vcf} \
#                --out {params.directory}
#        '''
#
