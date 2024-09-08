chroms = ['chr' + str(i) for i in range(1, 32)]
#chroms = ['chr1']

rule all:
    input:
        expand('combine_original_beagle/split/original/{chrom}/output.vcf.gz', chrom=chroms),
        expand('combine_original_beagle/split/original/{chrom}/reheader.vcf.gz', chrom=chroms),
        expand('combine_original_beagle/isec/{chrom}/0000.vcf.gz', chrom=chroms),
        expand('combine_original_beagle/isec/{chrom}/merge/output.vcf.gz', chrom=chroms),
        'combine_original_beagle/vcf_list.txt',
        'combine_original_beagle/final/original_and_beagle.vcf.gz',

rule split_original:
    input:
        original = 'original_cleaned.vcf.gz',
    output:
        split_original = 'combine_original_beagle/split/original/{chrom}/output.vcf.gz',
    resources:
        time    = 60,
        mem_mb    = 24000,
    shell:
        '''
            bcftools view \
            -r {wildcards.chrom} \
            -o {output.split_original} \
            {input.original}

            tabix {output.split_original}

        '''

rule reheader:
    input:
        split_original = 'combine_original_beagle/split/original/{chrom}/output.vcf.gz',
    output:
        reheadered_vcf = 'combine_original_beagle/split/original/{chrom}/reheader.vcf.gz',
    resources:
        time    = 60,
        mem_mb    = 24000,
    params:
        contig_file = 'EquCab3_nucl_.fasta.truncated.fai'
    shell:
        '''
            bcftools reheader -f {params.contig_file} -o {output.reheadered_vcf} {input.split_original} 

            tabix {output.reheadered_vcf}

        '''

rule split_imputed:
    input:
        imputed = 'results/imputation_only_results/FINAL/imputed.reheadered.vcf.gz',
    output:
        split_imputed = 'combine_original_beagle/split/imputed/{chrom}/output.vcf.gz',
    resources:
        time    = 60,
        mem_mb    = 24000,
    shell:
        '''
            bcftools view \
            -r {wildcards.chrom} \
            -o {output.split_imputed} \
            {input.imputed}

            tabix {output.split_imputed}

        '''
rule isec:
    input:
        reheadered_vcf = 'combine_original_beagle/split/original/{chrom}/reheader.vcf.gz',
        split_imputed = 'combine_original_beagle/split/imputed/{chrom}/output.vcf.gz',
    output:
        isec_vcf1 = 'combine_original_beagle/isec/{chrom}/0000.vcf.gz',
        isec_vcf2 = 'combine_original_beagle/isec/{chrom}/0003.vcf.gz',
    resources:
        time    = 60,
        mem_mb  = 24000,
    shell:
        '''
            bcftools isec \
            -p combine_original_beagle/isec/{wildcards.chrom} \
            {input.reheadered_vcf} \
            {input.split_imputed} \
            -Oz

        '''

rule merge:
    input:
        shared_beagle = 'combine_original_beagle/isec/{chrom}/0003.vcf.gz',
        unique_original = 'combine_original_beagle/isec/{chrom}/0000.vcf.gz',
    output:
        merged_vcf = 'combine_original_beagle/isec/{chrom}/merge/output.vcf.gz',
    resources:
        time    = 60,
        mem_mb  = 24000,
    shell:
        '''
            java -jar /panfs/jay/groups/27/mccuem/dimml002/picard.jar MergeVcfs \
            -I {input.shared_beagle} \
            -I {input.unique_original} \
            -O {output.merged_vcf}

        '''

rule vcf_list:
    input:
        merged_vcfs = expand('combine_original_beagle/isec/{chrom}/merge/output.vcf.gz', chrom=chroms)
    output:
        vcf_list = 'combine_original_beagle/vcf_list.txt',
    resources:
        time    = 60,
        mem_mb  = 24000,
    run:
        outfile = open('combine_original_beagle/vcf_list.txt', 'wt')
        for chrom in range(1,32):
            print('combine_original_beagle/isec/chr' + str(chrom) + '/merge/output.vcf.gz',
                    file = outfile)

rule merge_final:
    input:
        vcf_list = 'combine_original_beagle/vcf_list.txt',
    output:
        final_vcf = 'combine_original_beagle/final/original_and_beagle.vcf.gz',
    resources:
        time    = 60,
        mem_mb  = 24000,
    shell:
        '''
            bcftools concat \
            -Oz -o {output.final_vcf} \
            -f {input.vcf_list}

            tabix {output.final_vcf}

        '''
