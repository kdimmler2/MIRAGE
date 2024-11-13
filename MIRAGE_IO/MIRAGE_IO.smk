import gzip
import random

chroms = [str(i) for i in range(1,32)]
#chroms = ['1']

rule all:
    input:
#        expand('imputation_only_results/ref/split/chr{chrom}/chr{chrom}.vcf.gz', chrom=[str(i) for i in range(1,31)] + ['X']),
#        expand('imputation_only_results/ref/split/chr{chrom}/chr{chrom}.vcf.gz.tbi', chrom=[str(i) for i in range(1,31)] + ['X']),
        expand('imputation_only_results/split/chr{chrom}/chr{chrom}.vcf.gz', chrom=chroms),
        expand('imputation_only_results/split/chr{chrom}/chr{chrom}.vcf.gz.tbi', chrom=chroms),
        expand('imputation_only_results/split/prepped/chr{chrom}/chr{chrom}.vcf.gz', chrom=chroms),
        expand('imputation_only_results/imputed/chr{chrom}/imputed.chr{chrom}.vcf.gz', chrom=chroms),
        expand('imputation_only_results/imputed/chr{chrom}/imputed.chr{chrom}.vcf.gz.tbi', chrom=chroms),
        'imputation_only_results/imputed/imputed.list',
        'imputation_only_results/FINAL/imputed.vcf.gz',
        'imputation_only_results/FINAL/imputed.vcf.gz.tbi',


rule split_phased:
    input:
        reference_panel = config['reference_panel'],
        reference_panel_tbi = config['reference_panel_tbi'] 
    output:
        split_vcf = 'imputation_only_results/ref/split/chr{chrom}/chr{chrom}.vcf.gz',
        split_tbi = 'imputation_only_results/ref/split/chr{chrom}/chr{chrom}.vcf.gz.tbi'
    threads: 4
    resources:
        time    = 1440,
        mem_mb  = 60000,
        cpus    = 4
    shell:
        '''
            bcftools view \
            -r chr{wildcards.chrom} \
            -Oz -o {output.split_vcf} \
            {input.reference_panel}

            gatk IndexFeatureFile -I {output.split_vcf}
        '''

rule split_original:
    input:
        original_vcf = config['original_vcf'],
        original_vcf_tbi = config['original_vcf_tbi'],
    output:
        original_split_vcf = 'imputation_only_results/split/chr{chrom}/chr{chrom}.vcf.gz',
        original_split_tbi = 'imputation_only_results/split/chr{chrom}/chr{chrom}.vcf.gz.tbi',
    threads: 4
    resources:
        time    = 1440,
        mem_mb  = 60000,
        cpus    = 4
    shell:
        '''
            bcftools view \
            -r chr{wildcards.chrom} \
            -Oz -o {output.original_split_vcf} \
            {input.original_vcf}

            gatk IndexFeatureFile -I {output.original_split_vcf}
        '''
rule prep_original:
    input:
        original_vcf_split = 'imputation_only_results/split/chr{chrom}/chr{chrom}.vcf.gz', 
        original_vcf_split_tbi = 'imputation_only_results/split/chr{chrom}/chr{chrom}.vcf.gz.tbi',
    output:
        original_vcf_prepped = 'imputation_only_results/split/prepped/chr{chrom}/chr{chrom}.vcf.gz',
        original_vcf_prepped_tbi = 'imputation_only_results/split/prepped/chr{chrom}/chr{chrom}.vcf.gz.tbi',
    resources:
        time = 30,
        mem_mb = 24000,
        cpus = 4
    params:
        directory = 'imputation_only_results/split/prepped'
    shell:
        r'''
            bcftools norm -m -both {input.original_vcf_split} -Oz -o {params.directory}/chr{wildcards.chrom}.temp.vcf.gz

            bcftools view --min-ac 1 -Oz -o {params.directory}/chr{wildcards.chrom}.temp2.vcf.gz {params.directory}/chr{wildcards.chrom}.temp.vcf.gz

            bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' {params.directory}/chr{wildcards.chrom}.temp2.vcf.gz -o {output.original_vcf_prepped}

            gatk IndexFeatureFile -I {output.original_vcf_prepped}

            rm -f {params.directory}/chr{wildcards.chrom}.temp.vcf.gz
            rm -f {params.directory}/chr{wildcards.chrom}.temp2.vcf.gz
        '''


rule beagle40_impute:
    input:
        original_vcf_prepped = 'imputation_only_results/split/prepped/chr{chrom}/chr{chrom}.vcf.gz',
        original_vcf_prepped_tbi = 'imputation_only_results/split/prepped/chr{chrom}/chr{chrom}.vcf.gz.tbi',
        phased_vcf = '0.05_masking_results/split/ref/chr{chrom}/chr{chrom}.vcf.gz',
        phased_tbi = '0.05_masking_results/split/ref/chr{chrom}/chr{chrom}.vcf.gz.tbi',
    output:
        imputed_vcf = 'imputation_only_results/imputed/chr{chrom}/imputed.chr{chrom}.vcf.gz',
        imputed_tbi = 'imputation_only_results/imputed/chr{chrom}/imputed.chr{chrom}.vcf.gz.tbi'
    params:
        prefix = lambda wildcards, output: output.imputed_vcf.rsplit('.',2)[0],
        breed = config['breed'],
        breed_ab = config['breed_ab'],
    threads: 24
    resources:
        time    = 1440,
        mem_mb  = 60000,
        cpus    = 4
    shell:
        '''
            java -Xmx60g -jar beagle.27Jan18.7e1.jar \
                gtgl={input.original_vcf_prepped} \
                chrom=chr{wildcards.chrom} \
                ref={input.phased_vcf} \
                map=recombination_maps/{params.breed}/beagle_maps/BEAGLE_{params.breed_ab}_ECA{wildcards.chrom}_map.txt \
                impute=false \
                window=5000000 \
                overlap=20000 \
                ne=100 \
                err=0.05 \
                gprobs=true \
                out={params.prefix}

            gatk IndexFeatureFile -I {output.imputed_vcf}

        '''

rule imputed_list:
    input:
        imputed_vcfs = sorted(expand(
        'imputation_only_results/imputed/chr{chrom}/imputed.chr{chrom}.vcf.{ext}',
        chrom=[str(i) for i in range(1,31)], # NO CHROM M CORRECT?
        ext=['gz','gz.tbi']
        ))
    output:
        sorted_list = 'imputation_only_results/imputed/imputed.list',
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    run:
        # drop indices from input
        outfile = open('imputation_only_results/imputed/imputed.list', 'wt')
        for chrom in range(1,31):
            print('imputation_only_results/imputed/chr' + str(chrom) + '/imputed.chr' + str(chrom) + '.vcf.gz',
                    file = outfile)

rule combine_imputed:
    input:
        sorted_list = 'imputation_only_results/imputed/imputed.list', 
        imputed_vcfs = sorted(expand(
        'imputation_only_results/imputed/chr{chrom}/imputed.chr{chrom}.vcf.{ext}', 
        chrom=[str(i) for i in range(1,31)], # NO CHROM M CORRECT?
        ext=['gz','gz.tbi']
        ))
    output:
        sorted_vcf = 'imputation_only_results/FINAL/imputed.vcf.gz',
        sorted_tbi = 'imputation_only_results/FINAL/imputed.vcf.gz.tbi' 
    threads: 4
    resources:
        time   = 720,
        mem_mb = 24000
    shell:
        '''
            bcftools concat \
            -Oz -o {output.sorted_vcf} \
            -f {input.sorted_list}

            gatk IndexFeatureFile -I {output.sorted_vcf}
        '''

