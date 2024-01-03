import gzip
import random

rule all:
    input:
        expand('results/ref/split/{chrom}/{chrom}.vcf.gz', chrom=[f'chr{i}' for i in range(1,32)] + ['chrX']),
        expand('results/ref/split/{chrom}/{chrom}.vcf.gz.tbi', chrom=[f'chr{i}' for i in range(1,32)] + ['chrX']),
        expand('results/imputation_only/split/{chrom}/{chrom}.vcf.gz', chrom=[f'chr{i}' for i in range(1,32)] + ['chrX']),
        expand('results/imputation_only/split/{chrom}/{chrom}.vcf.gz.tbi', chrom=[f'chr{i}' for i in range(1,32)] + ['chrX']),
        expand('results/imputation_only/imputed/{chrom}/imputed.{chrom}.vcf.gz', chrom=[f'chr{i}' for i in range(1,32)] + ['chrX']),
        expand('results/imputation_only/imputed/{chrom}/imputed.{chrom}.vcf.gz.tbi', chrom=[f'chr{i}' for i in range(1,32)] + ['chrX']),
        'results/imputation_only/imputed/imputed.list',
        'results/imputation_only/FINAL/imputed.vcf.gz',
        'results/imputation_only/FINAL/imputed.vcf.gz.tbi',

#rule intersect:
#    input:
#        reference_panel = config['reference_panel'],
#        reference_panel_tbi = config['reference_panel_tbi'],
#        sample_vcf = config['sample_vcf'],
#        sample_vcf_tbi = config['sample_vcf_tbi'],
#    output:
#        split_vcf = 'results/ref/split/{chrom}/{chrom}.vcf.gz',
#        split_tbi = 'results/ref/split/{chrom}/{chrom}.vcf.gz.tbi'
#    threads: 4
#    resources:
#        time    = 1440,
#        mem_mb  = 60000,
#        cpus    = 4
#    shell:
#            '''
#                bcftools view \
#                -r {wildcards.chrom} \
#                -Oz -o {output.split_vcf} \
#                {input.reference_panel}
#
#                gatk IndexFeatureFile -I {output.split_vcf}
#            '''

rule split_phased:
    input:
        reference_panel = config['reference_panel'],
        reference_panel_tbi = config['reference_panel_tbi'] 
    output:
        split_vcf = 'results/ref/split/{chrom}/{chrom}.vcf.gz',
        split_tbi = 'results/ref/split/{chrom}/{chrom}.vcf.gz.tbi'
    threads: 4
    resources:
        time    = 1440,
        mem_mb  = 60000,
        cpus    = 4
    shell:
            '''
                bcftools view \
                -r {wildcards.chrom} \
                -Oz -o {output.split_vcf} \
                {input.reference_panel}

                gatk IndexFeatureFile -I {output.split_vcf}
            '''

rule split_sample:
   input:
       sample_vcf = config['sample_vcf'],
       sample_vcf_tbi = config['sample_vcf_tbi'],
   output:
       split_vcf = 'results/imputation_only/split/{chrom}/{chrom}.vcf.gz',
       split_tbi = 'results/imputation_only/split/{chrom}/{chrom}.vcf.gz.tbi',
   threads: 4
   resources:
       time    = 1440,
       mem_mb  = 60000,
       cpus    = 4
   shell:
       '''
           bcftools view \
           -r {wildcards.chrom} \
           -Oz -o {output.split_vcf} \
           {input.sample_vcf}

           tabix -p vcf {output.split_vcf}
       '''

rule beagle40_impute:
    input:
        sample_vcf = 'results/imputation_only/split/{chrom}/{chrom}.vcf.gz', 
        sample_tbi = 'results/imputation_only/split/{chrom}/{chrom}.vcf.gz.tbi', 
        phased_vcf = 'results/ref/split/{chrom}/{chrom}.vcf.gz', 
        phased_tbi = 'results/ref/split/{chrom}/{chrom}.vcf.gz.tbi' ,
    output:
        imputed_vcf = 'results/imputation_only/imputed/{chrom}/imputed.{chrom}.vcf.gz',
        imputed_tbi = 'results/imputation_only/imputed/{chrom}/imputed.{chrom}.vcf.gz.tbi'
    params:
        prefix = lambda wildcards, output: output.imputed_vcf.rsplit('.',2)[0]
    threads: 24
    resources:
        time    = 360,
        mem_mb  = 60000,
        cpus    = 4
    shell:
        '''
            java -jar beagle.27Jan18.7e1.jar \
            gl={input.sample_vcf} \
            ref={input.phased_vcf} \
            nthreads={threads} \
            impute=false \
            gprobs=true \
            window=50000 \
            overlap=20000 \
            out={params.prefix}

            gatk IndexFeatureFile -I {output.imputed_vcf}

        '''

rule imputed_list:
    input:
        imputed_vcfs = sorted(expand(
        'results/imputation_only/imputed/{chrom}/imputed.{chrom}.vcf.{ext}',
        chrom=[f'chr{i}' for i in range(1,31)] + ['chrX'], # NO CHROM M CORRECT?
        ext=['gz','gz.tbi'],
        itr=[f'set{i}' for i in range(1,11)]
        ))
    output:
        sorted_list = 'results/imputation_only/imputed/imputed.list',
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    run:
        # drop indices from input
        for itr in range(1,11):
            outfile = open('results/imputation_only/imputed/imputed.list', 'wt')
            for chrom in range(1,31):
                print('results/imputation_only/imputed/chr' + str(chrom) + '/imputed.chr' + str(chrom) + '.vcf.gz',
                    file = outfile)

rule combine_imputed:
    input:
        sorted_list = 'results/imputation_only/imputed/imputed.list', 
        imputed_vcfs = sorted(expand(
        'results/imputation_only/imputed/{chrom}/imputed.{chrom}.vcf.{ext}', 
        chrom=[f'chr{i}' for i in range(1,31)] + ['chrX'], # NO CHROM M CORRECT?
        ext=['gz','gz.tbi'],
        itr=[f'set{i}' for i in range(1,11)]
        ))
    output:
        sorted_vcf = 'results/imputation_only/FINAL/imputed.vcf.gz',
        sorted_tbi = 'results/imputation_only/FINAL/imputed.vcf.gz.tbi' 
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

