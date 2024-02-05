import gzip
import random

rule all:
    input:
        expand('PaintedPatches/results/masked_vcfs/{itr}/masked_{itr}.vcf.gz', itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
        expand('PaintedPatches/results/masked_vcfs/{itr}/masked_{itr}.vcf.gz.tbi', itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
        expand('PaintedPatches/results/ref/split/chr{chrom}/chr{chrom}.vcf.gz', chrom=[str(i) for i in range(1,31)] + ['X']),
        expand('PaintedPatches/results/ref/split/chr{chrom}/chr{chrom}.vcf.gz.tbi', chrom=[str(i) for i in range(1,31)] + ['X']),
        expand('PaintedPatches/results/split/masked_vcfs/{itr}/chr{chrom}/masked.chr{chrom}.vcf.gz', chrom=[str(i) for i in range(1,31)] + ['X'], itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
        expand('PaintedPatches/results/split/masked_vcfs/{itr}/chr{chrom}/masked.chr{chrom}.vcf.gz.tbi', chrom=[str(i) for i in range(1,31)] + ['X'], itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
        expand('PaintedPatches/results/imputed/masked_{itr}/chr{chrom}/imputed.chr{chrom}.vcf.gz', chrom=[str(i) for i in range(1,31)] + ['X'], itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
        expand('PaintedPatches/results/imputed/masked_{itr}/chr{chrom}/imputed.chr{chrom}.vcf.gz.tbi', chrom=[str(i) for i in range(1,31)] + ['X'], itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
        expand('PaintedPatches/results/imputed/masked_{itr}/imputed.list', itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
        expand('PaintedPatches/results/FINAL/masked_{itr}/imputed.vcf.gz', itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
        expand('PaintedPatches/results/FINAL/masked_{itr}/imputed.vcf.gz.tbi', itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
        expand('PaintedPatches/results/mask_assess/{itr}/imputed.table', itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
        expand('PaintedPatches/results/mask_assess/{itr}/original.table', itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
        expand('PaintedPatches/results/mask_assess/{itr}/masked.table', itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]),
#        expand('PaintedPatches/results/mask_assess/{itr}/mismatched_genos.txt', itr=[f'set{i}' for i in range(1, config['iterations'] + 1)]),
#        expand('PaintedPatches/results/mask_assess/{itr}/remove_imputed.list', itr=[f'set{i}' for i in range(1, config['iterations'] + 1)]),
#        expand('PaintedPatches/results/mask_assess/{itr}/README.txt', itr=[f'set{i}' for i in range(1, config['iterations'] + 1)]),
#        expand('PaintedPatches/results/mask_assess/{itr}/no_imputed.vcf.gz', itr=[f'set{i}' for i in range(1, config['iterations'] + 1)]),
#        expand('PaintedPatches/results/mask_assess/{itr}/no_imputed.vcf.gz.tbi', itr=[f'set{i}' for i in range(1, config['iterations'] + 1)]),
#        expand('PaintedPatches/results/mask_assess/{itr}/full.vcf.gz', itr=[f'set{i}' for i in range(1, config['iterations'] + 1)]),
#        expand('PaintedPatches/results/mask_assess/{itr}/full.vcf.gz.tbi', itr=[f'set{i}' for i in range(1, config['iterations'] + 1)]),
#        expand('PaintedPatches/results/mask_assess/{itr}/plink.imiss', itr=[f'set{i}' for i in range(1, config['iterations'] + 1)]),
#        expand('PaintedPatches/results/mask_assess/{itr}/plink.lmiss', itr=[f'set{i}' for i in range(1, config['iterations'] + 1)]),
#        expand('PaintedPatches/results/mask_assess/{itr}/plink.nosex', itr=[f'set{i}' for i in range(1, config['iterations'] + 1)]),
#        expand('PaintedPatches/results/mask_assess/{itr}/plink.log', itr=[f'set{i}' for i in range(1, config['iterations'] + 1)]),

rule mask:
    input:
        vcf_to_mask = config['vcf_to_mask'],
        vcf_to_mask_tbi = config['vcf_to_mask_tbi'],
    output:
        masked_vcf = 'PaintedPatches/results/masked_vcfs/{itr}/masked_{itr}.vcf.gz',
        masked_tbi = 'PaintedPatches/results/masked_vcfs/{itr}/masked_{itr}.vcf.gz.tbi',
    threads: 4
    resources:
        time    = 30,
        mem_mb  = 24000,
        cpus    = 4,
#    conda:
#        'mask.yaml' #this command for setGT only works with bcftools 1.17
    params:
        seed = lambda wildcards, input: random.randint(1, 100),
        mask_proportion = config['mask_proportion'],
        bcftools_directory = config['bcftools_directory']
    shell:
            '''
                {params.bcftools_directory} +setGT -o {output.masked_vcf} {input.vcf_to_mask} -- -t r:{params.mask_proportion} -s {params.seed} -n .
                
                gatk IndexFeatureFile -I {output.masked_vcf}

            '''

rule split_phased:
    input:
        reference_panel = config['reference_panel'],
        reference_panel_tbi = config['reference_panel_tbi'] 
    output:
        split_vcf = 'PaintedPatches/results/ref/split/chr{chrom}/chr{chrom}.vcf.gz',
        split_tbi = 'PaintedPatches/results/ref/split/chr{chrom}/chr{chrom}.vcf.gz.tbi'
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

rule split_masked:
   input:
       masked_vcf = 'PaintedPatches/results/masked_vcfs/{itr}/masked_{itr}.vcf.gz',
       masked_tbi = 'PaintedPatches/results/masked_vcfs/{itr}/masked_{itr}.vcf.gz.tbi',
   output:
       split_vcf = 'PaintedPatches/results/split/masked_vcfs/{itr}/chr{chrom}/masked.chr{chrom}.vcf.gz',
       split_tbi = 'PaintedPatches/results/split/masked_vcfs/{itr}/chr{chrom}/masked.chr{chrom}.vcf.gz.tbi',
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
           {input.masked_vcf}

           tabix -p vcf {output.split_vcf}
       '''

rule beagle40_impute:
    input:
        sample_vcf = 'PaintedPatches/results/split/masked_vcfs/{itr}/chr{chrom}/masked.chr{chrom}.vcf.gz', 
        sample_tbi = 'PaintedPatches/results/split/masked_vcfs/{itr}/chr{chrom}/masked.chr{chrom}.vcf.gz.tbi', 
        phased_vcf = 'PaintedPatches/results/ref/split/chr{chrom}/chr{chrom}.vcf.gz', 
        phased_tbi = 'PaintedPatches/results/ref/split/chr{chrom}/chr{chrom}.vcf.gz.tbi',
    output:
        imputed_vcf = 'PaintedPatches/results/imputed/masked_{itr}/chr{chrom}/imputed.chr{chrom}.vcf.gz',
        imputed_tbi = 'PaintedPatches/results/imputed/masked_{itr}/chr{chrom}/imputed.chr{chrom}.vcf.gz.tbi'
    params:
        prefix = lambda wildcards, output: output.imputed_vcf.rsplit('.',2)[0],
        breed = config['breed'],
        breed_ab = config['breed_ab'],
    threads: 24
    resources:
        time    = 360,
        mem_mb  = 60000,
        cpus    = 4
    shell:
        '''
            java -jar beagle.27Jan18.7e1.jar \
            gtgl={input.sample_vcf} \
            chrom=chr{wildcards.chrom} \
            map=~/NuGEN/Imputation/reference_panel/beagle_54_all_breeds/recombination_maps/beagle_maps/BEAGLE_{params.breed}_ECA{wildcards.chrom}_map.txt \
            ref={input.phased_vcf} \
			window=5000000 \
			overlap=20000 \
			ne=100 \
			err=0.05 \
			gprobs=true \
            nthreads={threads} \
            impute=true \
            out={params.prefix}

            gatk IndexFeatureFile -I {output.imputed_vcf} 
        '''

rule imputed_list:
    input:
        imputed_vcfs = sorted(expand(
        'PaintedPatches/results/imputed/masked_{itr}/{chrom}/imputed.{chrom}.vcf.{ext}',
        chrom=['chr' + str(i) for i in range(1,31)] + ['chrX'], # NO CHROM M CORRECT?
        ext=['gz','gz.tbi'],
        itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]
        ))
    output:
        sorted_list = 'PaintedPatches/results/imputed/masked_{itr}/imputed.list',
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    run:
        # drop indices from input
        for itr in range(1, config['iterations'] + 1):
            outfile = open('PaintedPatches/results/imputed/masked_set' + str(itr) + '/imputed.list', 'wt')
            for chrom in range(1,32):
                print('PaintedPatches/results/imputed/masked_set' + str(itr) + '/chr' + str(chrom) + '/imputed.chr' + str(chrom) + '.vcf.gz',
                    file = outfile)

rule combine_imputed:
    input:
        sorted_list = 'PaintedPatches/results/imputed/masked_{itr}/imputed.list', 
        imputed_vcfs = sorted(expand(
        'PaintedPatches/results/imputed/masked_{itr}/{chrom}/imputed.{chrom}.vcf.{ext}', 
        chrom=['chr' + str(i) for i in range(1,31)] + ['chrX'], # NO CHROM M CORRECT?
        ext=['gz','gz.tbi'],
        itr=['set' + str(i) for i in range(1, config['iterations'] + 1)]
        ))
    output:
        sorted_vcf = 'PaintedPatches/results/FINAL/masked_{itr}/imputed.vcf.gz',
        sorted_tbi = 'PaintedPatches/results/FINAL/masked_{itr}/imputed.vcf.gz.tbi' 
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

rule to_table:
    input:
        masked_vcf = 'PaintedPatches/results/masked_vcfs/{itr}/masked_{itr}.vcf.gz',
        original_vcf = config['vcf_to_mask'],
        imputed_vcf = 'PaintedPatches/results/FINAL/masked_{itr}/imputed.vcf.gz',
    output:
        imputed_table = 'PaintedPatches/results/mask_assess/{itr}/imputed.table',
        original_table = 'PaintedPatches/results/mask_assess/{itr}/original.table',
        masked_table = 'PaintedPatches/results/mask_assess/{itr}/masked.table',
    threads: 4
    resources:
        time   = 30,
        mem_mb = 24000
    shell:
        '''
             bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%DR2\t%QUAL\t%FILTER\t[%DS\t]\n' -H {input.imputed_vcf}  > {output.imputed_table}

             bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%AC\t%QUAL\t%FILTER\t[%GT\t]\n' -H {input.original_vcf}  > {output.original_table}

             bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\t%AC\t%QUAL\t%FILTER\t[%GT\t]\n' -H {input.masked_vcf}  > {output.masked_table}
    
        '''

rule compare_genos:
    input:
        imputed_table = 'PaintedPatches/results/mask_assess/{itr}/imputed.table',
        original_table = 'PaintedPatches/results/mask_assess/{itr}/original.table',
        masked_table = 'PaintedPatches/results/mask_assess/{itr}/masked.table',
    output:
        mismatched_genos = 'PaintedPatches/results/mask_assess/{itr}/mismatched_genos.txt',
        remove_imputed = 'PaintedPatches/results/mask_assess/{itr}/remove_imputed.list',
        readme = 'PaintedPatches/results/mask_assess/{itr}/README.txt',
    threads: 4
    resources:
        time   = 30,
        mem_mb = 24000
    run:
        infile1 = open(input.original_table, 'rt')
        infile2 = open(input.masked_table, 'rt')
        infile3 = open(input.imputed_table, 'rt')

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
        remove_imputed = 'PaintedPatches/results/mask_assess/{itr}/remove_imputed.list',
    output:
        no_imputed = 'PaintedPatches/results/mask_assess/{itr}/no_imputed.vcf.gz',
        no_imputed_tbi = 'PaintedPatches/results/mask_assess/{itr}/no_imputed.vcf.gz.tbi',
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
        no_imputed = 'PaintedPatches/results/mask_assess/{itr}/no_imputed.vcf.gz',
        imputed_vcf = 'PaintedPatches/results/FINAL/masked_{itr}/imputed.vcf.gz',
    output:
        full_vcf = 'PaintedPatches/results/mask_assess/{itr}/full.vcf.gz',
        full_tbi = 'PaintedPatches/results/mask_assess/{itr}/full.vcf.gz.tbi',
    threads: 4
    resources:
        time   = 30,
        mem_mb = 24000
    shell:
        '''
            bcftools concat \
                -Oz \
                -o {output.full_vcf} \
                --allow-overlaps \
                {input.imputed_vcf} {input.no_imputed}

            gatk IndexFeatureFile -I {output.full_vcf}
 
        '''

rule plink_missing:
    input:
        full_vcf = 'PaintedPatches/results/mask_assess/{itr}/full.vcf.gz', 
    output:
        imiss = 'PaintedPatches/results/mask_assess/{itr}/plink.imiss',
        lmiss = 'PaintedPatches/results/mask_assess/{itr}/plink.lmiss',
        nosex = 'PaintedPatches/results/mask_assess/{itr}/plink.nosex',
        logfile = 'PaintedPatches/results/mask_assess/{itr}/plink.log',
    threads: 4
    resources:
        time   = 30,
        mem_mb = 24000
    conda:
        'plink.yaml'
    params:
        directory = 'PaintedPatches/results/mask_assess/{itr}/plink'
    shell:
        '''
            plink \
                --missing \
                --horse \
                --vcf {input.full_vcf} \
                --out {params.directory}
        '''

