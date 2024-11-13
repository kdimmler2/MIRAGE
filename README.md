# MIRAGE
## Mixed Imputation for Reference and Alternative Genotype Estimation

The MIRAGE pipeline is designed to seamlessly integrate two complementary imputation methods, maximizing genotype completion for datasets with missing data. Many genomic analyses require complete datasets and are sensitive to missing values, which are common in sequencing data. Furthermore, reference panels may lack some of the key variants needed for specific studies. MIRAGE addresses these challenges by combining Beagle 4.1 (reference-panel-based imputation) with STITCH (reference-panel-free imputation), leveraging the strengths of both approaches.

In the MIRAGE workflow, Beagle 4.1 first uses a reference panel to impute as many genotypes as possible. The partially imputed VCF is then merged with the original data, creating a comprehensive VCF file that integrates both imputed and remaining un-imputed data. This enriched dataset is subsequently processed by STITCH, which completes additional genotypes without requiring a reference panel. The result is a robust, high-coverage dataset suitable for rigorous downstream analysis, even in the presence of complex variant structures.

There are two usages. 1) To assess imputation accuracy in your specific dataset through a masking approach or 2) To fully impute missing data

In both options, imputation is done in parallel by chromosome.

### Option 1
#### Assessing Imputation Through Masking

#### Inputs:
- A VCF to be imputed
- A reference panel (must not share any samples with the input VCF)
- Recombination Maps (must be in the format Beagle wants)
- The proportion of missingness to be evaluated

These are specificed in config.yaml
Be aware:
bcftools_directory needs to be changed to the path where your bcftools is installed
breed and breed_ab are for use when breed-specific recombination maps are available. These can be commented out and you will need to change the beagle40_impute rule to remove the use of this parameter.

#### Workflow:
##### MIRAGE.smk:
The VCFs are split by chromosome. The input VCF is cleaned, and any variants with missing data are removed. Each VCF is then divided into bins by minor allele frequency (MAF). From each bin, a random 30% of all variants are selected to have genotypes masked (set to missing) at the specified missingness proportion. All intermediate VCFs are combined and then imputed using Beagle 4.1. MAF annotations are added to the final VCF.

##### combine.smk:
The imputed VCF and the original VCF are combined so that all imputed variants are present, along with the variants from the original dataset that were not found in the reference panel and therefore not available for imputation.

#### Output:
A VCF that mirrors the original dataset (after removing variants with missing data) but with a specified proportion of genotypes masked and subsequently imputed. This output allows for accuracy and concordance assessment of imputed genotypes within your dataset.

*Evaluation is under construction.*

### Option 2:
#### Full Imputation of Missing Data:

#### Inputs:
- A VCF to impute (should only include variants that pass filters)
- A reference panel

These are specificed in config_IO.yaml
Be aware:
breed and breed_ab are for use when breed-specific recombination maps are available. These can be commented out and you will need to change the beagle40_impute rule to remove the use of this parameter.

#### Workflow:
##### MIRAGE_IO.smk:
The VCF files are split by chromosome, and the input VCF undergoes cleaning, including the separation of multi-allelic sites. Each chromosome-specific VCF is then imputed using Beagle 4.1, and the chromosome VCFs are combined into a fully imputed dataset.

#### Output:
A VCF with imputed genotypes for missing data. Note that this VCF will only contain variants that are present in the reference panel.

*This section is currently under development. In future versions, the resulting VCF will retain all variants from the original dataset, with imputed genotypes added wherever possible.*
