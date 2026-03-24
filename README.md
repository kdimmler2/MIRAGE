# MIRAGE  
## Mixed Imputation for Reference and Alternative Genotype Estimation

The MIRAGE pipeline integrates two complementary imputation approaches to maximize genotype completion in datasets with missing data. Many genomic analyses require complete datasets and are sensitive to missing values, which are common in sequencing data. In addition, reference panels may lack variants relevant to specific studies.

MIRAGE addresses these challenges by combining:
- **Beagle 4.1** (reference-panel-based imputation)  
- **STITCH** (reference-panel-free imputation)  

to leverage the strengths of both approaches.

In this workflow, Beagle 4.1 first imputes genotypes using a reference panel. The partially imputed VCF is then merged with the original data, preserving variants not present in the reference panel. This combined dataset is subsequently processed by STITCH, which imputes remaining genotypes without requiring a reference panel. The result is a high-coverage dataset suitable for downstream analyses, even with complex variant structures.

---

## Usage Overview

MIRAGE supports two primary use cases:

1. **Assessing imputation accuracy via masking**  
2. **Full imputation of missing data**

In both cases, imputation is performed in parallel by chromosome.

---

## Option 1: Assessing Imputation Accuracy (Masking)

### Inputs

- VCF to be imputed  
- Reference panel (must not share samples with input VCF)  
- Recombination maps (Beagle-compatible format)  
- Desired proportion of missingness  

These parameters are specified in `config.yaml`.

**Notes:**
- `bcftools_directory` must be set to your local installation  
- `breed` and `breed_ab` are optional (used for breed-specific recombination maps)  

---

### Workflow

#### `MIRAGE.smk`

- Split VCF by chromosome  
- Remove variants with missing data  
- Bin variants by minor allele frequency (MAF)  
- Randomly select ~30% of variants per bin  
- Mask genotypes at specified missingness proportion  
- Recombine intermediate VCFs  
- Impute using Beagle 4.1  
- Annotate final VCF with MAF  

#### `combine.smk`

- Merge imputed VCF with original VCF  
- Retain:
  - imputed variants  
  - variants absent from the reference panel  

---

### Output

A VCF mirroring the original dataset (after filtering), with masked genotypes re-imputed.  
This enables evaluation of imputation accuracy and concordance.

*Evaluation is currently under development.*

---

## Option 2: Full Imputation of Missing Data

### Inputs

- Filtered VCF (variants passing QC)  
- Reference panel  

Specified in `config_IO.yaml`.

**Notes:**
- `breed` and `breed_ab` are optional and may require adjusting the `beagle40_impute` rule  

---

### Workflow

#### `MIRAGE_IO.smk`

- Split VCF by chromosome  
- Clean input VCF (including separation of multi-allelic sites)  
- Perform imputation with Beagle 4.1 per chromosome  
- Combine chromosome-level VCFs  

---

### Output

A VCF containing imputed genotypes for missing data.

**Note:**  
Currently, only variants present in the reference panel are retained.

*Future versions will retain all original variants with imputed genotypes added where possible.*

---

## Notes

This repository reflects a research workflow and may require adaptation depending on the dataset and analysis goals. Scripts are modular and were developed iteratively rather than as a single automated pipeline.
