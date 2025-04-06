# OmniPRS
![](https://github.com/biostatShao/OmniPRS/blob/main/image.png)
===========================================================================
# 1 Overview

With the rapid development and maturation of high-throughput omics technologies, we face new opportunities and challenges in integrating multi-level and multi-omics data for computational modeling and disease risk prediction. Polygenic Risk Scores (PRS) play a crucial role in this systems epidemiology, effectively predicting disease risk at both individual and population levels, optimizing screening programs, and enabling precision prevention.

To enhance the accuracy of polygenic prediction, we propose a polygenic risk score R package based on a comprehensive genetic model, named **OmniPRS**. This package utilizes various functional annotation categories, showcasing high flexibility in accommodating different types of annotations. Notably, the OmniPRS-R framework integrates 11 tissue-specific annotation models to estimate the enrichment of GWAS signals by combining these models with summary statistics from genome-wide association studies (GWAS). 

We use stratified LD score regression (S-LDSC) to assess the contribution of each annotation model to polygenic heritability. OmniPRS models polygenic risk by using functional annotation-based heritability and covariance matrices, along with marginal summary statistics and LD matrices from reference panels, within a mixed model framework to infer the best linear unbiased predictors (BLUPs) for each SNP, i.e., the joint effects.

Additionally, we explore three alternative strategies for integrating tissue-specific annotation PRS models: Equal Weight (EW) models, Bayesian Model Averaging (BMA), and Least Absolute Shrinkage and Selection Operator (LASSO) models. These strategies assist in evaluating the relative contribution of each functional category to the overall genetic risk.

Ultimately, the OmniPRS package provides a comprehensive tool and metrics for assessing the overall genetic risk of diseases or traits by aggregating individual polygenic scores from multiple functional categories.

# 2 Installation Instructions

The OmniPRS package does not require installation and can be directly invoked. Before use, set the working directory to the location of the software:
```r
setwd("/your/path/")
require(data.table)
require(magrittr)
require(rhdf5)
require(BAS)
require(glmnet)
```
And the following Python packages:
```python
import getopt
import sys
import os
import math
import timeit
import time
import h5py
import scipy as sp
from scipy import stats
from scipy import linalg
from plinkio import plinkfile
import itertools as it
import gzip
import glob
import re
import traceback
```
Then, in the R environment, load the OmniPRS package by entering:
```r
source("/your/path/OmniPRS.R")
```
## 3 Usage Instructions

### 3.1 Input Files

#### 3.1.1 Plink Files

Each chromosome corresponds to a set of Plink files (binary PED format). These files are used as the LD reference panel and are later used to compute PRS.

#### 3.1.2 Annotation File

##### 3.1.2.1

First, you need to infer the heritability of each SNP using the S-LDSC model (Baseline-LD annotation can be downloaded [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/)) and use summary statistics that will later be used for training. When running S-LDSC, make sure to use the `--print-coefficients` flag to obtain the regression coefficients.

##### 3.1.2.2

After running S-LDSC:

a. Obtain the `h2g` estimate from the `*.log` file.

b. Retrieve the regression coefficients from the `*.results` file (8th column). Divide the regression coefficients by `h2g` and define it as `T=tau/h2g`. This results in a vector of dimension `C*1`, where `C` is the total number of annotations.

c. From the Baseline-LD annotation downloaded [here](https://data.broadinstitute.org/alkesgroup/LDSCORE/), read the annotation file `baselineLD.*.annot.gz` and retain only the annotation columns (i.e., remove the first 4 columns). This matrix is denoted as `X` with dimensions `M*C`, where `M` is the number of SNPs and `C` is the total number of annotations.

d. Define the expected heritability of each SNP as an `M*1` vector, which is the result of matrix `X` multiplied by `T`.

The final file format is:
- Column 1: SNP ID
- Column 2: Heritability of each SNP

#### 3.1.3 GWAS Summary Statistics File

The file must adhere strictly to the following format:
- `CHR`: Chromosome
- `SNP`: SNP ID
- `BP`: Physical position (base pairs)
- `A1`: Minor allele name
- `A2`: Major allele name
- `P`: GWAS p-value
- `BETA`: Effect size
- `Z`: Z-statistic

#### 3.1.4 Phenotype File

Format:
- Column 1: FID
- Column 2: Phenotype
- Column 3 (optional): Covariates

### 3.2 Input Parameters

- `traits`: Trait name (e.g., Height). The prefix of the summary data file should be consistent (e.g., Height_sums.txt). Directory name for output files.
- `chr`: Chromosome number (e.g., 1)
- `N`: The GWAS sample size
- `h2`: Estimated SNP Heritability (pre-compute this using your favorite method).
- `sums_p`: Absolute path to the summary data (e.g., `/your/path/sum_file/`)
- `base_p`: Prefix of Plink format LD data files, without the chromosome number (e.g., for chromosome 1, the file is `/your/path/plink/eur_hm3_chr1`, so enter `/your/path/plink/eur_hm3_chr`)
- `target_p`: Prefix of Plink format test set data files, without the chromosome number (e.g., for chromosome 1, the file is `/your/path/plink/ukb22828_eur_hm3_chr1`, so enter `/your/path/plink/ukb22828_eur_hm3_chr`)
- `pheno`: Phenotype file and its absolute path. If covariates are included, they should be in this file (e.g., `/your/path/ukbb.phen`)
- `phe_trait`: Column name of the outcome in the phenotype file (e.g., Height)
- `out`: Output path for result files (e.g., `/your/path/out/`)
- `temp`: Output path for temporary files (e.g., `/your/path/temp/`)
- `cova`: Covariates to consider; if none, enter `NULL` (e.g., `c("BaseAge","Sexgenetic")`)
- `bina`: Whether the outcome is binary data (e.g., `T`)
- `ct_result`: Whether to output 11 tissue-specific PRS prediction levels (e.g., `T`)
- `software_path`: The directory to OmniPRS (e.g., `/your/path/`)
- `plink_path`: The directory to PLINK1.9 and PLINK2.0

## 3.3 Running

Below is an example of running the OmniPRS package:

```r
# R version
R 4.1.2

# Set parameters
base_p <- "/data2/projects/bioinfo/1KG/"
base_f <- "eur_hm3_chr"
target_p <- "/data2/projects/bioinfo/UKB/"
target_f <- "ukb22828_eur_hm3_chr"
pheno <- "/data2/projects/bioinfo/ukbb.phen"
temp <- "/data2/projects/bioinfo/temp_ukb/"
cova <- c("BaseAge", "Sexgenetic", paste0("PC", 1:10))
traits <- "Height"
bina <- FALSE
sums_p <- "/data2/projects/bioinfo/sums/"
out <- "/data2/projects/bioinfo/out/"
phe_trait <- traits
print(traits)

# Set GWAS summary file path
sums_p <- paste0(sums_p, traits, "/")

# Run OmniPRS for each chromosome
for (chr in 1:22) {
  OmniPRS(traits, chr, sums_p,
    base_p, base_f, cova, target_p, target_f,
    pheno, phe_trait, out, temp, bina)
}
```

### 3.4 Output Files

#### 3.4.1 Posterior Estimates of Effect Sizes

- **Column 1**: SNP ID
- **Column 2**: Minor Allele
- **Columns 3-13**: Posterior estimates of tissue-specific effect sizes for each tissue

#### 3.4.2 Polygenic Risk Scores for Validation Set Individuals

- **Column 1**: SNP ID
- **Column 2**: Minor Allele
- **Columns 3-13**: Tissue-specific polygenic risk scores for each tissue
- **Column 14**: Polygenic risk score using the EW algorithm
- **Column 15**: Polygenic risk score using the LASSO algorithm
- **Column 16**: Polygenic risk score using the BMA algorithm

#### 3.4.3 Polygenic Risk Score Prediction Levels

Prediction levels are output as filenames (e.g., if the prediction level R² is 0.45, the output will be a file named `OmniPRS-0.45`).

## Acknowledgments
We refer to the preprocessor of GWAS data in LDpred.funct software developed by Carla Márquez-Luna for this software.

## Citations
Zhonghe Shao, Wangxia Tang, Hongji Wu, Yifan Kong, Xingjie Hao* (2024). lncorporating multiple functional annotations to improve polygenic risk prediction accuracy, Cell Genomics (2025), https://doi.org/10.1016/i.xgen.2025.100850

## Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Zhonghe Shao](https://github.com/biostatShao) via zhonghe@hust.edu.cn.

## Update
2024-06-17 OmniPRS version 1.0.0
