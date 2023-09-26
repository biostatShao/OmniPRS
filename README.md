# OmniPRS
===========================================================================
## Background
Here we provide a new PRS model using conventional algorithms BLUP that estimate individuals' total genetic effects as random effects in mixed linear models.

**[OmniPRS](https://github.com/biostatShao/OmniPRS)** is implemented in R statistical environment.

## Example
```ruby
library(data.table)
library(magrittr)
source("~/szh_code.R")
source("~/GRS.pipeline.R")

traits = "Height"
bina = F
phe_trait = traits
sums_p = "/your/path/test/sums/"
out = "/your/path/test/out/"

base_p = "/data2/projects/bioinfo/zhshao/hm3eur/"
base_f = "eur_hm3_chr"
target_p = "/data2/projects/bioinfo/zhshao/UKB_hm3eur/"
target_f = "ukb22828_eur_hm3_chr"
pheno = "/data2/projects/bioinfo/zhshao/GRS/ukbb.phen"
temp = "/data2/projects/bioinfo/zhshao/GRS/tem/"

sum.pro(paste0("/data2/projects/bioinfo/zhshao/GWAS.summary/raw/GLGC/",traits,".raw"), 
        traits, sums_p)
gen.inp(traits, sums_p, pheno, phe_trait)
sums_p = paste0(sums_p,traits,"/")

for(chr in 1:22){
  GRS.OmniPRS(GRS.input, traits, chr, sums_p,
              base_p, base_f, target_p, target_f,
              pheno, phe_trait, out, temp,bina)
}
```

## Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Zhonghe Shao](https://github.com/biostatShao) via zhonghe@hust.edu.cn.

## Update
2023-01-06 omniPRS version 1.0.
