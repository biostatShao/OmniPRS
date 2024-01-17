# OmniPRS
===========================================================================
## Background
We introduced OmniPRS, a framework that enhances the genetic risk prediction for complex diseases and traits using the GWAS summary statistics. Within the framework, we re-estimated the effect sizes of SNPs with the mixed model, where the genetic variance components are based on their functional annotations and can be tissue-specific. We constructed multiple tissue-specific PRSs and integrated these scores as the final OmniPRS. 

## Example
```ruby
library(data.table)
library(magrittr)
source("~/function.R")

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
## Citations
Zhonghe Shao, Wangxia Tang, Yifan Kong, Si Li, Yunlong Guan, Hongji Wu, Minghui Jiang, Xi Cao, Xingjie Hao* (2024). OmniPRS: Incorporating multiple functional annotations to improve polygenic risk prediction accuracy. Am. J. Hum. Genet. submitted.


## Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Zhonghe Shao](https://github.com/biostatShao) via zhonghe@hust.edu.cn.

## Update
2023-01-06 omniPRS version 1.0.
