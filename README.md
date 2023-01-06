# omniPRS
===========================================================================
## Background
Here we provide a new PRS model using conventional algorithms BLUP that estimate individuals' total genetic effects as random effects in mixed linear models.

In addition, we provided a set of common PRS models to use to generate the basic results. They will be combined later.

**[omniPRS](https://github.com/biostatShao/omniPRS)** is implemented in R statistical environment.

## Example
```ruby
args=as.numeric(commandArgs(TRUE))
taa = args[1]
tbb = args[2]

library(data.table)
library(magrittr)

if(T){
  traits = "height"[taa]
  chr = tbb
  raw_s <- "/data2/projects/bioinfo/zhshao/GWAS.summary/raw/height.txt"
  sums_p = "/data2/projects/bioinfo/zhshao/GWAS.summary/sums/"
  base_p = "/data2/projects/bioinfo/zhshao/hm3eur/"
  base_f = "eur_hm3_chr"
  target_p = "/data2/projects/bioinfo/zhshao/UKB_hm3eur/"
  target_f = "ukb22828_eur_hm3_chr"
  pheno = "/data2/projects/bioinfo/zhshao/GRS/ukbb.phen"
  phe_trait = "Height"
  out = "/home/zhshao/project/GRS/111_result/"
  temp = "/data2/projects/bioinfo/zhshao/GRS/temp/"
  
  source("/home/zhshao/project/szh_code.R")
  source("/home/zhshao/project/GRS/GRS.pipeline.R")
}

if(F){
  sum.pro(raw_s, traits, sums_p)
  gen.inp(raw_s, traits, sums_p, pheno, phe_trait)
  
  system(paste0("rm /home/zhshao/out/temp/*"))
  system(paste0("rm /home/zhshao/out/PRS/*"))
  
  GRS.LDpred2(raw_s, traits,
              sums_p, base_p, base_f,
              target_p, target_f,
              pheno, phe_trait,
              out, temp)
  
  res.pro(out)
}

for(GRS.model in c("CT","DBSLMM", "PRS.CS","lassosum",
                   "LDpred.funct","AnnoPred","SDPR","omniPRS")){
  GRS.input <- NA; attr(GRS.input,"class") <- GRS.model
  GRS(GRS.input, traits, chr, paste0(sums_p,traits,"/"), 
      base_p, base_f, target_p, target_f,
      pheno, phe_trait, out, temp)
}

```

## Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Zhonghe Shao](https://github.com/biostatShao) via zhonghe@hust.edu.cn.

## Update
2023-01-06 omniPRS version 1.0.
