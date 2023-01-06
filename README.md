# omniPRS
===============================================================================
## Background
As an effective supplementary strategy of single-marker analysis, multilocus methods have been increasingly applied. Multilocus analysis often jointly examines a set of SNPs that are pre-defined within a functional unit such as gene to evaluate the overall association evidence at the gene level; it is thus also referred to as SNP-set or gene-based approach. Due to the usefulness, distinct SNP-set methods have been recently developed, many of which can be implemented with only GWAS summary association statistics, greatly generalizing their applicability due to the widespread availability of summary-level data. With distinct SNP-set approaches for multilocus studies, one naturally wonders which one should be chosen in practice. Moreover, existing SNP-set methods are not used without deficiencies, potential limitations include insufficient power, inability to provide statistically valid tests under certain parameter settings, and reliance on permutation sampling. Unfortunately, despite the importance of multilocus analysis in GWAS and the vast number of SNP-set methods, few comprehensive comparison studies have been performed to evaluate their effectiveness. Subsequently, due to the lack of consensus on the most suitable SNP-set method, the realization of the above advantages and benefits is to some extent currently hindered.

In the present work we sought to fill this knowledge gap by conducting a comprehensive comparison for 22 commonly-used summary-statistics based SNP-set methods in the hope that our results could serve as an important guidance for practitioners on how to choose appropriate methods for SNP-set analysis in post-GWAS era. Including: MLR: multiple linear regression; FLM: functional multiple linear regression model; HC: higher criticism test; GHC: generalized higher criticism test; BJ: Berk-Jones test; GBJ: generalized Berk-Jones test; DOT: decorrelation by orthogonal transformation; BT: burden test; SKATO: optimal sequence kernel association test; SKAT: sequence kernel association test; Simes: Simes’s test; FCP: fisher combined probability; TPM: truncated product method; RTP: rank truncated product; ART: augmented rank truncation; ART-A: adaptive augmented rank truncation; GM: gamma method; GATES: gene-based association test that uses extended Simes procedure; HMP: The harmonic mean P value test; ACAT: aggregated Cauchy association test.

**[MCA](https://github.com/biostatpzeng/MCA)** is implemented in R statistical environment.

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
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatShao) via zhonghe@hust.edu.cn.

## Update
2023-01-06 omniPRS version 1.0.
