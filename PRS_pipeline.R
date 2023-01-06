# summary data preparation ------------------------------------------------

sum.pro <- function(summs, trait.name = "Height", 
                    out = "/data2/projects/bioinfo/zhshao/GWAS.summary/sums/"){
  require(data.table)
  require(magrittr)
  
  out_p <- paste0(out, trait.name,"/")
  dir.check(out_p); setwd(out_p)
  
  cat("Start processing trait: ",trait.name,"\n")
  data <- fread(summs) #[1] 1373020 
  hm3 <- fread("/data2/projects/bioinfo/zhshao/GWAS.summary/w_hm3.snplist") #[1] 1217311 
  cat("A total of ",nrow(data),"SNPs are loaded \n")
  
  data <- merge(hm3, data, by.x = "SNP",by.y = "RSID") %>% na.omit() #[1] 1101852 
  cat(nrow(data),"belong to HapMap3 SNPs \n")
  data$P <- as.numeric(data$P)
  attach(data)
  index1 <- A1 == EFFECT_ALLELE & A2 == OTHER_ALLELE
  index2 <- A2 == EFFECT_ALLELE & A1 == OTHER_ALLELE
  index3 <- N > 0.67 * quantile(N,0.90)
  index4 <- P > 0 & P <= 1
  BETA[index2] <- -BETA[index2]
  EFFECT_ALLELE_FREQ[index2] <- 1 - EFFECT_ALLELE_FREQ[index2]
  detach(data)
  data <- data[(index1 | index2) & index3 & index4, 
               c("SNP","CHR","POS","A1","A2","BETA","SE","P","N","EFFECT_ALLELE_FREQ")] #[1] 1072395 
  cat("A total of ",nrow(data),"SNPs passed the QC \n")
  # fwrite(data,paste0(out_p,trait.name,".txt"), sep="\t",quote=F, row.names=F)
  fwrite(as.data.frame(max(data$N)),paste0(out_p, trait.name,"-N.txt"), sep="\t",quote=F, row.names=F)
  
  # data <- fread("/data2/projects/bioinfo/zhshao/GWAS.summary/raw/GIANT_height/Height.txt")
  data$Z <- data$BETA / data$SE
  data <- data[,c(1,4,5,2,3,9,6,7,11,8,10)]
  colnames(data) <- c("SNP","A1","A2","CHR","BP","N","Beta","SE","Z","P","MAF")
  fwrite(data,paste0(out_p, trait.name, ".sumstats"), sep="\t",quote=F, row.names=F)
  
  # munge_sumstats <- paste0("python /home/opt/software/anaconda3/envs/ldsc/bin/munge_sumstats.py",
  #                          " --merge-alleles /data2/projects/bioinfo/zhshao/GWAS.summary/w_hm3.snplist",
  #                          " --sumstats ", out_p, trait.name, ".txt",
  #                          " --N ", max(data$N),
  #                          " --out ", out_p, trait.name)
  # system(munge_sumstats)
  # 
  # system(paste0("gzip -d ", out_p, trait.name,".sumstats.gz"))
  
  ldsc_p <- paste0(out_p, "ldsc/")
  dir.check(ldsc_p)
  cat("Start ldsc \n")
  ldsc <- paste0("python /home/opt/software/anaconda3/envs/ldsc/bin/ldsc.py",
                 " --ref-ld-chr /data2/projects/bioinfo/xingjie/src/ldsc/eur_w_ld_chr/",
                 " --w-ld-chr /data2/projects/bioinfo/zhshao/LDSC/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.",
                 " --h2 ", out_p, trait.name, ".sumstats",
                 " --out ", ldsc_p, trait.name,
                 " > ", ldsc_p, trait.name, ".ldsc.temp")
  system(ldsc)
  cat("Start baseline ldsc \n")
  baseline_ldsc <- paste0("python /home/opt/software/anaconda3/envs/ldsc/bin/ldsc.py",
                          " --ref-ld-chr /data2/projects/bioinfo/xingjie/src/ldsc/baselineLD_v1.1/baselineLD.",
                          " --w-ld-chr /data2/projects/bioinfo/xingjie/src/ldsc/eur_w_ld_chr/",
                          " --not-M-5-50",
                          " --overlap-annot",
                          " --print-coefficients",
                          " --h2 ", out_p, trait.name, ".sumstats",
                          " --out ", ldsc_p, trait.name, "_baselineLD",
                          " > ", ldsc_p, trait.name, ".baselineLD.temp")
  system(baseline_ldsc)
  
  for(ct_index in 1:10){
    cat("ldsc cell type :", ct_index ," \n")
    celltype_ldsc <- paste0("python /home/opt/software/anaconda3/envs/ldsc/bin/ldsc.py",
                            " --ref-ld-chr /data2/projects/bioinfo/xingjie/src/ldsc/baselineLD_v1.1/baselineLD.",
                            ",/data2/projects/bioinfo/xingjie/src/ldsc/1000G_Phase3_cell_type_groups/cell_type_group.", ct_index,".",
                            " --w-ld-chr /data2/projects/bioinfo/xingjie/src/ldsc/eur_w_ld_chr/",
                            " --not-M-5-50",
                            " --overlap-annot",
                            " --print-coefficients",
                            " --h2 ", out_p, trait.name, ".sumstats",
                            " --out ", ldsc_p, trait.name,"_baselineLD_cell_type_group.", ct_index,
                            " > ", ldsc_p, trait.name, ct_index, ".baselineLD.temp")
    system(celltype_ldsc)
  }
  
  if(T){
    setwd(ldsc_p)
    h2_est <- fread(paste0(trait.name,".log"),fill= T) %>% paste
    h2_est <- strsplit(h2_est,split = "Total Observed scale h2: ")[[1]][2] 
    h2_est <- strsplit(h2_est,"[ (]")[[1]][1] %>% as.numeric()
    cat("Total Observed scale h2 of ", trait.name," is ",h2_est,"\n")
    fwrite(as.data.frame(h2_est),paste0(out_p, trait.name,"-h2.txt"), sep="\t",quote=F, row.names=F)
    
    ctg_names <- read.table("/data2/projects/bioinfo/zhshao/GWAS.summary/ctnames", header = T)
    hm3snp <- data.frame(fread("/data2/projects/bioinfo/zhshao/GWAS.summary/snpinfo_mult_1kg_hm3.snp", header = F))
    for(chr_index in 1:22) {
      baseline_ann_dat <- data.frame(fread(paste0("/data2/projects/bioinfo/xingjie/src/ldsc/baselineLD_v1.1/baselineLD.",
                                                  chr_index, ".annot.gz"), header = F))
      hm3snp_index <- which(baseline_ann_dat[, 3] %in% hm3snp[, 1])
      baseline_ann_dat <- baseline_ann_dat[hm3snp_index, ]
      ##
      est_snp_vg_mat <- NULL
      ##
      baseline_coef <- read.table(paste0(trait.name, "_baselineLD.results"), header = T)
      est_snp_vg <- as.matrix(baseline_ann_dat[, -c(1:4)]) %*% 
        matrix(baseline_coef$Coefficient/ h2_est, 
               length(baseline_coef$Coefficient), 1)
      est_snp_vg_mat <- cbind(est_snp_vg_mat, est_snp_vg)
      
      for(ctg_index in 1:10) {
        ctg_dann_dat <- data.frame(fread(paste0("/data2/projects/bioinfo/xingjie/src/ldsc/1000G_Phase3_cell_type_groups/cell_type_group.",
                                                ctg_index, ".", chr_index, ".annot.gz"), header = T))
        ctg_dann_dat <- ctg_dann_dat[hm3snp_index, ]
        #
        baseline_ctg_coef <- read.table(paste0(trait.name, "_baselineLD_cell_type_group.", 
                                               ctg_index, ".results"), header = T)
        #
        est_snp_vg <- as.matrix(cbind(baseline_ann_dat[, -c(1:4)], ctg_dann_dat[, -c(1:4)])) %*% 
          matrix(baseline_ctg_coef$Coefficient / h2_est, length(baseline_ctg_coef$Coefficient), 1)
        est_snp_vg_mat <- cbind(est_snp_vg_mat, est_snp_vg)
        # cat("The per snp h2/genetic variance has been calculated for Chr ", 
        #     chr_index, " and cell-type group", ctg_names[ctg_index, 2], fill = T)
      }
      colnames(est_snp_vg_mat) <- c("baseline", ctg_names[, 2])
      est_snp_vg_mat <- cbind(ctg_dann_dat[, 1:3], est_snp_vg_mat)
      ## output
      fwrite(est_snp_vg_mat, file = paste0(out_p, trait.name, ".snpvg"), sep = "\t", append = T)
      # cat("The per snp h2/genetic variance has been calculated for Chr ", chr_index, fill = T)
      # cat("########", fill = T)
    }
  }# calculate SNP Vg
  
  cat("Summary data preparation has done!")
}

gen.inp <- function(summs, trait.name = "Height",
                    out = "/data2/projects/bioinfo/zhshao/GWAS.summary/sums/",
                    pheno, phe_trait){
  require(data.table)
  require(magrittr)
  
  out_p <- paste0(out, trait.name,"/")
  dir.check(out_p); setwd(out_p)
  
  traits <- trait.name
  data <- fread(paste0(out_p,traits,".sumstats"))
  
  if(F){
    cat("generating DBSLMM input data\n")
    sumstats <- data
    sumstats$n_mis <- max(sumstats$N)-sumstats$N
    N <- max(sumstats$N)
    sumstats <- sumstats[,c("CHR","SNP","BP","n_mis","N","A1","A2","MAF","Beta","SE","P")]
    fwrite(sumstats,paste0(out_p, traits, ".DBSLMM")
           ,sep="\t",quote=F,row.names=F)
    
    cat("generating PRS.CS input data\n")
    sumstats <- data
    sumstats <- data.frame(sumstats)[,c("SNP","A1","A2","Beta","P")]
    colnames(sumstats) <- c("SNP","A1","A2","BETA","P")
    fwrite(sumstats,paste0(out_p, traits, ".PRS.CS")
           ,sep="\t",quote=F,row.names=F)
    
    cat("generating lassosum input data\n")
    require(lassosum)
    ss <- data
    ss <- ss[!P == 0]
    cor <- p2cor(p = as.numeric(ss$P),
                 n = as.numeric(ss$N),
                 sign = as.numeric(ss$Beta))
    cor[which(is.na(cor))] = 0
    save(ss,file=paste0(out_p, traits, "ss.lassosum.RData"))
    save(cor,file=paste0(out_p, traits, "cor.lassosum.RData"))
    
    cat("generating LDpred.funct input data\n")
    p = fread(pheno) %>% data.frame()
    p = p[,c(1,which(colnames(p) == phe_trait))]
    fwrite(p,paste0(out_p,traits,"_trait.txt"),
           sep="\t",quote=F, row.names=F, col.names=F)
    h2 <- fread(paste0(out_p, traits, ".snpvg"))
    bim <- fread("/data2/projects/bioinfo/zhshao/UKB_all/EUR_all.bim")
    h2_est <- h2[h2$SNP %in% bim$V2, c(3,4)]
    fwrite(h2_est,paste0(out_p,traits,"_fc.txt"),
           sep="\t",quote=F, row.names=F)
    ss = as.data.frame(data)
    ss = ss[,c("CHR","SNP","BP","A1","A2","P","Beta","Z")] # colnames(ss) %in% 
    colnames(ss)[7] = "BETA"
    fwrite(ss,paste0(out_p,traits,"_sums.txt"),
           sep="\t", quote=F, row.names=F)
    
    cat("generating AnnoPred input data\n")
    sums <- data
    sums <- sums[,c("CHR","SNP","A1","A2","BP","Beta","P")]
    sums$CHR <- paste0("chr",sums$CHR)
    sums$Beta <- exp(sums$Beta)
    colnames(sums) <- c("hg19chrc","snpid","a1","a2","bp","or","p")
    fwrite(sums,paste0(out_p, traits, ".AnnoPred")
           ,sep="\t",quote=F,row.names=F)
    
    cat("generating SDPR input data\n")
    sums <- data
    sums <- sums[,c("SNP","A1","A2","Beta","P")]
    colnames(sums)[4] <- "BETA"
    fwrite(sums,paste0(out_p, traits, ".SDPR")
           ,sep="\t",quote=F,row.names=F)
  }
  
  if(F){
    for(chr in 1:22){
      SDPR.ref <- paste0("/home/zhshao/project/GRS/SDPR/software/SDPR",
                         " -make_ref",
                         " -ref_prefix /data2/projects/bioinfo/zhshao/hm3eur/eur_hm3_chr",chr,
                         " -chr ",chr,
                         " -ref_dir /data2/projects/bioinfo/zhshao/GRS/SDPR/ref/")
      system(SDPR.ref)
    }
  }
  
  cat("All done\nNow step into PRS procedure\n")
}

# auxiliary function ------------------------------------------------------
plink.R2 <- function(rf = "internal_pred_.profile",
                     GRS.input,traits,chr,temp_p,
                     pheno,phe_trait,out_p){
  system(paste0("touch /home/zhshao/out/temp/",class(GRS.input),".done.",traits,"_",chr))
  cat(class(GRS.input)," in chr",chr,"done\n")
  
  files <- list.files(path = "/home/zhshao/out/temp/",pattern = paste0(class(GRS.input),".done."))
  if(length(files) == 22){
    cat("Begin merge 22 chromosome result\n")
    options(scipen=100)
    
    res <- NULL
    for(c in 1:22){
      d1 = paste0(temp_p, "/", class(GRS.input),"-",c,rf) %>% fread(.) %>%
        .[,c("FID", "SCORESUM")]
      if(c == 1){
        res <- d1[,c("FID","SCORESUM")]
      }else{
        dd = merge(res, d1[,c("FID","SCORESUM")], by="FID")
        dd$SCORESUM <- dd$SCORESUM.x + dd$SCORESUM.y
        res <- dd[,c("FID","SCORESUM")]
      }
    }
    
    cat("Calculating R2\n")
    phen <- fread(pheno)
    data <- merge(phen, res, by.x = "UDI", by.y = "FID") %>% data.frame()
    R2 <- cor(data[,which(colnames(data) == phe_trait)], data$SCORESUM, method='pearson')**2
    fwrite(as.data.frame(1),paste0(out_p,"/", class(GRS.input),"-", R2))
    cat(class(GRS.input)," method is done!\n")
  }
  cat("~\n")
  cat("~~\n")
  cat("~~~\n")
  cat("~~\n")
  cat("~\n")
}

score.block <- function(traits, chr, p_sums, 
                        p_ref, target_p, target_f, clu, temp_path, out_path){
  cat("start analysis. trait:", traits," CHR", chr,"\n")
  out_path_ef <- paste0(out_path,"effect/")
  if(!dir.exists(out_path_ef)) dir.create(out_path_ef)
  out_path_prs <- paste0(out_path,"score/")
  if(!dir.exists(out_path_prs)) dir.create(out_path_prs)
  
  sums <- fread(paste0(p_sums, traits, ".sumstats"))
  n <- max(sums$N)
  # block <- data.frame(fread(paste0(p_block,"chr",chr,".bed")))
  Vg <- fread(paste0(p_sums, traits, ".snpvg"))
  load(paste0(p_ref, "chr", chr, ".RData"))
  sums <- sums[,c("SNP","A1","CHR","BP","N","Beta")]
  data1 <- sums[CHR == chr] %>% 
    merge(., Vg, by = c("SNP", "CHR", "BP"))
  
  cat("start estimate joint effect. A total of",length(LD_list),"LD blocks are calculated in chr", chr,"\n")
  for(i in 1:length(LD_list))try({
    LD_R <- cbind(1:length(LD_list[[i]][[2]]),
                  LD_list[[i]][[2]],
                  LD_list[[i]][[1]]) %>% data.frame
    colnames(LD_R)[1:2] <- c("index", "SNP")
    index <- which(LD_R$SNP %in% data1$SNP)
    LD_R <- LD_R[index, c(1,2,index+2)]

    data <- merge(data1, LD_R, by = "SNP") %>%
      .[order(as.numeric(.$index))]

    joint_est(data, n, out_path_ef ,clu)
    cat("Block",i,"finish\n")
  })
  cat("Joint effect estimation done\n",
      "Start calculate score by using Plink1.9\n")
  
  anno <- c("baseline", "AdrenalPancreas", "Cardiovascular",
            "CNS", "ConnectiveBone", "GI", "Hematopoietic",
            "Kidney", "Liver", "Other", "SkeletalMuscle")
  for(method in c("BLUP.main")){
    d1 = fread(paste0(out_path_ef, 
                      method, "_chr",chr,".txt"))
    d1 = d1[!duplicated(d1$SNP),]
    fwrite(d1,paste0(out_path_ef, 
                     method, "_chr",chr,".txt"),sep = "\t", quote=F, row.names=F)
    for(i in 1:length(anno)){
      print(paste(anno[i],method))
      
      ii = 4 + i
      prs <- paste0("/home/opt/software/plink_linux_x86_64_20210606/plink",
                    " --bfile ", target_p, target_f, chr,
                    " --score ", out_path_ef, 
                    method, "_chr",chr,".txt 1 4 ", ii," header sum",
                    " --out ", out_path_prs, 
                    method,"_",traits,"_chr",chr,"_anno_",anno[i],
                    "> ",temp_path, traits, method, anno[i], chr,".log")
      system(prs)
    }
  }
  cat("Polygenic score for baseline and tissue-type models are finished\ndone\n")
}

ETM_beta <- function(f_data, n){
  if(!file.exists(f_data)) stop("make sure file exists!")
  load(f_data)
  if(nrow(data) == 1){
    return(NULL)
  }else{
    b_hat <- data$Beta %>% as.matrix() %>% apply(.,2,as.numeric)
    index1 <- which(colnames(data) == "baseline")
    index2 <- which(colnames(data) == "SkeletalMuscle")
    R_square <- data[,-(1:(index2+1))] %>% 
      as.matrix() %>% apply(.,2,as.numeric)
    h2 <- data[,index1:index2] %>% apply(.,2,function(x){
      # x[which(x <= 0)] = 0;x # approaxmite 0
      c <- mean(x) / mean(x[x > 0])
    }) %>% data.frame()
    res.BLUP <- parSapply(clu, h2, BLUP.main, R_square, n, b_hat) %>% 
      cbind(data[,c("SNP","CHR","BP","A1")],.) %>% data.frame
    res.Bayes.Ana <- parSapply(clu, h2, Bayes.Ana, R_square, n, b_hat) %>% 
      cbind(data[,c("SNP","CHR","BP","A1")],.) %>% data.frame
    res.Bayes.MCMC <- parSapply(clu, h2, Bayes.MCMC, R_square, n, b_hat) %>% 
      cbind(data[,c("SNP","CHR","BP","A1")],.) %>% data.frame
    return(list(res.BLUP, res.Bayes.Ana, res.Bayes.MCMC))
  }
}

joint_est <- function(data, n, out_path_ef, clu){
  if(nrow(data) == 1){
    return(NULL)
  }else{
    b_hat <- data$Beta %>% as.matrix() %>% apply(.,2,as.numeric)
    index1 <- which(colnames(data) == "baseline")
    index2 <- which(colnames(data) == "SkeletalMuscle")
    R_square <- data[,-(1:(index2+1))] %>% 
      as.matrix() %>% apply(.,2,as.numeric)
    h2 <- data[,index1:index2] %>% apply(.,2,function(x){
      # x[which(x <= 0)] = 0;x # approaxmite 0
      c <- sum(x) / sum(x[x > 0])
      x[which(x <= 0)] = 0
      x <- x * c
    }) %>% data.frame()
    for(method in c("BLUP.main")){
      # res <- paste0("parApply(clu, h2,2,", method,", R_square, n, b_hat)") %>% 
      #   parse(text=.) %>% eval() %>% cbind(data[,c("SNP","CHR","BP","A1")],.) %>% data.frame
      res <- paste0("apply(h2,2,", method,", R_square, n, b_hat)") %>% 
        parse(text=.) %>% eval() %>% cbind(data[,c("SNP","CHR","BP","A1")],.) %>% data.frame
      file.name <- paste0(out_path_ef, method, "_chr",chr,".txt")
      if(!file.exists(file.name)) rr1 <- res else rr1 <- fread(file.name) %>% rbind(.,res)
      rr1 <- data.frame(rr1)
      write.table(rr1, file.name, quote=F, row.names=F)
    }
  }
}

BLUP.main <- function(h2, R_square, n, b_hat){
  # Dg <- h2 %>% as.matrix %>% as.numeric() %>% diag %>% 
  #   as.matrix() %>% apply(.,2,as.numeric)
  index <- (h2 != 0)
  BETA <- rep(0,length(h2))
  
  h2 <- h2[index]
  R_square <- R_square[index, index]
  b_hat <- b_hat[index]
  
  Dg <- apply(as.matrix(diag(as.numeric(as.matrix(h2)))),2,as.numeric)
  Dg_inv <- solve(Dg)
  eta <- solve(Dg_inv + n * R_square)
  beta_est <- n * Dg %*% (diag(nrow(Dg)) - n * R_square %*% eta) %*% b_hat
  
  BETA[index] <- beta_est
  return(BETA)
}

Bayes.Ana <- function(h2, R_square, n, b_hat){
  index <- (h2 != 0)
  BETA <- rep(0,length(h2))
  
  h2 <- mean(h2[index])
  R_square <- R_square[index, index]
  b_hat <- b_hat[index]
  
  A = n * R_square + diag(1/h2, sum(index))
  beta_est = (solve(A) * n) %*% b_hat  # Adjust the beta_hats
  
  BETA[index] <- beta_est
  return(BETA)
}

Bayes.MCMC <- function(h2, R_square, n, b_hat, 
                       Pi = 0.1, ld_radius=100, num_iter=60, burn_in=5, zero_jump_prob=0.05){
  index <- (h2 != 0)
  BETA <- rep(0,length(h2))
  
  h2 <- h2[index]
  R_square <- R_square[index, index]
  b_hat <- b_hat[index]
  
  A = n * R_square + diag(1/mean(h2), sum(index))
  updated_betas = (solve(A) * n) %*% b_hat  # Adjust the beta_hats
  
  m = length(b_hat)
  curr_betas = updated_betas
  curr_post_means = rep(0,m)
  avg_betas = rep(0,m)
  
  for(k in 1:num_iter){
    h2_est = max(0.00001, sum(curr_betas ** 2))
    alpha = min(1-zero_jump_prob, 1.0 / h2_est, (sum(h2) + 1 / sqrt(n)) / h2_est)
    rand_ps = runif(m)
    for(snp_i in 1:m){
      hdmp = (h2[snp_i]/Pi)#(sum(h2) / Mp)
      hdmpn = hdmp + 1/n
      hdmp_hdmpn = (hdmp / hdmpn)
      c_const = (Pi / sqrt(hdmpn))
      d_const = (1 - Pi) / (sqrt(1/n))
      
      #Local LD matrix
      D_i = R_square[,snp_i]
      
      #Local (most recently updated) effect estimates
      local_betas = curr_betas
      
      #Calculate the local posterior mean, used when sampling.
      res_beta_hat_i = b_hat[snp_i] - D_i %*% local_betas
      b2 = res_beta_hat_i ** 2
      
      d_const_b2_exp = d_const * exp(-b2 / (2.0*1/n))
      numerator = c_const * exp(-b2 / (2.0 * hdmpn))
      if(numerator == 0){
        postp = 0
      }else{
        postp = numerator / (numerator + d_const_b2_exp)
      }
      curr_post_means[snp_i] = hdmp_hdmpn * postp * res_beta_hat_i
      
      if(rand_ps[snp_i] < postp * alpha){
        #Sample from the posterior Gaussian dist.
        proposed_beta = rnorm(1, 0, (hdmp_hdmpn) * 1/n) + hdmp_hdmpn * res_beta_hat_i
      }else{
        proposed_beta = 0
      }
      
      curr_betas[snp_i] = proposed_beta  #UPDATE BETA
    }
  }
  if(k >= burn_in) avg_betas = avg_betas + curr_post_means 
  #Averaging over the posterior means instead of samples.
  
  avg_betas = avg_betas/(num_iter-burn_in)
  
  BETA[index] <- avg_betas
  return(BETA)
}

lasso <- function(data, lambdas=exp(seq(log(0.001), log(0.1), length.out = 20))){
  time_start <- Sys.time()
  
  Y <- data[,1]; X <-  as.matrix(data[,-1])
  # lambdas=exp(seq(log(0.001), log(0.1), length.out = 20))
  lasso_model <- cv.glmnet(X,Y,alpha = 1,lambda = lambdas, nfolds = 10)
  lambda.min <- lasso_model$lambda.min
  lasso_best <- glmnet(X,Y,alpha = 1,lambda = lambda.min)
  beta <- as.matrix(coef(lasso_best))
  PRS <- cbind(1, X) %*% beta
  
  exc_time <- difftime(Sys.time(),time_start,units = 'mins')
  return(list(beta = beta, PRS = PRS, time = exc_time))
}

elasticnet <- function(data, alphas = seq(0,1,by=0.1),
                       lambdas = exp(seq(log(0.001), log(0.1), length.out = 20))){
  time_start <- Sys.time()
  
  Y <- data[,1]; X <-  as.matrix(data[,-1])
  alphas = seq(0,1,by=0.1)
  lambdas = exp(seq(log(0.001), log(0.1), length.out = 20))
  control <- caret::trainControl(method = "repeatedcv",
                                 number = 5, repeats = 5,
                                 search = "random", verboseIter = TRUE)
  elastic_model <- caret::train(y ~ ., data = data,
                                method = "glmnet",
                                preProcess = c("center", "scale"),
                                tuneLength = 10,
                                trControl = control,
                                tuneGrid = expand.grid(alpha = alphas,
                                                       lambda = lambdas))
  elastic_best <- glmnet(X,Y,alpha = elastic_model$bestTune$alpha, 
                         lambda = elastic_model$bestTune$lambda)
  beta <- as.matrix(coef(elastic_best))
  PRS <- cbind(1, X) %*% beta
  
  exc_time <- difftime(Sys.time(),time_start,units = 'mins')
  return(list(beta = beta, PRS = PRS, time = exc_time))
}

BMA <- function(data){
  time_start <- Sys.time()
  
  fit <- bas.lm(y~., data = data,
                method = "MCMC", MCMC.iterations = 10000,
                prior = "g-prior", alpha = sqrt(nrow(data)))
  sum.fit <- summary(fit)
  beta <- sum.fit[1:12, 1] * sum.fit[1:12, 2]
  PRS <- cbind(1, as.matrix(data[,-1])) %*% as.matrix(beta)
  
  exc_time <- difftime(Sys.time(),time_start,units = 'mins')
  return(list(beta = beta, PRS = PRS, time = exc_time))
}

Additive <- function(data){
  beta = rep(1,ncol(data)-1)
  return(list(beta = beta))
}

fold_fun <- function(fold, sa, method, data){
  index_sa <- (1+(fold - 1)*sa):min(fold*sa, nrow(data))
  train <- data[-index_sa,]
  test <- data[index_sa,]
  beta.train <- paste0(method, "(data = train)$beta") %>% parse(text=.) %>% eval()
  if(method == "Additive"){
    R2 <- cor(test$y, as.matrix(test[,-1]) %*% as.matrix(beta.train), method = "pearson")
  }else{
    R2 <- cor(test$y, cbind(1,as.matrix(test[,-1])) %*% as.matrix(beta.train), method = "pearson")
  }
  R2
}

crossV <- function(fold = 5, data, method = "LAsso"){
  time_start <- Sys.time()
  cat(paste0("start ", fold,"-fold cross validation using method ", method), '\n')
  sa <- round(1/fold*nrow(data), digits = 0)
  # res <- parSapply(clu, 1:fold, fold_fun, sa, method, data)
  resu <- NULL
  for(f in 1:fold){
    res <- fold_fun(f, sa, method, data)
    resu <- c(resu, res)
    cat(paste0("fold ", f," done !"), '\n')
  }
  exc_time <- difftime(Sys.time(),time_start,units = 'mins')
  return(list(result = resu, time = exc_time))
}

res.pro <- function(out){
  require(magrittr)
  res <- NULL
  traits <- list.files(out)
  for(i in 1:length(traits)){
    files <- list.files(paste0(out,traits[i]))
    files <- files[-which(files %in% c("effect","score","RR.csv"))]
    files <- strsplit(files,split = "-") %>% unlist() %>% matrix(.,ncol=2,byrow = T)
    res <- cbind(res,files)
  }
  fwrite(as.data.frame(res),paste0(out,"result.csv"))
}

# GRS model ---------------------------------------------------------------

if(F){
  library(data.table)
  library(magrittr)
  
  traits = "height"
  GRS.input <- NA; attr(GRS.input,"class") <- "omniPRS"
  chr = 1
  sums_p = "/data2/projects/bioinfo/zhshao/GWAS.summary/sums/height/"
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

GRS <- function(x, ...) UseMethod("GRS")

GRS.CT <- function(GRS.input, traits, chr, sums_p, base_p, base_f,
                   target_p, target_f, pheno, phe_trait, out, temp){
  cat("Start C+T method \n")
  out_p <- paste0(out,traits); dir.check(out_p)
  temp_p <- paste0(temp,traits); dir.check(temp_p)

  setwd("/data2/projects/bioinfo/zhshao/GRS/P+T")
  clumping <- paste0("/home/opt/software/plink_linux_x86_64_20210606/plink",
                     " --bfile ",base_p, base_f, chr,
                     " --clump-p1 1",
                     " --clump-r2 0.1",
                     " --clump-kb 1000",
                     " --clump ", sums_p, traits,".sumstats",
                     " --clump-snp-field SNP",
                     " --clump-field P",
                     " --out /data2/projects/bioinfo/zhshao/GRS/P+T/",traits,chr,
                     " > ", temp_p, "/", class(GRS.input),"-",chr)
  system(clumping)
  
  va.snp <- paste0("awk 'NR!=1{print $3}' /data2/projects/bioinfo/zhshao/GRS/P+T/",
                   traits, chr,".clumped >  ", traits, chr,".valid.snp; ",
                   "awk '{print $1,$10}' ", sums_p, traits, ".sumstats"," > ",
                   traits, chr,"SNP.pvalue")
  system(va.snp)
  
  prs <- paste0("/home/opt/software/plink_linux_x86_64_20210606/plink",
                " --bfile ", target_p, target_f, chr,
                " --score ", sums_p, traits,".sumstats 1 2 7 header sum",
                " --q-score-range /data2/projects/bioinfo/zhshao/GRS/range_list ",
                traits, chr, "SNP.pvalue",
                " --extract ", traits, chr,".valid.snp",
                " --out ", temp_p, "/", class(GRS.input),"-",chr,
                " > ", temp_p, "/", class(GRS.input),"-",chr)
  system(prs)
  
  # file.remove(paste0(traits,chr,c(".clumped",".log",".nopred",".valid.snp","SNP.pvalue",".nosex")))
  
  options(scipen=100)
  plink.R2(".0.00000005.profile",GRS.input,traits,chr,temp_p,pheno,phe_trait,out_p)
}

GRS.DBSLMM <- function(GRS.input, traits, chr, sums_p, base_p, base_f,
                       target_p, target_f, pheno, phe_trait, out, temp){
  cat("Start ",class(GRS.input)," method \n")
  out_p <- paste0(out,traits); dir.check(out_p)
  temp_p <- paste0(temp,traits); dir.check(temp_p)

  library(data.table)
  library(magrittr)
  source("/home/zhshao/software/DBSLMM/DBSLMM_szh.R")
  
  if(!dir.exists(paste0("/home/zhshao/project/GRS/DBSLMM/result/", traits))) 
    dir.create(paste0("/home/zhshao/project/GRS/DBSLMM/result/", traits))
  
  N <- fread(paste0(sums_p,traits,"-N.txt"))
  h2 <- fread(paste0(sums_p,traits,"-h2.txt"))
  
  DBSLMM(summary = paste0(sums_p, traits, ".",class(GRS.input)),
         ref = paste0(base_p,base_f,chr),
         test = paste0(target_p,target_f,chr),
         block = paste0("/home/zhshao/software/DBSLMM/block_data/EUR/chr",chr,".bed"),
         outPath = paste0(temp_p, "/", class(GRS.input),"-",chr),
         n = N, h2f = h2, GRS.input = GRS.input, temp_p = temp_p,
         plink = "/home/opt/software/plink_linux_x86_64_20210606/plink",
         dbslmm = "/home/zhshao/software/DBSLMM/software/dbslmm",
         model = "DBSLMM",type = "auto", mafMax = "0.2",thread = "5")
  
  plink.R2("internal_pred_.profile",GRS.input,traits,chr,temp_p,pheno,phe_trait,out_p)
}

GRS.PRS.CS <- function(GRS.input, traits, chr, sums_p, base_p, base_f,
                       target_p, target_f, pheno, phe_trait, out, temp){
  cat("Start ",class(GRS.input)," method \n")
  out_p <- paste0(out,traits); dir.check(out_p)
  temp_p <- paste0(temp,traits); dir.check(temp_p)

  library(data.table)
  library(magrittr)
  
  N <- fread(paste0(sums_p,traits,"-N.txt"))
  
  PRS.CS <- paste0("python /home/zhshao/software/PRScs/PRScs.py",
                   " --ref_dir=/data2/projects/bioinfo/zhshao/REF/1000G/ldblk_1kg_eur",
                   # " --ref_dir=/data2/projects/bioinfo/zhshao/REF/UKB/ldblk_ukbb_eur",
                   " --bim_prefix=",target_p, target_f, chr,
                   " --sst_file=",sums_p, traits, ".",class(GRS.input),
                   " --n_gwas=",N,
                   " --chrom=",chr,
                   " --phi=1e-2",
                   " --out_dir=/home/zhshao/project/GRS/PRS-CS/out/",traits)
  system(PRS.CS)
  
  d1 = fread(paste0("/home/zhshao/project/GRS/PRS-CS/out/",
                    traits,"_pst_eff_a1_b0.5_phi1e-02_chr",chr,".txt"),fill = T)
  d1 = na.omit(d1) %>% .[!duplicated(.$V2),]
  fwrite(d1,paste0("/home/zhshao/project/GRS/PRS-CS/out/",
                   traits,"_pst_eff_a1_b0.5_phi1e-02_chr",chr,".txt"),sep="\t",quote=F,row.names=F)
  
  prs <- paste0("/home/opt/software/plink_linux_x86_64_20210606/plink",
                " --bfile ",target_p, target_f, chr,
                " --score /home/zhshao/project/GRS/PRS-CS/out/",
                traits,"_pst_eff_a1_b0.5_phi1e-02_chr",chr,".txt 2 4 6 header sum",
                " --out ", temp_p, "/", class(GRS.input),"-",chr,
                " > ", temp_p, "/", class(GRS.input),"-",chr)
  system(prs)
  
  plink.R2(".profile",GRS.input,traits,chr,temp_p,pheno,phe_trait,out_p)
}

GRS.LDpred2 <- function(GRS.input, traits, chr, sums_p, base_p, base_f,
                        target_p, target_f, pheno, phe_trait, out, temp){
  cat("Start LDpred2 method \n")
  out_p <- paste0(out,traits); dir.check(out_p)
  temp_p <- paste0(temp,traits); dir.check(temp_p)
  
  require(data.table)
  require(magrittr)
  require(bigsnpr)
  options(bigstatsr.check.parallel.blas = FALSE)
  options(default.nproc.blas = NULL)
  NCORES <- nb_cores()
  
  if(F){
    cat("generating LDpred2 input \n")
    sumstats <- bigreadr::fread2(raw_s)
    sumstats <- sumstats[,-9]
    colnames(sumstats)[1:9] <- c("rsid", "a1", "a0",
                                 "chr", "pos", "n_eff",
                                 "beta", "beta_se", "p")
    setwd("/data2/projects/bioinfo/zhshao/GRS/LDpred2")
    NCORES <- nb_cores()
    tmp <- tempfile(tmpdir = "tmp-data")
    on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
    corr <- NULL
    ld <- NULL
    info_snp <- NULL
    fam.order <- NULL
    for (chr in 1:22) {
      # snp_readBed(paste0("/data2/projects/bioinfo/zhshao/hm3eur/eur_hm3_chr",chr,".bed"))
      obj.bigSNP <- snp_attach(paste0(base_p,base_f,chr,".rds"))
      map <- obj.bigSNP$map[-3]
      names(map) <- c("chr", "rsid", "pos", "a1", "a0")
      tmp_snp <- snp_match(sumstats[sumstats$chr==chr,], map)
      info_snp <- rbind(info_snp, tmp_snp)
      genotype <- obj.bigSNP$genotypes
      CHR <- map$chr
      POS <- map$pos
      POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
      ind.chr <- which(tmp_snp$chr == chr)
      ind.chr2 <- tmp_snp$`_NUM_ID_`[ind.chr]
      corr0 <- snp_cor(
        genotype,
        ind.col = ind.chr2,
        ncores = NCORES,
        infos.pos = POS2[ind.chr2],
        size = 3 / 1000
      )
      if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, tmp)
      } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
      }
      cat("chr",chr,"\n")
    }
    
    df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
    ldsc <- snp_ldsc(   ld, 
                        length(ld), 
                        chi2 = (df_beta$beta / df_beta$beta_se)^2,
                        sample_size = df_beta$n_eff, 
                        blocks = NULL)
    h2_est <- ldsc[["h2"]]
    
    save(info_snp,file = paste0(sums_p,traits,"/info_snp.RData"))
    save(ld,file = paste0(sums_p,traits,"/ld.RData"))
    save(corr,file = paste0(sums_p,traits,"/corr.RData"))
    fwrite(as.data.frame(h2_est),paste0(sums_p,traits,"/h2_est"))
  }
  
  load(paste0(sums_p,traits,"/info_snp.RData"))
  load(paste0(sums_p,traits,"/ld.RData"))
  load(paste0(sums_p,traits,"/corr.RData"))
  h2_est <- fread(paste0(sums_p,traits,"/h2_est")) %>% as.numeric()
  df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
  
  if(T){
    cat("Start inf model\n")
    beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
    data <- cbind(info_snp[,5:4], beta_inf)
    write.table(data, paste0("/home/zhshao/project/GRS/LDpred2/result/",traits,"beta_inf.txt"), 
                quote=F, row.names=F)
    cat("inf beta done\n")
    
    prs <- paste0("/home/opt/software/plink_linux_x86_64_20210606/plink",
                  " --bfile /data2/projects/bioinfo/zhshao/UKB_all/EUR_all",
                  " --score /home/zhshao/project/GRS/LDpred2/result/",traits,"beta_inf.txt 1 2 3 header sum",
                  " --out /data2/projects/bioinfo/zhshao/GRS/LDpred2/score/",traits,"inf")
    system(prs)
    
    res <- fread(paste0("/data2/projects/bioinfo/zhshao/GRS/LDpred2/score/",traits,"inf.profile"))
    phen <- fread(pheno)
    data <- merge(phen, res, by.x = "UDI", by.y = "FID") %>% data.frame()
    R2 <- cor(data[,which(colnames(data) == phe_trait)], data$SCORESUM, method='pearson')**2
    fwrite(as.data.frame(1),paste0(out_p,"/LDpred.inf-", R2))
    cat("LDpred inf model is done!\n")
  }
  
  if(T){
    cat("Start auto model\n")
    multi_auto <- snp_ldpred2_auto(
      corr,
      df_beta,
      h2_init = h2_est,
      vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
      ncores = NCORES
    )
    beta_auto <- sapply(multi_auto, function(auto)
      auto$beta_est)
    data <- cbind(info_snp[,5:4], beta_auto)
    write.table(data, paste0("/home/zhshao/project/GRS/LDpred2/result/",traits,"beta_auto.txt"), 
                quote=F, row.names=F)
    cat("auto beta done\n")
    
    prs <- paste0("/home/opt/software/plink_linux_x86_64_20210606/plink",
                  " --bfile /data2/projects/bioinfo/zhshao/UKB_all/EUR_all",
                  " --score /home/zhshao/project/GRS/LDpred2/result/",traits,"beta_auto.txt 1 2 3 header sum",
                  " --out /data2/projects/bioinfo/zhshao/GRS/LDpred2/score/",traits,"auto")
    system(prs)
    
    res <- fread(paste0("/data2/projects/bioinfo/zhshao/GRS/LDpred2/score/",traits,"auto.profile"))
    phen <- fread(pheno)
    data <- merge(phen, res, by.x = "UDI", by.y = "FID") %>% data.frame()
    R2 <- cor(data[,which(colnames(data) == phe_trait)], data$SCORESUM, method='pearson')**2
    fwrite(as.data.frame(1),paste0(out_p,"/LDpred.auto-", R2))
    cat("LDpred auto model is done!\n")
  }
}

GRS.lassosum <- function(GRS.input, traits, chr, sums_p, base_p, base_f,
                         target_p, target_f, pheno, phe_trait, out, temp){
  cat("Start ",class(GRS.input)," method \n")
  out_p <- paste0(out,traits); dir.check(out_p)
  temp_p <- paste0(temp,traits); dir.check(temp_p)
  
  library(lassosum)
  library(data.table)
  library(methods)
  library(magrittr)
  library(parallel)
  cl <- makeCluster(detectCores())
  
  load(paste0(sums_p, traits, "ss.lassosum.RData"))
  load(paste0(sums_p, traits, "cor.lassosum.RData"))
  
  out <- lassosum.pipeline(
    cor = cor, chr = ss$CHR, pos = ss$BP,
    A1 = ss$A1, A2 = ss$A2,
    ref.bfile = paste0(base_p, base_f, chr),
    test.bfile = paste0(target_p, target_f, chr),
    LDblocks = "EUR.hg19", cluster=cl)
  save(out,file = paste0(temp_p,"/lassosum.",chr,".RData"))
  
  system(paste0("touch /home/zhshao/out/temp/",class(GRS.input),".done.",traits,"_",chr))
  cat(class(GRS.input)," in chr",chr,"done\n")
  files <- list.files(path = "/home/zhshao/out/temp/",pattern = paste0(class(GRS.input),".done."))
  if(length(files) == 22){
    
    load(paste0(temp_p,"/lassosum.1.RData"))
    lasso_out <- out
    for(chr in 2:22){
      load(paste0(temp_p,"/lassosum.",chr,".RData"))
      lasso_out <- merge(lasso_out,out)
    }
    target.pheno <- fread(pheno)
    target.pheno$FID <- target.pheno$IID <- target.pheno$UDI
    target.pheno <- data.frame(target.pheno)[,c("FID", "IID", phe_trait)]
    target.res <- validate(lasso_out, pheno = as.data.frame(target.pheno),plot=F)
    fwrite(as.data.frame(1),paste0(out_p,"/", class(GRS.input),"-", max(target.res$validation.table$value)))
    cat(class(GRS.input)," method is done!\n")
  }
  cat("~\n")
  cat("~~\n")
  cat("~~~\n")
  cat("~~\n")
  cat("~\n")
}

GRS.LDpred.funct <- function(GRS.input, traits, chr, sums_p, base_p, base_f,
                             target_p, target_f, pheno, phe_trait, out, temp){
  cat("Start ",class(GRS.input)," method \n")
  out_p <- paste0(out,traits); dir.check(out_p)
  temp_p <- paste0(temp,traits); dir.check(temp_p)
  
  N <- fread(paste0(sums_p,traits,"-N.txt"))
  h2 <- fread(paste0(sums_p,traits,"-h2.txt"))
  
  outCoord = paste0("/data2/projects/bioinfo/zhshao/GRS/LDpred_funct/coord",traits,chr)
  if(file.exists(outCoord)) file.remove(outCoord)
  setwd("/home/zhshao/software/LDpred-funct")
  LDpf <- paste0("python /home/zhshao/software/LDpred-funct/ldpredfunct.py",
                 " --gf=",base_p,base_f,chr,
                 " --pf=",sums_p,traits,"_trait.txt",
                 " --FUNCT_FILE=",sums_p,traits,"_fc.txt",
                 " --coord=",outCoord,
                 " --ssf=",sums_p,traits,"_sums.txt",
                 " --N=",N,
                 " --posterior_means=",temp_p,"/",chr,"_PM",
                 " --H2=",h2,
                 " --out=",temp_p,"/",chr,"_PRS",
                 " > ",outCoord,".log")
  system(LDpf)
  
  prs <- paste0("/home/opt/software/plink_linux_x86_64_20210606/plink",
                " --bfile ",target_p, target_f, chr,
                " --score ", temp_p,"/",chr,"_PM_LDpred-inf-ldscore.txt 3 5 7 header sum",
                " --out ", temp_p, "/", class(GRS.input),"-",chr,
                " > ", temp_p, "/", class(GRS.input),"-",chr)
  system(prs)
  
  plink.R2(".profile",GRS.input,traits,chr,temp_p,pheno,phe_trait,out_p)
}

GRS.AnnoPred <- function(GRS.input, traits, chr, sums_p, base_p, base_f,
                         target_p, target_f, pheno, phe_trait, out, temp){
  cat("Start ",class(GRS.input)," method \n")
  out_p <- paste0(out,traits); dir.check(out_p)
  temp_p <- paste0(temp,traits); dir.check(temp_p)
  
  N <- fread(paste0(sums_p,traits,"-N.txt"))
  
  outCoord = paste0("/data2/projects/bioinfo/zhshao/GRS/AnnoPred/out/coord",traits,chr)
  if(file.exists(outCoord)) file.remove(outCoord)
  setwd("/home/zhshao/software/AnnoPred")
  anno_cmd <- paste0("python /home/zhshao/software/AnnoPred/AnnoPred.py",
                     " --sumstats=", sums_p, traits, ".AnnoPred",
                     " --ref_gt=",base_p,base_f,chr,
                     " --val_gt=",target_p,target_f,chr,
                     " --coord_out=",outCoord,
                     " --N_sample=",N,
                     " --annotation_flag=tier3",
                     " --P=0.1",
                     " --local_ld=/data2/projects/bioinfo/zhshao/GRS/AnnoPred/out/ld_",traits,chr,
                     " --out=/home/zhshao/project/GRS/AnnoPred/out/",traits,chr,
                     " --temp_dir=/data2/projects/bioinfo/zhshao/GRS/AnnoPred/temp")
  system(anno_cmd)
  print(chr)
  
  # plink.R2(".profile",GRS.input,traits,chr,temp_p,pheno,phe_trait,out_p)
}

GRS.SDPR <- function(GRS.input, traits, chr, sums_p, base_p, base_f,
                     target_p, target_f, pheno, phe_trait, out, temp){
  cat("Start ",class(GRS.input)," method \n")
  out_p <- paste0(out,traits); dir.check(out_p)
  temp_p <- paste0(temp,traits); dir.check(temp_p)
  
  N <- fread(paste0(sums_p,traits,"-N.txt"))
  
  SDPR.mcmc <- paste0("/home/zhshao/project/GRS/SDPR/software/SDPR",
                      " -mcmc"," -n_threads 3",
                      " -ref_dir /data2/projects/bioinfo/zhshao/GRS/SDPR/ref/",
                      " -ss ", sums_p, traits, ".SDPR",
                      " -N ", N,
                      " -chr ", chr,
                      " -out ", temp_p,"/SDPR", chr)
  system(SDPR.mcmc)
  
  prs <- paste0("/home/opt/software/plink_linux_x86_64_20210606/plink",
                " --bfile ",target_p, target_f, chr,
                " --score ", temp_p,"/SDPR", chr," 1 2 3 header sum",
                " --out ", temp_p, "/", class(GRS.input),"-",chr,
                " > ", temp_p, "/", class(GRS.input),"-",chr)
  system(prs)
  
  plink.R2(".profile",GRS.input,traits,chr,temp_p,pheno,phe_trait,out_p)
}

GRS.omniPRS <- function(GRS.input, traits, chr, sums_p, base_p, base_f,
                        target_p, target_f, pheno, phe_trait, out, temp){
  cat("Start ",class(GRS.input)," method \n")
  out_path <- paste0(out,traits,"/"); dir.check(out_path)
  temp_path <- paste0(temp,traits,"/"); dir.check(temp_path)
  
  library(parallel)
  clu <- makeCluster(detectCores())
  
  p_ref <- "/data2/projects/bioinfo/xingjie/1000g/hm3eurldmat/"
  p_sums <- sums_p
  
  score.block(traits, chr, p_sums, 
              p_ref, target_p, target_f, clu, temp_path, out_path)
  
  stopCluster(clu)
  system(paste0("touch /home/zhshao/out/temp/",class(GRS.input),".done.",traits,"_",chr))
  cat(class(GRS.input)," in chr",chr,"done\n")
  
  files <- list.files(path = "/home/zhshao/out/temp/",pattern = paste0(class(GRS.input),".done."))
  if(length(files) == 22){
    cat("Begin merge 22 chromosome result\n")
    options(scipen=100)
    
    anno <- c("baseline", "AdrenalPancreas", "Cardiovascular",
              "CNS", "ConnectiveBone", "GI", "Hematopoietic",
              "Kidney", "Liver", "Other", "SkeletalMuscle")
    t = traits
    for(m in c("BLUP.main")){
      ress <- NULL
      index <- 1
      for(a in anno){
        res <- NULL
        for(c in 1:22){
          d1 = paste0(out_path,"score/", m,"_",t,"_chr",c,"_anno_",
                      a,".profile") %>% fread(.) %>% 
            .[,c("FID", "CNT","SCORESUM")]
          if(c == 1){
            res <- d1[,c("FID","SCORESUM")]
          }else{
            dd = merge(res, d1[,c("FID","SCORESUM")], by="FID")
            dd$SCORESUM <- dd$SCORESUM.x + dd$SCORESUM.y
            res <- dd[,c("FID","SCORESUM")]
          }
        }
        if(a == "baseline"){
          ress <- res
        }else{
          ress <- merge(ress, res, by = "FID")
        }
        colnames(ress)[index+1] <- a
        index <- index + 1
        fwrite(ress,paste0("/data2/projects/bioinfo/zhshao/GRS/data/",
                           m,"_",t), quote=F, row.names=F)
      }
    }
    
    cat("Calculating R2\n")
    RR <- NULL
    phen <- fread(pheno)
    for(m in c("BLUP.main")){
      score <- paste0("/data2/projects/bioinfo/zhshao/GRS/data/",m,"_",traits) %>% fread(.)
      data <- merge(phen, score, by.x = "UDI", by.y = "FID") %>% data.frame()
      data_adj <- cbind(y=data[,which(colnames(data) == phe_trait)], 
                        data[,-(2:which(colnames(data) == "baseline")-1)]) %>% data.frame
      for(c in colnames(data_adj)[-1]){
        R2 <- paste0("cor(data_adj$y,data_adj$",c,",method='pearson')") %>% parse(text=.) %>% eval()
        RR <- rbind(RR,c(traits,m,c,R2**2))
      }
      
      require(BAS)
      require(glmnet)
      options(scipen=100)
      
      # baseline + tissue-type model
      Additive.re <- crossV(10, data_adj, "Additive")
      Lasso.re <- crossV(10, data_adj, "lasso")
      # Elasticnet.re <- crossV(10, data_adj, "elasticnet")
      BMA.re <- crossV(10, data_adj, "BMA")
      
      fwrite(as.data.frame(1),paste0(out_path, class(GRS.input),".A-", max(Additive.re$result)**2))
      fwrite(as.data.frame(1),paste0(out_path, class(GRS.input),".L-", max(Lasso.re$result)**2))
      fwrite(as.data.frame(1),paste0(out_path, class(GRS.input),".B-", max(BMA.re$result)**2))
      
      print(paste(traits,m))
    }
    fwrite(RR,paste0(out_path,"RR.csv"))
    cat(class(GRS.input)," method is done!\n")
  }
}
  
# phenotype data procession -----------------------------------------------

if(F){
  library(data.table)
  data <- fread("/data2/projects/bioinfo/xingjie/biobank/ukb/phen/ukbb_select_adj_lipids.phen")
  # turn categorical variable into factor
  data$Sex[data$Sex == "Male"] = 1; data$Sex[data$Sex == "Female"] = 0; data$Sex <- factor(data$Sex)
  data$Smoking[data$Smoking == "Never"] = 0; data$Smoking[data$Smoking == "Prefer not to answer"] = 0
  data$Smoking[data$Smoking == "Previous"] = 1; data$Smoking[data$Smoking == "Current"] = 2
  data$Smoking <- factor(data$Smoking)
  data$Alcohol[data$Alcohol == "Never"] = 0; data$Alcohol[data$Alcohol == "Prefer not to answer"] = 0
  data$Alcohol[data$Alcohol == "Previous"] = 1; data$Alcohol[data$Alcohol == "Current"] = 2
  data$Alcohol <- factor(data$Alcohol)
  
  # QC
  data <- data.frame(data)
  for(i in 2:ncol(data))try({
    x <- data[,i]
    x[is.na(x)] = mean(x,na.rm=T)
    data[,i] <- scale(x)
    # print(i)
  })
  fwrite(data,"/data2/projects/bioinfo/zhshao/GRS/ukbb.phen",
         sep="\t",quote=F, row.names=F)
} # phenotype data pro
