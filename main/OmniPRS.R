#  traits: Trait name (e.g., Height). The prefix of the summary data file should be consistent (e.g., Height_sums.txt). Directory name for output files.
#  chr: Chromosome number (e.g., 1)
#  N: The GWAS sample size
#  h2: Estimated SNP Heritability (precompute this using your favorite method).
#  sums_p: Absolute path to the summary data (e.g., /your/path/sum_file/)
#  base_p: Prefix of Plink format LD data files, without the chromosome number (e.g., for chromosome 1, the file is /your/path/plink/eur_hm3_chr1, so enter /your/path/plink/eur_hm3_chr)
#  target_p: Prefix of Plink format test set data files, without the chromosome number (e.g., for chromosome 1, the file is /your/path/plink/ukb22828_eur_hm3_chr1, so enter /your/path/plink/ukb22828_eur_hm3_chr)
#  pheno: Phenotype file and its absolute path. If covariates are included, they should be in this file (e.g., /your/path/ukbb.phen)
#  phe_trait: Column name of the outcome in the phenotype file (e.g., Height)
#  out: Output path for result files (e.g., /your/path/out/)
#  temp: Output path for temporary files (e.g., /your/path/temp/)
#  cova: Covariates to consider; if none, enter NULL (e.g., c("BaseAge","Sexgenetic"))
#  bina: Whether the outcome is binary data (e.g., T)
#  ct_result: Whether to output 11 tissuespecific PRS prediction levels (e.g., T)
#  software_path: The directory to OmniPRS (e.g., /your/path/)
#  plink_path: The directory to PLINK1.9 and PLINK2.0

GRS.OmniPRS <- function(traits, chr, N, h2 sums_p, base_p,
                        target_p, pheno, phe_trait, out, temp, cova,
                        bina = F,software_path, plink_path){
  t1 = Sys.time()
  cat("Start OmniPRS method \n")
  out_p <- paste0(out,traits); dir.create(out_p)
  temp_p <- paste0(temp,traits); dir.create(temp_p)
  
  snpvg = fread(paste0(sums_p,"snpvg.txt"))
  funct_temp = cbind(snpvg[,3],1) %>% as.data.frame()
  fwrite(funct_temp,paste0(sums_p,"funct.txt"),sep="\t",quote=F, row.names=F)
  phen = fread(pheno)
  phen = phen[,c("UDI",phe_trait)]
  fwrite(phen,paste0(sums_p,"pheno.txt"),sep="\t",quote=F, row.names=F)
  
  rere <- NULL
  outCoord = paste0(temp_p,"/",traits,"_",chr)
  if(file.exists(outCoord)) file.remove(outCoord)
  
  omni <- paste0("python2 ",software_path,"/sum_cord.py",
                 " --gf=",base_p,chr,
                 " --pf=",sums_p,"pheno.txt",
                 " --FUNCT_FILE=",sums_p,"funct.txt",
                 " --coord=",outCoord,
                 " --ssf=",sums_p,"GWAS_sums.txt",
                 " --N=",N,
                 " --H2=",h2,
                 " --out=",temp_p,"/",chr,"_PRS")
  system(omni)
  
  require(rhdf5)
  require(parallel)
  hd5 <- H5Fopen(outCoord)
  hd5_cor <- paste0("hd5$cord_data$chrom_",chr) %>% parse(text=.) %>% eval()
  H5Fclose(hd5)
  t2 = Sys.time()
  
  raw_snps = apply(hd5_cor$raw_snps_ref,2 ,scale) %>% t
  for(ct in 1:11){
    cat("cell type ",ct,"\n")
    sid = hd5_cor$sid; sid = cbind(sid,1:length(sid)) %>% data.frame
    snp_h2 = snpvg[,c(3,3+ct)]
    
    colnames(snp_h2)[2] = "tiss"
    snp_h2 = snp_h2[snp_h2$tiss > 0, ]
    snp_h2 = snp_h2[na.omit(match(sid$sid,snp_h2$SNP)),]
    sid = sid[na.omit(match(snp_h2$SNP,sid$sid)),]
    snp_stds = hd5_cor$snp_stds_ref
    ok_snps_filter = intersect(which(snp_stds > 0),sid$V2) %>% as.numeric()
    snp_stds = snp_stds[ok_snps_filter]
    pos = hd5_cor$positions[ok_snps_filter]
    
    snps = raw_snps[ok_snps_filter,]
    pval_derived_betas = hd5_cor$betas
    pval_derived_betas = pval_derived_betas[ok_snps_filter]
    m = length(pval_derived_betas)
    ld_window=ifelse(round((0.15/100)*m)<100, 100,
                     (round((0.15/100)*m)- round((0.15/100)*m)%%100))
    snp_h2 = snp_h2[,2] %>% as.matrix()
    tot_snp_sum = sum(snp_h2)
    Cval = h2 / tot_snp_sum
    beta_est = rep(NA,m)
    for(wi in seq(1, m, (ld_window))){
      start_i = wi
      stop_i = min(m, wi + ld_window)
      curr_window_size = stop_i - start_i
      if(curr_window_size == 0){
        beta_est[start_i] <- pval_derived_betas[start_i]
      }else{
        X = snps[start_i: stop_i,]
        num_indivs = ncol(X)
        R = X %*% t(X) / num_indivs
        
        Dg_inv <- diag(1/(Cval * snp_h2[start_i: stop_i]))
        beta_est[start_i: stop_i] <- solve(R + Dg_inv/N)%*%
          pval_derived_betas[start_i: stop_i]
      }
    }
    beta_est = beta_est / (snp_stds)
    result <- cbind(hd5_cor$sid[ok_snps_filter],t(hd5_cor$nts)[ok_snps_filter,],beta_est) %>% data.frame()
    fwrite(result,paste0(temp_p,"/",chr,"_",ct,"_PM_Omni.txt"),sep="\t",quote=F,row.names=F)
    
    prs <- paste0(plink_path,
                  " --bfile ",target_p, chr,
                  " --score ", temp_p,"/",chr,"_",ct,"_PM_Omni.txt 1 2 4 header sum",
                  " --out ", temp_p, "/OmniPRS-",chr,"_",ct,
                  " > ", temp_p, "/OmniPRS-",chr,"_",ct)
    system(prs)
    
    system(paste0("touch ",temp_p,"/OmniPRS.done.",traits,"_",ct,"_chr_",chr))
    
    files <- list.files(path = temp_p,
                        pattern = paste0("OmniPRS",".done.",traits,"_",ct,"_chr_"))
    cat(length(files), " chromosome are completed in cell type",ct,"\n")
    if(length(files) == 22){
      cat("Begin merge 22 chromosome result\n")
      options(scipen=100)
      
      res <- NULL
      for(c in 1:22){
        d1 = paste0(temp_p, "/OmniPRS-",c,"_",ct,".profile") %>% fread(.) %>%
          .[,c("FID", "SCORESUM")]
        if(c == 1){
          res <- d1[,c("FID","SCORESUM")]
        }else{
          dd = merge(res, d1[,c("FID","SCORESUM")], by="FID")
          dd$SCORESUM <- dd$SCORESUM.x + dd$SCORESUM.y
          res <- dd[,c("FID","SCORESUM")]
        }
      }
      colnames(res)[2] <- c("baseline","AdrenalPancreas","Cardiovascular","CNS","ConnectiveBone",
                            "GI","Hematopoietic","Kidney","Liver","Other","SkeletalMuscle")[ct]
      if(ct == 1){
        rere <- res
      }else{
        rere <- merge(rere,res,by = "FID")
      }
    }
  }
  system(paste0("touch ",temp_p,"/OmniPRS.done.",traits,"_chr_",chr))
  cat("OmniPRS in chr",chr,"done\n")
  
  files <- list.files(path = temp_p,
                      pattern = paste0("OmniPRS",".done.",traits,"_chr_"))
  cat("A total of ",length(files), " chromosome files has created\n")
  if(length(files) == 22){
    require(BAS)
    require(glmnet)
    options(scipen=100)
    
    # system(paste0("rm ",temp_p,"/*"))
    
    bet <- NULL
    for(c in 1:22){
      bet_tem <- fread(paste0(temp_p,"/1_",ct,"_PM_Omni.txt"))[,c(1,4)]
      colnames(bet_tem) = c("RSID","baseline")
      for(ct in 2:11){
        d2 = fread(paste0(temp_p,"/",c,"_",ct,"_PM_Omni.txt"))[,c(1,4)]
        colnames(d2) = c("RSID",c("baseline","AdrenalPancreas","Cardiovascular","CNS","ConnectiveBone",
                                  "GI","Hematopoietic","Kidney","Liver","Other","SkeletalMuscle")[ct])
        bet_tem = merge(bet_tem,d2,by = "RSID")
      }
      bet = rbind(bet,bet_tem)
    }
    fwrite(bet,paste0(out_p,"/OmniPRS_beta"))
    
    phen <- fread(pheno)
    data <- merge(phen, rere, by.x = "UDI", by.y = "FID") %>% data.frame()
    data_adj <- cbind(y=data[,which(colnames(data) == phe_trait)], 
                      data[,-(2:which(colnames(data) == "baseline")-1)]) %>% data.frame
    
    # baseline + tissue-type model
    EW.re <- crossV(10, data_adj, "Additive", bina = bina)
    Lasso.re <- crossV(10, data_adj, "lasso", bina = bina)
    # Elasticnet.re <- crossV(10, data_adj, "elasticnet")
    BMA.re <- crossV(10, data_adj, "BMA", bina = bina)
    
    final_result = as.data.frame(cbind(rere,EW.re$score,Lasso.re$score,BMA.re$score))
    fwrite(final_result,paste0(out_p,"/OmniPRS_score"))
  }
  
  lasso <- function(data, bina = F, lambdas=exp(seq(log(0.001), log(0.1), length.out = 20))){
    time_start <- Sys.time()
    
    Y <- data[,1]; X <-  as.matrix(data[,-1])
    # lambdas=exp(seq(log(0.001), log(0.1), length.out = 20))
    if(bina == F){
      lasso_model <- cv.glmnet(X,Y,alpha = 1,lambda = lambdas, nfolds = 10)
      lambda.min <- lasso_model$lambda.min
      lasso_best <- glmnet(X,Y,alpha = 1,lambda = lambda.min)
    }else{
      lasso_model <- cv.glmnet(X,Y,alpha = 1,lambda = lambdas, nfolds = 10,family = "binomial")
      lambda.min <- lasso_model$lambda.min
      lasso_best <- glmnet(X,Y,alpha = 1,lambda = lambda.min,family = "binomial")
    }
    beta <- as.matrix(coef(lasso_best))
    PRS <- cbind(1, X) %*% beta
    
    exc_time <- difftime(Sys.time(),time_start,units = 'mins')
    return(list(beta = beta, PRS = PRS, time = exc_time))
  }
  
  BMA <- function(data, bina = F){
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
  
  Additive <- function(data, bina = F){
    beta = rep(1,ncol(data)-1)
    return(list(beta = beta))
  }
  
  fold_fun <- function(fold, sa, method, data, bina = F){
    index_sa <- (1+(fold - 1)*sa):min(fold*sa, nrow(data))
    train <- data[-index_sa,]
    test <- data[index_sa,]
    beta.train <- paste0(method, "(data = train[,1:12], bina = bina)$beta") %>% parse(text=.) %>% eval()
    if(method == "Additive"){
      if(bina == F){
        R2 <- cor(test$y, as.matrix(test[,2:12]) %*% 
                    as.matrix(beta.train), method = "pearson")**2
      }else{
        R2 <- auc(roc(test$y, as.matrix(test[,2:12]) %*% 
                        as.matrix(beta.train))) %>% as.numeric()
      }
    }else{
      if(bina == F){
        R2 <- cor(test$y, cbind(1,as.matrix(test[,2:12])) %*% 
                    as.matrix(beta.train), method = "pearson")**2
      }else{
        R2 <- auc(roc(test$y, cbind(1,as.matrix(test[,2:12])) %*% 
                        as.matrix(beta.train))) %>% as.numeric()
      }
    }
    return(list(R2=R2,beta.train=beta.train))
  }
  
  crossV <- function(fold = 5, data, method = "LAsso",bina =F,print_b = T){
    # time_start <- Sys.time()
    cat(paste0("start ", fold,"-fold cross validation using method ", method), '\n')
    test_data <- data[1:(0.5*nrow(data)),]
    train_data <- data[-(1:(0.5*nrow(data))),]
    sa <- round(1/fold*nrow(train_data), digits = 0)
    # res <- parSapply(clu, 1:fold, fold_fun, sa, method, data)
    resu_R <- NULL; resu_b <- NULL
    for(f in 1:fold){
      res <- fold_fun(f, sa, method, train_data, bina = bina)
      resu_R <- c(resu_R, res$R2)
      resu_b <- rbind(resu_b, res$beta.train)
      cat(paste0("fold ", f," done !"), '\n')
    }
    # exc_time <- difftime(Sys.time(),time_start,units = 'mins')
    # if(sum(resu_R == max(resu_R)) != 10){
    if(method == "Additive"){
      if(bina == F){
        f1 = as.formula(paste0("y~",paste("Sex","BaseAge",
                                          paste0("PC",1:10,collapse = "+"),sep = "+")))
        eta <- summary(lm(f1,test_data))$residuals
        score <- as.matrix(test_data[,2:12]) %*% 
          as.matrix(resu_b[resu_R == max(resu_R),])
        R2 <- cor(eta, as.matrix(test_data[,2:12]) %*% 
                    as.matrix(resu_b[resu_R == max(resu_R),]), method='pearson')**2
      }else{
        score <- as.matrix(test_data[,2:12]) %*% 
          as.matrix(resu_b[resu_R == max(resu_R),])
        test_data = cbind(test_data,x = as.matrix(test_data[,2:12]) %*% 
                            as.matrix(resu_b[resu_R == max(resu_R),]))
        r1 = rms::lrm(as.formula(paste0("y~x+",
                                        paste("Sex","BaseAge", 
                                              paste0("PC",1:10,collapse = "+"),sep = "+"))), 
                      data = test_data)
        R2 <- r1$stats[6]
      }
    }else{
      if(bina == F){
        f1 = as.formula(paste0("y~",paste("Sex","BaseAge",
                                          paste0("PC",1:10,collapse = "+"),sep = "+")))
        eta <- summary(lm(f1,test_data))$residuals
        score <- cbind(1,as.matrix(test_data[,2:12])) %*% 
          as.matrix(resu_b[resu_R == max(resu_R),])
        R2 <- cor(eta, cbind(1,as.matrix(test_data[,2:12])) %*% 
                    as.matrix(resu_b[resu_R == max(resu_R),]), method='pearson')**2
      }else{
        score <- cbind(1,as.matrix(test_data[,2:12])) %*% 
          as.matrix(resu_b[resu_R == max(resu_R),])
        test_data = cbind(test_data,x = cbind(1,as.matrix(test_data[,2:12])) %*% 
                            as.matrix(resu_b[resu_R == max(resu_R),]))
        r1 = rms::lrm(as.formula(paste0("y~x+",
                                        paste("Sex","BaseAge", 
                                              paste0("PC",1:10,collapse = "+"),sep = "+"))), 
                      data = test_data)
        R2 <- r1$stats[6]
      }
    }
    # }else{
    #   R2 <- 0
    # }
    if(print_b) return(list(score = score, R2=R2,beta=as.matrix(resu_b[resu_R == max(resu_R),]))) else return(R2)
  }
  
  sum.pro <- function(summs, trait.name = "Height", 
                      out = "/data2/projects/bioinfo/zhshao/GWAS.summary/sums/"){
    require(data.table)
    require(magrittr)
    
    out_p <- paste0(out, trait.name,"/")
    dir.check(out_p); setwd(out_p)
    ldsc_p <- paste0(out_p, "ldsc/")
    dir.check(ldsc_p)
    
    cat("Start processing trait: ",trait.name,"\n")
    data <- fread(summs) #[1] 1373020 
    hm3 <- fread("/data2/projects/bioinfo/zhshao/GWAS.summary/w_hm3.snplist") #[1] 1217311 
    cat("A total of ",nrow(data),"SNPs are loaded \n")
    # data1 <- data[, c("RSID","CHR","POS","EFFECT_ALLELE","OTHER_ALLELE","BETA","SE","P","N","EFFECT_ALLELE_FREQ")]
    # data1$Z <- data1$BETA / data1$SE
    # data1 <- data1[,c(1,4,5,2,3,9,6,7,11,8,10)]
    # colnames(data1) <- c("SNP","A1","A2","CHR","BP","N","Beta","SE","Z","P","MAF")
    # fwrite(data1,paste0(out_p, trait.name, "_CT.sumstats"), sep="\t",quote=F, row.names=F)
    # 
    data <- merge(hm3, data,by.x = "SNP",by.y = "RSID") %>% na.omit() #[1] 1101852 
    cat(nrow(data),"belong to HapMap3 SNPs \n")
    data$P <- as.numeric(data$P)
    data$P[data$P == 0] <- 1e-323
    attach(data)
    index1 <- A1 == EFFECT_ALLELE & A2 == OTHER_ALLELE
    index2 <- A2 == EFFECT_ALLELE & A1 == OTHER_ALLELE
    index3 <- N > 0.67 * quantile(N,0.90)
    index4 <- P > 0 & P <= 1
    detach(data)
    data$BETA[index2] <- (-1) * data$BETA[index2]
    try({data$EFFECT_ALLELE_FREQ[index2] <- 1 - data$EFFECT_ALLELE_FREQ[index2]})
    
    data <- data[(index1 | index2) & index3 & index4, 
                 c("SNP","CHR","POS","A1","A2","BETA","SE","P","N","EFFECT_ALLELE_FREQ")] #[1] 1072395 
    cat("A total of ",nrow(data),"SNPs passed the QC \n")
    # fwrite(data,paste0(out_p,trait.name,".txt"), sep="\t",quote=F, row.names=F)
    fwrite(as.data.frame(max(data$N)),paste0(out_p, trait.name,"-N.txt"), sep="\t",quote=F, row.names=F)
    
    data <- data[, c("SNP","CHR","POS","A1","A2","BETA","SE","P","N","EFFECT_ALLELE_FREQ")]
    # data <- data[, c("RSID","CHR","POS","EFFECT_ALLELE","OTHER_ALLELE","BETA","SE","P","N","EFFECT_ALLELE_FREQ")]
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
    
    cat("Start ldsc \n")
    ldsc <- paste0("python /home/opt/software/anaconda3/envs/ldsc/bin/ldsc.py",
                   " --ref-ld-chr /data2/projects/bioinfo/xingjie/src/ldsc/eur_w_ld_chr/",
                   " --w-ld-chr /data2/projects/bioinfo/xingjie/src/ldsc/eur_w_ld_chr/",
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
    
    re_ct = trait.name
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
      data_ct = fread(paste0(ldsc_p, trait.name,"_baselineLD_cell_type_group.", ct_index,".results"))
      re_ct = c(re_ct, 1-pnorm(data_ct$`Coefficient_z-score`[nrow(data_ct)]))
    }
    fwrite(as.data.frame(re_ct),paste0(out_p, trait.name,"-cell_type.txt"), sep="\t",quote=F, row.names=F)
    
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
        if(length(hm3snp_index) == 0) next
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
    
    cat("Summary data preparation has done!\n")
  }
  
  output.check <- function(output_path){ #/home/zhshao/out/PRS
    setwd(output_path)
    files <- list.files()
    res <- sapply(files,FUN = function(x){
      d1 <- fread(x,fill=T)
      ifelse(any(grep(paste(d1),pattern="done")), 1, 0)
    })
    file.remove(names(res[res == 1]))
    if(any(res != 1)){
      return(res[which(res != 1)])
    }else{
      return("OK !")
    }
  }
  
  dir.check <- function(path){
    if(!dir.exists(path)) dir.create(path)
  }
  
  pseudo.R2 <- function (fit, null = NULL, restrictNobs = FALSE){
    TOGGLE = (class(fit)[1] == "lm" | class(fit)[1] == 
                "gls" | class(fit)[1] == "lme" | class(fit)[1] == 
                "glm" | class(fit)[1] == "negbin" | class(fit)[1] == 
                "zeroinfl" | class(fit)[1] == "clm" | class(fit)[1] == 
                "vglm" | class(fit)[1] == "betareg" | class(fit)[1] == 
                "rq")
    BOGGLE = (class(fit)[1] == "nls" | class(fit)[1] == 
                "lmerMod" | class(fit)[1] == "glmerMod" | 
                class(fit)[1] == "merModLmerTest" | class(fit)[1] == 
                "lmerModLmerTest" | class(fit)[1] == "clmm")
    SMOGGLE = (class(fit)[1] == "lmerMod" | class(fit)[1] == 
                 "glmerMod" | class(fit)[1] == "merModLmerTest" | 
                 class(fit)[1] == "lmerModLmerTest" | class(fit)[1] == 
                 "vglm")
    ZOGGLE = (class(fit)[1] == "zeroinfl")
    ZOGGLE2 = (class(fit)[1] == "rq")
    NOGGLE = is.null(null)
    ERROR = "Note: For models fit with REML, these statistics are based on refitting with ML"
    ERROR2 = "None"
    if (!restrictNobs & NOGGLE & TOGGLE) {
      print(1)
      null = update(fit, ~1)
    }
    if (restrictNobs & NOGGLE & TOGGLE) {
      print(2)
      null = update(fit, ~1, data = fit$model)
    }
    if (restrictNobs & !NOGGLE) {
      print(3)
      null = update(null, data = fit$model)
    }
    if (NOGGLE & BOGGLE) {
      print(4)
      ERROR = "You need to supply a null model for nls, lmer, glmer, or clmm"
    }
    if ((!TOGGLE) & (!BOGGLE)) {
      print(5)
      ERROR = "This function will work with lm, gls, lme, lmer, glmer, glm, negbin, zeroinfl, nls, clm, clmm, and vglm"
    }
    
    SMOGGLE2 = (class(null)[1] == "lmerMod" | class(null)[1] == 
                  "glmerMod" | class(null)[1] == "merModLmerTest" | 
                  class(null)[1] == "lmerModLmerTest" | class(null)[1] == 
                  "vglm")
    Y = matrix(rep(NA, 2), ncol = 1)
    colnames(Y) = ""
    rownames(Y) = c("Model:", "Null:")
    Z = matrix(rep(NA, 3), ncol = 1)
    colnames(Z) = c("Pseudo.R.squared")
    rownames(Z) = c("McFadden", "Cox and Snell (ML)", 
                    "Nagelkerke (Cragg and Uhler)")
    X = matrix(rep(NA, 4), ncol = 4)
    colnames(X) = c("Df.diff", "LogLik.diff", "Chisq", 
                    "p.value")
    rownames(X) = ""
    U = matrix(rep(NA, 2), ncol = 1)
    colnames(U) = ""
    rownames(U) = c("Model:", "Null:")
    if (TOGGLE | BOGGLE) {
      if (!SMOGGLE) {
        Y[1] = toString(fit$call)
      }
      if (SMOGGLE) {
        Y[1] = toString(fit@call)
      }
    }
    if (TOGGLE | (BOGGLE & !NOGGLE)) {
      if (!SMOGGLE2) {
        Y[2] = toString(null$call)
      }
      if (SMOGGLE2) {
        Y[2] = toString(null@call)
      }
      if (!ZOGGLE & !ZOGGLE2) {
        N = nobs(fit)
        U[1, 1] = nobs(fit)
        U[2, 1] = nobs(null)
      }
      if (!ZOGGLE & ZOGGLE2) {
        N = length(fit$y)
        U[1, 1] = length(fit$y)
        U[2, 1] = length(null$y)
      }
      if (ZOGGLE) {
        N = fit$n
        U[1, 1] = fit$n
        U[2, 1] = null$n
      }
      if (U[1, 1] != U[2, 1]) {
        ERROR2 = "WARNING: Fitted and null models have different numbers of observations"
      }
      m = suppressWarnings(logLik(fit, REML = FALSE))[1]
      n = suppressWarnings(logLik(null, REML = FALSE))[1]
      mf = 1 - m/n
      Z[1, ] = signif(mf, digits = 6)
      cs = 1 - exp(-2/N * (m - n))
      Z[2, ] = signif(cs, digits = 6)
      nk = cs/(1 - exp(2/N * n))
      Z[3, ] = signif(nk, digits = 6)
    }
    W = ERROR
    WW = ERROR2
    V = list(Y, Z, X, U, W, WW)
    names(V) = c("Models", "Pseudo.R.squared.for.model.vs.null", 
                 "Likelihood.ratio.test", "Number.of.observations", 
                 "Messages", "Warnings")
    return(Z[3])
  }
  
  r2_diff <- function (dat, v1, v2, nv){
    dat = scale(dat)
    omat = cor(dat)
    if (length(v1) == 1 & length(v2) == 1) {
      ord = c(1, (1 + v1), (1 + v2))
      m1 = lm(dat[, 1] ~ dat[, (1 + v1)])
      s1 = summary(m1)
      m2 = lm(dat[, 1] ~ dat[, (1 + v2)])
      s2 = summary(m2)
      R2 = s1$r.squared
      mv2 = 1
      t100 = (1/(nv) * (1 - R2)^2)
      lamda = R2/t100
      t100 = t100^2 * 2 * (mv2 + 2 * lamda)
      var1 = t100
      R2 = s2$r.squared
      mv2 = 1
      t100 = (1/(nv) * (1 - R2)^2)
      lamda = R2/t100
      t100 = t100^2 * 2 * (mv2 + 2 * lamda)
      var2 = t100
      dvr2 = s1$r.squared - s2$r.squared
      aoa = olkin1_2(omat[ord, ord], nv)
      chi_dum = dvr2^2/aoa
      p3 = pchisq(chi_dum, 1, lower.tail = F)
      uci = dvr2 + 1.96 * aoa^0.5
      lci = dvr2 - 1.96 * aoa^0.5
      z = list(rsq1 = s1$r.squared, rsq2 = s2$r.squared, var1 = var1, 
               var2 = var2, var_diff = aoa, r2_based_p = p3, r2_based_p_one_tail = p3/2, 
               mean_diff = dvr2, upper_diff = uci, lower_diff = lci)
      return(z)
    }
  }
  
  olkin1_2 <- function (omat, nv) {
    av = array(0, 3)
    av[1] = 2 * omat[2, 1]
    av[2] = -2 * omat[3, 1]
    av[3] = 0
    ov = matrix(0, 3, 3)
    ov[1, 1] = (1 - omat[2, 1]^2)^2/nv
    ov[2, 2] = (1 - omat[3, 1]^2)^2/nv
    ov[3, 3] = (1 - omat[3, 2]^2)^2/nv
    ov[2, 1] = (0.5 * (2 * omat[3, 2] - omat[2, 1] * omat[3, 
                                                          1]) * (1 - omat[3, 2]^2 - omat[2, 1]^2 - omat[3, 1]^2) + 
                  omat[3, 2]^3)/nv
    ov[1, 2] = ov[2, 1]
    ov[3, 1] = (0.5 * (2 * omat[3, 1] - omat[2, 1] * omat[3, 
                                                          2]) * (1 - omat[3, 2]^2 - omat[2, 1]^2 - omat[3, 1]^2) + 
                  omat[3, 1]^3)/nv
    ov[1, 3] = ov[3, 1]
    ov[3, 2] = (0.5 * (2 * omat[2, 1] - omat[3, 1] * omat[3, 
                                                          2]) * (1 - omat[3, 2]^2 - omat[2, 1]^2 - omat[3, 1]^2) + 
                  omat[2, 1]^3)/nv
    ov[2, 3] = ov[3, 2]
    aova = t(av) %*% ov %*% (av)
  }
  
  
