
  GRS.OmniPRS <- function(GRS.input, traits, chr, sums_p, base_p, base_f,
                          target_p, target_f, pheno, phe_trait, out, temp, bina = F){
    t1 = Sys.time()
    cat("Start OmniPRS method \n")
    out_p <- paste0(out,traits); dir.check(out_p)
    temp_p <- paste0(temp,traits); dir.check(temp_p)
    if(file.exists(paste0(temp_p, "/", "OmniPRS","-",chr,".profile"))) return(NULL)
    
    N <- fread(paste0(sums_p,traits,"-N.txt")) %>% as.numeric
    h2 <- fread(paste0(sums_p,traits,"-h2.txt")) %>% as.numeric
    
    rere <- NULL
    setwd("/home/zhshao/project/GRS/AnnoBLUP")
    
    outCoord = paste0(temp_p,"/",traits,"_",chr)
    if(file.exists(outCoord)) file.remove(outCoord)
    
    omni <- paste0("python /home/zhshao/project/GRS/AnnoBLUP/sum_cord.py",
                   " --gf=",base_p,base_f,chr,
                   " --pf=",sums_p,traits,"_trait.txt",
                   " --FUNCT_FILE=",sums_p,traits,"_fa.txt",
                   " --coord=",outCoord,
                   " --ssf=",sums_p,traits,"_sums.txt",
                   " --N=",N,
                   # " --posterior_means=",temp_p,"/",chr,"_",ct,"_PM",
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
      snp_h2 = fread(paste0(sums_p,traits,"_fc_",ct,".txt"))
      
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
      
      prs <- paste0("/home/opt/software/plink_linux_x86_64_20210606/plink",
                    " --bfile ",target_p, target_f, chr,
                    " --score ", temp_p,"/",chr,"_",ct,"_PM_Omni.txt 1 2 4 header sum",
                    " --out ", temp_p, "/", "OmniPRS","-",chr,"_",ct,
                    " > ", temp_p, "/", "OmniPRS","-",chr,"_",ct)
      system(prs)
      
      system(paste0("touch /home/zhshao/out/temp/OmniPRS.done.",traits,"_",ct,"_chr_",chr))
      
      files <- list.files(path = "/home/zhshao/out/temp/",
                          pattern = paste0("OmniPRS",".done.",traits,"_",ct,"_chr_"))
      cat(length(files), " chromosome are completed in cell type",ct,"\n")
      if(length(files) == 22){
        cat("Begin merge 22 chromosome result\n")
        options(scipen=100)
        
        res <- NULL
        for(c in 1:22){
          d1 = paste0(temp_p, "/", "OmniPRS","-",c,"_",ct,".profile") %>% fread(.) %>%
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
        
        cat("Calculating R2\n")
        phen <- fread(pheno)
        data <- merge(phen, res, by.x = "UDI", by.y = "FID") %>% data.frame()
        
        if(bina == F){
          R2 <- cor(data[,which(colnames(data) == phe_trait)],
                    data[,which(colnames(data) == colnames(res)[2])], method='pearson')**2
        }else{
          R2 <- auc(roc(data[,which(colnames(data) == phe_trait)], 
                        data[,which(colnames(data) == colnames(res)[2])])) %>% as.numeric()
          # model_data = cbind(data[,which(colnames(data) == phe_trait)], data$SCORESUM) %>% as.data.frame()
          # colnames(model_data) = c("y","x")
          # model = glm(y ~ x,
          #             data = model_data,
          #             family = binomial(link="logit"))
          # R2 <- pscl::pR2(model, type = "mcfadden")[6]
        }
        fwrite(as.data.frame(1),paste0(out_p,"/", "OmniPRS",ct,"-", R2))
      }
    }
    system(paste0("touch /home/zhshao/out/temp/OmniPRS.done.",traits,"_chr_",chr))
    cat("OmniPRS"," in chr",chr,"done\n")
    
    files <- list.files(path = "/home/zhshao/out/temp/",
                        pattern = paste0("OmniPRS",".done.",traits,"_chr_"))
    cat("A total of ",length(files), " chromosome files has created\n")
    if(length(files) == 22){
      require(BAS)
      require(glmnet)
      options(scipen=100)
      
      # system(paste0("rm ",temp_p,"/*"))
      fwrite(rere,paste0(sums_p,"/OmniPRS"))
      phen <- fread(pheno)
      data <- merge(phen, rere, by.x = "UDI", by.y = "FID") %>% data.frame()
      data_adj <- cbind(y=data[,which(colnames(data) == phe_trait)], 
                        data[,-(2:which(colnames(data) == "baseline")-1)]) %>% data.frame
      
      # baseline + tissue-type model
      Additive.re <- crossV(10, data_adj, "Additive", bina = bina)
      Lasso.re <- crossV(10, data_adj, "lasso", bina = bina)
      # Elasticnet.re <- crossV(10, data_adj, "elasticnet")
      BMA.re <- crossV(10, data_adj, "BMA", bina = bina)
      
      fwrite(as.data.frame(1),paste0(out_p,"/", "OmniPRS",".A-", max(Additive.re)))
      fwrite(as.data.frame(1),paste0(out_p,"/", "OmniPRS",".L-", max(Lasso.re)))
      fwrite(as.data.frame(1),paste0(out_p,"/", "OmniPRS",".B-", max(BMA.re)))
    }
    
    cat("~\n")
    cat("~~\n")
    cat("~~~\n")
    cat("~~\n")
    cat("~\n")
  }
  
