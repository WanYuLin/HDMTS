
# Author: Wan-Yu Lin, Institute of Health Data Analytics and Statistics, National Taiwan University
# If you use this program, please cite: High?dimensional mediation testing with stepwise regression: Exploring epigenome-wide mediation effects 


HDMTS <- function(X, M, Y, COV = NULL, FDRcut = 0.05, VIFcut = 5, nboot = 1000, bootCI = 0.95) {
if (!requireNamespace("HDMT", quietly=TRUE)) {
  stop("Please install the 'HDMT' package before running this code")
}

if (!requireNamespace("medScan", quietly=TRUE)) {
  stop("Please install the 'medScan' package before running this code")
}

if (!requireNamespace("car", quietly=TRUE)) {
  stop("Please install the 'car' package before running this code")
}

if (is.null(COV)) {
    datacombined <- data.frame(cbind(Y, X))
    colnames(datacombined)[1] <- "Y"
    fixed_vars <- c("X")
    fixed_vars1 <- c("m", "X")
} else {
    datacombined <- data.frame(cbind(Y, X, COV))
    colnames(datacombined)[1] <- "Y"    
    cov_names <- colnames(COV)
    fixed_vars <- c("X", cov_names)
    fixed_vars1 <- c("m", "X", cov_names)
}
    fixed_part <- paste(fixed_vars, collapse = " + ")
    fixed_part1 <- paste(fixed_vars1, collapse = " + ")

    talpha <- apply(M, 2, function(m){tmp <- summary(lm(as.formula(paste("m ~", fixed_part)), data=datacombined))$coef; tmp[which(rownames(tmp)=="X"),3]})
    tbeta <- apply(M, 2, function(m){tmp <- summary(lm(as.formula(paste("Y ~", fixed_part1)), data=datacombined))$coef; tmp[which(rownames(tmp)=="m"),3]})


obj4 = medScan::medScan(z.alpha = talpha, z.beta = tbeta, method = "HDMT")

fdr <- p.adjust(obj4$pvalue, method="BH")

if(sum(fdr < FDRcut) >= 1){
  m <- M[,which(fdr < FDRcut),drop=FALSE]
  if (is.null(COV)) {
    datacombined1 <- data.frame(cbind(Y, X, m))
    colnames(datacombined1)[1] <- "Y"
  } else {
    datacombined1 <- data.frame(cbind(Y, X, COV, m))
    colnames(datacombined1)[1] <- "Y"
    vif_value <- max(car::vif(lm(as.formula(paste("Y ~", fixed_part)), data = datacombined1)))
    if(vif_value > VIFcut){
       message("The variance inflation factor (VIF) exceeds VIFcut before putting mediators.")
       return(NULL)
    }
  }
  if(sum(fdr < FDRcut) > 1){          
    FitAll <- lm(Y ~ . , data = datacombined1) # Fit reg model with all variables
    FitStart <- lm(Y ~ 1, data = datacombined1) # Fit reg model with just intercept
    
    modelF <- step(FitAll, direction = "both"  , trace = 0, scope=list(lower = as.formula(paste("Y ~", fixed_part)))) # stepwise, "both"=forward&backward
    if(max(car::vif(modelF)) > VIFcut){
       datacombined1 <- datacombined1[,append("Y",names(coef(modelF))[-1])]
       modelF <- lm(Y ~ . , data = datacombined1)  
       done <- FALSE
       while(!done){
          coef_tab <- summary(modelF)$coef
          coef_tab <- coef_tab[rownames(coef_tab) %in% colnames(m), , drop=FALSE]
          remov <- rownames(coef_tab)[which.max(coef_tab[, "Pr(>|t|)"])]
          datacombined1 <- datacombined1[,!(names(datacombined1) %in% remov)]
          modelF <- lm(Y ~ . ,data = datacombined1)
          if(ncol(datacombined1) > 2){ 
            if(max(car::vif(modelF)) <= VIFcut){
               done <- TRUE
            }
          }else{
            message("No mediators can be found. Please check whether FDRcut and VIFcut are reasonable.")
            return(NULL)
          }
       }
    }
  }else{
    modelF <- lm(Y ~ . ,data = datacombined1) 
    if(max(car::vif(modelF)) <= VIFcut){
      done <- TRUE
    }
  }  
  final_mediators <- intersect(colnames(m), all.vars(formula(modelF))[-1])
  coef_tab <- summary(modelF)$coef
  direct_effect <- coef_tab[rownames(coef_tab) == "X", "Estimate"]
  mediator_beta <- coef_tab[rownames(coef_tab) %in% colnames(m), c("Estimate", "Pr(>|t|)"), drop = FALSE]

  alpha_list <- sapply(final_mediators, function(med) {
    fit <- summary(lm(as.formula(paste(med, "~", fixed_part)), data = datacombined1))$coef
    mediator_alpha <- fit[rownames(fit) == "X", c("Estimate", "Pr(>|t|)"), drop = FALSE] })
    
  pmatrix <- cbind(t(alpha_list)[,2], mediator_beta[,2])    
  result <- cbind(t(alpha_list), mediator_beta, apply(pmatrix,1,max,na.rm=T), t(alpha_list)[,1]*mediator_beta[,1])
       
  colnames(result) <- c("alpha", "p_alpha", "beta", "p_beta", "Pmax", "Mediation.Effect")
  result <- as.data.frame(result)
  result$Mediation.Prop <- (result[ , "Mediation.Effect"]/(sum(result[ , "Mediation.Effect"])+direct_effect))
  result$p_Bonferroni <- p.adjust(result$Pmax, method="bonferroni")
  result$p_FDR  <- p.adjust(result$Pmax, method="BH")

  boot_results <- sapply(final_mediators, function(m){
    replicate(nboot, {
      idx <- sample(nrow(datacombined1), replace=TRUE)
      d <- datacombined1[idx, ]
      a <- coef(lm(as.formula(paste(m, "~", fixed_part)), data = d))[2]
      b <- coef(lm(as.formula(paste("Y ~", m, "+", fixed_part)), data = d))[2]
      a * b
    })
  })

  CI <- apply(boot_results, 2, function(x){
    c(quantile(x, probs=c((1-bootCI)/2, 1-(1-bootCI)/2)), 2*min(sum(x <= 0), sum(x >= 0))/nboot )
  })

  result$CI_lower <- CI[1, ]  
  result$CI_upper <- CI[2, ]  
  result$p_Boot <- CI[3, ]  
  result$p_Boot_FDR <- p.adjust(CI[3, ], method="BH")  

  return(result)
    
}else{
  message("No mediators can be found.")
  }
}


