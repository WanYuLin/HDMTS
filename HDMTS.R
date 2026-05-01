
# Author: Wan-Yu Lin, Institute of Health Data Analytics and Statistics, National Taiwan University
# If you use this program, please cite: A Systematic Evaluation of High-Dimensional Mediation Methods in Epigenome-Wide Studies: Implications for Method Selection 


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

safe_vif <- function(model){
  if(length(coef(model)) <= 2) return(0)
  max(car::vif(model))
}

if(!is.numeric(X))
  stop("X must be numeric.")

if(nrow(M) != length(Y))
  stop("Sample size mismatch.")

if (is.null(COV)) {
    datacombined <- data.frame(Y = Y, X = X)
    fixed_vars <- c("X")
    fixed_vars1 <- c("m", "X")
} else {
    datacombined <- data.frame(Y = Y, X = X, COV)
    cov_names <- colnames(COV)
    fixed_vars <- c("X", cov_names)
    fixed_vars1 <- c("m", "X", cov_names)
}
    fixed_part <- paste(fixed_vars, collapse = " + ")
    fixed_part1 <- paste(fixed_vars1, collapse = " + ")

    
    form_alpha <- as.formula(paste("m ~", fixed_part))
    form_beta  <- as.formula(paste("Y ~", fixed_part1))
    
    talpha <- apply(M, 2, function(med){
      tmpdata <- datacombined
      tmpdata$m <- med
      summary(lm(form_alpha, data=tmpdata))$coef["X",3]
    })
    talpha <- as.numeric(talpha)
    names(talpha) <- colnames(M)
    
    tbeta <- apply(M, 2, function(med){
      tmpdata <- datacombined
      tmpdata$m <- med
      summary(lm(form_beta, data=tmpdata))$coef["m",3]
    })
    tbeta <- as.numeric(tbeta)
    names(tbeta) <- colnames(M)



obj4 = medScan::medScan(z.alpha = talpha, z.beta = tbeta, method = "HDMT")

fdr <- p.adjust(obj4$pvalue, method="BH")

if(sum(fdr < FDRcut) >= 1){
  m <- M[,which(fdr < FDRcut),drop=FALSE]
  mediator_names <- colnames(m)
  ## ---- HDMTS SAFE PATCH 1: remove constant mediators ----
  var_m <- apply(m, 2, var, na.rm = TRUE)
  m <- m[, var_m > 0, drop = FALSE]
  mediator_names <- colnames(m)

  if(length(mediator_names) == 0){
    message("All selected mediators are constant.")
    return(NULL)
  }
  
  if (is.null(COV)) {
    datacombined1 <- data.frame(Y=Y, X=X, m)
    ## ---- HDMTS SAFE PATCH 2: lock model frame ----
    datacombined1 <- droplevels(datacombined1)
    rownames(datacombined1) <- NULL
    colnames(datacombined1) <- make.names(colnames(datacombined1))
    colnames(m) <- make.names(colnames(m))
    mediator_names <- colnames(m)
  } else {
    datacombined1 <- data.frame(Y = Y, X = X, COV, m, check.names = FALSE)
    ## ---- HDMTS SAFE PATCH 2: lock model frame ----
    datacombined1 <- droplevels(datacombined1)
    rownames(datacombined1) <- NULL
    colnames(datacombined1) <- make.names(colnames(datacombined1))
    colnames(m) <- make.names(colnames(m))
    mediator_names <- colnames(m)
    fit_fixed <- lm(as.formula(paste("Y ~", fixed_part)), data = datacombined1)
    if(length(coef(fit_fixed)) > 2){
      vif_value <- max(car::vif(fit_fixed))
      if(vif_value > VIFcut){
        message("The variance inflation factor (VIF) exceeds VIFcut before putting mediators.")
        return(NULL)
      }
    }    
  }
  if(sum(fdr < FDRcut) > 1){           
    full_formula  <- as.formula(paste("Y ~", fixed_part, "+", paste(mediator_names, collapse="+")))
    start_formula <- as.formula(paste("Y ~", fixed_part))
    mf <- model.frame(full_formula, data = datacombined1, na.action = na.omit)
    FitAll   <- lm(full_formula,  data = mf)
    FitStart <- lm(start_formula, data = mf)
    if(length(mediator_names) + length(fixed_vars) >= nrow(datacombined1)-1){
      message("Too many predictors relative to sample size.")
      return(NULL)
    }
    oldopt <- options(na.action="na.omit")
    on.exit(options(oldopt), add=TRUE)
    modelF <- step(FitAll, direction = "backward", trace = 0)
    ## ---- HDMTS SAFE PATCH 3: remove aliased coefficients ----
    aliased <- is.na(coef(modelF))
    if(any(aliased)){
      keep_vars <- names(coef(modelF))[!aliased]
      keep_vars <- keep_vars[keep_vars != "(Intercept)"]
      datacombined1 <- datacombined1[, c("Y", keep_vars), drop=FALSE]
      modelF <- lm(Y ~ ., data = datacombined1)
    }         
    if(safe_vif(modelF) > VIFcut){
       keep_vars <- unique(c("Y","X", names(coef(modelF))[-1]))
       datacombined1 <- datacombined1[, keep_vars, drop = FALSE]
       modelF <- lm(Y ~ . , data = datacombined1)  
       done <- FALSE
       iter <- 0
       current_meds <- intersect(colnames(datacombined1), colnames(m))
       max_iter <- length(current_meds)
       while(!done && iter < max_iter){
         iter <- iter + 1
          coef_tab <- summary(modelF)$coef
          coef_tab <- coef_tab[rownames(coef_tab) %in% colnames(m), , drop=FALSE]
          pvals <- coef_tab[, "Pr(>|t|)"]
          pvals[is.na(pvals)] <- Inf
          remov <- rownames(coef_tab)[which.max(pvals)]
          mediator_names <- colnames(m)
          datacombined1 <- datacombined1[ ,!(names(datacombined1) %in% remov & names(datacombined1) %in% mediator_names)]
          modelF <- lm(Y ~ . ,data = datacombined1)
          if(ncol(datacombined1) > 2){ 
            if(safe_vif(modelF) <= VIFcut){
               done <- TRUE
            }
          }else{
            message("No mediators can be found. Please check whether FDRcut and VIFcut are reasonable.")
            return(NULL)
          }
       }
       if(iter == max_iter){
         message("VIF pruning reached iteration limit.")
       }
    }
  }else{
    modelF <- lm(Y ~ . ,data = datacombined1) 
    if(max(car::vif(modelF)) <= VIFcut){
      done <- TRUE
    }
  }  
  final_mediators <- intersect(colnames(m), all.vars(formula(modelF))[-1]) 
  if(length(final_mediators)==0){
    message("No mediators remain after stepwise/VIF pruning.")
    return(NULL)
  }
  coef_tab <- summary(modelF)$coef
  direct_effect <- coef_tab["X","Estimate"]
  if(is.na(direct_effect)) {direct_effect <- 0}
  mediator_beta <- coef_tab[rownames(coef_tab) %in% final_mediators, c("Estimate","Pr(>|t|)"), drop = FALSE]

  alpha_list <- sapply(final_mediators, function(med) {
    fit <- summary(lm(as.formula(paste(med, "~", fixed_part)), data = datacombined1))$coef
    mediator_alpha <- fit[rownames(fit) == "X", c("Estimate", "Pr(>|t|)"), drop = FALSE] })
    
  if(is.null(dim(alpha_list))){
    alpha_list <- matrix(alpha_list, nrow=2)
    rownames(alpha_list) <- c("Estimate","Pr(>|t|)")
  }
    
  pmatrix <- cbind(t(alpha_list)[,2], mediator_beta[,2])    
  result <- cbind(t(alpha_list), mediator_beta, apply(pmatrix,1,max,na.rm=T), t(alpha_list)[,1]*mediator_beta[,1])
       
  colnames(result) <- c("alpha", "p_alpha", "beta", "p_beta", "Pmax", "Mediation.Effect")
  result <- as.data.frame(result)
  Mediation.Prop.Denom <- (sum(result[ , "Mediation.Effect"])+direct_effect)
  if(abs(Mediation.Prop.Denom) < 1e-8) {Mediation.Prop.Denom <- NA}
  result$Mediation.Prop <- (result[ , "Mediation.Effect"])/Mediation.Prop.Denom
  result$p_Bonferroni <- p.adjust(result$Pmax, method="bonferroni")
  result$p_FDR  <- p.adjust(result$Pmax, method="BH")

  boot_results <- sapply(final_mediators, function(med){
    replicate(nboot, {
      idx <- sample(nrow(datacombined1), replace=TRUE)
      d <- datacombined1[idx, ]
      v <- var(d[[med]], na.rm = TRUE)
      if (is.na(v) || v == 0) return(NA)
      current_fixed <- intersect(fixed_vars, colnames(d))
      if(length(current_fixed) == 0) return(NA)
      form_a <- as.formula(paste(med, "~", paste(current_fixed, collapse = "+")))
      fit_a <- try(lm(form_a, data = d), silent = TRUE)
      if(inherits(fit_a,"try-error")) return(NA)
      coef_a <- coef(fit_a)
      if(!("X" %in% names(coef_a))) return(NA)
      a <- coef_a["X"]
      form_b <- as.formula(paste("Y ~", med, "+", paste(current_fixed, collapse = "+")))
      fit_b <- try(lm(form_b, data = d), silent=TRUE)
      if(inherits(fit_b,"try-error")) return(NA)
      coef_b <- coef(fit_b)
      if(!(med %in% names(coef_b))) return(NA)
      b <- coef_b[med]
      a * b
    })
  })
  
  if(is.null(dim(boot_results))){
    boot_results <- matrix(boot_results, ncol=1)
  }

  CI <- apply(boot_results, 2, function(x){
    x <- x[!is.na(x)]
    if(length(x) < 10) return(c(NA,NA,NA))
    ## percentile CI
    ci <- quantile(x, probs = c((1-bootCI)/2, 1-(1-bootCI)/2), type = 6)
    ## two-sided bootstrap p-value (correct)
    pval <- 2 * min(mean(x <= 0), mean(x >= 0))
    c(ci, pval)
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


