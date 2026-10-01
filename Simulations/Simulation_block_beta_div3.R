
# Author: Wan-Yu Lin, Institute of Health Data Analytics and Statistics, College of Public Health, National Taiwan University, Taipei, Taiwan
# If you use this program, please cite: Benchmarking High-Dimensional Mediation Methods in Epigenome-Wide Studies, Bioinformatics Advances 2026. 

library(MASS)
library(HIMA)
library(medScan)
library(car)

hdma_env <- new.env()

# utils.R was downloaded from https://github.com/YinanZheng/HIMA/blob/master/R/utils.R

sys.source("utils.R", envir = hdma_env)

# HDMA.R was downloaded from https://github.com/YuzhaoGao/High-dimensional-mediation-analysis-R/blob/master/HDMA.R

sys.source("HDMA.R",  envir = hdma_env)

# HDMTS.R was downloaded from https://github.com/WanYuLin/HDMTS/blob/main/HDMTS.R

source("HDMTS.R")  

N <- c(300,600)
gamma1 <- 0.5
delta <- c(0.3,0.3)
eta <- c(0.5,0.5)
p <- c(1000,5000)
n.rep <- 1000
rho <- c(0,0.25,0.5,0.75)


############################################################
## Block-structured AR(1) correlation matrix
############################################################

blockAR1.mat <- function(p,
                         block.size = 100,
                         rho = 0.5){

  Sigma <- matrix(0, p, p)

  n.block <- ceiling(p / block.size)

  for(b in 1:n.block){

    idx <- ((b-1)*block.size + 1):min(b*block.size, p)

    k <- length(idx)

    Sigma[idx, idx] <- outer(
      1:k,
      1:k,
      function(i, j) rho^abs(i-j)
    )
  }

  diag(Sigma) <- 1

  Sigma
}




for(Nloop in 1:length(N)){
for(ploop in 1:length(p)){
beta1 <- append(c(0.2,0.25,0.15,0.3,0.35,0.1,0),rep(0,p[ploop]-7))/3
alpha <- append(c(0.2,0.25,0.15,0.3,0.35,0,0.1),rep(0,p[ploop]-7))
ANS <- paste0("M",which(beta1*alpha!=0))
for(rholoop in 1:length(rho)){
for(k in 1:n.rep){



X <- rnorm(N[Nloop],0,sqrt(2))
Z1 <- rnorm(N[Nloop],0,sqrt(2))
Z2 <- rnorm(N[Nloop],0,sqrt(2))
Z <- cbind(Z1,Z2)
epsilon <- rnorm(N[Nloop],0,1)



Sigma <- blockAR1.mat(
    p = p[ploop],
    block.size = 100,
    rho = rho[rholoop]
)

e <- mvrnorm(
    N[Nloop],
    mu = rep(0, p[ploop]),
    Sigma = Sigma
)



M <- matrix(NA,N[Nloop],p[ploop])
for(i in 1:p[ploop]){
  M[,i] <- alpha[i]*X+delta[1]*Z1+delta[2]*Z2+e[,i]  
    }

Y <- matrix(gamma1*X,ncol=1)+M%*%beta1+Z%*%eta+matrix(epsilon,ncol=1)


x <- X
m <- M
colnames(m) <- paste0("M",seq(1,p[ploop],1))
y <- Y
Cov <- Z



time1 <- Sys.time()
hima_cla = hima_classic(
  X=x,
  M=m,
  Y=y,
  COV.XM = Cov,
  COV.MY = Cov,  
  Y.type = c("continuous"),
  M.type = c("gaussian"),
  penalty = c("MCP"),  
  topN = NULL,
  scale = FALSE,
  Bonfcut = 0.05,
  verbose = FALSE,
  parallel = FALSE,
  ncore = 1
)
timeHc <- as.numeric(difftime(Sys.time(), time1, units = "secs"))



time1 <- Sys.time()
hima_result = hima_dblasso(
  X=x,
  M=m,
  Y=y,
  COV = Cov,
  topN = NULL,
  scale = FALSE,
  FDRcut = 0.05,
  verbose = FALSE,
  parallel = FALSE,
  ncore = 1
)
timeH <- as.numeric(difftime(Sys.time(), time1, units = "secs"))


time1 <- Sys.time()
hdma_res <- hdma_env$hdma(X=x, Y=y, M=m, COV.XM = Cov, COV.MY = Cov, family = c("gaussian"), method = c("lasso"), topN = NULL,
		  parallel = FALSE, ncore = 1, verbose = FALSE)	
timeA <- as.numeric(difftime(Sys.time(), time1, units = "secs"))




time1 <- Sys.time()
HDMTSres <- HDMTS(X=x, M=m, Y=y, COV = Cov, FDRcut = 0.05, VIFcut = 5.0, nboot = 1000, bootCI = 0.95)
timeM <- as.numeric(difftime(Sys.time(), time1, units = "secs"))



timeALL <- c(timeHc, timeH, timeA, timeM)


hdma_res1 <- hdma_res[which(p.adjust(hdma_res$P.value, method = "BH") < 0.05), , drop = FALSE]

powerHIMAc <- sum(hima_cla[,1] %in% ANS)/length(ANS)
powerHIMA <- sum(hima_result[,1] %in% ANS)/length(ANS)
powerA <- sum(rownames(hdma_res1) %in% ANS)/length(ANS)
powermy <- sum(rownames(HDMTSres) %in% ANS)/length(ANS)

zero_discovery <- c(length(hima_cla$Index) == 0, length(hima_result$Index) == 0, nrow(hdma_res1) == 0, nrow(HDMTSres) == 0)

fdpHIMAc <- ifelse(length(hima_cla[,1])>0, (length(hima_cla[,1])-sum(hima_cla[,1] %in% ANS))/length(hima_cla[,1]), 0)
fdpHIMA <- ifelse(length(hima_result[,1])>0, (length(hima_result[,1])-sum(hima_result[,1] %in% ANS))/length(hima_result[,1]), 0)
fdpA <- ifelse(length(rownames(hdma_res1))>0, (length(rownames(hdma_res1))-sum(rownames(hdma_res1) %in% ANS))/length(rownames(hdma_res1)), 0)
fdpmy <- ifelse(length(rownames(HDMTSres))>0, (length(rownames(HDMTSres))-sum(rownames(HDMTSres) %in% ANS))/length(rownames(HDMTSres)), 0)


ME <- matrix(0,p[ploop],4)
findHc <- hima_cla$Index
findHc1 <- substring(findHc,2)
ME[as.double(findHc1),1] <- hima_cla$IDE

findH <- hima_result$Index
findH1 <- substring(findH,2)
ME[as.double(findH1),2] <- hima_result$IDE

findA <- rownames(hdma_res1)
findH1 <- substring(findA,2)
ME[as.double(findH1),3] <- hdma_res1[,4]

findmy1 <- substring(rownames(HDMTSres),2)
ME[as.double(findmy1),4] <- HDMTSres$Mediation.Effect

MAE <- c()
for(mm in 1:4){
    MAE[mm] <- sum(abs(ME[,mm]-alpha*beta1))/length(alpha*beta1)
    } 
    
    
result <- c(powerHIMAc,powerHIMA,powerA,powermy,fdpHIMAc,fdpHIMA,fdpA,fdpmy,MAE,timeALL,N[Nloop],p[ploop],rho[rholoop],zero_discovery)

  write.table(matrix(result,nrow=1), "1_block_beta_div3.csv",
              col.names = FALSE, row.names = FALSE, sep = ",", append = TRUE)
              
ME2 <- rbind(ME[1:7,], apply(ME[8:nrow(ME),],2,mean,na.rm=T))
Mediation.Effect <- append(matrix(ME2,nrow=1),c(N[Nloop],p[ploop],rho[rholoop]))

  write.table(matrix(Mediation.Effect,nrow=1), "1ME_block_beta_div3.csv",
              col.names = FALSE, row.names = FALSE, sep = ",", append = TRUE)

}
}
}
}







