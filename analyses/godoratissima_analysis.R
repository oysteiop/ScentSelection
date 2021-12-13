##################################################################################
##### Using reduced-rank-regression to measure selection on floral volatiles #####
##################################################################################

# Gymnadenia odoratissima (Gross et al. 2016 Plos ONE)

rm(list=ls())
setwd("C:/data/SelectionRRR/")

library(Hmsc)
library(abind)
library(mgcv)
library(knitr)
library(reshape2)

indat = read.table("data/godoratissima.txt", header=T)
indat$study = as.factor(paste(indat$Population, indat$Year, sep="_"))
table(indat$study)

studies = levels(indat$study)

#### Fit models and compute selection gradients ####
MFList = list()
resList = list()
betaList = list()
betaList_MR = list()

# Set sampling parameters
samples = 1000
thin = 100
transient = .5*(thin*samples)
adaptNf = 0.4*(thin*samples)
nChains = 1

a = Sys.time()

for(s in 1:length(studies)){

indatpop = indat[indat$study==studies[s],]

Y = cbind(indatpop$NrFruits)
XData = indatpop[,c(7, 8, 9)] #3 other traits
XRRRData = indatpop[,10:31] #volatiles

XFormula = as.formula(paste("~", paste0(colnames(XData), collapse=" + ")))
XRRRFormula = as.formula(paste("~", " -1 + ", paste0(colnames(XRRRData), collapse=" + ")))

XRRRDataScaled = data.frame(scale(XRRRData))

m = Hmsc(Y = as.matrix(Y/mean(Y)), XData=XData, XFormula=XFormula,
         distr="normal", 
         XRRRData=XRRRDataScaled, XRRRFormula=XRRRFormula, ncRRR=1, XRRRScale = TRUE)

XData2 = data.frame(XData, XRRRData)
XFormula2 = as.formula(paste("~", paste0(colnames(XData2), collapse=" + ")))
XFormula2

m2 = Hmsc(Y = as.matrix(Y/mean(Y)), XData=XData2, XFormula=XFormula2,
          distr="normal")

m = sampleMcmc(m, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
               transient = transient, nChains = nChains)

m2 = sampleMcmc(m2, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
               transient = transient, nChains = nChains)

filename = paste0("analyses/godoratissima/model_", studies[s], "_thin_", thin, ".RData")
save(m, file=filename)

filename2 = paste0("analyses/godoratissima/model2_", studies[s], "_thin_", thin, ".RData")
save(m2, file=filename2)

# Assess MCMC convergence and evaluate model fit
post = convertToCodaObject(m) #wRRR needs to be added here

# Posterior of wRRR
bind0 = function(...) {abind(..., along = 0)}
postList = poolMcmcChains(m$postList, chainIndex = 1:length(m$postList), 
                          start = 1, thin = 1)
valList = lapply(postList, function(a) a[["wRRR"]])
wRRR_post = do.call(bind0, valList)

predY = computePredictedValues(m)
predYm = apply(simplify2array(predY), 1:2, mean)
MF = evaluateModelFit(m, predY)

predY2 = computePredictedValues(m2)
predYm2 = apply(simplify2array(predY2), 1:2, mean)
MF2 = evaluateModelFit(m2, predY2)

# Using 5-fold cross-validation to evaluate the predictive power of the two models
partition = createPartition(m, nfolds=5)

predY_CV = computePredictedValues(m, partition)
predYm_CV = apply(simplify2array(predY_CV), 1:2, mean)
MF_CV = evaluateModelFit(m, predY_CV)

predY_CV2 = computePredictedValues(m2, partition)
predYm_CV2 = apply(simplify2array(predY_CV2), 1:2, mean)
MF_CV2 = evaluateModelFit(m2, predY_CV2)

MFmat = data.frame(matrix(c(MF, MF2, MF_CV, MF_CV2), nrow=2))
rownames(MFmat) = c("RMSE", "R2")
colnames(MFmat) = c("RRR_E", "MR_E", "RRR_P", "MR_P")

MFList[[s]] = MFmat

# Compute selection gradients
post = convertToCodaObject(m)
beta_raw_post = post$Beta[,2:4][[1]] #posterior for betas

beta_var_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=ncol(beta_raw_post))
for(i in 1:ncol(beta_raw_post)){
  beta_var_post[,i] = beta_raw_post[,i]*sd(XData[,i])
}

beta_mean_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=ncol(beta_raw_post))
for(i in 1:ncol(beta_raw_post)){
  beta_mean_post[,i] = beta_raw_post[,i]*mean(XData[,i])
}

beta_raw = apply(beta_raw_post, 2, mean)
beta_raw_lower = apply(beta_raw_post, 2, quantile, c(0.025))
beta_raw_upper = apply(beta_raw_post, 2, quantile, c(0.975))

beta_var = apply(beta_var_post, 2, mean)
beta_var_lower = apply(beta_var_post, 2, quantile, c(0.025))
beta_var_upper = apply(beta_var_post, 2, quantile, c(0.975))

beta_mean = apply(beta_mean_post, 2, mean)
beta_mean_lower = apply(beta_mean_post, 2, quantile, c(0.025))
beta_mean_upper = apply(beta_mean_post, 2, quantile, c(0.975))

resmat_RRR = data.frame(beta_raw = round(beta_raw, 3),
                        L = round(beta_raw_lower, 3),
                        U = round(beta_raw_upper, 3),
                        beta_var = round(beta_var, 3),
                        L = round(beta_var_lower, 3),
                        U = round(beta_var_upper, 3),
                        beta_mean = round(beta_mean*100, 3),
                        L = round(beta_mean_lower, 3),
                        U = round(beta_mean_upper, 3))

betaRRR = post$Beta[,5:(m$ncRRR+4)][[1]] #posterior for beta_RRR

# The composite scent trait
comp1_post = matrix(NA, nrow=length(betaRRR), ncol=m$ny)
for(i in 1:length(betaRRR)){
  comp1_post[i,]=m$XRRR %*% as.matrix((wRRR_post[i,1,]))
}

beta_comp1_post=NULL
for(i in 1:samples){
  beta_comp1_post[i]=betaRRR[i]*sd(comp1_post[i,])
}

resmat_RRR[nrow(resmat_RRR)+1,]=c("NA","NA","NA", 
                                  round(mean(beta_comp1_post), 3),
                                  round(quantile(beta_comp1_post, c(0.025)), 3),
                                  round(quantile(beta_comp1_post, c(0.975)), 3),
                                  "NA","NA","NA")
rownames(resmat_RRR)[nrow(resmat_RRR)]="Scent selection axis 1"
kable(resmat_RRR) #Note the sign of beta_RRR is arbitrary, only the magnitude (and statistical support) can be interpreted. 

resList[[s]] = resmat_RRR

# Projecting selection on scent back to the original measurements
betapost = post$Beta[,5:(m$ncRRR+4)][[1]] #posterior for beta_RRR

# Compute beta = wA for each posterior sample
beta_raw_post = matrix(NA, nrow=length(betapost), ncol=m$ncORRR)
for(i in 1:length(betapost)){
  beta_raw_post[i,] = t(as.matrix(wRRR_post[i,,])) * betapost[i] #Need matrix * for more than one RRR
}

beta_var_post = matrix(NA, nrow=length(betapost), ncol=m$ncORRR)
for(i in 1:m$ncORRR){
  beta_var_post[,i] = beta_raw_post[,i]*sd(m$XRRRData[,i])
}

#beta_mean_post = matrix(NA, nrow=length(betapost), ncol=m$ncORRR)
#for(i in 1:m$ncORRR){
#  beta_mean_post[,i] = beta_raw_post[,i]*mean(XRRRData[,i])
#}

# Alternative if XRRRData scaled before model fitting
beta_mean_post = matrix(NA, nrow=length(betapost), ncol=m$ncORRR)
for(i in 1:m$ncORRR){
  beta_mean_post[,i] = beta_raw_post[,i]*(mean(XRRRData[,i])/sd(XRRRData[,i]))
}

beta_raw = apply(beta_raw_post, 2, mean)
beta_raw_lower = apply(beta_raw_post, 2, quantile, c(0.025))
beta_raw_upper = apply(beta_raw_post, 2, quantile, c(0.975))

beta_var = apply(beta_var_post, 2, mean)
beta_var_lower = apply(beta_var_post, 2, quantile, c(0.025))
beta_var_upper = apply(beta_var_post, 2, quantile, c(0.975))

beta_mean = apply(beta_mean_post, 2, mean)
beta_mean_lower = apply(beta_mean_post, 2, quantile, c(0.025))
beta_mean_upper = apply(beta_mean_post, 2, quantile, c(0.975))
beta_mean_SE2 = apply(beta_mean_post, 2, var)

betamat_RRR = data.frame(beta_raw = round(beta_raw, 3),
                         lower = round(beta_raw_lower, 3),
                         upper = round(beta_raw_upper, 3),
                         beta_var = round(beta_var, 3),
                         beta_mean = round(beta_mean*100, 3),
                         L = round(beta_mean_lower*100, 3),
                         U = round(beta_mean_upper*100, 3),
                         SE2 = round(beta_mean_SE2*100, 3))
rownames(betamat_RRR) = colnames(m$XRRRData)

betaList[[s]] = betamat_RRR

# Selection gradients from the standard model
post = convertToCodaObject(m2)

beta_raw_post = post$Beta[,5:ncol(m2$X)][[1]] #posterior for betas

beta_var_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=ncol(beta_raw_post))
for(i in 1:ncol(beta_raw_post)){
  beta_var_post[,i] = beta_raw_post[,i]*sd(XData2[,i+3])
}

beta_mean_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=ncol(beta_raw_post))
for(i in 1:ncol(beta_raw_post)){
  beta_mean_post[,i] = beta_raw_post[,i]*mean(XData2[,i+3])
}

beta_raw = apply(beta_raw_post, 2, mean)
beta_raw_lower = apply(beta_raw_post, 2, quantile, c(0.025))
beta_raw_upper = apply(beta_raw_post, 2, quantile, c(0.975))

beta_var = apply(beta_var_post, 2, mean)
beta_var_lower = apply(beta_var_post, 2, quantile, c(0.025))
beta_var_upper = apply(beta_var_post, 2, quantile, c(0.975))

beta_mean = apply(beta_mean_post, 2, mean)
beta_mean_lower = apply(beta_mean_post, 2, quantile, c(0.025))
beta_mean_upper = apply(beta_mean_post, 2, quantile, c(0.975))

betamat_MR = data.frame(beta_raw = round(beta_raw, 3),
                       L = round(beta_raw_lower, 3),
                       U = round(beta_raw_upper, 3),
                       beta_var = round(beta_var, 3),
                       L = round(beta_var_lower, 3),
                       U = round(beta_var_upper, 3),
                       beta_mean = round(beta_mean*100, 3),
                       L = round(beta_mean_lower, 3),
                       U = round(beta_mean_upper, 3))

betaList_MR[[s]] = betamat_MR

}
Sys.time()-a

save(resList, file ="analyses/godoratissima/resList.RData")
save(MFList, file ="analyses/godoratissima/MFList.RData")
save(betaList, file ="analyses/godoratissima/betaList.RData")
save(betaList_MR, file ="analyses/godoratissima/betaList_MR.RData")

# Extracting selection gradients from standard models ####
load(file ="analyses/godoratissima/betaList.RData")
betaList_MR = list()
xcor = NULL

for(s in 1:13){
  thin = 100
  study = studies[s]
  
  indatpop = indat[indat$study==studies[s],]
  
  Y = cbind(indatpop$NrFruits, indatpop$NrFruits)
  XData = indatpop[,c(7, 8, 9)] #3 other traits
  XRRRData = indatpop[,10:31] #volatiles
  
  XRRRDataScaled = data.frame(scale(XRRRData))
  
  filename = paste0("analyses/godoratissima/model2_",studies[s],"_thin_",thin,".RData")
  load(filename)
  
  xcor[[s]] = median(cor(m2$XData))
  
  # Selection gradients from the standard model
  post = convertToCodaObject(m2)
  
  beta_raw_post = post$Beta[,5:ncol(m2$X)][[1]] #posterior for betas
  
  beta_var_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=ncol(beta_raw_post))
  for(i in 1:ncol(beta_raw_post)){
    beta_var_post[,i] = beta_raw_post[,i]*sd(m2$XData[,i+3])
  }
  
  beta_mean_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=ncol(beta_raw_post))
  for(i in 1:ncol(beta_raw_post)){
    beta_mean_post[,i] = beta_raw_post[,i]*mean(m2$XData[,i+3])
  }
  
  beta_raw = apply(beta_raw_post, 2, mean)
  beta_raw_lower = apply(beta_raw_post, 2, quantile, c(0.025))
  beta_raw_upper = apply(beta_raw_post, 2, quantile, c(0.975))
  
  beta_var = apply(beta_var_post, 2, mean)
  beta_var_lower = apply(beta_var_post, 2, quantile, c(0.025))
  beta_var_upper = apply(beta_var_post, 2, quantile, c(0.975))
  
  beta_mean = apply(beta_mean_post, 2, mean)
  beta_mean_lower = apply(beta_mean_post, 2, quantile, c(0.025))
  beta_mean_upper = apply(beta_mean_post, 2, quantile, c(0.975))
  
  betamat_MR = data.frame(beta_raw = round(beta_raw, 3),
                          L = round(beta_raw_lower, 3),
                          U = round(beta_raw_upper, 3),
                          beta_var = round(beta_var, 3),
                          L = round(beta_var_lower, 3),
                          U = round(beta_var_upper, 3),
                          beta_mean = round(beta_mean*100, 3),
                          L = round(beta_mean_lower, 3),
                          U = round(beta_mean_upper, 3))
  
  betaList_MR[[s]] = betamat_MR
  
}

# Compare selection gradients via RRR and MR ####
corvals = NULL
for(s in 1:13){
  MR_betas = as.numeric(betaList_MR[[s]]$beta_mean)
  RRR_betas = as.numeric(betaList[[s]]$beta_mean)
  corvals[s] = cor(MR_betas, RRR_betas)
  }

# Load results ####
load(file ="analyses/godoratissima/resList.RData")
load(file ="analyses/godoratissima/MFList.RData")
load(file ="analyses/godoratissima/betaList.RData")

# Compile model fit results ####
MFvals = unlist(lapply(MFList, function(x) x[2,1]))
hist(MFvals)
mean(MFvals)
range(MFvals)

MF2vals = unlist(lapply(MFList, function(x) x[2,2]))
MFCVvals = unlist(lapply(MFList, function(x) x[2,3]))
MFCV2vals = unlist(lapply(MFList, function(x) x[2,4]))

# Model fit per study
mfvals = data.frame(studies, apply(data.frame(100*MFvals, 100*MF2vals, 100*MFCVvals, 100*MFCV2vals), 2, round, 1))

mfvals$cor_MR_RRR = round(corvals, 2) # From further down

save(mfvals, file="results/mf_godo.RData")

# Summary of betas ####
summary_mat = matrix(NA, nrow = 4, ncol = length(studies))
for(i in 1:length(studies)){
  summary_mat[,i] = round(as.numeric(resList[[i]]$beta_var), 3)
}

colnames(summary_mat) = studies
rownames(summary_mat) = rownames(resList[[1]])
summary_mat = t(summary_mat)
summary_mat

upper_mat = matrix(NA, nrow = 4, ncol = length(studies))
for(i in 1:length(studies)){
  upper_mat[,i] = round(as.numeric(resList[[i]]$U.1), 3)
}

colnames(upper_mat) = studies
rownames(upper_mat) = rownames(resList[[1]])
upper_mat = t(upper_mat)

lower_mat = matrix(NA, nrow = 4, ncol = length(studies))
for(i in 1:length(studies)){
  lower_mat[,i] = round(as.numeric(resList[[i]]$L.1), 3)
}

colnames(lower_mat) = studies
rownames(lower_mat) = rownames(resList[[1]])
lower_mat = t(lower_mat)

for(i in 1:13){
  if(summary_mat[i,4]<0){
    summary_mat[i,4] = summary_mat[i,4]*-1 
    lower_mat[i,4] = lower_mat[i,4]*-1
    upper_mat[i,4] = upper_mat[i,4]*-1
      }
}

bdf = data.frame(summary_mat[,4], lower_mat[,4], upper_mat[,4])
bdf[,2:3] = t(apply(bdf[,2:3], 1, sort))
bdf

save(bdf, file="results/bdf_godo.RData")

# Figure 1 ####
load(file ="analyses/godoratissima/resList.RData")

summary_mat = matrix(NA, nrow = 4, ncol = length(studies))
for(i in 1:length(studies)){
  summary_mat[,i] = round(as.numeric(resList[[i]]$beta_var),2)
}

colnames(summary_mat) = studies
rownames(summary_mat) = rownames(resList[[1]])
summary_mat = t(summary_mat)
summary_mat = as.data.frame(summary_mat)

summary_mat$altitude = factor(c("Mountain2010", "Mountain2011","Lowland2010","Lowland2011","Lowland2010","Lowland2011",
                                "Mountain2010","Mountain2011","Lowland2010","Lowland2011","Lowland2011","Mountain2010","Mountain2011"))
head(summary_mat)
names(summary_mat) = c("Beta_PlantHeight", "Beta_InfLength", "Beta_NrFlowers", "Beta_Scent", "Altitude")
str(summary_mat)

library(fBasics)
cols = divPalette(4, "Spectral")

x11(height=5, width=5.5)
par(mfrow=c(1,1))
plot(summary_mat$Altitude, summary_mat$Beta_PlantHeight, at=1:4, xlab="", xaxt="n", las=1, xlim=c(0,20), ylim=c(-.2, 1), col=cols, ylab="")
mtext("Selection intensity", 2, line=2.5)
par(new=T)
plot(summary_mat$Altitude, summary_mat$Beta_InfLength, at=6:9, xlab="", xaxt="n", yaxt="n",las=1,xlim=c(0,20), ylim=c(-.2, 1), col=cols, ylab="")
par(new=T)
plot(summary_mat$Altitude, summary_mat$Beta_NrFlowers, at=11:14, xlab="", xaxt="n", yaxt="n",las=1,xlim=c(0,20), ylim=c(-.2, 1), col=cols, ylab="")
par(new=T)
plot(summary_mat$Altitude, abs(summary_mat$Beta_Scent), at=16:19, xlab="", xaxt="n", yaxt="n", las=1,xlim=c(0,20), ylim=c(-.2, 1), col=cols, ylab="")
axis(1, at=c(2.5, 7.5, 12.5, 17.5), c("Plant height", "Inf. length", "Flowers", "Scent"))

legend("topleft", legend=c("Lowland 2010", "Lowland 2011", "Mountains 2010", "Mountains 2011"), 
       pch=15, pt.cex = 1.5, col=cols, bty="n")

# Summary of selection gradients per compound ####
XRRRData = indat[,10:31]

beta_summary_mat = matrix(NA, nrow = ncol(XRRRData), ncol = length(studies))
for(i in 1:length(studies)){
  beta_summary_mat[,i] = round(betaList[[i]]$beta_mean, 2)
}
colnames(beta_summary_mat) = studies
rownames(beta_summary_mat) = colnames(XRRRData)

#View(beta_summary_mat)

beta_SE2_mat = matrix(NA, nrow = ncol(XRRRData), ncol = length(studies))
for(i in 1:length(studies)){
  beta_SE2_mat[,i] = round(betaList[[i]]$SE2, 2)
}
colnames(beta_SE2_mat) = studies
rownames(beta_SE2_mat) = colnames(XRRRData)

#View(beta_SE2_mat)

# Error-corrected among-study SD
rawVar = apply(beta_summary_mat, 1, var)
meanSV = apply(beta_SE2_mat, 1, mean)
sigmaC = sqrt(rawVar-meanSV)

meltmeans = melt(beta_summary_mat)
head(meltmeans)
#meltmeans$value = abs(meltmeans$value)

summary_mat = as.data.frame(summary_mat)

altvals=NULL
for(i in 1:nrow(meltmeans)){
  altvals[i] = as.character(summary_mat$Altitude[which(rownames(summary_mat)==meltmeans$Var2[i])])
}

meltmeans$altitude = altvals

x11(height=8, width=7)
par(mfrow=c(2,1), mar=c(0,4,2,2))

means = tapply(meltmeans$value, meltmeans$Var1, var)
meltmeans$Var1 = factor(meltmeans$Var1, levels=names(sort(means, decreasing=T)))
sigmaC = sigmaC[match(levels(meltmeans$Var1), names(sigmaC))]
names(sigmaC)==levels(meltmeans$Var1)

plot(1:10, 1:10, col = "white", las=1, xaxt="n", ylab="", xlab="", xlim=c(0,22), ylim=c(-20, 50))
points(meltmeans$Var1, meltmeans$value, pch=16,
       col = divPalette(4, "Spectral")[as.numeric(as.factor(meltmeans$altitude))], 
       las=1, xaxt="n", ylab="", xlab="")
mtext("Selection gradient (%)", 2, line=2.5, cex=1)

points(1:22, tapply(meltmeans$value, meltmeans$Var1, mean), pch=3)

abline(h=0)

legend("topright", pch=16, legend=levels(as.factor(meltmeans$altitude)), col=divPalette(4, "Spectral"))
axis(1, 1:22, labels=F)
#text(1:22, par("usr")[3] - 2, srt = 45, adj = 1,cex=.7,
#     labels = gsub("Z","",gsub("_ngPerL","",levels(meltmeans$Var1))), xpd = TRUE)

par(mar=c(6,4,2,2))

plot(1:22, sigmaC, pch=16, ylab="", xlab="", xaxt="n", las=1, ylim=c(-1, 15))
mtext("Error-corrected SD (%)", 2, line=2.5, cex=1)
axis(1, 1:22, labels=F)
text(1:22, par("usr")[3] - 1.3, srt = 45, adj = 1,cex=.7,
     labels = gsub("Z","",gsub("_ngPerL","",levels(meltmeans$Var1))), xpd = TRUE)

# Posterior support per model ####
post_supp = NULL
for(s in 1:13){
thin = 100
study = studies[s]

indatpop = indat[indat$study==studies[s],]

Y = cbind(indatpop$NrFruits, indatpop$NrFruits)
XData = indatpop[,c(7, 8, 9)] #3 other traits
XRRRData = indatpop[,10:31] #volatiles

XRRRDataScaled = data.frame(scale(XRRRData))

filename = paste0("analyses/godoratissima/model_",studies[s],"_thin_",thin,".RData")
load(filename)

# Posterior support for beta_scent
post_supp[s] = getPostEstimate(m, "Beta")$support[5,]
}
post_supp[which(post_supp<0.5)] = 1 - post_supp[which(post_supp<0.5)]
post_supp

