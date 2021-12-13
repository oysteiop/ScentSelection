##################################################################################
##### Using reduced-rank-regression to measure selection on floral volatiles #####
##################################################################################

# Gymnadenia coposea (Chapurlat et al. 2019 New Phytologist) 

rm(list=ls())

# Load packages ####
library(Hmsc)
library(abind)
library(mgcv)
library(knitr)
library(reshape2)

# Read data ####
indat = read.table("data/gconopsea.txt", header=T)
levels(indat$treat)

# Open-pollinated treatment ####
indat = indat[indat$treat=="C",]

Y = cbind(indat$fitness)
head(Y)
hist(Y)

names(indat)[1:20]
XData = indat[,c(3, 8, 9, 6, 7)] #5 other traits
head(XData)

XRRRData = indat[,19:130] #volatiles
names(XRRRData)
XRRRData = XRRRData[,-c(55,56,111,112)]
XRRRData = XRRRData[,which(colSums(XRRRData)>0)]

# Volatiles included in original analysis
sel = c("phenylacetaldehyde_D",
        "X2_phenylethanol_D",
        "X2_phenylethanol_N",
        "X2_phenylethylacetate_N",
        "indole_D",
        "indole_N",
        "X2_aminobenzaldehyde_D",
        "elemicin_D",
        "elemicin_N",
        "methyleugenol_N",
        "benzyl_alcohol_D",
        "benzyl_alcohol_N",
        "p_cresol_D",
        "Z_5_dodecenyl_acetate_D")
  
XRRRData = XRRRData[, colnames(XRRRData) %in% sel]
XRRRData = XRRRData[, match(sel, colnames(XRRRData))]
head(XRRRData)
sort(colSums(XRRRData))

XFormula = as.formula(paste("~", paste0(colnames(XData), collapse=" + ")))
XFormula

XRRRFormula = as.formula(paste("~", " -1 + ", paste0(colnames(XRRRData), collapse=" + ")))
XRRRFormula

XRRRDataScaled = data.frame(scale(XRRRData))

dim(XRRRData)

m = Hmsc(Y = as.matrix(Y/mean(Y)), XData=XData, XFormula=XFormula,
         distr="normal", 
         XRRRData=XRRRDataScaled, XRRRFormula=XRRRFormula, ncRRR=1, XRRRScale = TRUE)

#m$a2RRR = 5

XData2 = data.frame(XData, XRRRData)
XFormula2 = as.formula(paste("~", paste0(colnames(XData2), collapse=" + ")))
XFormula2

m2 = Hmsc(Y = as.matrix(Y/mean(Y)), XData=XData2, XFormula=XFormula2,
         distr="normal")

# Set sampling parameters
samples = 1000
thin = 100
transient = .5*(thin*samples)
adaptNf = 0.4*(thin*samples)
nChains = 1

# Run MCMC samling ####
a1 = Sys.time()
m = sampleMcmc(m, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
                transient = transient, nChains = nChains)
m2 = sampleMcmc(m2, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
               transient = transient, nChains = nChains)
b1 = Sys.time()
b1-a1

save(m, file="analyses/gconopsea/RRR1mod_thin100.RData")
save(m2, file="analyses/gconopsea/RRR1mod2_thin100.RData")

# Assess MCMC convergence and evaluate model fit ####
load(file="analyses/gconopsea/RRR1mod_thin100.RData")
load(file="analyses/gconopsea/RRR1mod2_thin100.RData")

getPostEstimate(m, "Beta") # Support for betas, including beta_scent

post2 = convertToCodaObject(m2)
summary(effectiveSize(post2$Beta))

post = convertToCodaObject(m)
str(post$Beta[[1]])
plot(post$Beta[,7])
effectiveSize(post$Beta[,7]) # ESS for beta_scent

# Posterior of wRRR
bind0 = function(...) {abind(..., along = 0)}
postList = poolMcmcChains(m$postList, chainIndex = 1:length(m$postList), 
                          start = 1, thin = 1)
valList = lapply(postList, function(a) a[["wRRR"]])
wRRR_post = do.call(bind0, valList)
dim(wRRR_post)

summary(effectiveSize(wRRR_post[,1,])) #ESS for weights

x11()
par(mfrow=c(4,4))
for(i in 1:14){
plot(wRRR_post[,1,i], type="l")
}

# Model fit
predY = computePredictedValues(m)
predYm = apply(simplify2array(predY), 1:2, mean)
MF = evaluateModelFit(m, predY)
MF

predY2 = computePredictedValues(m2)
predYm2 = apply(simplify2array(predY2), 1:2, mean)
MF2 = evaluateModelFit(m2, predY2)
MF2

# Using 5-fold cross-validation to evaluate the predictive power of the two models
partition = createPartition(m, nfolds=5)

predY_CV = computePredictedValues(m, partition)
predYm_CV = apply(simplify2array(predY_CV), 1:2, mean)
MF_CV = evaluateModelFit(m, predY_CV)
save(MF_CV, file="analyses/gconopsea/MF_CV.RData")

predY_CV2 = computePredictedValues(m2, partition)
predYm_CV2 = apply(simplify2array(predY_CV2), 1:2, mean)
MF_CV2 = evaluateModelFit(m2, predY_CV2)
save(MF_CV2, file="analyses/gconopsea/MF_CV2.RData")

load(file="analyses/gconopsea/MF_CV.RData")
load(file="analyses/gconopsea/MF_CV2.RData")

MF_CV
MF_CV2

# Compute selection gradients
nX = ncol(m$XData)

post = convertToCodaObject(m)
beta_raw_post = post$Beta[,2:(nX+1)][[1]] #posterior for betas
head(beta_raw_post)

beta_var_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=nX)
for(i in 1:nX){
  beta_var_post[,i] = beta_raw_post[,i]*sd(XData[,i])
}

beta_mean_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=nX)
for(i in 1:nX){
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

betaRRR = post$Beta[,(nX+2):(m$ncRRR+(nX+1))][[1]] #posterior for beta_RRR
length(betaRRR)

# The composite scent trait
dim(m$XRRR)

comp1_post = matrix(NA, nrow=length(betaRRR), ncol=m$ny)
for(i in 1:length(betaRRR)){
  comp1_post[i,]=m$XRRR %*% as.matrix((wRRR_post[i,1,]))
}

beta_comp1_post=NULL
for(i in 1:1000){
  beta_comp1_post[i]=betaRRR[i]*sd(comp1_post[i,])
}

resmat_RRR[nrow(resmat_RRR)+1,]=c("NA","NA","NA", 
                 round(mean(beta_comp1_post), 3),
                 round(quantile(beta_comp1_post, c(0.025)), 3),
                 round(quantile(beta_comp1_post, c(0.975)), 3),
                 "NA","NA","NA")
kable(resmat_RRR) # Note the sign of beta_RRR is arbitrary, only the magnitude (and statistical support) can be interpreted. 

bdf = as.numeric(c(resmat_RRR[6, 4:6]))
bdf = bdf*(-1) #Because of  negative beta
bdf[2:3] = bdf[3:2]
bdf = t(as.data.frame(bdf))
bdf

save(bdf, file="results/bdf_gcon.RData")

# Projecting selection on scent back to the original measurements
betapost = post$Beta[,(nX+2):(m$ncRRR+(nX+1))][[1]] #posterior for beta_RRR
head(betapost)

# Compute beta = wA for each posterior sample
beta_raw_post = matrix(NA, nrow=length(betapost), ncol=m$ncORRR)
for(i in 1:length(betapost)){
  beta_raw_post[i,] = t(as.matrix(wRRR_post[i,,])) * betapost[i] #Need matrix * for more than one RRR
}

beta_var_post = matrix(NA, nrow=length(betapost), ncol=m$ncORRR)
for(i in 1:m$ncORRR){
  beta_var_post[,i] = beta_raw_post[,i]*sd(m$XRRRData[,i])
}

#Because XRRRData was scaled before model fitting
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

betamat_RRR = data.frame(beta_raw = round(beta_raw, 3),
                        lower = round(beta_raw_lower, 3),
                        upper = round(beta_raw_upper, 3),
                        beta_var = round(beta_var, 3),
                        beta_mean = round(beta_mean*100, 3),
                        L = round(beta_mean_lower*100, 3),
                        U = round(beta_mean_upper*100, 3))
rownames(betamat_RRR) = colnames(m$XRRRData)

# Selection gradients from the standard model
post = convertToCodaObject(m2)
beta_raw_post = post$Beta[,(nX+2):ncol(m2$X)][[1]] #posterior for betas

beta_var_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=ncol(beta_raw_post))
for(i in 1:ncol(beta_raw_post)){
  beta_var_post[,i] = beta_raw_post[,i]*sd(XData2[,i+nX])
}

beta_mean_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=ncol(beta_raw_post))
for(i in 1:ncol(beta_raw_post)){
  beta_mean_post[,i] = beta_raw_post[,i]*mean(XData2[,i+nX])
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

resmat_MR = data.frame(beta_raw = round(beta_raw, 3),
                        L = round(beta_raw_lower, 3),
                        U = round(beta_raw_upper, 3),
                        beta_var = round(beta_var, 3),
                        L = round(beta_var_lower, 3),
                        U = round(beta_var_upper, 3),
                        beta_mean = round(beta_mean*100, 3),
                        L = round(beta_mean_lower, 3),
                        U = round(beta_mean_upper, 3))

plot(resmat_MR$beta_mean, betamat_RRR$beta_mean, pch=16)
lines(-10:10, -10:10)
cor(resmat_MR$beta_mean, betamat_RRR$beta_mean)

# Save model-fit values
mfvals = c(studies="G_conopsea", apply(data.frame(100*MF$R2, 100*MF2$R2, 100*MF_CV$R2, 100*MF_CV2$R2), 2, round, 1))

mfvals[6] = cor(resmat_MR$beta_mean, betamat_RRR$beta_mean) # From further down

save(mfvals, file="results/mf_gcon.RData")

# Plotting the estimated fitness surface
wRRR = getPostEstimate(m, "wRRR")
wRRR_C = wRRR$mean
rrr = as.matrix(m$XRRRData) %*% t(as.matrix(wRRR$mean))

gamdat = data.frame(W=m$Y[,1], r1=rrr[,1], m$XData)
head(gamdat)

gm = gam(W ~ r1 + date_flow + pl_height + nb_tot_fl + CA + SL, family=gaussian, data=gamdat)
#gm = gam(W ~ s(r1) + s(date_flow) + s(pl_height) + s(nb_tot_fl) + s(CA) + s(SL), family=gaussian, data=gamdat)
x11(height=4, width=8)
par(mfrow=c(1,2))
par(mar=c(4,4,2,2))
res = vis.gam(gm, c("nb_tot_fl","r1"), plot.type="contour", type="response",
              xlab="Flower number",
              ylab="Scent selection axis 1",
              main="", color="bw")

points(gamdat$nb_tot_fl, gamdat$r1, col="grey", pch=16, cex=1*gamdat$W)

par(mar=c(4,6,2,2), xpd=T)
out = barplot(sort(betamat_RRR$beta_mean, dec=F), xlim=c(-.7, 2), col="royalblue",
            xlab="Mean-scaled selection gradient (%)", las=1, horiz=T)
text(y=out, par("usr")[3] + -.00, srt = 0, adj = 1,cex=.7,
     labels = gsub("X","",colnames(XRRRData)[order(betamat_RRR$beta_mean, decreasing=F)]), xpd = TRUE)


