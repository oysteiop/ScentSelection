#######################################################################################
##### Using HMSC reduced-rank-regression to measure selection on floral volatiles #####
#######################################################################################

# Chapurlat et al. 2019 New Phytologist

rm(list=ls())
#setwd("C:/data/reduced_rank_selection/")

library(Hmsc)
library(abind)
library(mgcv)
library(gsg)
library(knitr)

indat = read.table("data/chapurlat.txt", header=T)

levels(indat$treat)
indat = indat[indat$treat=="C",]

Y = cbind(indat$fitness)
head(Y)
hist(Y)

names(indat)[1:20]
XData=indat[,c(3, 8, 9, 6, 7)] #5 other traits
head(XData)

XRRRData=indat[,19:130] #volatiles
names(XRRRData)

XRRRData=XRRRData[,-c(55,56,111,112)]

XRRRData=XRRRData[,which(colSums(XRRRData)>0)]


sel=c("phenylacetaldehyde_D",
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

studyDesign = data.frame(Plant=as.factor(indat$plant_ID))
rL = HmscRandomLevel(units=unique(studyDesign$Plant))
rL

XFormula = as.formula(paste("~", paste0(colnames(XData), collapse=" + ")))
XFormula

XRRRFormula = as.formula(paste("~", " -1 + ", paste0(colnames(XRRRData), collapse=" + ")))
XRRRFormula

XRRRDataScaled=data.frame(scale(XRRRData))
dim(XRRRData)

m = Hmsc(Y = as.matrix(Y/mean(Y)), XData=XData, XFormula=XFormula,
         distr="normal", 
         XRRRData=XRRRDataScaled, XRRRFormula=XRRRFormula, ncRRR=1, XRRRScale = TRUE,
         studyDesign=studyDesign, ranLevels=list(Plant = rL))

#Set sampling parameters
samples = 1000
thin = 10
transient = .5*(thin*samples)
adaptNf = 0.4*(thin*samples)
nChains = 1

a1 = Sys.time()
m = sampleMcmc(m, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
                transient = transient, nChains = nChains)
b1 = Sys.time()
b1-a1

save(m, file="analyses/chapurlat/RRR1mod_thin10.RData")

#### - Assess MCMC convergence and evaluate model fit - ####
load(file="analyses/chapurlat/RRR1mod_thin10.RData")

post = convertToCodaObject(m) #wRRR needs to be added here
str(post$Beta[[1]])
#plot(post$Beta[,7])

#Temp. posterior of wRRR
bind0 = function(...) {abind(..., along = 0)}
postList = poolMcmcChains(m$postList, chainIndex = 1:length(m$postList), 
                          start = 1, thin = 1)
valList = lapply(postList, function(a) a[["wRRR"]])
wRRR_post = do.call(bind0, valList)

dim(wRRR_post)

x11()
par(mfrow=c(4,4))
for(i in 1:14){
plot(wRRR_post[,1,i], type="l")
}

predY = computePredictedValues(m)
predYm = apply(simplify2array(predY), 1:2, mean)
MF = evaluateModelFit(m, predY)
MF

partition = createPartition(m, nfolds=2)
#predY_CV = computePredictedValues(m, partition)
#predYm_CV = apply(simplify2array(predY_CV), 1:2, mean)
MF_CV = evaluateModelFit(m, predY_CV)
MF_CV

par(mfrow=c(1,2))
plot(m$Y[,1], predYm[,1])
lines(-10:10, -10:10)

plot(m$Y[,1], predYm_CV[,1])
lines(-10:10, -10:10)

# Compute selection gradients
post = convertToCodaObject(m)
beta_raw_post = post$Beta[,2:6][[1]] #posterior for betas
#head(beta_raw_post)

beta_var_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=5)
for(i in 1:5){
  beta_var_post[,i] = beta_raw_post[,i]*sd(XData[,i])
}

beta_mean_post = matrix(NA, nrow=nrow(beta_raw_post), ncol=5)
for(i in 1:5){
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

betaRRR = post$Beta[,7:(m$ncRRR+6)][[1]] #posterior for beta_RRR
length(betaRRR)

#The first composite scent trait
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
kable(resmat_RRR) #Note the sign of beta_RRR is arbitrary, only the magnitude (and statistical support) can be interpreted. 

# Projecting selection on scent back to the original measurements
betapost = post$Beta[,7:(m$ncRRR+6)][[1]] #posterior for beta_RRR
head(betapost)

#Compute beta = wA for each posterior sample
beta_raw_post = matrix(NA, nrow=length(betapost), ncol=m$ncORRR)
for(i in 1:length(betapost)){
  beta_raw_post[i,] = t(as.matrix(wRRR_post[i,,])) * betapost[i] #Need matrix * for more than one RRR
}

beta_var_post = matrix(NA, nrow=length(betapost), ncol=m$ncORRR)
for(i in 1:m$ncORRR){
  beta_var_post[,i] = beta_raw_post[,i]*sd(m$XRRRData[,i])
}

beta_mean_post = matrix(NA, nrow=length(betapost), ncol=m$ncORRR)
for(i in 1:m$ncORRR){
  beta_mean_post[,i] = beta_raw_post[,i]*mean(XRRRData[,i])
}

#Alternative if XRRRData scaled before model fitting
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
rownames(betamat_RRR)=colnames(m$XRRRData)

#Plotting the estimated fitness surface
wRRR = getPostEstimate(m, "wRRR")
rrr = as.matrix(m$XRRRData) %*% t(as.matrix(wRRR$mean))

gamdat = data.frame(W=m$Y[,1], r1=rrr[,1], m$XData)
head(gamdat)

gm = gam(W ~ r1 + date_flow + pl_height + nb_tot_fl + CA + SL, family=gaussian, data=gamdat)
#gm = gam(W ~ s(r1) + s(date_flow) + s(pl_height) + s(nb_tot_fl) + s(CA) + s(SL), family=gaussian, data=gamdat)
x11(height=6.8, width=6.8)
par(mfrow=c(2,2))
par(mar=c(4,4,2,2))
res = vis.gam(gm, c("nb_tot_fl","r1"), plot.type="contour", type="response",
              xlab="Flower number",
              ylab="Scent selection axis 1",
              main="", color="bw")

points(gamdat$nb_tot_fl, gamdat$r1, col="grey", pch=16, cex=1*gamdat$W)

par(mar=c(4,6,2,2), xpd=T)
out = barplot(sort(betamat_RRR$beta_mean, dec=F), xlim=c(-2,5), col="royalblue",
            xlab="Mean-scaled selection gradient (%)", las=1, horiz=T)
text(y=out, par("usr")[3] + -1, srt = 0, adj = 1,cex=.7,
     labels = gsub("X","",colnames(XRRRData)[order(betamat_RRR$beta_mean, decreasing=F)]), xpd = TRUE)


#arrows(out, betamat_RRR$L, out, betamat_RRR$U, angle=90, length=0.05, code=3)

lm_res=summary(lm(W~scale(date_flow)+scale(pl_height)+scale(nb_tot_fl)+scale(CA)+scale(SL)+scale(r1), data=gamdat))
lm_res

#### Multiple-regression ####
regdat = data.frame(XData, XRRRData)

formula=as.formula(paste("Y[,1]/mean(Y[,1]) ~ ", paste0(colnames(regdat), collapse=" + ")))
formula

linmod = lm(formula=formula, data=regdat)

beta_raw = summary(linmod)$coef[-1,1]
se = summary(linmod)$coef[-1,2]
beta_var = beta_raw * apply(regdat, 2, sd)
beta_mean = beta_raw * colMeans(regdat)
resmat_mv = data.frame(beta_raw = round(beta_raw, 2),
                       se = round(se, 2),
                       beta_var = round(beta_var, 3),
                       beta_mean = round(beta_mean*100, 3))

rownames(resmat_mv) = colnames(regdat)

par(mfrow=c(1,1), xpd=F)
plot(resmat_mv[-c(1:5),4], betamat_RRR[,5], las=1, pch=16,
     xlab="Selection gradient via multiple-regression (%)",
     ylab="Selection gradient via reduced-rank regression (%)")
abline(h=0, lty=2)
abline(v=0, lty=2)
lines(-100:100, -100:100)

#### HP treatment ####
indat = read.table("data/chapurlat.txt", header=T)
indat = indat[indat$treat=="HP",]

Y = cbind(indat$fitness)

XData = indat[,c(3, 8, 9, 6, 7)] #5 other traits

XRRRData = indat[,19:130] #volatiles
XRRRData = XRRRData[,-c(55,56,111,112)]

sel=c("phenylacetaldehyde_D",
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

studyDesign = data.frame(Plant=as.factor(indat$plant_ID))
rL = HmscRandomLevel(units=unique(studyDesign$Plant))
rL

XFormula = as.formula(paste("~", paste0(colnames(XData), collapse=" + ")))
XFormula

XRRRFormula = as.formula(paste("~", " -1 + ", paste0(colnames(XRRRData), collapse=" + ")))
XRRRFormula

XRRRDataScaled=data.frame(scale(XRRRData))
dim(XRRRData)

#Some probable typos, need to pass >1 Y variable
m = Hmsc(Y = as.matrix(Y/mean(Y)), XData=XData, XFormula=XFormula,
         distr="normal", 
         XRRRData=XRRRDataScaled, XRRRFormula=XRRRFormula, ncRRR=1, XRRRScale = TRUE,
         studyDesign=studyDesign, ranLevels=list(Plant = rL))

a1 = Sys.time()
m = sampleMcmc(m, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
               transient = transient, nChains = nChains)
b1 = Sys.time()
b1-a1


save(m, file="RRR1mod_HP_thin10.RData")

#### - Assess MCMC convergence and evaluate model fit - ####
load(file="analyses/chapurlat/RRR1mod_HP_thin10.RData")

post_HP = convertToCodaObject(m) #wRRR needs to be added here
str(post_HP$Beta)
str(post_HP$Beta[[1]])
plot(post_HP$Beta[,7])

#Temp. posterior of wRRR
bind0 = function(...) {abind(..., along = 0)}
postList = poolMcmcChains(m$postList, chainIndex = 1:length(m$postList), 
                          start = 1, thin = 1)
valList = lapply(postList, function(a) a[["wRRR"]])
wRRR_post_HP = do.call(bind0, valList)

dim(wRRR_post_HP)

x11()
par(mfrow=c(4,4))
for(i in 1:14){
  plot(wRRR_post_HP[,1,i], type="l")
}

predY = computePredictedValues(m)
predYm = apply(simplify2array(predY), 1:2, mean)
MF = evaluateModelFit(m, predY)
MF

partition = createPartition(m, nfolds=2)
#predY_CV = computePredictedValues(m, partition)
predYm_CV = apply(simplify2array(predY_CV), 1:2, mean)
MF_CV = evaluateModelFit(m, predY_CV)
save(MF_CV, file="MF_CV.RData")
MF_CV

#Selection gradients
beta_raw_post_HP = post_HP$Beta[,2:6][[1]] #posterior for betas
head(beta_raw_post_HP)

beta_var_post_HP = matrix(NA, nrow=nrow(beta_raw_post_HP), ncol=5)
for(i in 1:5){
  beta_var_post_HP[,i] = beta_raw_post_HP[,i]*sd(XData[,i])
}

beta_mean_post_HP = matrix(NA, nrow=nrow(beta_raw_post_HP), ncol=5)
for(i in 1:5){
  beta_mean_post_HP[,i] = beta_raw_post_HP[,i]*mean(XData[,i])
}

beta_raw_HP = apply(beta_raw_post_HP, 2, mean)
beta_raw_lower_HP = apply(beta_raw_post_HP, 2, quantile, c(0.025))
beta_raw_upper_HP = apply(beta_raw_post_HP, 2, quantile, c(0.975))

beta_var_HP = apply(beta_var_post_HP, 2, mean)
beta_var_lower_HP = apply(beta_var_post_HP, 2, quantile, c(0.025))
beta_var_upper_HP = apply(beta_var_post_HP, 2, quantile, c(0.975))

beta_mean_HP = apply(beta_mean_post_HP, 2, mean)
beta_mean_lower_HP = apply(beta_mean_post_HP, 2, quantile, c(0.025))
beta_mean_upper_HP = apply(beta_mean_post_HP, 2, quantile, c(0.975))

resmat_RRR_HP = data.frame(beta_raw = round(beta_raw_HP, 3),
                        L = round(beta_raw_lower_HP, 3),
                        U = round(beta_raw_upper_HP, 3),
                        beta_var = round(beta_var_HP, 3),
                        L = round(beta_var_lower_HP, 3),
                        U = round(beta_var_upper_HP, 3),
                        beta_mean = round(beta_mean_HP*100, 3),
                        L = round(beta_mean_lower_HP*100, 3),
                        U = round(beta_mean_upper_HP*100, 3))

betaRRR_HP = post_HP$Beta[,7:(m$ncRRR+6)][[1]] #posterior for beta_RRR
length(betaRRR_HP)

#The first composite scent trait
dim(m$XRRR)

comp1_post_HP = matrix(NA, nrow=length(betaRRR_HP), ncol=m$ny)
for(i in 1:length(betaRRR_HP)){
  comp1_post_HP[i,]=m$XRRR %*% as.matrix((wRRR_post_HP[i,1,]))
}

beta_comp1_post_HP=NULL
for(i in 1:1000){
  beta_comp1_post_HP[i]=betaRRR_HP[i]*sd(comp1_post_HP[i,])
}

resmat_RRR_HP[nrow(resmat_RRR_HP)+1,]=c("NA","NA","NA", 
                 round(mean(beta_comp1_post_HP), 3),
                 round(quantile(beta_comp1_post_HP, c(0.025)), 3),
                 round(quantile(beta_comp1_post_HP, c(0.975)), 3),
                 "NA", "NA", "NA")
resmat_RRR_HP
resmat_RRR

# Projecting selection on scent back to the original measurements
betapost_HP = post_HP$Beta[,7:(m$ncRRR+6)][[1]] #posterior for beta_RRR

#Compute beta = wA for each posterior sample
beta_raw_post_HP = matrix(NA, nrow=length(betapost_HP), ncol=m$ncORRR)
for(i in 1:length(betapost_HP)){
  beta_raw_post_HP[i,] = t(as.matrix(wRRR_post_HP[i,,])) * betapost_HP[i] #Need matrix * for more than one RRR
}

beta_var_post_HP = matrix(NA, nrow=length(betapost_HP), ncol=m$ncORRR)
for(i in 1:m$ncORRR){
  beta_var_post_HP[,i] = beta_raw_post_HP[,i]*sd(m$XRRRData[,i])
}

beta_mean_post_HP = matrix(NA, nrow=length(betapost_HP), ncol=m$ncORRR)
for(i in 1:m$ncORRR){
  beta_mean_post_HP[,i] = beta_raw_post_HP[,i]*mean(XRRRData[,i])
}

#Alternative if XRRRData scaled before model fitting
beta_mean_post_HP = matrix(NA, nrow=length(betapost_HP), ncol=m$ncORRR)
for(i in 1:m$ncORRR){
  beta_mean_post_HP[,i] = beta_raw_post_HP[,i]*(mean(XRRRData[,i])/sd(XRRRData[,i]))
}

beta_raw_HP = apply(beta_raw_post_HP, 2, mean)
beta_raw_lower_HP = apply(beta_raw_post_HP, 2, quantile, c(0.025))
beta_raw_upper_HP = apply(beta_raw_post_HP, 2, quantile, c(0.975))

beta_var_HP = apply(beta_var_post_HP, 2, mean)
beta_var_lower_HP = apply(beta_var_post_HP, 2, quantile, c(0.025))
beta_var_upper_HP = apply(beta_var_post_HP, 2, quantile, c(0.975))

beta_mean_HP = apply(beta_mean_post_HP, 2, mean)
beta_mean_lower_HP = apply(beta_mean_post_HP, 2, quantile, c(0.025))
beta_mean_upper_HP = apply(beta_mean_post_HP, 2, quantile, c(0.975))

betamat_RRR_HP = data.frame(beta_raw = round(beta_raw_HP, 3),
                         lower = round(beta_raw_lower_HP, 3),
                         upper = round(beta_raw_upper_HP, 3),
                         beta_var = round(beta_var_HP, 3),
                         beta_mean = round(beta_mean_HP*100, 3),
                         L = round(beta_mean_lower_HP*100, 3),
                         U = round(beta_mean_upper_HP*100, 3))
rownames(betamat_RRR_HP)=colnames(m$XRRRData)

#Plotting the estimated fitness surface
wRRR_HP = getPostEstimate(m, "wRRR")
rrr_HP = as.matrix(m$XRRRData) %*% t(as.matrix(wRRR_HP$mean))

gamdat_HP = data.frame(W=m$Y[,1], r1=rrr_HP[,1], m$XData)
head(gamdat_HP)

gm_HP = gam(W ~ r1 + date_flow + pl_height + nb_tot_fl + CA + SL, family=gaussian, data=gamdat_HP)
#gm = gam(W ~ s(r1) + s(date_flow) + s(pl_height) + s(nb_tot_fl) + s(CA) + s(SL), family=gaussian, data=gamdat)

#par(mfrow=c(1,2))
par(mar=c(4,4,2,2))
res = vis.gam(gm_HP, c("nb_tot_fl","r1"), plot.type="contour", type="response",
              xlab="Flower number",
              ylab="Scent selection axis 1",
              main="", color="bw")

points(gamdat_HP$nb_tot_fl, gamdat_HP$r1, col="grey", pch=16, cex=1*gamdat_HP$W)

par(mar=c(4,6,2,2), xpd=T)
out = barplot(sort(betamat_RRR_HP$beta_mean, dec=F), xlim=c(-2,5), col="royalblue",
              xlab="Mean-scaled selection gradient (%)", las=1, horiz=T)
text(y=out, par("usr")[3] + -1, srt = 0, adj = 1,cex=.7,
     labels = gsub("X","",colnames(XRRRData)[order(betamat_RRR_HP$beta_mean, decreasing=F)]), xpd = TRUE)


lm_res_HP=summary(lm(W~scale(date_flow)+scale(pl_height)+scale(nb_tot_fl)+scale(CA)+scale(SL)+scale(r1), data=gamdat_HP))

#### Compare control and HP ####

lm_res$coef
lm_res_HP

resmat_RRR
resmat_RRR_HP

plot(lm_res$coef[-1,1], resmat_RRR$beta_var)
lines(-100:100, -100:100)

plot(lm_res_HP$coef[-1,1], resmat_RRR_HP$beta_var)
lines(-100:100, -100:100)

betamat_RRR_HP
betamat_RRR

par(mfrow=c(1,2))
plot(rrr_HP, gamdat_HP$W)
plot(rrr, gamdat$W)

par(mfrow=c(1,1))
plot(betamat_RRR$beta_mean, betamat_RRR_HP$beta_mean,
     las=1, pch=16,
     xlab="Selection gradient control treatment (%)",
     ylab="Selection gradient HP treatment (%)")
cor(betamat_RRR$beta_mean, betamat_RRR_HP$beta_mean)^2
lines(-100:100, -100:100)
abline(h=0, lty=2)
abline(v=0, lty=2)


beta_poll_post=beta_mean_post - beta_mean_post_HP

beta_poll = apply(beta_poll_post, 2, mean)*100
beta_poll_lower = apply(beta_poll_post, 2, quantile, 0.025)*100
beta_poll_upper = apply(beta_poll_post, 2, quantile, 0.975)*100


beta_poll_mat = data.frame(betamat_RRR[,5:7], betamat_RRR_HP[,5:7], 
                           round(beta_poll,3), round(beta_poll_lower, 3), round(beta_poll_upper, 3))
colnames(beta_poll_mat)=c("Control", "L", "U", "HP", "L", "U","Delta", "L", "U")
kable(beta_poll_mat)

plot(1:14, beta_poll, ylim=c(-5,10))
segments(1:14, beta_poll_lower, 1:14, beta_poll_upper)
abline(h=0)

library(reshape2)
vals=melt(rbind(beta_poll_mat$Control, beta_poll_mat$HP, beta_poll_mat$Delta))$value
lowervals=melt(rbind(beta_poll_mat[,2], beta_poll_mat[,5], beta_poll_mat[,8]))$value
uppervals=melt(rbind(beta_poll_mat[,3], beta_poll_mat[,6], beta_poll_mat[,9]))$value

cols=topo.colors(3)
out=barplot(vals, col=rep(cols, 14), 
            space=c(.1, rep(c(.1,.1,1),13), .1, .1), 
            ylim=c(-3,5),
            las=1)
legend("top", col=cols, pch=15, pt.cex=1.5,bty="n", 
       legend=c(expression(paste(beta[net])),
                expression(paste(beta[HP])),
                expression(paste(Delta,beta[pollinators]))))
mtext("Selection gradient (%)", 2, line=2.5, cex=1)

text(x=out[seq(2,41,3)], par("usr")[3]+1, srt = 45, adj = 1,cex=.7,
     labels = gsub("X","",colnames(XRRRData)), xpd = TRUE)

segments(out,lowervals, out, uppervals, col="lightgrey")


14*3
plot(seq(1, 40,3), betamat_RRR$beta_mean, pch=16, xlim=c(1,42))
points(seq(2, 41,3), betamat_RRR_HP$beta_mean, pch=1, xlim=c(1,16))
points(seq(3, 42,3), beta_poll, pch=16,col="red")
abline(v=seq(3.5, 39.5, 3), lty=2)
abline(h=0)


wRRR_unit=wRRR$mean/sqrt(sum(wRRR$mean^2))
wRRR_HP_unit=(wRRR_HP$mean/sqrt(sum(wRRR_HP$mean^2))) #Minus because slope in opposite direction 

cbind(c(wRRR_unit), c(wRRR_HP_unit))
plot(wRRR_unit, wRRR_HP_unit)

acos(wRRR_unit %*% t(wRRR_HP_unit))*(180/pi)



betapollpost=beta_comp1_post-beta_comp1_post_HP
hist(betapollpost)
mean(betapollpost)
quantile(betapollpost, c(0.025, .5, .975))


#plot(XRRRData$phenylacetaldehyde_N, m$Y[,1])
#summary(lm(m$Y[,1]~XRRRData$phenylacetaldehyde_N))$coef[2,1]*mean(XRRRData$phenylacetaldehyde_N)*100

#### Random directions ####
library(evolvability)
randombetas=randomBeta(1000, 14)

#XRRRscaled=apply(XRRRData, 2, scale)
XRRRscaled=apply(XRRRData, 2, function(x) x/mean(x))

out=NULL
for(b in 1:1000){
  z = as.matrix(XRRRData)%*%randombetas[,b]
  out[b]=summary(lm(Y[,1]/mean(Y)~scale(date_flow) + scale(pl_height) + scale(nb_tot_fl)
                    + scale(CA) + scale(SL) + scale(z), data = XData))$coef[7,1]
}

mean(abs(out))
max(abs(out))
hist(abs(out))

z_max = as.matrix(XRRRData)%*%randombetas[,which.max(out)]
plot(z_max, apply(comp1_post, 2, mean))

#### PC regression ####
pca = princomp(XRRRData, cor=FALSE)
str(pca)
summary(pca)

loadings = pca$loadings[,1:14]

pc1=pca$scores[,1]
pc2=pca$scores[,2]
pc3=pca$scores[,3]
pc4=pca$scores[,4]
pc5=pca$scores[,5]
pc6=pca$scores[,6]
pc7=pca$scores[,7]
pc8=pca$scores[,8]
pc9=pca$scores[,9]
pc10=pca$scores[,10]
pc11=pca$scores[,11]
pc12=pca$scores[,12]
pc13=pca$scores[,13]
pc14=pca$scores[,14]


linmod = lm(Y[,1]/mean(Y[,1]) ~ date_flow + pl_height + nb_tot_fl + CA + SL + 
              pc1+pc2+pc3+pc4+pc5+pc6+pc7+pc8+pc9+pc10+pc11+pc12+pc13+pc14, data = XData)

betaPC = as.matrix(summary(linmod)$coef[-c(1:6), 1])

betaPC
dim(betaPC)
dim(loadings[,1:8])

beta_raw = c(loadings[,1:14] %*% betaPC)
beta_var = beta_raw * apply(XRRRData, 2, sd)
beta_mean = beta_raw * colMeans(XRRRData)

resmat_PC = data.frame(beta_raw = round(beta_raw, 3),
                        #se = round(se, 3),
                        beta_var = round(beta_var, 3),
                        beta_mean = round(beta_mean, 3))

resmat_PC


#### Univariate selection gradients ####

resmat_uv = matrix(NA, nrow = ncol(XRRRData), ncol = 4)

for(i in 1:ncol(XRRRData)){
linmod = lm(Y[,1]/mean(Y[,1]) ~  date_flow + pl_height + nb_tot_fl + CA + SL + XRRRData[,i], data=XData)
beta_raw = summary(linmod)$coef[7,1]
se = summary(linmod)$coef[7,2]
beta_var = beta_raw * sd(XRRRData[,i])
beta_mean = beta_raw * mean(XRRRData[,i])
resmat_uv[i, ] = round(c(beta_raw, se, beta_var, beta_mean), 3)
}

rownames(resmat_uv) = colnames(XRRRData)
resmat_uv = as.data.frame(resmat_uv)
colnames(resmat_uv) = c("beta", "se", "beta_var", "beta_mean")

resmat_uv

#### Multivariate selection gradients ####

XFormula = as.formula(paste("~", paste0(c(colnames(XData), colnames(XRRRData)), collapse=" + ")))
XFormula

regdat = data.frame(XData, XRRRData)

linmod = lm(Y[,1]/mean(Y[,1]) ~ date_flow + pl_height + nb_tot_fl + CA + SL + 
              phenylacetaldehyde_D + 
              X2_phenylethanol_D + 
              X2_phenylethanol_N + 
              X2_phenylethylacetate_N + 
              indole_D + 
              indole_N + 
              X2_aminobenzaldehyde_D + 
              elemicin_D + 
              elemicin_N + 
              methyleugenol_N + 
              benzyl_alcohol_D + 
              benzyl_alcohol_N + 
              p_cresol_D + 
              Z_5_dodecenyl_acetate_D,
              data=regdat)
  
beta_raw = summary(linmod)$coef[-1,1]
se = summary(linmod)$coef[-1,2]
beta_var = beta_raw * apply(regdat, 2, sd)
beta_mean = beta_raw * colMeans(regdat)
resmat_mv = data.frame(beta_raw = round(beta_raw, 2),
                    se = round(se, 2),
                    beta_var = round(beta_var, 3),
                    beta_mean = round(beta_mean, 3))

rownames(resmat_mv) = colnames(regdat)

resmat_mv



#### Projection pursuit regression ####

library(gsg)

regdat = data.frame(XData, XRRRData)
regdat$fitness=Y[,1]
head(regdat)

?ppr
ppr_mod = ppr(Y[,1]/mean(Y[,1]) ~ 
                #date_flow + pl_height + nb_tot_fl + CA + SL + 
                phenylacetaldehyde_D + 
                X2_phenylethanol_D + 
                X2_phenylethanol_N + 
                X2_phenylethylacetate_N + 
                indole_D + 
                indole_N + 
                X2_aminobenzaldehyde_D + 
                elemicin_D + 
                elemicin_N + 
                methyleugenol_N + 
                benzyl_alcohol_D + 
                benzyl_alcohol_N + 
                p_cresol_D + 
                Z_5_dodecenyl_acetate_D,
                data=regdat,
                optlevel=3,
                nterms=1)
ppr_mod
summary(ppr_mod)
ppr_mod$alpha
ppr_mod$beta


comp=as.matrix(XRRRData) %*% as.matrix(ppr_mod$alpha)
comp=c(comp)
hist(comp)
plot(comp[,1], predict(ppr_mod))
plot(comp[,2], predict(ppr_mod))

gm=gam(Y[,1]/mean(Y[,1]) ~ date_flow + pl_height + nb_tot_fl + CA + SL + s(comp), data=regdat)
summary(gm)
plot(comp[,2], predict(gm))

res = vis.gam(gm, c("nb_tot_fl", "comp"), plot.type="contour", type="response",
              xlab="Flower number",
              ylab="Major selection axis 1",
              main="", color="heat")
points(regdat$nb_tot_fl, comp, col="grey", pch=16, cex=1*Y[,1]/mean(Y[,1]))

?gam.gradients
res=gam.gradients(gm, phenotype="comp", covariates=c("date_flow", "pl_height", "nb_tot_fl", "CA", "SL"), standardized=T, n.boot=10)
res

plot(predict(ppr_mod), Y[,1]/mean(Y[,1]))
lines(-1:5, -1:5)

getPostEstimate(m, "Beta")$mean[2,1]
wRRR = getPostEstimate(m, "wRRR")$mean #rows are RRR axes, cols are original variables
wRRR

plot(ppr_mod$alpha, wRRR)
abline(h=0, lty=2)
abline(v=0, lty=2)
lines(-10:10, -10:10)
plot(ppr_mod$alpha[,2], wRRR$mean[2,])


summary(ppr_mod)


ppr1=
W=predict(ppr_mod)

tmp1=as.matrix(regdat[,6:19])
tmp2=as.matrix(ppr_mod$alpha)
dim(tmp1)
dim(tmp2)
ppr1=as.matrix(regdat[,6:19]) %*% as.matrix(ppr_mod$alpha)
ppr1

head(regdat)

gamdat = data.frame(W=W, r1=ppr1[,1], r2=ppr1[,2])
head(gamdat)

gm = gam(W ~ r1 + r2, family=gaussian, data=gamdat)
summary(gm)
res = vis.gam(gm, c("r1","r2"), plot.type="contour", type="response",
              xlab="Major selection axis 1",
              ylab="Major selection axis 2",
              main="", color="heat")
points(gamdat$r1, gamdat$r2, col="grey", pch=16, cex=1*gamdat$W)




colMeans(regdat)

regdat_scaled=as.data.frame(apply(regdat[,1:19], 2, scale))
regdat_scaled$fitness = regdat$fitness

head(regdat_scaled)

ppr_mod = gppr(y = "fitness",
               xterms = c("phenylacetaldehyde_D", 
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
                          "Z_5_dodecenyl_acetate_D"),
               family = "poisson",
               data = regdat_scaled,
               nterms = 1,
               maxit = 50)


ppr_mod = gppr(y = "fitness",
               xterms = c("nb_tot_fl"), 
               family = "poisson",
               data = regdat_scaled,
               nterms = 1,
               tol = 1,
               maxit = 50)

summary(ppr_mod)
ppr_mod$ppr

gppr.gradients
out = gppr.gradients(ppr_mod, phenotype=c("pl_height","nb_tot_fl"),
                     n.boot=10)
str(out)
out$ests$estimates



# simulated data (two traits, stabilizing selection on trait 1)
n<-250
z<-cbind(rnorm(n,0,1),rnorm(n,0,1))
W<-rpois(n,exp(2-0.6*z[,1]^2))
d<-as.data.frame(cbind(W,z))
names(d)<-c("W","z1","z2")

fit.func<-gppr(y="W", xterms=c("z1","z2"), data=d, family="poisson",
               nterms=2, max.terms=2)
fit.func$ppr
gppr.gradients(mod= fit.func,phenotype=c("z1","z2"),se.method='n',standardize=FALSE)











?gppr.gradients
out = gppr.gradients(ppr_mod, phenotype=c("phenylacetaldehyde_D", 
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
                              "Z_5_dodecenyl_acetate_D"),
                              n.boot=10)
str(out)
out$ests$estimates

#### Compare results ####
resmat_uv
resmat_mv
resmat_RRR
resmat_PC


plot(resmat_uv[,1], resmat_mv[-c(1:5),1])
abline(h=0, lty=2)
abline(v=0, lty=2)
lines(-100:100, -100:100)

plot(resmat_uv[,1], resmat_RRR[,1])
abline(h=0, lty=2)
abline(v=0, lty=2)
lines(-100:100, -100:100)

plot(resmat_uv[,3], resmat_RRR[,1]) #Variance-scaled uv vs. raw RRR
abline(h=0, lty=2)
abline(v=0, lty=2)
lines(-100:100, -100:100)

plot(resmat_mv[-c(1:5),1], resmat_RRR[,2]) #Raw
abline(h=0, lty=2)
abline(v=0, lty=2)
lines(-100:100, -100:100)

plot(resmat_mv[-c(1:5),3], resmat_RRR[,1]) #Variance-scaled
abline(h=0, lty=2)
abline(v=0, lty=2)
lines(-100:100, -100:100)

plot(resmat_mv[-c(1:5),3], resmat_PC[,2])
abline(h=0, lty=2)
abline(v=0, lty=2)
lines(-100:100, -100:100)

plot(resmat_RRR[,2], resmat_PC[,2], ylim=c(-.1, .1), xlim=c(-.1, .1))
abline(h=0, lty=2)
abline(v=0, lty=2)
lines(-100:100, -100:100)








y ="fitness"
xterms =c("pl_height","nb_tot_fl")
data  = regdat
nterms = 1
tol = 0.001 
gcvpen = 1
maxit = 50
family = "gaussian"
max.terms = 2

gppr2=function (y, xterms, data, nterms = 1, tol = 0.001, gcvpen = 1, 
                maxit = 50, family = "binomial", max.terms = 2) 
{
  fam <- get(family, mode = "function", envir = parent.frame())
  fam <- fam()
  g <- fam$linkfun
  g.inv <- fam$linkinv
  V <- fam$variance
  g.prime <- function(x) {
    1/V(x)
  }
  response <- data[, as.character(y)]
  mu <- rep(mean(response), length(response))
  mu <- rep(0.5, length(response))
  deltaMu <- 1
  iter <- 0
  gppr <- NULL
  f <- paste("Z", "~", paste(xterms, collapse = "+"), collapse = "")
  while (deltaMu >= tol & iter < maxit) {
    Z <- g(mu) + g.prime(mu) * (response - mu)
    cur.weights <- 1/(g.prime(mu)^2 * V(mu))
    data$Z <- Z
    gppr <- ppr(as.formula(f), weights = cur.weights, sm.method = "gcvspline", 
                gcvpen = gcvpen, nterms = nterms, data = data, max.terms = max.terms, 
                optlevel = 3)
    eta <- predict(gppr, type = "raw")
    muPrime <- mu
    mu <- g.inv(eta)
    if (family == "binomial") {
      boundary.tol <- 10^(-5)
      mu[which(mu < boundary.tol)] <- boundary.tol
      mu[which(mu > (1 - boundary.tol))] <- (1 - boundary.tol)
    }
    deltaMu <- sum(abs(muPrime - mu))/sum(abs(muPrime))
    iter <- iter + 1
  }
  output <- list(ppr = gppr, family = fam, iterations = iter, 
                 data = data, nterms = nterms, tol = tol, gcvpen = gcvpen, 
                 maxit = maxit, max.terms = max.terms, y = y, xterms = xterms, 
                 formula = f)
  class(output) <- c("gppr")
  return(output)
}




#### RRR package ####
library(rrr)
?rrr

head(Y)
rrr_fit=rrr(x = scale(as.matrix(XRRRData)),
            y = as.matrix(Y[,1]/mean(Y)),
            type = "identity",
            rank = 1)
rrr_fit

as.matrix(XRRRData)            
