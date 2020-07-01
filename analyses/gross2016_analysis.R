#######################################################################################
##### Using HMSC reduced-rank-regression to measure selection on floral volatiles #####
#######################################################################################

# Gross et al. 2016 Plos ONE

rm(list=ls())
#setwd("C:/data/reduced_rank_selection/gross")

library(Hmsc)
library(abind)
library(mgcv)
library(gsg)
library(knitr)

indat = read.table("gross.txt", header=T)

indat$study = as.factor(paste(indat$Population, indat$Year, sep="_"))

table(indat$study)

#Add for loop here
studies = levels(indat$study)

#### Fit models and compute selection gradients ####
MFList = list()
resList = list()
betaList = list()

#Set sampling parameters
samples = 1000
thin = 10
transient = .5*(thin*samples)
adaptNf = 0.4*(thin*samples)
nChains = 1

a = Sys.time()

for(s in 1:length(studies)){

indatpop = indat[indat$study==studies[s],]

Y = cbind(indatpop$NrFruits)
XData = indatpop[,c(7, 8, 9)] #3 other traits
XRRRData = indatpop[,10:31] #volatiles

studyDesign = data.frame(Plant=as.factor(indatpop$PlantID))
studyDesign$Plant = factor(studyDesign$Plant)
rL = HmscRandomLevel(units=unique(studyDesign$Plant))

XFormula = as.formula(paste("~", paste0(colnames(XData), collapse=" + ")))
XRRRFormula = as.formula(paste("~", " -1 + ", paste0(colnames(XRRRData), collapse=" + ")))

XRRRDataScaled = data.frame(scale(XRRRData))

#Some probable typos, need to pass >1 Y variable
m = Hmsc(Y = as.matrix(Y/mean(Y)), XData=XData, XFormula=XFormula,
         distr="normal", 
         XRRRData=XRRRDataScaled, XRRRFormula=XRRRFormula, ncRRR=1, XRRRScale = TRUE,
         studyDesign=studyDesign, ranLevels=list(Plant = rL))

m = sampleMcmc(m, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
               transient = transient, nChains = nChains)

filename = paste0("model_",studies[s],"_thin_",thin,".RData")
save(m, file=filename)

#### - Assess MCMC convergence and evaluate model fit - ####

post = convertToCodaObject(m) #wRRR needs to be added here

#Temp. posterior of wRRR
bind0 = function(...) {abind(..., along = 0)}
postList = poolMcmcChains(m$postList, chainIndex = 1:length(m$postList), 
                          start = 1, thin = 1)
valList = lapply(postList, function(a) a[["wRRR"]])
wRRR_post = do.call(bind0, valList)

predY = computePredictedValues(m)
predYm = apply(simplify2array(predY), 1:2, mean)
MF = evaluateModelFit(m, predY)
MFList[[s]] = MF

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

#The first composite scent trait
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

#Compute beta = wA for each posterior sample
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

betaList[[s]] = betamat_RRR

}
Sys.time()-a
#thin 10 1.11 hours

save(resList, file ="resList.RData")
save(MFList, file ="MFList.RData")
save(betaList, file ="betaList.RData")

#### Load results ####
load(file ="analyses/gross/resList.RData")
load(file ="analyses/gross/MFList.RData")
load(file ="analyses/gross/betaList.RData")

MFvals = unlist(lapply(MFList, function(x) x$R2[1]))
hist(MFvals)
mean(MFvals)
range(MFvals)

summary_mat = matrix(NA, nrow = 4, ncol = length(studies))
for(i in 1:length(studies)){
  summary_mat[,i] = round(as.numeric(resList[[i]]$beta_var),2)
}

colnames(summary_mat) = studies
rownames(summary_mat) = rownames(resList[[1]])
summary_mat = t(summary_mat)
summary_mat

upper_mat = matrix(NA, nrow = 4, ncol = length(studies))
for(i in 1:length(studies)){
  upper_mat[,i] = round(as.numeric(resList[[i]]$U.1),2)
}

colnames(upper_mat) = studies
rownames(upper_mat) = rownames(resList[[1]])
upper_mat = t(upper_mat)

lower_mat = matrix(NA, nrow = 4, ncol = length(studies))
for(i in 1:length(studies)){
  lower_mat[,i] = round(as.numeric(resList[[i]]$L.1),2)
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

par(mar=c(6,4,2,2))
plot(1:13, summary_mat[,4], pch=16, ylim=c(-.5, 1), las=1, 
     xaxt="n", xlab="", ylab="Selection intensity on scent")
segments(1:13, lower_mat[,4], 1:13, upper_mat[,4])
abline(h=0, lty=2)
axis(1, 1:13, labels=F)
text(1:13, par("usr")[3] - .1, srt = 45, adj = 1,cex=.7,
     labels = rownames(summary_mat), xpd = TRUE)



me = melt(data.frame(indat[,7:9], indat$study))
popmeans = tapply(me$value, list(me$indat.study, me$variable), mean)
dmat = cov(log(popmeans))

B = cov(summary_mat[,1:3])

plot(apply(summary_mat[,1:3], 2, sd), diag(dmat))

XRRRData = indat[,10:31]
beta_summary_mat = matrix(NA, nrow = ncol(XRRRData), ncol = length(studies))
for(i in 1:length(studies)){
  beta_summary_mat[,i] = round(betaList[[i]]$beta_mean,1)
}
colnames(beta_summary_mat) = studies
rownames(beta_summary_mat) = colnames(XRRRData)

#View(beta_summary_mat)

meltmeans=melt(beta_summary_mat)
head(meltmeans)
#meltmeans$value = abs(meltmeans$value)

means= tapply(meltmeans$value, meltmeans$Var1, median)
meltmeans$Var1=factor(meltmeans$Var1, levels=names(sort(means,decreasing=T)))

par(mfrow=c(1,1), mar=c(6,4,2,2))
plot(meltmeans$Var1, meltmeans$value, las=1, xaxt="n", ylab="Selection gradient (%)")
abline(h=0)

axis(1, 1:22, labels=F)
text(1:22, par("usr")[3] - 8, srt = 45, adj = 1,cex=.7,
     labels = gsub("Z","",gsub("_ngPerL","",levels(meltmeans$Var1))), xpd = TRUE)



me=melt(data.frame(indat[,10:31], indat$study))
popmeans = tapply(me$value, list(me$indat.study, me$variable), mean)


#popmeans
dmat = cov(log(popmeans))

B=sqrt(cov(t(beta_summary_mat)))

plot(apply(beta_summary_mat, 1, sd), diag(dmat))
plot(B, dmat)
cor(c(B), c(dmat), "pairwise")

library(evolvability)



#### Model output ####
load(file ="analyses/gross/resList.RData")

summary_mat = matrix(NA, nrow = 4, ncol = length(studies))
for(i in 1:length(studies)){
  summary_mat[,i] = round(as.numeric(resList[[i]]$beta_var),2)
}

colnames(summary_mat) = studies
rownames(summary_mat) = rownames(resList[[1]])
summary_mat = t(summary_mat)
summary_mat=as.data.frame(summary_mat)

summary_mat$altitude=factor(c("Mountain2010", "Mountain2011","Lowland2010","Lowland2011","Lowland2010","Lowland2011",
                       "Mountain2010","Mountain2011","Lowland2010","Lowland2011","Lowland2011","Mountain2010","Mountain2011"))
head(summary_mat)

names(summary_mat) = c("Beta_PlantHeight", "Beta_InfLength", "Beta_NrFlowers", "Beta_Scent", "Altitude")
str(summary_mat)

cols=topo.colors(4)

x11(height=5, width=5.5)
par(mfrow=c(1,1))
plot(summary_mat$Altitude, summary_mat$Beta_PlantHeight, at=1:4, xaxt="n", las=1, xlim=c(0,20), ylim=c(-.2, 1), col=cols)
mtext("Selection intensity", 2, line=2.5)
par(new=T)
plot(summary_mat$Altitude, summary_mat$Beta_InfLength, at=6:9, xaxt="n", yaxt="n",las=1,xlim=c(0,20), ylim=c(-.2, 1), col=cols)
par(new=T)
plot(summary_mat$Altitude, summary_mat$Beta_NrFlowers, at=11:14, xaxt="n", yaxt="n",las=1,xlim=c(0,20), ylim=c(-.2, 1), col=cols)
par(new=T)
plot(summary_mat$Altitude, abs(summary_mat$Beta_Scent), at=16:19, xaxt="n", yaxt="n", las=1,xlim=c(0,20), ylim=c(-.2, 1), col=cols)
axis(1, at=c(2.5, 7.5, 12.5, 17.5), c("Plant height", "Inf. length", "Flowers", "Scent"))

legend("topleft", legend=c("Lowland 2010", "Lowland 2011", "Mountains 2010", "Mountains 2011"), 
       pch=15, pt.cex = 1.5, col=cols, bty="n")


load(file ="analyses/gross/betaList.RData")

beta_summary_mat = matrix(NA, nrow = ncol(XRRRData), ncol = length(studies))
for(i in 1:length(studies)){
  beta_summary_mat[,i] = round(betaList[[i]]$beta_mean,1)
}
colnames(beta_summary_mat) = studies
rownames(beta_summary_mat) = colnames(XRRRData)

beta_summary_mat = as.data.frame(t(beta_summary_mat))

beta_summary_mat$altitude=factor(c("Mountain2010", "Mountain2011","Lowland2010","Lowland2011","Lowland2010","Lowland2011",
                              "Mountain2010","Mountain2011","Lowland2010","Lowland2011","Lowland2011","Mountain2010","Mountain2011"))
head(beta_summary_mat)

str(beta_summary_mat)

sort(colMeans(beta_summary_mat[,-23]))
sort(colMeans(abs(beta_summary_mat[,-23])))

pdf("beta_compound_plots.pdf")
for(i in 1:22){
plot(beta_summary_mat$altitude, beta_summary_mat[,i], las=1, main=colnames(beta_summary_mat[i]))
abline(h=0, lwd=2)
}
dev.off()

#### Output per model ####
studies
s=10
thin=10
study=studies[s]

indatpop = indat[indat$study==studies[s],]
  
Y = cbind(indatpop$NrFruits, indatpop$NrFruits)
XData = indatpop[,c(7, 8, 9)] #3 other traits
XRRRData = indatpop[,10:31] #volatiles

XRRRDataScaled=data.frame(scale(XRRRData))

filename = paste0("model_",studies[s],"_thin_",thin,".RData")
load(filename)

#Model fit
post = convertToCodaObject(m) #wRRR needs to be added here
plot(post$Beta[,5])

#Temp. posterior of wRRR
bind0 = function(...) {abind(..., along = 0)}
postList = poolMcmcChains(m$postList, chainIndex = 1:length(m$postList), 
                          start = 1, thin = 1)
valList = lapply(postList, function(a) a[["wRRR"]])
wRRR_post = do.call(bind0, valList)

x11()
par(mfrow=c(4,5))
for(i in 1:20){
  plot(wRRR_post[,1,i], type="l")
}

predY = computePredictedValues(m)
predYm = apply(simplify2array(predY), 1:2, mean)
MF = evaluateModelFit(m, predY)
MF
plot(m$Y[,1], predYm[,1])
lines(-10:10, -10:10)

#Plotting the estimated fitness surface
wRRR = getPostEstimate(m, "wRRR")
rrr = as.matrix(m$XRRRData) %*% t(as.matrix(wRRR$mean))

gamdat = data.frame(W=m$Y[,1], r1=rrr[,1], m$XData)
gm = gam(W ~ r1 + PlantHeight_cm + InflorescenceLength_cm + NrFlowers, family=gaussian, data=gamdat)

gam.gradients(gm, phenotype=c("r1"), covariates = c("PlantHeight_cm", "InflorescenceLength_cm", "NrFlowers"), standardized=T, se.method="n")

x11(height = 3.5, width = 7)
par(mfrow=c(1,2))
par(mar=c(4,4,2,2))
res = vis.gam(gm, c("NrFlowers","r1"), plot.type="contour", type="response",
              xlab="Flower number",
              ylab="Scent selection axis 1",
              main="", color="bw")

points(gamdat$NrFlowers, gamdat$r1, col="grey", pch=16, cex=1*gamdat$W)

par(mar=c(4,6,2,2), xpd=T)
out = barplot(sort(beta_summary_mat[,s], dec=F), xlim=c(-50,50), col="royalblue", yaxt="n",
              xlab="Mean-scaled selection gradient (%)", las=1, horiz=T)
text(y=out, par("usr")[3] + -40, srt = 0, adj = 1,cex=.7,
     labels = gsub("X","",colnames(XRRRData)[order(beta_summary_mat[,s], decreasing=F)]), xpd = TRUE)

#arrows(out, betamat_RRR$L, out, betamat_RRR$U, angle=90, length=0.05, code=3)

lm_res=summary(lm(W~scale(NrFlowers)+scale(PlantHeight_cm)+scale(InflorescenceLength_cm)+scale(r1), data=gamdat))
lm_res

#### Multiple-regression ####
regdat = data.frame(XData, XRRRData)
names(regdat)

formula=as.formula(paste("Y[,1]/mean(Y[,1]) ~ ", paste0(colnames(regdat), collapse=" + ")))
#formula

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
plot(resmat_mv[-c(1:3),4], beta_summary_mat[,s], las=1, pch=16,
     xlab="Selection gradient via multiple-regression (%)",
     ylab="Selection gradient via reduced-rank regression (%)")
abline(h=0, lty=2)
abline(v=0, lty=2)
lines(-100:100, -100:100)

