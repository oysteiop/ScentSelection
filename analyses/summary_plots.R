rm(list=ls())

# Betas
load(file="results/bdf_godo.RData")
bdat = bdf
names(bdat) = c("Mean", "Lower", "Upper")

load(file="results/bdf_anac.RData")
names(bdf) = c("Mean", "Lower", "Upper")
bdat = rbind(bdat, bdf)

load(file="results/bdf_pens.RData")
bdat[nrow(bdat)+1, ] = bdf

load(file="results/bdf_gcon.RData")
bdat[nrow(bdat)+1, ] = bdf

rownames(bdat)[21:22] = c("Penstemon", "G_conopsea")

bdat

mean(bdat$Mean)
median(bdat$Mean)
range(bdat$Mean)

# Plot
par(mfrow=c(1,1))
plot(1:22, bdat$Mean, pch=16, ylim=c(-0.3, 1))
segments(1:22, bdat$Lower, 1:22, bdat$Upper)
abline(h=0, lty=2)

hist(bdat$Mean, breaks=20, las=1, main="", xlab="Selection intensity")

po = rev(order(bdat$Mean))

plot(1:22, bdat$Mean[po], pch=16, ylim=c(-0.3, 1),
     ylab="Selection intensity (95% CI)")
segments(1:22, bdat$Lower[po], 1:22, bdat$Upper[po])
abline(h=0, lty=2)



#####

load(file="results/mf_godo.RData")
mfdat = mfvals

load(file="results/mf_anac.RData")
mfdat = rbind(mfdat, mfvals)
mfdat$studies = as.character(mfdat$studies)

load(file="results/mf_pens.RData")
mfdat[nrow(mfdat)+1, ] = t(as.data.frame(mfvals))

load(file="results/mf_gcon.RData")
mfdat[nrow(mfdat)+1, ] = t(as.data.frame(mfvals))

mfdat[,2:6] = apply(mfdat[2:6], 2, as.numeric)
names(mfdat) = c("Study", "MFvals", "MF2vals", "MFCVvals", "MFCV2vals", "cor_MR_RRR")

mfdat
hist(mfdat$cor)
median(mfdat$cor_MR_RRR)
mean(mfdat$cor_MR_RRR)
range(mfdat$cor_MR_RRR)

apply(mfdat[,-1], 2, mean)
apply(mfdat[,-1], 2, median)

# Plot

x11(height=4.5, width=9)
par(mfrow=c(1,2))
plot(mfdat$MFvals, mfdat$MF2vals, xlim=c(0,100), ylim=c(0,100), las=1,
     xlab="Reduced-rank regression", ylab="Multiple regression")
points(mfdat$MFCVvals, mfdat$MFCV2vals, pch=16)
lines(-10:1000, -10:1000)
legend("topleft", pch=c(1,16), c("Explanatory R^2", "Predictive R^2"))

plot(mfdat$MFvals, mfdat$MFCVvals, xlim=c(0,100), ylim=c(0,100), las=1,
     xlab="Explanatory R^2", ylab="Predictive R^2")
points(mfdat$MF2vals, mfdat$MFCV2vals, pch=16)
lines(-10:1000, -10:1000)
legend("topleft", pch=c(1,16), c("Reduced-rank regression", "Multiple regression"))

# Combine all ####
bdat$Study = rownames(bdat)

alldat = merge(bdat, mfdat, by="Study")
alldat

alldat$Study
alldat$n = c(54, 48, 55, 44, 60, 54, 54, 69, 82, 73, 92, 169, 92, 92, 96, 94, 88, 88, 56, 72, 47, 75)
alldat$nvol = c(32, 27, 31, 29, 26, 30, 32, 22, 22, 22, 22, 14, 22, 22, 22 ,22, 23, 22, 22, 22, 22, 22)

names(alldat)
par(mfrow=c(1,1))
plot(alldat$Mean, alldat$cor_MR_RRR)
plot(alldat$n, alldat$cor_MR_RRR)
plot(alldat$nvol, alldat$cor_MR_RRR)

plot(alldat$MFvals, alldat$cor_MR_RRR)
points(alldat$MF2vals, alldat$cor_MR_RRR, pch=16)

plot(alldat$MFCVvals, alldat$cor_MR_RRR)
points(alldat$MFCV2vals, alldat$cor_MR_RRR, pch=16)

plot(alldat$MFCVvals-alldat$MFCV2vals, alldat$cor_MR_RRR)

summary(lm(cor_MR_RRR~Mean + n + nvol + MFvals, data=alldat))
