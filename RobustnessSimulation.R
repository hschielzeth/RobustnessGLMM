rm(list=ls())

basepath = "TOBEADDED"
sim = FALSE
plotYN = FALSE
nrep = 10000

require(lme4)
require(lmerTest)
require(PearsonDS)
require(arm)

for(set in c("Basic", "SmallSample", "Balanced", "FewGroups", "LowRep", "HighRep", "Poisson", "Proportion")) {

####################################################################################
## Setting  ########################################################################
####################################################################################

if(set == "FewGroups") {
path = paste0(basepath, "FewGroups\\")
nobs = 120
ngroups = 6
nhilevelgroups = 3
nlolevelgroups = 12
ncrossedgroups = 6
unbalancedYN = TRUE
vargroups = 0.5
varres = 0.5
intercept = 1
GLMfamily = "Gaussian"
xlimresp=c(-2,5)
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
breaksVarComp=seq(0,2,by=0.025)
breaksSlopes=seq(0,2,by=0.025)
breaksSlopesLast=seq(0,2,by=0.025)
SEfactor=1.98
}


if(set == "Basic") {
path = paste0(basepath, "Basic\\")
nobs = 120
ngroups = 30
nhilevelgroups = 10
nlolevelgroups = 60
ncrossedgroups = 30
unbalancedYN = TRUE
vargroups = 0.5
varres = 0.5
intercept = 1
GLMfamily = "Gaussian"
xlimresp=c(-2,5)
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
breaksVarComp=seq(0,2,by=0.025)
breaksSlopes=seq(0,2,by=0.025)
breaksSlopesLast=seq(0,2,by=0.025)
SEfactor=1.98
}

if(set == "SmallSample") {
path = paste0(basepath, "SmallSample\\")
nobs = 40
ngroups = 10
nhilevelgroups = 5
nlolevelgroups = 20
ncrossedgroups = 10
unbalancedYN = TRUE
vargroups = 0.5
varres = 0.5
intercept = 1
GLMfamily = "Gaussian"
xlimresp=c(-2,5)
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
SEfactor=2.03
}

if(set == "LowRep") {
path = paste0(basepath, "LowRep\\")
nobs = 120
ngroups = 30
nhilevelgroups = 10
nlolevelgroups = 60
ncrossedgroups = 30
unbalancedYN = TRUE
vargroups = 0.1
varres = 0.9
intercept = 1
GLMfamily = "Gaussian"
xlimresp=c(-2,5)
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
}

if(set == "HighRep") {
path = paste0(basepath, "HighRep\\")
nobs = 120
ngroups = 30
nhilevelgroups = 10
nlolevelgroups = 60
ncrossedgroups = 30
unbalancedYN = TRUE
vargroups = 0.9
varres = 0.1
intercept = 1
GLMfamily = "Gaussian"
xlimresp=c(-2,5)
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
SEfactor=1.98
}

if(set == "Balanced") {
path = paste0(basepath, "Balanced\\")
nobs = 120
ngroups = 30
nhilevelgroups = 10
nlolevelgroups = 60
ncrossedgroups = 30
unbalancedYN = FALSE 
vargroups = 0.5
varres = 0.5
intercept = 1
GLMfamily = "Gaussian"
xlimresp=c(-2,5)
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
SEfactor=1.98
}

if(set == "Poisson") {
path = paste0(basepath, "Poisson\\")
nobs = 120
ngroups = 30
nhilevelgroups = 10
nlolevelgroups = 60
ncrossedgroups = 30
unbalancedYN = TRUE
vargroups = 0.5
varres = 0.5
intercept = 1
GLMfamily = "Poisson"
xlimresp=c(0,20)
ylimBias1 = c(-0.10,0.10)
ylimBias2 = c(-0.3,0.40)
ylimBias3 = c(-0.2,0.3)
ylimRMSE = c(1.2, 1.5, 1.2)
SEfactor=1.98
# Note at latent values are divided by 2 and variance components post-multiplied by 4 and slopes by 2
}

if(set == "Proportion") {
path = paste0(basepath, "Proportion\\")
nobs = 120
ngroups = 30
nhilevelgroups = 10
nlolevelgroups = 60
ncrossedgroups = 30
unbalancedYN = TRUE
vargroups = 0.5
varres = 0.5
intercept = 0
GLMfamily = "Proportion"
xlimresp=c(0,1)
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.10,0.20)
ylimBias3 = c(-0.10,0.10)
ylimRMSE = c(0.9, 1.3, 0.7)
SEfactor=1.98
}

####################################################################################
## Functions  ######################################################################
####################################################################################

datsim = function(vargroups = 0.3, varres = 0.7, slope1 = 0.5, slope2 = 0.5, ngroups = 10, nobs = 50, intercept = 1, nrep = 100,
					skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor = 0, datlevelpred=c(0,0), grouplevelpred=c(0,0),
					varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, 
					nhilevelgroups = 5, nlolevelgroups = 40, ncrossedgroups = 20,
					heteroscedfactor = 0, unbalancedYN=TRUE, GLMfamily="Gaussian", PoissonScaler = 2) {
	if(nhilevelgroups > ngroups) stop("nhilevelgroups needs to be smaller or equal to ngroups") 
	if(nlolevelgroups < ngroups) stop("nlolevelgroups needs to be smaller or equal to ngroups") 
	res = data.frame(Intercept=rep(NA,nrep), slope1=NA, slope2=NA, slope1SE=NA, slope2SE=NA, vargroups=NA, varres=NA, SimSlope=slope1, SimGroups=vargroups, SimRes=varres, xCor=xCor)
	for(i in 1:nrep) {
		# Covariates (level 1)
		x = MASS::mvrnorm(nobs,c(0,0),matrix(c(1,xCor,xCor,1),ncol=2))
		x1 = x[,1]; x2 = x[,2]
		# Group labels
		if(unbalancedYN) {
			logroups = c(1:nlolevelgroups, sample(1:nlolevelgroups, nobs-nlolevelgroups, replace=TRUE))
			groups = c(1:ngroups, sample(1:ngroups, nlolevelgroups-ngroups, replace=TRUE))[logroups]
			higroups = c(1:nhilevelgroups, sample(1:nhilevelgroups, ngroups-nhilevelgroups, replace=TRUE))[groups]
			crgroups = sample(c(1:ncrossedgroups, sample(1:ncrossedgroups, nobs-ncrossedgroups, replace=TRUE)), nobs, replace=FALSE)
		}
		if(!unbalancedYN) {
			if(nobs%%nlolevelgroups!=0) stop("the number of observations has to be a multiple of the number of low-level groups with balanced sampling") 
			if(nlolevelgroups%%ngroups!=0) stop("the number of low-groups has to be a multiple of the number of groups with balanced sampling") 
			if(ngroups%%nhilevelgroups!=0) stop("the number of groups has to be a multiple of the number of high-level groups with balanced sampling") 
			if(ngroups!=ncrossedgroups) stop("the number of crossed level groups has to be a multiple of the number of groups with balanced sampling") 
			logroups = rep(1:nlolevelgroups, nobs/nlolevelgroups)
			groups = rep(1:ngroups, nlolevelgroups/ngroups)[logroups]
			higroups = rep(1:nhilevelgroups, ngroups/nhilevelgroups)[groups]
			crgroups = rep(1:ncrossedgroups, each=nobs/ncrossedgroups)
		}	
		# Group means
		grmeans = PearsonDS::rpearson(ngroups, moments = c(mean = 0,variance = vargroups, skewness = skewGroups, kurtosis = kurtGroups))
		if(!all(grouplevelpred==0)) {
			sdgrmeans = sd(grmeans)
			x4 = sample(grouplevelpred, ngroups, replace=TRUE)[groups]
			grmeans = grmeans + x4
			grmeans = grmeans/sd(grmeans)*sdgrmeans #rnorm(1, sqrt(vargroups), sd(replicate(1000, sd(rnorm(ngroups)))))
		}
		crgrmeans = rnorm(ncrossedgroups, 0, sqrt(varcrossedgroups))
		logrmeans = rnorm(nlolevelgroups, 0, sqrt(varlolevelgroups))
		higrmeans = rnorm(nhilevelgroups, 0, sqrt(varhilevelgroups))
		# Residuals
		if(heteroscedfactor!=0) {
			sdresdev = sd(rnorm(nobs, 0, sqrt(varres)))
			resdev = rnorm(nobs, 0, sqrt(varres) + heteroscedfactor*(x1-min(x1)) ) 
			resdev = resdev/sd(resdev)*sdresdev #rnorm(1, sqrt(varres), sd(replicate(1000, sd(rnorm(nobs)))))
		}
		else resdev = PearsonDS::rpearson(nobs, moments = c(mean = 0,variance = varres, skewness = skewRes, kurtosis = kurtRes))
		if(!all(datlevelpred==0)) {
			sdresdev = sd(resdev)
			x3 = sample(datlevelpred, nobs, replace=TRUE)
			resdev = resdev + x3
			resdev = resdev/sd(resdev)*sdresdev # rnorm(1, sqrt(varres), sd(replicate(1000, sd(rnorm(nobs)))))
		}
		# Data generating model 
		resp =  intercept + slope1*x1 + slope2*x2 + grmeans[groups] + logrmeans[logroups] + higrmeans[higroups] + crgrmeans[crgroups] + resdev
		if(GLMfamily=="Poisson") resp = rpois(nobs, exp(resp/PoissonScaler))
		if(GLMfamily=="Proportion") resp = rbinom(nobs, 20, invlogit(resp))
		# Model fitting
		if(GLMfamily == "Poisson" | GLMfamily == "Proportion") obsid = as.factor(1:nobs)
		if(GLMfamily=="Gaussian") modraw = lmer(resp ~ x1 + x2 + (1|groups))
		if(GLMfamily=="Poisson") modraw = glmer(resp ~ x1 + x2 + (1|groups) + (1|obsid), family="poisson")
		if(GLMfamily=="Proportion") modraw = glmer(cbind(resp,20-resp) ~ x1 + x2 + (1|groups) + (1|obsid), family="binomial")
		mod = summary(modraw)
		res[i,c("Intercept", "slope1", "slope2")] = coefficients(mod)[,1]
		res[i,c("slope1SE", "slope2SE")] = coefficients(mod)[2:3,2]
		if(GLMfamily == "Gaussian") {
			res[i,"vargroups"] = mod$varcor$groups
			res[i,"varres"] = attr(mod$varcor,"sc")^2
		}
		if(GLMfamily == "Poisson") {
			res[i,"vargroups"] = mod$varcor$groups[[1]] * PoissonScaler^2
			res[i,"varres"] = mod$varcor$obsid[[1]]	* PoissonScaler^2		
			res[i,c("Intercept", "slope1", "slope2")] = coefficients(mod)[,1] * PoissonScaler
			res[i,c("slope1SE", "slope2SE")] = coefficients(mod)[2:3,2] * PoissonScaler
		}
		if(GLMfamily == "Proportion") {
			res[i,"vargroups"] = mod$varcor$groups[[1]]
			res[i,"varres"] = mod$varcor$obsid[[1]]			
		}
	}
	lastSim = data.frame(grmeans=grmeans, resdev=resdev, resp=resp)
	res = list(res=res, lastSim=lastSim, lastModraw=modraw, lastMod=mod)
	return(res)
}

plotLastSimHist = function(lastSim, histcol=rep("grey", 4), xlimgrmeans=c(-2.2,2.2), xlimresdev=c(-2.2,2.2), xlimresp=c(-2,5), prop=FALSE) {
	nclass=30
	breaksGrMeans = c(-999, seq(xlimgrmeans[1], xlimgrmeans[2], length.out=nclass+1), 999)
	breaksResDev = c(-999, seq(xlimresdev[1], xlimresdev[2], length.out=nclass+1),999)
	breaksResp = c(-999, seq(xlimresp[1], xlimresp[2], length.out=nclass+1), 999)
	if(prop) lastSim$resp = lastSim$resp/20
	hist(lastSim$grmeans, breaks=breaksGrMeans, xlab="", main="", col=histcol[1], las=1, yaxt="n", xaxt="n", xlim=xlimgrmeans) 
		axis(1, at=c(-2,0,2))
		title(xlab="Simulated\ngroup means", line=3) 
	hist(lastSim$resdev, breaks=breaksResDev, xlab="", main="", col=histcol[2], las=1, yaxt="n", xaxt="n", xlim=xlimresdev)
		axis(1, at=c(-2,0,2))
		title(xlab="Simulated\nresiduals", line=3) 
	hist(lastSim$resp, breaks=breaksResp, xlab="", main="", col=histcol[3], las=1, yaxt="n", xaxt="n", xlim=xlimresp)
		axis(1, at=c(-2,0,0.5,1,2,4))
		title(xlab="Phenotypic values", line=2.5) 
}

plotRes = function(res, plotHist=TRUE, plotBLUPs=FALSE, histcol=rep("grey", 4), xlimSlopeRes=c(0,0.6), xlimGrmeansRes=c(0,2.5), xlimResidRes=c(0,1.5), ...) {
	nclass=30
	breaksSlopeRes = c(-999, seq(xlimSlopeRes[1], xlimSlopeRes[2], length.out=nclass+1), 999)
	breaksGrmeansRes = c(-999, seq(xlimGrmeansRes[1], xlimGrmeansRes[2], length.out=nclass+1),999)
	breaksResidRes = c(-999, seq(xlimResidRes[1], xlimResidRes[2], length.out=nclass+1), 999)
	if(plotHist) plotLastSimHist(res$lastSim, histcol=histcol, ...)
	if(plotBLUPs) plotBLUPs(res$lastSim, res$lastModraw, res$lastMod, histcol=histcol)
	res = res$res
	h <- hist(res$slope1, breaks=breaksSlopeRes, col="grey", xlab="", main="", las=1, yaxt="n", xaxt="n", xlim=xlimSlopeRes); 
		axis(1, at=c(0,0.2,0.4,0.6))
		title(xlab="Estimate\nslope", line=3)
		points(res$SimSlope[1], 0, pch=24, cex=2, col="red", bg="red")
		segments(x0=mean(res$slope1), x1=mean(res$slope1), y0=0, y1=max(h$counts), col="red", lwd=2); 
		segments(x0=quantile(res$slope1, c(0.025, 0.975)), x1=quantile(res$slope1, c(0.025, 0.975)), y0=0, y1=max(h$counts), col="red", lty=2, lwd=2)
	h <- hist(res$vargroups, breaks=breaksGrmeansRes, col="grey", xlab="", main="", las=1, yaxt="n", xaxt="n", xlim=xlimGrmeansRes); 
		axis(1, at=c(0,0.5,1,1.5,2))
		title(xlab="Estimate\ngroup variance", line=3)
		points(res$SimGroups[1], 0, pch=24, cex=2, col="red", bg="red")
		segments(x0=mean(res$vargroups), x1=mean(res$vargroups), y0=0, y1=max(h$counts), col="red", lwd=2); 
		segments(x0=quantile(res$vargroups, c(0.025, 0.975)), x1=quantile(res$vargroups, c(0.025, 0.975)), y0=0, y1=max(h$counts), col="red", lty=2, lwd=2)
	h <- hist(res$varres, breaks=breaksResidRes, col="grey", xlab="", main="", las=1, yaxt="n", xaxt="n", xlim=xlimResidRes); 
		axis(1, at=c(0,0.5,1,1.5))
		title(xlab="Estimate\nresidual variance", line=3)
		points(res$SimRes[1], 0, pch=24, cex=2, col="red", bg="red")
		segments(x0=mean(res$varres), x1=mean(res$varres), y0=0, y1=max(h$counts), col="red", lwd=2); 
		segments(x0=quantile(res$varres, c(0.025, 0.975)), x1=quantile(res$varres, c(0.025, 0.975)), y0=0, y1=max(h$counts), col="red", lty=2, lwd=2)
}

plotResCor = function(res, meanSEnorm, xlimSlopeSE=c(0.05,0.15), xlimSlopeEst=c(-0.2,0.8), xlimVarGroups=c(0,1.5), xlimVarRes=c(0,1)) {
	nclass=30
	breaksSlopeSE = c(-999, seq(xlimSlopeSE[1], xlimSlopeSE[2], length.out=nclass+1), 999)
	breaksSlopeEst = c(-999, seq(xlimSlopeEst[1], xlimSlopeEst[2], length.out=nclass+1), 999)
	breaksVarGroups = c(-999, seq(xlimVarGroups[1], xlimVarGroups[2], length.out=nclass+1),999)
	breaksVarRes = c(-999, seq(xlimVarRes[1], xlimVarRes[2], length.out=nclass+1), 999)
	h <- hist(res$slope1SE, breaks=breaksSlopeSE, col="grey", xlab="", main="", las=1, yaxt="n", xaxt="n", xlim=xlimSlopeSE); 
		axis(1, at=c(0.05, 0.1,0.15, 0.2))
		title(xlab="Standard error\nslope", line=3)
		points(meanSEnorm, 0, pch=24, cex=2, col="red", bg="red")
		segments(x0=mean(res$slope1SE), x1=mean(res$slope1SE), y0=0, y1=max(h$counts), col="red", lwd=2); 
		segments(x0=quantile(res$slope1SE, c(0.025, 0.975)), x1=quantile(res$slope1SE, c(0.025, 0.975)), y0=0, y1=max(h$counts), col="red", lty=2, lwd=2)
	h <- hist(res$slope1, breaks=breaksSlopeEst, col="grey", xlab="", main="", las=1, yaxt="n", xaxt="n", xlim=xlimSlopeEst); 
		axis(1, at=c(-0.2,0.2,0.6))
		title(xlab="Estimate\nslope", line=3)
		points(res$SimSlope[1], 0, pch=24, cex=2, col="red", bg="red")
		segments(x0=mean(res$slope1), x1=mean(res$slope1), y0=0, y1=max(h$counts), col="red", lwd=2); 
		segments(x0=quantile(res$slope1, c(0.025, 0.975)), x1=quantile(res$slope1, c(0.025, 0.975)), y0=0, y1=max(h$counts), col="red", lty=2, lwd=2)
	h <- hist(res$vargroups, breaks=breaksVarGroups, col="grey", xlab="", main="", las=1, yaxt="n", xaxt="n", xlim=xlimVarGroups); 
		axis(1, at=c(0,0.5,1,1.5))
		title(xlab="Estimate\ngroup variance", line=3)
		points(res$SimGroups[1], 0, pch=24, cex=2, col="red", bg="red")
		segments(x0=mean(res$vargroups), x1=mean(res$vargroups), y0=0, y1=max(h$counts), col="red", lwd=2); 
		segments(x0=quantile(res$vargroups, c(0.025, 0.975)), x1=quantile(res$vargroups, c(0.025, 0.975)), y0=0, y1=max(h$counts), col="red", lty=2, lwd=2)
	h <- hist(res$varres, breaks=breaksVarRes, col="grey", xlab="", main="", las=1, yaxt="n", xaxt="n", xlim=xlimVarRes); 
		axis(1, at=c(0,0.5,1))
		title(xlab="Estimate\nresidual variance", line=3)
		points(res$SimRes[1], 0, pch=24, cex=2, col="red", bg="red")
		segments(x0=mean(res$varres), x1=mean(res$varres), y0=0, y1=max(h$counts), col="red", lwd=2); 
		segments(x0=quantile(res$varres, c(0.025, 0.975)), x1=quantile(res$varres, c(0.025, 0.975)), y0=0, y1=max(h$counts), col="red", lty=2, lwd=2)
}


tableRes = function(res, SEfactor = 1.98) {
		SimSlope = res$SimSlope[1]
		vargroupsSim = res$SimGroups[1]
		varresSim =res$SimRes[1]
		slope1SE = res$slope1SE
	print("hallo")
		resTabRow = c(nrep=nrep, 
				   SlopeSim= SimSlope,
				   SlopeMean=mean(res$slope1), 
				   SlopeSD=sd(res$slope1), 
				   SlopeBias=mean(res$slope1-SimSlope), 
				   SlopeRMSE=sqrt(mean((res$slope1-SimSlope)^2)), 
				   SlopeRMSEsd=sd(sqrt((res$slope1-SimSlope)^2)), 
				   vargroupsSim=vargroupsSim, 
				   vargroupsMean=mean(res$vargroups), 
				   vargroupsSD=sd(res$vargroups), 
				   vargroupsBias=mean(res$vargroups-vargroupsSim), 
				   vargroupsRMSE=sqrt(mean((res$vargroups-res$SimGroups[1])^2)), 
				   vargroupsRMSEsd=sd(sqrt((res$vargroups-res$SimGroups[1])^2)), 
				   varresSim = varresSim, 
				   varresMean=mean(res$varres), 
				   varresSD=sd(res$varres), 
				   varresBias=mean(res$varres-varresSim), 
				   varresRMSE=sqrt(mean((res$varres-varresSim)^2)),
				   varresRMSEsd=sd(sqrt((res$varres-varresSim)^2)),
				   Slope1SEmean = mean(slope1SE),
				   Slope1SEsd = sd(slope1SE),
				   SlopeCICoverage = mean(res$slope1-SEfactor*slope1SE<=SimSlope & res$slope1+SEfactor*slope1SE>=SimSlope))  
	print("hallo")
		return(resTabRow)
}

##########################################################################################
## Simulation code  ######################################################################
##########################################################################################

biasMSEtable = data.frame(Figure=NA, Type=rep(NA,4*5), nrep=nrep, 
					SlopeSim=-9999, SlopeMean=-9999, SlopeSD=-9999, SlopeBias=-9999, SlopeRMSE=-9999, SlopeRMSEsd=-9999, 
					vargroupsSim=-9999, vargroupsMean=-9999, vargroupsSD=-9999, vargroupsBias=-9999, vargroupsRMSE=-9999, vargroupsRMSEsd=-9999,  
					varresSim =-9999, varresMean =-9999, varresSD =-9999, varresBias=-9999, varresRMSE=-9999, varresRMSEsd=-9999,
					Slope1SEmean=-9999, Slope1SEsd=-9999, SlopeCICoverage=-9999)
				   
#######################################
# Skew: Figure 1  #####################
#######################################

if(plotYN) windows(14,10)
par(mfrow=c(4,6), mar=c(5,2,1,1)) 
if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=0,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
	save(res, file=paste0(path, "Figure1_Row1_Reference.RData"))}
load(file=paste0(path, "Figure1_Row1_Reference.RData"))
biasMSEtable[1,] = c("Figure1", "Control", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey", "grey", "grey"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 3, skewGroups = 0, kurtRes = 20, kurtGroups = 3, xCor=0, heteroscedfactor=0,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure1_Row2_ResSkew.RData"))}
load(file=paste0(path, "Figure1_Row2_ResSkew.RData"))
biasMSEtable[2,] = c("Figure1", "ResSkew", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey", "grey40", "grey40"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 3, kurtRes = 3, kurtGroups = 20, xCor=0, heteroscedfactor=0,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure1_Row3_GroupsSkew.RData")) }
load(file=paste0(path, "Figure1_Row3_GroupsSkew.RData")) 
biasMSEtable[3,] = c("Figure1", "GroupsSkew", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey40", "grey", "grey40"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 3, skewGroups = 3, kurtRes = 20, kurtGroups = 20, xCor=0, heteroscedfactor=0,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure1_Row4_ResGroupsSkew.RData")) }
load(file=paste0(path, "Figure1_Row4_ResGroupsSkew.RData"))	
biasMSEtable[4,] = c("Figure1", "ResGroupsSkew", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey40", "grey40", "grey40"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

########################################################## 
# Bimodal: Figure S1  ####################################
########################################################## 

if(plotYN) windows(14,10)
par(mfrow=c(4,6), mar=c(5,2,1,1)) 
if(sim) {
	res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=0,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure2_Row1_Reference.RData")) }
load(file=paste0(path, "Figure2_Row1_Reference.RData"))
biasMSEtable[4+1,] = c("Figure2", "Control", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey", "grey", "grey"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=0,
		datlevelpred=c(-1.5,1.5), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure2_Row2_ResBimod.RData")) }
load(file=paste0(path, "Figure2_Row2_ResBimod.RData"))
biasMSEtable[4+2,] = c("Figure2", "ResBimod", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey", "grey40", "grey40"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=0,
		datlevelpred=c(0,0), grouplevelpred=c(-1.5,1.5), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure2_Row3_GroupsBimod.RData")) }
load(file=paste0(path, "Figure2_Row3_GroupsBimod.RData"))
biasMSEtable[4+3,] = c("Figure2", "GroupsBimod", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey40", "grey", "grey40"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=0,
		datlevelpred=c(-1.5,1.5), grouplevelpred=c(-1.5,1.5), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily) 
save(res, file=paste0(path, "Figure2_Row4_ResGroupsBimod.RData")) }
load(file=paste0(path, "Figure2_Row4_ResGroupsBimod.RData"))
biasMSEtable[4+4,] = c("Figure2", "ResGroupsBimod", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey40", "grey40", "grey40"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

##################################################################### 
# Heteroskedasticity: Figure S2  ####################################
##################################################################### 

if(plotYN) windows(14,10)
par(mfrow=c(4,6), mar=c(5,2,1,1))
if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=0,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure3_Row1_Reference.RData")) }
load(file=paste0(path, "Figure3_Row1_Reference.RData"))
biasMSEtable[8+1,] = c("Figure3", "Control", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey", "grey", "grey"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
			skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=2,
			datlevelpred=c(0,0), grouplevelpred=c(0,0), 
			nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
			varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure3_Row2_lowHet.RData")) }
load(file=paste0(path, "Figure3_Row2_lowHet.RData"))
biasMSEtable[8+2,] = c("Figure3", "lowHet", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey", "grey40", "grey40"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=4,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure3_Row3_modHet.RData")) }
load(file=paste0(path, "Figure3_Row3_modHet.RData"))
biasMSEtable[8+3,] = c("Figure3", "modHet", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey", "grey40", "grey40"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=8,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure3_Row4_hiHet.RData")) }
load(file=paste0(path, "Figure3_Row4_hiHet.RData"))
biasMSEtable[8+4,] = c("Figure3", "hiHet", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, histcol=c("grey", "grey40", "grey40"), xlimresp=xlimresp, prop=GLMfamily=="Proportion")

######################################################################## 
# Missing random effects: Figure 3  ####################################
######################################################################## 

if(plotYN) windows(14,10)
par(mfrow=c(4,6), mar=c(5,2,1,1)) 
if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=0,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure4_Row1_Reference.RData")) }
load(file=paste0(path, "Figure4_Row1_Reference.RData"))
biasMSEtable[12+1,] = c("Figure4", "Control", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=0,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0.5, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure4_Row2_VarHiMissing.RData")) }
load(file=paste0(path, "Figure4_Row2_VarHiMissing.RData"))
biasMSEtable[12+2,] = c("Figure4", "VarHiMissing", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=0,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0.5, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure4_Row3_VarLoMissing.RData")) }
load(file=paste0(path, "Figure4_Row3_VarLoMissing.RData"))
biasMSEtable[12+3,] = c("Figure4", "VarLoMissing", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, xlimresp=xlimresp, prop=GLMfamily=="Proportion")

if(sim) {
res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
		skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=0,
		datlevelpred=c(0,0), grouplevelpred=c(0,0), 
		nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
		varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0.5, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
save(res, file=paste0(path, "Figure4_Row4_VarCrossMissing.RData")) }
load(file=paste0(path, "Figure4_Row4_VarCrossMissing.RData"))
biasMSEtable[12+4,] = c("Figure4", "VarCrossMissing", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotRes(res, xlimresp=xlimresp, prop=GLMfamily=="Proportion")

####################################################################### 
# Correlated predictors: Figure S3  ###################################
####################################################################### 

if(plotYN) windows(14/6*4,10)
par(mfrow=c(4,4), mar=c(5,2,1,1)) 
if(sim) {
	res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
			skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0, heteroscedfactor=0,
			datlevelpred=c(0,0), grouplevelpred=c(0,0), 
			nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
			varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
	save(res, file=paste0(path, "Figure5_Row1_Reference.RData"))
}
load(file=paste0(path, "Figure5_Row1_Reference.RData"))
meanSEnorm = mean(res$res$slope1SE)
biasMSEtable[16+1,] = c("Figure5", "Control", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotResCor(res$res, meanSEnorm)

if(sim) {
	res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
			skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0.2, heteroscedfactor=0,
			datlevelpred=c(0,0), grouplevelpred=c(0,0),  
			nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
			varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
	save(res, file=paste0(path, "Figure5_Row2_LowCorr.RData")) }
load(file=paste0(path, "Figure5_Row1_Reference.RData"))
meanSEnorm = mean(res$res$slope1SE)
load(file=paste0(path, "Figure5_Row2_LowCorr.RData"))
biasMSEtable[16+2,] = c("Figure5", "LowCorr", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotResCor(res$res, meanSEnorm)

if(sim) {
	res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
			skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0.5, heteroscedfactor=0,
			datlevelpred=c(0,0), grouplevelpred=c(0,0), 
			nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
			varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
	save(res, file=paste0(path, "Figure5_Row3_ModCorr.RData")) }
load(file=paste0(path, "Figure5_Row1_Reference.RData"))
meanSEnorm = mean(res$res$slope1SE)
load(file=paste0(path, "Figure5_Row3_ModCorr.RData"))
biasMSEtable[16+3,] = c("Figure5", "ModCorr", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotResCor(res$res, meanSEnorm)

if(sim) {
	res = datsim(vargroups = vargroups, varres = varres, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = nrep,
			skewRes = 0, skewGroups = 0, kurtRes = 3, kurtGroups = 3, xCor=0.8, heteroscedfactor=0,
			datlevelpred=c(0,0), grouplevelpred=c(0,0), 
			nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
			varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
	save(res, file=paste0(path, "Figure5_Row4_HiCorr.RData")) }
load(file=paste0(path, "Figure5_Row1_Reference.RData"))
meanSEnorm = mean(res$res$slope1SE)
load(file=paste0(path, "Figure5_Row4_HiCorr.RData"))
biasMSEtable[16+4,] = c("Figure5", "HiCorr", tableRes(res$res, SEfactor=SEfactor))
if(plotYN) plotResCor(res$res, meanSEnorm)

##################################################### 
# Bias and MSE plots: Figure 4 and S5-S10 ###########
##################################################### 

for(i in 3:24) biasMSEtable[,i] = as.numeric(biasMSEtable[,i])
save(biasMSEtable, file=paste0(path, "Table1.RData"))

basepath = "C:\\Users\\Schielzeth\\Documents\\A_PaperDrafts\\D_CurrentProjects\\_Holger\\HS_LMMRobustness\\Simulation\\"

if(set == "FewGroups") {
path = paste0(basepath, "FewGroups\\")
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
}

if(set == "Basic") {
path = paste0(basepath, "Basic\\")
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
}

if(set == "SmallSample") {
path = paste0(basepath, "SmallSample\\")
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
}

if(set == "LowRep") {
path = paste0(basepath, "LowRep\\")
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
}

if(set == "HighRep") {
path = paste0(basepath, "HighRep\\")
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
}

if(set == "Balanced") {
path = paste0(basepath, "Balanced\\")
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.02,0.18)
ylimBias3 = c(-0.02,0.05)
ylimRMSE = c(0.9, 1.3, 0.7)
}

if(set == "Poisson") {
path = paste0(basepath, "Poisson\\")
ylimBias1 = c(-0.10,0.10)
ylimBias2 = c(-0.3,0.40)
ylimBias3 = c(-0.2,0.3)
ylimRMSE = c(1.2, 1.5, 1.2)
}

if(set == "Proportion") {
path = paste0(basepath, "Proportion\\")
ylimBias1 = c(-0.02,0.02)
ylimBias2 = c(-0.3,0.30)
ylimBias3 = c(-0.3,0.30)
ylimRMSE = c(0.9, 1.3, 0.7)
}

load(paste0(path, "Table1.RData"))

biasMSEtable$slopeBiasPerc = biasMSEtable$SlopeBias / biasMSEtable$SlopeSim *100
biasMSEtable$slopeRMSEPerc = biasMSEtable$SlopeRMSE / biasMSEtable$SlopeSim *100
biasMSEtable$slopeSEPerc = (biasMSEtable$SlopeSD / sqrt(biasMSEtable$nrep)) / biasMSEtable$SlopeSim *100
biasMSEtable$slopeRMSEPercSE = (biasMSEtable$SlopeRMSEsd / sqrt(biasMSEtable$nrep)) / biasMSEtable$SlopeSim *100
biasMSEtable$groupvarBiasPerc = biasMSEtable$vargroupsBias / biasMSEtable$vargroupsSim *100
biasMSEtable$groupvarRMSEPerc = biasMSEtable$vargroupsRMSE / biasMSEtable$vargroupsSim *100
biasMSEtable$groupvarSEPerc = (biasMSEtable$vargroupsSD / sqrt(biasMSEtable$nrep)) / biasMSEtable$vargroupsSim *100
biasMSEtable$groupvarRMSEPercSE = (biasMSEtable$vargroupsRMSEsd / sqrt(biasMSEtable$nrep)) / biasMSEtable$vargroupsSim *100
biasMSEtable$resvarBiasPerc = biasMSEtable$varresBias / biasMSEtable$varresSim *100
biasMSEtable$resvarRMSEPerc = biasMSEtable$varresRMSE / biasMSEtable$varresSim *100
biasMSEtable$resvarSEPerc = (biasMSEtable$varresSD / sqrt(biasMSEtable$nrep)) / biasMSEtable$varresSim *100
biasMSEtable$resvarRMSEPercSE = (biasMSEtable$varresRMSEsd / sqrt(biasMSEtable$nrep)) / biasMSEtable$varresSim *100

windows(7, 7)
par(mfrow=c(3,2), mar=c(3.5,4.5,2,1))
plot(0, ylim=ylimBias1*100, xlim=c(0.5,20.5), las=1, xaxt="n", xlab="", ylab="Bias [%]", main="Bias in slopes", type="n")
	abline(v=seq(0,20,by=4)+0.5, h=0, lwd=c(1,2,1,1,2,2,2))
	axis(1,at=1:20, labels=paste0(rep(LETTERS[1:5],each=4), rep(1:4,5)), las=2)
	rect(xleft=1:20-0.3, xright=1:20+0.3, ybottom=0, ytop=biasMSEtable$slopeBiasPerc, col=rep(c("white", "grey70", "grey50", "grey30"),4))
	segments(x0=1:20, x1=1:20, y0=biasMSEtable$slopeBiasPerc-biasMSEtable$slopeSEPerc, y1=biasMSEtable$slopeBiasPerc+biasMSEtable$slopeSEPerc, lwd=1)
	text(1:20-0.3, ylimBias1[2]*100, ifelse(biasMSEtable$slopeBiasPerc>ylimBias1[2]*100,round(biasMSEtable$slopeBiasPerc,0),""), pos=1, col="white",srt=90, cex=0.8)
plot(0, ylim=c(0,ylimRMSE[1])*100, xlim=c(0.5,20.5), las=1, xaxt="n", xlab="", ylab="RMSE [%]", main="Prediction error", type="n")
	abline(v=seq(0,20,by=4)+0.5, h=0, lwd=c(1,2,1,1,2,2,2))
	axis(1,at=1:20, labels=paste0(rep(LETTERS[1:5],each=4), rep(1:4,5)), las=2)
	rect(xleft=1:20-0.3, xright=1:20+0.3, ybottom=0, ytop=biasMSEtable$slopeRMSEPerc, col=rep(c("white", "grey70", "grey50", "grey30"),4))
	segments(x0=1:20, x1=1:20, y0=biasMSEtable$slopeRMSEPerc-biasMSEtable$slopeRMSEPercSE, y1=biasMSEtable$slopeRMSEPerc+biasMSEtable$slopeRMSEPercSE, lwd=1)
	text(1:20-0.3, ylimRMSE[1]*100, ifelse(biasMSEtable$slopeRMSEPerc>ylimRMSE[1]*100,round(biasMSEtable$slopeRMSEPerc,0),""), pos=1, col="white",srt=90, cex=0.8)
plot(0, ylim=ylimBias2*100, xlim=c(0.5,20.5), las=1, xaxt="n", xlab="", ylab="Bias [%]", main="Bias in group variance", type="n")
	abline(v=seq(0,20,by=4)+0.5, h=0, lwd=c(1,2,1,1,2,2,2))
	axis(1,at=1:20, labels=paste0(rep(LETTERS[1:5],each=4), rep(1:4,5)), las=2)
	rect(xleft=1:20-0.3, xright=1:20+0.3, ybottom=0, ytop=biasMSEtable$groupvarBiasPerc, col=rep(c("white", "grey70", "grey50", "grey30"),4))
	segments(x0=1:20, x1=1:20, y0=biasMSEtable$groupvarBiasPerc-biasMSEtable$groupvarSEPerc, y1=biasMSEtable$groupvarBiasPerc+biasMSEtable$groupvarSEPerc, lwd=1)
	text(1:20-0.3, ylimBias2[2]*100, ifelse(biasMSEtable$groupvarBiasPerc>ylimBias2[2]*100,round(biasMSEtable$groupvarBiasPerc,0),""), pos=1, col="white",srt=90, cex=0.8)
plot(0, ylim=c(0,ylimRMSE[2])*100, xlim=c(0.5,20.5), las=1, xaxt="n", xlab="", ylab="RMSE [%]", main="Prediction error", type="n")
	abline(v=seq(0,20,by=4)+0.5, h=0, lwd=c(1,2,1,1,2,2,2))
	axis(1,at=1:20, labels=paste0(rep(LETTERS[1:5],each=4), rep(1:4,5)), las=2)
	rect(xleft=1:20-0.3, xright=1:20+0.3, ybottom=0, ytop=biasMSEtable$groupvarRMSEPerc, col=rep(c("white", "grey70", "grey50", "grey30"),4))
	segments(x0=1:20, x1=1:20, y0=biasMSEtable$groupvarRMSEPerc-biasMSEtable$groupvarRMSEPercSE, y1=biasMSEtable$groupvarRMSEPerc+biasMSEtable$groupvarRMSEPercSE, lwd=1)
	text(1:20-0.3, ylimRMSE[2]*100, ifelse(biasMSEtable$groupvarRMSEPerc>ylimRMSE[2]*100, round(biasMSEtable$groupvarRMSEPerc,0), ""), pos=1, col="white",srt=90, cex=0.8)
plot(0, ylim=ylimBias3*100, xlim=c(0.5,20.5), las=1, xaxt="n", xlab="", ylab="Bias [%]", main="Bias in error variance", type="n")
	abline(v=seq(0,20,by=4)+0.5, h=0, lwd=c(1,2,1,1,2,2,2))
	axis(1,at=1:20, labels=paste0(rep(LETTERS[1:5],each=4), rep(1:4,5)), las=2)
	rect(xleft=1:20-0.3, xright=1:20+0.3, ybottom=0, ytop=biasMSEtable$resvarBiasPerc, col=rep(c("white", "grey70", "grey50", "grey30"),4))
	segments(x0=1:20, x1=1:20, y0=biasMSEtable$resvarBiasPerc-biasMSEtable$resvarSEPerc, y1=biasMSEtable$resvarBiasPerc+biasMSEtable$resvarSEPerc, lwd=1)
	text(1:20-0.3, ylimBias3[2]*100, ifelse(biasMSEtable$resvarBiasPerc>ylimBias3[2]*100,round(biasMSEtable$resvarBiasPerc,0),""), pos=1, col="white",srt=90, cex=0.8)
plot(0, ylim=c(0,ylimRMSE[3])*100, xlim=c(0.5,20.5), las=1, xaxt="n", xlab="", ylab="RMSE [%]", main="Prediction error", type="n")
	abline(v=seq(0,20,by=4)+0.5, h=0, lwd=c(1,2,1,1,2,2,2))
	axis(1,at=1:20, labels=paste0(rep(LETTERS[1:5],each=4), rep(1:4,5)), las=2)
	rect(xleft=1:20-0.3, xright=1:20+0.3, ybottom=0, ytop=biasMSEtable$resvarRMSEPerc, col=rep(c("white", "grey70", "grey50", "grey30"),4))
	segments(x0=1:20, x1=1:20, y0=biasMSEtable$resvarRMSEPerc-biasMSEtable$resvarRMSEPercSE, y1=biasMSEtable$resvarRMSEPerc+biasMSEtable$resvarRMSEPercSE, lwd=1)
	text(1:20-0.3, ylimRMSE[3]*100, ifelse(biasMSEtable$resvarRMSEPerc>ylimRMSE[3]*100,round(biasMSEtable$resvarRMSEPerc,0),""), pos=1, col="white",srt=90, cex=0.8)


}

############################################# 
# SE plot: Figures S10 and S11 ##############
############################################# 

plotResSE = function(res, meanSEnorm, ylimSlopeSE=c(0.05,0.15), ylimRMSE=c(0,1), main="") {
	plot(0, ylim=ylimSlopeSE, xlim=c(0.5,20.5), las=1, xaxt="n", xlab="", ylab="Slope standard error", main=main, type="n")
		abline(v=seq(0,20,by=4)+0.5, h=0, lwd=c(1,2,1,1,2,2,2))
		axis(1,at=1:20, labels=paste0(rep(LETTERS[1:5],each=4), rep(1:4,5)), las=2)
		rect(xleft=1:20-0.3, xright=1:20+0.3, ybottom=0, ytop=biasMSEtable$Slope1SEmean, col=rep(c("white", "grey70", "grey50", "grey30"),4))
		segments(x0=1:20, x1=1:20, y0=biasMSEtable$Slope1SEmean-biasMSEtable$Slope1SEsd, y1=biasMSEtable$Slope1SEmean+biasMSEtable$Slope1SEsd, lwd=1)
		#text(1:20-0.3, ylimSlopeSE[2], ifelse(biasMSEtable$Slope1SEmean>ylimSlopeSE[2], round(biasMSEtable$Slope1SEmean,0),""), pos=1, col="white",srt=90, cex=0.8)
		text(1:20-0.3, ylimSlopeSE[2], format(round(1-biasMSEtable$SlopeCICoverage, 3), nsmall=3), pos=1, col="black",srt=90, cex=0.8)
	plot(0, ylim=ylimRMSE, xlim=c(0.5,20.5), las=1, xaxt="n", xlab="", ylab="RMSE (abolute)", main=main, type="n")
		abline(v=seq(0,20,by=4)+0.5, h=0, lwd=c(1,2,1,1,2,2,2))
		axis(1,at=1:20, labels=paste0(rep(LETTERS[1:5],each=4), rep(1:4,5)), las=2)
		rect(xleft=1:20-0.3, xright=1:20+0.3, ybottom=0, ytop=biasMSEtable$SlopeRMSE, col=rep(c("white", "grey70", "grey50", "grey30"),4))
		segments(x0=1:20, x1=1:20, y0=biasMSEtable$SlopeRMSE-biasMSEtable$SlopeRMSEsd, y1=biasMSEtable$SlopeRMSE+biasMSEtable$SlopeRMSEsd, lwd=1)
		text(1:20-0.3, ylimRMSE[2], ifelse(biasMSEtable$SlopeRMSE>ylimRMSE[2], round(biasMSEtable$SlopeRMSE,2),""), pos=1, col="white",srt=90, cex=0.8)
}

set = c("Basic", "Balanced", "SmallSample", "FewGroups")
windows(7,10)
par(mfrow=c(4,2), mar=c(3,4,2,1)) 
for(i in 1:length(set)) {
	if(set[i] == "Basic") path = paste0(basepath, "Basic\\")
	if(set[i] == "Balanced") path = paste0(basepath, "Balanced\\")
	if(set[i] == "SmallSample") path = paste0(basepath, "SmallSample\\")
	if(set[i] == "FewGroups") path = paste0(basepath, "FewGroups\\")
	load(paste0(path, "Table1.RData"))
	plotResSE(biasMSEtable, ylimSlopeSE=c(0,0.3), ylimRMSE=c(0,0.3), main=set[i])
}
set = c("LowRep", "HighRep", "Poisson", "Proportion")
windows(7,10)
par(mfrow=c(4,2), mar=c(3,4,2,1)) 
for(i in 1:length(set)) {
	if(set[i] == "LowRep") path = paste0(basepath, "LowRep\\")
	if(set[i] == "HighRep") path = paste0(basepath, "HighRep\\")
	if(set[i] == "Poisson") path = paste0(basepath, "Poisson\\")
	if(set[i] == "Proportion") path = paste0(basepath, "Proportion\\")
	load(paste0(path, "Table1.RData"))
	plotResSE(biasMSEtable, ylimSlopeSE=c(0,0.3), ylimRMSE=c(0,0.3), main=set[i])
}

############################################# 
# BLUP plot: Figure 2  ######################
############################################# 

plotBLUPs = function(res, histcol=rep("grey", 4), xlimGrmeansSim=c(-1,2.5), xlimGrmeansBLUPS=c(-1,2), xlimResdevSim=c(-1,3), xlimResdevResid=c(-2,3)) {
	nclass=30
	breaksGrmeansSim = c(-999, seq(xlimGrmeansSim[1], xlimGrmeansSim[2], length.out=nclass+1), 999)
	breaksGrmeansBLUPS = c(-999, seq(xlimGrmeansBLUPS[1], xlimGrmeansBLUPS[2], length.out=nclass+1),999)
	breaksResdevSim = c(-999, seq(xlimResdevSim[1], xlimResdevSim[2], length.out=nclass+1), 999)	
	breaksResdevResid = c(-999, seq(xlimResdevResid[1], xlimResdevResid[2], length.out=nclass+1), 999)	
	hist(res$lastSim$grmeans, breaks=breaksGrmeansSim, xlab="", main="", col=histcol[1], las=1, yaxt="n", xaxt="n", xlim=xlimGrmeansSim) 
		axis(1, at=c(-2,0,2))
		title(xlab="Simulated\ngroup means", line=3) 
	hist(ranef(res$lastModraw)$'groups'[,1], breaks=breaksGrmeansBLUPS, xlab="", main="", col=histcol[2], las=1, yaxt="n", xaxt="n", xlim=xlimGrmeansBLUPS)
		axis(1, at=c(-2,0,2))
		title(xlab="Estimated\ngroup means", line=3) 
	hist(res$lastSim$resdev, breaks=breaksResdevSim, xlab="", main="", col=histcol[3], las=1, yaxt="n", xaxt="n", xlim=xlimResdevSim)
		axis(1, at=c(-2,0,2))
		title(xlab="Simulated\nresiduals", line=3) 
	hist(res$lastMod$residuals, breaks=breaksResdevResid, xlab="", main="", col=histcol[4], las=1, yaxt="n", xaxt="n", xlim=xlimResdevResid)
		axis(1, at=c(-2,0,2))
		title(xlab="Estimated\nresiduals", line=3) 
}

windows(10/4*3,14/6*4)
par(mar=c(5,2,1,1))
layout(matrix(1:12, ncol=3, byrow=FALSE)) 
res = datsim(vargroups = 0.1, varres = 0.9, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = 1,
			skewRes = 3, skewGroups = 3, kurtRes = 20, kurtGroups = 20, xCor=0, heteroscedfactor=0,
			datlevelpred=c(0,0), grouplevelpred=c(0,0), 
			nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
			varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
	plotBLUPs(res)
res = datsim(vargroups = 0.5, varres = 0.5, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = 1,
			skewRes = 3, skewGroups = 3, kurtRes = 20, kurtGroups = 20, xCor=0, heteroscedfactor=0,
			datlevelpred=c(0,0), grouplevelpred=c(0,0), 
			nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
			varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
	plotBLUPs(res)
res = datsim(vargroups = 0.9, varres = 0.1, slope1 = 0.2, slope2 = -0.2, ngroups = ngroups, nobs = nobs, intercept = intercept, nrep = 1,
			skewRes = 3, skewGroups = 3, kurtRes = 20, kurtGroups = 20, xCor=0, heteroscedfactor=0,
			datlevelpred=c(0,0), grouplevelpred=c(0,0), 
			nhilevelgroups = nhilevelgroups, nlolevelgroups = nlolevelgroups, ncrossedgroups = ncrossedgroups,
			varhilevelgroups = 0, varlolevelgroups = 0, varcrossedgroups = 0, unbalancedYN=unbalancedYN, GLMfamily=GLMfamily)
	plotBLUPs(res)

