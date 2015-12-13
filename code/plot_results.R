## Plot results from stochastic dynamic optimization and forward simulations for 
## "Adaptive information use during growth explains long-term effects of early-life experiences" 
## English et al., American Naturalist (in press)

## last updated 11 December 2015

## contact: sineadenglish@cantab.net if any questions or suggestions


library(fields)
library(plyr)

## EXPLAIN THESE TWO FILES ###
source("~/Documents/work/modelling/info-use-dev/github/R/visualise_sims.r")

# FUNCTIONS

# need to make new function to return NA if iqr na
iqr1 = function(x)
{
	y = try(quantile(x)[2])
	return(ifelse(class(y)=="try-error",NA,y))
}
iqr2 = function(x)
{
	y = try(quantile(x)[4])
	return(ifelse(class(y)=="try-error",NA,y))
}

# find fitness at maturity in good environment
RV_mat = function(size, env, parno)
{
	phi = pars[[parno]][1]
	const = c(pars[[parno]][4],pars[[parno]][5])
	kappa = c(pars[[parno]][2],pars[[parno]][3])
	
	return (1/(const[env+1] + exp((-1.0*phi)*(size -kappa[env+1]))));      
}    


setwd("~/Documents/work/modelling/info-use-dev/github/data/") 

# read in data
pars = list()
fitList = list()
optactList = list()
infoList = list()

for(parno in 1:3) {
	pars[[parno]] <- as.numeric(read.table(file=paste("Params",parno,".txt",sep="")))
	fitList[[parno]] <- read.table(paste("FitnessValue",parno,".txt",sep=""))
	optactList[[parno]] <- read.table(paste("Decision",parno,".txt",sep="")); 
	infoList[[parno]] <- read.table(paste("InfoVal",parno,".txt",sep=""))
	}
	

checkpars = 1:3
# create lists to store simulation results

# name: "simdat_[res]_[parno]_[mort]_[prior]_[pfood]_[age]" 

# res = c("pop","ind") ## whether plotting population means or individual trajectories
# parno = 1:3 ## pred-specific mortality or not
# mort = c(0,1) ## sort of redundant for res: population means include mortality, trajectories don't
# prior = c(0.5, 0.1, 0.9)
# pfood = c(0.0, 0.1, 0.9) # experimental change in p(food discovery) [pfood = 0.0 is baseline]
# age = c(0,10,20,30,40,50) # experimental age at which start consistent experience (fixed to 0 in baseline)


for(res in c("pop","ind"))
{
	for(mort in c(0,1))
	{
		# baseline simulations & manipulate priors
		pfood = "0.0"
		age = 0
	
		for(prior in c(0.5,0.1,0.9)) 
		{
			for (parno in checkpars) 
			{
				assign(paste("simdat","pop",parno,mort,prior,pfood,age,sep="_"),read.table(file=paste("SimOutputAll",parno,mort,prior,pfood,age,"txt",sep="."),header=TRUE))
				assign(paste("simdat","ind",parno,mort,prior,pfood,age,sep="_"),read.table(file=paste("SimOutputInd",parno,mort,prior,pfood,age,"txt",sep="."),header=TRUE))
			}
		}
		
		# manipulate experiences at stages of development
		prior = 0.5
		for(parno in checkpars)
		{
			for(pfood in c(0.1,0.9))
			{
				for(age in seq(0,50,10)) 
				{
					assign(paste("simdat","pop",parno,mort,prior,pfood,age,sep="_"),read.table(file=paste("SimOutputAll",parno,mort,prior,pfood,age,"txt",sep="."),header=TRUE))
					assign(paste("simdat","ind",parno,mort,prior,pfood,age,sep="_"),read.table(file=paste("SimOutputInd",parno,mort,prior,pfood,age,"txt",sep="."),header=TRUE))
				}
			}
		}
	}
}

# summary data for plots of experiences
allfit <- c()

# population means: 
res = "pop"
mort = 1

# baseline simulations & manipulate priors
pfood = "0.0" # numeric not factor
age = 0

for(prior in c(0.1,0.5,0.9)) 
{
	for (parno in checkpars) 
	{
		thisdat = get(paste("simdat",res,parno,mort,prior,pfood,age,sep="_"))
				
		thisdat$fit <- RV_mat(thisdat$mass, thisdat$env, parno=parno)
		thisdat$fit[thisdat$state==0] <- 0
		n<-length(thisdat$fit[thisdat$env==0])
			
		allfit_i <- c(parno,prior,0,age,
			mean(thisdat$fit[thisdat$env==0]),
			mean(thisdat$fit[thisdat$env==1]),
			length(thisdat$state[thisdat$state==2&thisdat$env==0]),
			length(thisdat$state[thisdat$state==2&thisdat$env==1]),
			mean(thisdat$mass[thisdat$state==2&thisdat$env==0]),
			mean(thisdat$mass[thisdat$state==2&thisdat$env==1]),
			mean(thisdat$age[thisdat$state==2&thisdat$env==0]),
			mean(thisdat$age[thisdat$state==2&thisdat$env==1]),
			sd(thisdat$fit[thisdat$env==0])/sqrt(length(thisdat$state[thisdat$state==2&thisdat$env==0])),
			sd(thisdat$fit[thisdat$env==1])/sqrt(length(thisdat$state[thisdat$state==2&thisdat$env==1])),
			sd(thisdat$mass[thisdat$state==2&thisdat$env==0])/sqrt(n),
			sd(thisdat$mass[thisdat$state==2&thisdat$env==1])/sqrt(n),
			sd(thisdat$age[thisdat$state==2&thisdat$env==0])/sqrt(n),
			sd(thisdat$age[thisdat$state==2&thisdat$env==1])/sqrt(n),
			n)
		allfit <- rbind(allfit, allfit_i)
		
	}
}

# manipulate experiences at stages of development
prior = 0.5
for(parno in checkpars)
{
	for(pfood in c(0.1,0.9))
	{
		for(age in seq(0,50,10)) 
		{
			thisdat = get(paste("simdat",res,parno,mort,prior,pfood,age,sep="_"))

			thisdat$fit <- RV_mat(thisdat$mass, thisdat$env, parno=parno)
			thisdat$fit[thisdat$state==0] <- 0
			n<-length(thisdat$fit[thisdat$env==0])
			
				
			allfit_i <- c(parno,prior, pfood,age,
				mean(thisdat$fit[thisdat$env==0]),
				mean(thisdat$fit[thisdat$env==1]),
				length(thisdat$state[thisdat$state==2&thisdat$env==0]),
				length(thisdat$state[thisdat$state==2&thisdat$env==1]),
				mean(thisdat$mass[thisdat$state==2&thisdat$env==0]),
				mean(thisdat$mass[thisdat$state==2&thisdat$env==1]),
				mean(thisdat$age[thisdat$state==2&thisdat$env==0]),
				mean(thisdat$age[thisdat$state==2&thisdat$env==1]),
				sd(thisdat$fit[thisdat$env==0])/sqrt(length(thisdat$state[thisdat$state==2&thisdat$env==0])),
				sd(thisdat$fit[thisdat$env==1])/sqrt(length(thisdat$state[thisdat$state==2&thisdat$env==1])),
				sd(thisdat$mass[thisdat$state==2&thisdat$env==0])/sqrt(n),
				sd(thisdat$mass[thisdat$state==2&thisdat$env==1])/sqrt(n),
				sd(thisdat$age[thisdat$state==2&thisdat$env==0])/sqrt(n),
				sd(thisdat$age[thisdat$state==2&thisdat$env==1])/sqrt(n),
				n)
			allfit <- rbind(allfit, allfit_i)
		}
	}
}


allfit = data.frame(allfit)

colnames(allfit) = c("parno","prior","pfood","age","fit1","fit2","surv1","surv2","size1","size2","age1","age2","sdfit1","sdfit2","sdsize1","sdsize2","sdage1","sdage2","n")

allfit$propsurv1 = allfit$surv1/allfit$n
allfit$propsurv2 = allfit$surv2/allfit$n


# add unique id for each simulation type
allfit$simtype = paste(allfit$parno, allfit$prior, allfit$pfood, allfit$age, sep="-")
allfit$simno = as.numeric(as.factor(allfit$simtype))


meandat <- allfit  # because of old code

baseline_simno = meandat$simno[meandat$prior==0.5 & meandat$pfood==0.0 & meandat$age==0]



## DO PLOTS

setwd("~/Documents/work/modelling/info-use-dev/plots/")

bg.col = grey(0.75)
try.col = rev(heat.colors(501)) # gray.colors(20, start=0.3, end=0.9)

checkpars = 1:3
col.good.shade = rgb(0,0,1,0.65)
col.bad.shade = rgb(0.7,0.3,0,0.3)
col.good.line  = "darkblue"
col.bad.line = "orange"

plot.cex = 1.2
xlim.age = c(0,75) # 2 units more than earliest age at maturity in env0 (73)
xlim.age2 = c(0,65) # for mean change in foraging effort: to earliest age at maturity across both environments (if average across both, min=65)
ylim.info = c(0,1.5*1e-6) #c(0,1.5*max(median2))
ylim.freqSize = c(0,2000)
ylim.freqAge = c(0,310)
xlim.freqAge = c(60,141)
xlim.freqSize = c(30,48)
ylim.alpha = c(0,0.025) #c(0,max(plotdat1_summ$forage.mean))	


# plot all optact then val info, mean change behav


for(parno in 2) 
# sort out data
	# calculate change in information value with age
	simdat1 <- get(paste("simdat","ind",parno,0,0.5,"0.0",0,sep="_"))
	simdat1$sim_it <- paste(simdat1$simno,simdat1$it,simdat1$env,sep="-")
	summsim<-ddply(simdat1,.(age,env),summarize,
	mass_mean=mean(mass),belief_mean=mean(belief),info_mean=mean(na.omit(info)),
	mass_median=median(mass),belief_median=median(belief),info_median= median(na.omit(info)),
	mass_sd=sd(mass)/sqrt(length(mass)),belief_sd=sd(belief)/sqrt(length(mass)),info_sd=sd(na.omit(info))/sqrt(length(na.omit(info))),
	mass_iqr1=iqr1(mass),belief_iqr1= iqr1(belief),info_iqr1= iqr1(na.omit(info)),
	mass_iqr2=iqr2(mass),belief_iqr2= iqr2(belief),info_iqr2= iqr2(na.omit(info)),
	n=length(sim_it))

	# avoid bias of small sample size
	summsim <- summsim[summsim$age>1&summsim$n>=100,]  # CHECK OK TO DO THIS? 

	median1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("info","median",sep="_"))])
	median2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("info","median",sep="_"))])
	iqrL1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("info","iqr1",sep="_"))])
	iqrL2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("info","iqr1",sep="_"))])
	iqrU1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("info","iqr2",sep="_"))])
	iqrU2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("info","iqr2",sep="_"))])

	# get difference in alpha as a result of experience
	plotdat0 <- simdat1[simdat1$state==1 & simdat1$age>1,] # only works while still foraging and after first time-step (alpha = 0)
	
	# limit to age until first maturity in each environment
	plotdat <- plotdat0[(plotdat0$env==0 & plotdat0$age<tapply(summsim$age,summsim$env,max)[1])|(plotdat0$env==1 & plotdat0$age<tapply(summsim$age,summsim$env,max)[2]),]
	
	
	# calculate differences
	plotdat1 <- cbind(plotdat[2:length(plotdat[,1]),],plotdat[-length(plotdat[,1]),])
	plotdat1$simno.diff = plotdat$simno[2:length(plotdat[,1])] - plotdat$simno[-length(plotdat[,1])]
	plotdat1$it.diff = plotdat$it[2:length(plotdat[,1])] - plotdat$it[-length(plotdat[,1])]
	plotdat1$forage.diff = abs(plotdat$forage[2:length(plotdat[,1])] - plotdat$forage[-length(plotdat[,1])])
	plotdat1$mass.diff = abs(plotdat$mass[2:length(plotdat[,1])] - plotdat$mass[-length(plotdat[,1])])
	
	plotdat1 <- plotdat1[plotdat1$simno.diff==0 & plotdat1$it.diff==0,]
	
	# get mean delta alpha with age for successess and failures
	plotdat1_summ <- ddply(plotdat1,.(age,mass.diff),summarize,forage.mean=mean(forage.diff),forage.se=sd(forage.diff)/sqrt(length(forage.diff)),forage.n=length(forage.diff))
	
	## as averaged across both environments, should limit to age until first maturity across both environments
	min.age.mat = min(tapply(summsim$age,summsim$env,max)[1],tapply(summsim$age,summsim$env,max)[2]) # 65 (env0)
	plotdat1_summ <- plotdat1_summ[plotdat1_summ$age<min.age.mat,]

# do plot
plot.cex=1.8
pdf(file=paste("Fig3.pdf",sep=""), width=9, height=3.8, family="Times", pointsize=9, useDingbat=FALSE)

par(mfrow=c(1,3),mai=c(0.6,0.6,0.3,0.4),oma=c(0.2,0.2,0.2,0.2))

counter=0

{
	counter=counter+1
	optact = optactList[[parno]][1:60,]
	try.zlim = c(0.6,1) # so get clearer gradation of color
	image(t(matrix(0,60,501)),axes=F,cex.main=1.5,
		xlim=c(0,1),ylim=c(0,1),zlim=c(0,0),col=bg.col,
		ylab=expression(Size~italic(S)),xlab=expression(Belief~italic(P)),cex.lab=plot.cex)
		
	par(new=TRUE)
	image(t(optact),axes=F,cex.main=1.5,xlim=c(0,1),ylim=c(0,1),zlim=try.zlim,col=try.col,xlab="",ylab="",cex.lab=plot.cex)
			
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1.2, adj=-0.25,cex=plot.cex*0.8)
	
	box()
	axis(1,seq(0,1,length.out=6),labels=seq(0,1,length.out=6),cex.axis=plot.cex) 
	axis(2, at=seq(0,1,length.out=7), labels=seq(0,60,length.out=7),cex.axis=plot.cex)

	#par(new=TRUE,mai=c(0.6,0.6,0.2,1.2),oma=c(1,1,1.5,1.5))
	image.plot(t(optact),zlim=try.zlim,col=try.col,axes=F,cex.lab=4,legend.only=T,legend.width=2.5,
		axis.args=list(cex.axis=1.5))

	counter=counter+1
	plot(1:length(median2),type="n",xlab="Age",ylab=expression(Value~of~information~(x10^-5)),cex.axis=plot.cex,cex.lab=plot.cex,axes=F,ylim=ylim.info,xlim=xlim.age) #
	polygon(c(1:length(median1),rev(1:length(median1))), c(iqrU1, rev(iqrL1)), col=col.good.shade, border= col.good.shade)
	polygon(c(1:length(median2),rev(1:length(median2))), c(iqrU2, rev(iqrL2)), col= col.bad.shade,border= col.bad.shade)
	lines(1:length(median1),median1,lwd=2,col=col.good.line)
	lines(1:length(median2),median2,lwd=2,lty=2,col=col.bad.line)
	axis(1,seq(xlim.age[1],xlim.age[2],length.out=6),cex.axis= plot.cex) 
	axis(2,at=seq(ylim.info[1],ylim.info[2],length.out=4),lab=seq(0,1.5,length.out=4),cex.axis= plot.cex)
	box()
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1.2, adj=-0.25,cex= plot.cex*0.8)
		
	counter=counter+1
	plot(0,0,type="n",ylim= ylim.alpha,xlim=xlim.age2,xlab="Age",ylab="Change in foraging effort",cex=plot.cex,cex.lab=plot.cex,cex.axis=plot.cex,axes=F)
	lines(plotdat1_summ$age[plotdat1_summ$mass.diff==0], plotdat1_summ$forage.mean[plotdat1_summ$mass.diff==0],lwd=2, col="red",lty=2)
	lines(plotdat1_summ$age[plotdat1_summ$mass.diff==1], plotdat1_summ$forage.mean[plotdat1_summ$mass.diff==1],lwd=2)
	axis(1,seq(0,60,10),cex.axis=plot.cex) 
	axis(2,at=seq(ylim.alpha[1], ylim.alpha[2],length.out=6),lab=seq(0,ylim.alpha[2],length.out=6),cex.axis= plot.cex)
	box()
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1.2, adj=-0.25,cex= plot.cex*0.8)
	
}		

dev.off()




### PLOT EFFECT OF PRIORS AND EXPERIENCES


use.adj=-0.35

# PLOT MATERNAL PRIORS
# get baselines sims: 

col.good = c("dodgerblue2","dodgerblue3","dodgerblue4") #c("grey50","grey40","grey30") #c("lightskyblue","grey60","royalblue")
col.bad = c("orange3","orange2","orange1") #c("grey70","grey80","grey90") #c("sienna3","grey60","lightsalmon")

# YLIM PRIORS
ylim.age.bad <- c(90,110)
ylim.age.good <- c(70,90)

ylim.size.bad <- c(25,40)
ylim.size.good <- c(40,55)

ylim.fit.bad <- c(0,0.1)
ylim.fit.good <- c(0.1,0.2)

	
# data
for(parno in 2)
{
	plot_simno <- meandat$simno[meandat$parno==parno & meandat$age==0 & meandat$pfood==0]
	
	age1dat <- rbind(meandat$age1[meandat$simno %in% plot_simno],meandat$sdage1[meandat$simno %in% plot_simno])
	age2dat <- rbind(meandat$age2[meandat$simno %in% plot_simno],meandat$sdage2[meandat$simno %in% plot_simno])
	
	size1dat <- rbind(meandat$size1[meandat$simno %in% plot_simno],meandat$sdsize1[meandat$simno %in% plot_simno])
	size2dat <- rbind(meandat$size2[meandat$simno %in% plot_simno],meandat$sdsize2[meandat$simno %in% plot_simno])

	optact <- optactList[[parno]]
	opt.size.mat <- c(min(which(is.na(optact[,1]))),min(which(is.na(optact[,501]))))


	plot1dat <- rbind(meandat$fit1[meandat$simno %in% plot_simno],meandat$sdfit1[meandat$simno %in% plot_simno])
	plot2dat <- rbind(meandat$fit2[meandat$simno %in% plot_simno],meandat$sdfit2[meandat$simno %in% plot_simno])
}	



pdf(file="Fig4.pdf", width=7, height=4.5, family="Times", pointsize=9, useDingbat=FALSE)

par(mfrow=c(2,3),mai=c(0.6,0.6,0.2,0.2),oma=c(1,1,1.5,1.5))
	
counter=0
	
for(parno in 2)
{
	
	y.lim=ylim.age.bad
	
	counter=counter+1
	bpage2 <- barplot(age2dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Age at maturity",cex.lab=plot.cex,cex.names= plot.cex, col=col.bad,cex.axis=plot.cex)
	abline(h=0)
	arrows(bpage2, age2dat[1,]+ age2dat[2,]-y.lim[1], bpage2, age2dat[1,]-age2dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	axis(2,at=seq(0,y.lim[2]-y.lim[1],10),labels=seq(y.lim[1],y.lim[2],10),cex.axis=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1.5, adj=-0.25,cex=plot.cex*0.8)
	
	
	y.lim=ylim.size.bad

	counter=counter+1	
	bpsize2 <- barplot(size2dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Size at maturity",cex.lab= plot.cex,cex.names= plot.cex, col=col.bad,cex.axis=plot.cex)
	abline(h=opt.size.mat[1]-y.lim[1],lty=2)
	abline(h=0)
	arrows(bpsize2, size2dat[1,]+ size2dat[2,]-y.lim[1], bpsize2, size2dat[1,]-size2dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis= plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)

	y.lim=ylim.fit.bad

	counter=counter+1	
	bpplot2 <- barplot(plot2dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Fitness",cex.lab= plot.cex,cex.names= plot.cex, col=col.bad,cex.axis=plot.cex)
	arrows(bpplot2, plot2dat[1,]+ plot2dat[2,]-y.lim[1], bpplot2, plot2dat[1,]-plot2dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	axis(2,at=seq(0,y.lim[2]-y.lim[1],0.05),labels=seq(y.lim[1],y.lim[2],0.05),cex.axis= plot.cex)
	abline(h=0)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)
	

	y.lim=ylim.age.good
	
	counter=counter+1
	bpage1 <- barplot(age1dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Age at maturity",cex.lab= plot.cex,cex.names= plot.cex, col=col.good,cex.axis=plot.cex)
	abline(h=0)
	arrows(bpage1, age1dat[1,]+ age1dat[2,]-y.lim[1], bpage1, age1dat[1,]-age1dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	axis(2,at=seq(0,y.lim[2]-y.lim[1],10),labels=seq(y.lim[1],y.lim[2],10),cex.axis= plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)
	
	y.lim=ylim.size.good
	
	counter=counter+1	
	bpsize1 <- barplot(size1dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Size at maturity",cex.lab= plot.cex,cex.names= plot.cex, col=col.good,cex.axis=plot.cex)
	abline(h=opt.size.mat[2]-y.lim[1],lty=2)
	abline(h=0)
	arrows(bpsize1, size1dat[1,]+ size1dat[2,]-y.lim[1], bpsize1, size1dat[1,]-size1dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis= plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)

	y.lim=ylim.fit.good

	counter=counter+1	
	bpplot1 <- barplot(plot1dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Fitness",cex.lab= plot.cex,cex.names= plot.cex, col=col.good,cex.axis=plot.cex)
	arrows(bpplot1, plot1dat[1,]+ plot1dat[2,]-y.lim[1], bpplot1, plot1dat[1,]-plot1dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	abline(h=0)
	axis(2,at=seq(0,y.lim[2]-y.lim[1],0.05),labels=seq(y.lim[1],y.lim[2],0.05),cex.axis= plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)
}		
mtext(expression(Starting~belief~italic(P)),side=1,outer=T,adj=0.5,line=-0.8,cex=plot.cex*0.8)

dev.off()


	
### PLOT FORAGING EXPERIENCES

# Plot poor/rich conditions = PLOT AGE, SIZE, FITNESS
plot.cex=2
plot.col=c("grey50","orange1","darkblue")
use.lwd=1.5
use.9

# YLIM EXPERIENCES
ylim.age.bad <- c(80,110)
ylim.age.good <- c(70,100)

ylim.size.bad <- c(25,40)
ylim.size.good <- c(40,55)

ylim.fit.bad <- c(0,0.1)
ylim.fit.good <- c(0.1,0.2)

# data
parno = 2

plot_simno <- c(baseline_simno[parno],meandat$simno[meandat$parno==parno & meandat$pfood!=0])
	
plotage1 <- rbind(meandat$age1[meandat$simno %in% plot_simno],meandat$sdage1[meandat$simno %in% plot_simno])
plotage2 <- rbind(meandat$age2[meandat$simno %in% plot_simno],meandat$sdage2[meandat$simno %in% plot_simno])

plotsize1 <- rbind(meandat$size1[meandat$simno %in% plot_simno],meandat$sdsize1[meandat$simno %in% plot_simno])
plotsize2 <- rbind(meandat$size2[meandat$simno %in% plot_simno],meandat$sdsize2[meandat$simno %in% plot_simno])

plotfit1 <- rbind(meandat$fit1[meandat$simno %in% plot_simno],meandat$sdfit1[meandat$simno %in% plot_simno])
plotfit2 <- rbind(meandat$fit2[meandat$simno %in% plot_simno],meandat$sdfit2[meandat$simno %in% plot_simno])

	
pdf(file="Fig5.pdf", width=7, height=4.5, family="Times", pointsize=9, useDingbat=FALSE)

par(mfrow=c(2,3),mai=c(0.6,0.6,0.2,0.2),oma=c(1,1,1.5,1.5))
		
counter=0
	
## LOW-RESOURCE ENVIRONMENT

# AGE

y.lim=ylim.age.bad
plot2dat <- plotage2

counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Age at maturity",cex.lab=plot.cex)

# baseline
points(1, plot2dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot2dat[1,1]+ plot2dat[2,1]-y.lim[1],1, plot2dat[1,1]-plot2dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot2dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot2dat[1,2:7]-y.lim[1],type="o",pch=1,lty=5,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot2dat[1,2:7]+ plot2dat[2,2:7]-y.lim[1],2:7, plot2dat[1,2:7]-plot2dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=1,col=plot.col[2])	
# rich conditions
points((2:7)+0.1, plot2dat[1,8:13]-y.lim[1],type="o",pch=16,lty=1,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot2dat[1,8:13]+ plot2dat[2,8:13]-y.lim[1],(2:7)+0.1, plot2dat[1,8:13]-plot2dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[2])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)



# SIZE

y.lim=ylim.size.bad
plot2dat <- plotsize2

counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Size at maturity",cex.lab=plot.cex)

# baseline
points(1, plot2dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot2dat[1,1]+ plot2dat[2,1]-y.lim[1],1, plot2dat[1,1]-plot2dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot2dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot2dat[1,2:7]-y.lim[1],type="o",pch=1,lty=5,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot2dat[1,2:7]+ plot2dat[2,2:7]-y.lim[1],2:7, plot2dat[1,2:7]-plot2dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=5,col=plot.col[2])	
# rich conditions
points((2:7)+0.1, plot2dat[1,8:13]-y.lim[1],type="o",pch=16,lty=1,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot2dat[1,8:13]+ plot2dat[2,8:13]-y.lim[1],(2:7)+0.1, plot2dat[1,8:13]-plot2dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[2])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)



# FITNESS

y.lim=ylim.fit.bad
plot2dat <- plotfit2

## LOW-RESOURCE ENVIRONMENT
counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Fitness",cex.lab=plot.cex)

# baseline
points(1, plot2dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot2dat[1,1]+ plot2dat[2,1]-y.lim[1],1, plot2dat[1,1]-plot2dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot2dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot2dat[1,2:7]-y.lim[1],type="o",pch=1,lty=5,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot2dat[1,2:7]+ plot2dat[2,2:7]-y.lim[1],2:7, plot2dat[1,2:7]-plot2dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=5,col=plot.col[2])	
# rich conditions
points((2:7)+0.1, plot2dat[1,8:13]-y.lim[1],type="o",pch=16,lty=1,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot2dat[1,8:13]+ plot2dat[2,8:13]-y.lim[1],(2:7)+0.1, plot2dat[1,8:13]-plot2dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[2])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],0.05),labels=seq(y.lim[1],y.lim[2],0.05),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)



## HIGH-RESOURCE ENVIRONMENT

# AGE

y.lim=ylim.age.good
plot1dat <- plotage1

counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Age at maturity",cex.lab=plot.cex)

# baseline
points(1, plot1dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot1dat[1,1]+ plot1dat[2,1]-y.lim[1],1, plot1dat[1,1]-plot1dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot1dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot1dat[1,2:7]-y.lim[1],type="o",pch=1,lty=5,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot1dat[1,2:7]+ plot1dat[2,2:7]-y.lim[1],2:7, plot1dat[1,2:7]-plot1dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=1,col=plot.col[3])	
# rich conditions
points((2:7)+0.1, plot1dat[1,8:13]-y.lim[1],type="o",pch=16,lty=1,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot1dat[1,8:13]+ plot1dat[2,8:13]-y.lim[1],(2:7)+0.1, plot1dat[1,8:13]-plot1dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[3])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)



# SIZE

y.lim=ylim.size.good
plot1dat <- plotsize1

counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Size at maturity",cex.lab=plot.cex)

# baseline
points(1, plot1dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot1dat[1,1]+ plot1dat[2,1]-y.lim[1],1, plot1dat[1,1]-plot1dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot1dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot1dat[1,2:7]-y.lim[1],type="o",pch=1,lty=5,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot1dat[1,2:7]+ plot1dat[2,2:7]-y.lim[1],2:7, plot1dat[1,2:7]-plot1dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=5,col=plot.col[3])	
# rich conditions
points((2:7)+0.1, plot1dat[1,8:13]-y.lim[1],type="o",pch=16,lty=1,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot1dat[1,8:13]+ plot1dat[2,8:13]-y.lim[1],(2:7)+0.1, plot1dat[1,8:13]-plot1dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[3])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)



# FITNESS

y.lim=ylim.fit.good
plot1dat <- plotfit1

counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Fitness",cex.lab=plot.cex)

# baseline
points(1, plot1dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot1dat[1,1]+ plot1dat[2,1]-y.lim[1],1, plot1dat[1,1]-plot1dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot1dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot1dat[1,2:7]-y.lim[1],type="o",pch=1,lty=5,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot1dat[1,2:7]+ plot1dat[2,2:7]-y.lim[1],2:7, plot1dat[1,2:7]-plot1dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=5,col=plot.col[3])	
# rich conditions
points((2:7)+0.1, plot1dat[1,8:13]-y.lim[1],type="o",pch=16,1,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot1dat[1,8:13]+ plot1dat[2,8:13]-y.lim[1],(2:7)+0.1, plot1dat[1,8:13]-plot1dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[3])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],0.05),labels=seq(y.lim[1],y.lim[2],0.05),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)


# add outer margin text
mtext("Age at onset of 10-step manipulated foraging experiences",side=1,outer=T,adj=0.5,line=-1,cex=plot.cex*0.7)


dev.off()








### FOR SUPPLEMENTARY MATERIAL - APPENDIX A
 
 # add plots of info arrays and value of info etc for par 1 and 3
 
 # do plot
pdf(file=paste("FigA1.pdf",sep=""), width=9, height=5.5, family="Times", pointsize=9, useDingbat=FALSE)

par(mfrow=c(2,3),mai=c(0.6,0.6,0.2,0.2),oma=c(1,1,1.5,1.5))

plot.cex=2

for(parno in 1) 
# sort out data
	# calculate change in information value with age
	simdat1 <- get(paste("simdat","ind",parno,0,0.5,"0.0",0,sep="_"))
	simdat1$sim_it <- paste(simdat1$simno,simdat1$it,simdat1$env,sep="-")
	summsim<-ddply(simdat1,.(age,env),summarize,
	mass_mean=mean(mass),belief_mean=mean(belief),info_mean=mean(na.omit(info)),
	mass_median=median(mass),belief_median=median(belief),info_median= median(na.omit(info)),
	mass_sd=sd(mass)/sqrt(length(mass)),belief_sd=sd(belief)/sqrt(length(mass)),info_sd=sd(na.omit(info))/sqrt(length(na.omit(info))),
	mass_iqr1=iqr1(mass),belief_iqr1= iqr1(belief),info_iqr1= iqr1(na.omit(info)),
	mass_iqr2=iqr2(mass),belief_iqr2= iqr2(belief),info_iqr2= iqr2(na.omit(info)),
	n=length(sim_it))

	# avoid bias of small sample size
	summsim <- summsim[summsim$age>1&summsim$n>=100,]  # CHECK OK TO DO THIS? 

	median1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("info","median",sep="_"))])
	median2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("info","median",sep="_"))])
	iqrL1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("info","iqr1",sep="_"))])
	iqrL2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("info","iqr1",sep="_"))])
	iqrU1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("info","iqr2",sep="_"))])
	iqrU2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("info","iqr2",sep="_"))])

	# get difference in alpha as a result of experience
	plotdat0 <- simdat1[simdat1$state==1 & simdat1$age>1,] # only works while still foraging and after first time-step (alpha = 0)
	
	# ## limit to age until first maturity in each environment
	# plotdat <- plotdat0[(plotdat0$env==0 & plotdat0$age<tapply(summsim$age,summsim$env,max)[1])|(plotdat0$env==1 & plotdat0$age<tapply(summsim$age,summsim$env,max)[2]),]
	
	## limit to age until first maturity in BOTH environment
	plotdat <- plotdat0[plotdat0$age<min(c(tapply(summsim$age,summsim$env,max)[1],tapply(summsim$age,summsim$env,max)[2])),]

		
	# calculate differences
	plotdat1 <- cbind(plotdat[2:length(plotdat[,1]),],plotdat[-length(plotdat[,1]),])
	plotdat1$simno.diff = plotdat$simno[2:length(plotdat[,1])] - plotdat$simno[-length(plotdat[,1])]
	plotdat1$it.diff = plotdat$it[2:length(plotdat[,1])] - plotdat$it[-length(plotdat[,1])]
	plotdat1$forage.diff = abs(plotdat$forage[2:length(plotdat[,1])] - plotdat$forage[-length(plotdat[,1])])
	plotdat1$mass.diff = abs(plotdat$mass[2:length(plotdat[,1])] - plotdat$mass[-length(plotdat[,1])])
	
	plotdat1 <- plotdat1[plotdat1$simno.diff==0 & plotdat1$it.diff==0,]
	
	# get mean delta alpha with age for successess and failures
	plotdat1_summ <- ddply(plotdat1,.(age,mass.diff),summarize,forage.mean=mean(forage.diff),forage.se=sd(forage.diff)/sqrt(length(forage.diff)),forage.n=length(forage.diff))
	

counter=0

{
	counter=counter+1
	optact = optactList[[parno]][1:60,]
	try.zlim = c(0.6,1) # so get clearer gradation of color
	image(t(matrix(0,60,501)),axes=F,cex.main=1.5,
		xlim=c(0,1),ylim=c(0,1),zlim=c(0,0),col=bg.col,
		ylab=expression(Size~italic(S)),xlab=expression(Belief~italic(P)),cex.lab=plot.cex)
		
	par(new=TRUE)
	image(t(optact),axes=F,cex.main=1.5,xlim=c(0,1),ylim=c(0,1),zlim=try.zlim,col=try.col,xlab="",ylab="",cex.lab=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1, adj=-0.25,cex=plot.cex*0.8)
	
	box()
	axis(1,seq(0,1,length.out=6),labels=seq(0,1,length.out=6),cex.axis=plot.cex) 
	axis(2, at=seq(0,1,length.out=7), labels=seq(0,60,length.out=7),cex.axis=plot.cex)

	counter=counter+1
	plot(1:length(median2),type="n",xlab="Age",ylab=expression(Value~of~information~(x10^-5)),cex.axis=plot.cex,cex.lab=plot.cex,axes=F,ylim=ylim.info,xlim=xlim.age) #
	polygon(c(1:length(median1),rev(1:length(median1))), c(iqrU1, rev(iqrL1)), col=col.good.shade, border= col.good.shade)
	polygon(c(1:length(median2),rev(1:length(median2))), c(iqrU2, rev(iqrL2)), col= col.bad.shade,border= col.bad.shade)
	lines(1:length(median1),median1,lwd=2,col=col.good.line)
	lines(1:length(median2),median2,lwd=2,lty=2,col=col.bad.line)
	axis(1,seq(xlim.age[1],xlim.age[2],length.out=6),cex.axis= plot.cex) 
	axis(2,at=seq(ylim.info[1],ylim.info[2],length.out=4),lab=seq(0,1.5,length.out=4),cex.axis= plot.cex)
	box()
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1, adj=-0.4,cex=plot.cex*0.8)
		
	counter=counter+1
	plot(0,0,type="n",ylim= ylim.alpha,xlim=xlim.age,xlab="Age",ylab="Change in foraging effort",cex=plot.cex,cex.lab=plot.cex,cex.axis=plot.cex)
	lines(plotdat1_summ$age[plotdat1_summ$mass.diff==0], plotdat1_summ$forage.mean[plotdat1_summ$mass.diff==0],lwd=2, col="red",lty=2)
	lines(plotdat1_summ$age[plotdat1_summ$mass.diff==1], plotdat1_summ$forage.mean[plotdat1_summ$mass.diff==1],lwd=2)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1, adj=-0.4,cex=plot.cex*0.8)
	
}		




for(parno in 3) 
# sort out data
	# calculate change in information value with age
	simdat1 <- get(paste("simdat","ind",parno,0,0.5,"0.0",0,sep="_"))
	simdat1$sim_it <- paste(simdat1$simno,simdat1$it,simdat1$env,sep="-")
	summsim<-ddply(simdat1,.(age,env),summarize,
	mass_mean=mean(mass),belief_mean=mean(belief),info_mean=mean(na.omit(info)),
	mass_median=median(mass),belief_median=median(belief),info_median= median(na.omit(info)),
	mass_sd=sd(mass)/sqrt(length(mass)),belief_sd=sd(belief)/sqrt(length(mass)),info_sd=sd(na.omit(info))/sqrt(length(na.omit(info))),
	mass_iqr1=iqr1(mass),belief_iqr1= iqr1(belief),info_iqr1= iqr1(na.omit(info)),
	mass_iqr2=iqr2(mass),belief_iqr2= iqr2(belief),info_iqr2= iqr2(na.omit(info)),
	n=length(sim_it))

	# avoid bias of small sample size
	summsim <- summsim[summsim$age>1&summsim$n>=100,]  # CHECK OK TO DO THIS? 

	median1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("info","median",sep="_"))])
	median2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("info","median",sep="_"))])
	iqrL1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("info","iqr1",sep="_"))])
	iqrL2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("info","iqr1",sep="_"))])
	iqrU1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("info","iqr2",sep="_"))])
	iqrU2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("info","iqr2",sep="_"))])

	# get difference in alpha as a result of experience
	plotdat0 <- simdat1[simdat1$state==1 & simdat1$age>1,] # only works while still foraging and after first time-step (alpha = 0)
	
	# # limit to age until first maturity in each environment
	# plotdat <- plotdat0[(plotdat0$env==0 & plotdat0$age<tapply(summsim$age,summsim$env,max)[1])|(plotdat0$env==1 & plotdat0$age<tapply(summsim$age,summsim$env,max)[2]),]
	
		## limit to age until first maturity in BOTH environment
	plotdat <- plotdat0[plotdat0$age<min(c(tapply(summsim$age,summsim$env,max)[1],tapply(summsim$age,summsim$env,max)[2])),]

	# calculate differences
	plotdat1 <- cbind(plotdat[2:length(plotdat[,1]),],plotdat[-length(plotdat[,1]),])
	plotdat1$simno.diff = plotdat$simno[2:length(plotdat[,1])] - plotdat$simno[-length(plotdat[,1])]
	plotdat1$it.diff = plotdat$it[2:length(plotdat[,1])] - plotdat$it[-length(plotdat[,1])]
	plotdat1$forage.diff = abs(plotdat$forage[2:length(plotdat[,1])] - plotdat$forage[-length(plotdat[,1])])
	plotdat1$mass.diff = abs(plotdat$mass[2:length(plotdat[,1])] - plotdat$mass[-length(plotdat[,1])])
	
	plotdat1 <- plotdat1[plotdat1$simno.diff==0 & plotdat1$it.diff==0,]
	
	# get mean delta alpha with age for successess and failures
	plotdat1_summ <- ddply(plotdat1,.(age,mass.diff),summarize,forage.mean=mean(forage.diff),forage.se=sd(forage.diff)/sqrt(length(forage.diff)),forage.n=length(forage.diff))
	

# do plot

{
	counter=counter+1
	optact = optactList[[parno]][1:60,]
	try.zlim = c(0.6,1) # so get clearer gradation of color
	image(t(matrix(0,60,501)),axes=F,cex.main=1.5,
		xlim=c(0,1),ylim=c(0,1),zlim=c(0,0),col=bg.col,
		ylab=expression(Size~italic(S)),xlab=expression(Belief~italic(P)),cex.lab=plot.cex)
		
	par(new=TRUE)
	image(t(optact),axes=F,cex.main=1.5,xlim=c(0,1),ylim=c(0,1),zlim=try.zlim,col=try.col,xlab="",ylab="",cex.lab=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1, adj=-0.25,cex=plot.cex*0.8)
	
	box()
	axis(1,seq(0,1,length.out=6),labels=seq(0,1,length.out=6),cex.axis=plot.cex) 
	axis(2, at=seq(0,1,length.out=7), labels=seq(0,60,length.out=7),cex.axis=plot.cex)

	counter=counter+1
	plot(1:length(median2),type="n",xlab="Age",ylab=expression(Value~of~information~(x10^-5)),cex.axis=plot.cex,cex.lab=plot.cex,axes=F,ylim=ylim.info,xlim=xlim.age) #
	polygon(c(1:length(median1),rev(1:length(median1))), c(iqrU1, rev(iqrL1)), col=col.good.shade, border= col.good.shade)
	polygon(c(1:length(median2),rev(1:length(median2))), c(iqrU2, rev(iqrL2)), col= col.bad.shade,border= col.bad.shade)
	lines(1:length(median1),median1,lwd=2,col=col.good.line)
	lines(1:length(median2),median2,lwd=2,lty=2,col=col.bad.line)
	axis(1,seq(xlim.age[1],xlim.age[2],length.out=6),cex.axis= plot.cex) 
	axis(2,at=seq(ylim.info[1],ylim.info[2],length.out=4),lab=seq(0,1.5,length.out=4),cex.axis= plot.cex)
	box()
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1, adj=-0.4,cex=plot.cex*0.8)
		
	counter=counter+1
	plot(0,0,type="n",ylim= ylim.alpha,xlim=xlim.age,xlab="Age",ylab="Change in foraging effort",cex=plot.cex,cex.lab=plot.cex,cex.axis=plot.cex)
	lines(plotdat1_summ$age[plotdat1_summ$mass.diff==0], plotdat1_summ$forage.mean[plotdat1_summ$mass.diff==0],lwd=2, col="red",lty=2)
	lines(plotdat1_summ$age[plotdat1_summ$mass.diff==1], plotdat1_summ$forage.mean[plotdat1_summ$mass.diff==1],lwd=2)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1, adj=-0.4,cex=plot.cex*0.8)
	
}		

dev.off()




#### PAPER FIGURES PHENOTYPIC OUTCOMES: change in belief with age, and distributions of size and age at maturity 

pdf(file="FigA2.pdf", width=6.7, height=6.7, family="Times", pointsize=9, useDingbat=FALSE)

par(mfrow=c(3,3),mar=c(5,6,3,3))
plot.cex=1.5			
counter = 0

for(parno in 1:3)
{	
	
	simdat1 <- get(paste("simdat","ind",parno,0,0.5,"0.0",0,sep="_"))
	simdat1$sim_it <- paste(simdat1$simno,simdat1$it,simdat1$env,sep="-")
	summsim<-ddply(simdat1,.(age,env),summarize,
		mass_mean=mean(mass),belief_mean=mean(belief),info_mean=mean(na.omit(info)),
		mass_median=median(mass),belief_median=median(belief),info_median= median(na.omit(info)),
		mass_sd=sd(mass)/sqrt(length(mass)),belief_sd=sd(belief)/sqrt(length(mass)),info_sd=sd(na.omit(info))/sqrt(length(na.omit(info))),
		mass_iqr1=iqr1(mass),belief_iqr1= iqr1(belief),info_iqr1= iqr1(na.omit(info)),
		mass_iqr2=iqr2(mass),belief_iqr2= iqr2(belief),info_iqr2= iqr2(na.omit(info)),
		n=length(sim_it))

	# avoid bias of small sample size
	summsim <- summsim[summsim$n>=100,]  

	mean1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("belief","mean",sep="_"))])
	mean2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("belief","mean",sep="_"))])
	sd1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("belief","sd",sep="_"))])
	sd2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("belief","sd",sep="_"))])	
	median1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("belief","median",sep="_"))])
	median2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("belief","median",sep="_"))])
	iqrL1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("belief","iqr1",sep="_"))])
	iqrL2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("belief","iqr1",sep="_"))])
	iqrU1 = na.omit(summsim[summsim$env==0,which(names(summsim) %in% paste("belief","iqr2",sep="_"))])
	iqrU2 = na.omit(summsim[summsim$env==1,which(names(summsim) %in% paste("belief","iqr2",sep="_"))])

	
	# get average across ALL simulations of size and age at maturation
	simpop1 <- get(paste("simdat","pop",parno,1,0.5,"0.0",0,sep="_"))
	
	agesizemat <- simpop1[simpop1$state==2,]
	summsize <- ddply(agesizemat,.(env,mass),summarise,n=length(mass))
	summage <- ddply(agesizemat,.(env,age),summarise,n=length(age))
		
	summsize1 <- data.frame(mass=seq(xlim.freqSize[1],xlim.freqSize[2],1)) 
	summage1 <- data.frame(age=seq(xlim.freqAge[1],xlim.freqAge[2],1)) 
		
	summsize1 <- merge(summsize1, summsize[summsize$env==0,],by="mass",all.x=TRUE)
	summsize1 <- merge(summsize1, summsize[summsize$env==1,],by="mass",all.x=TRUE)
	summage1 <- merge(summage1, summage[summage$env==0,],by="age",all.x=TRUE)
	summage1 <- merge(summage1, summage[summage$env==1,],by="age",all.x=TRUE)
		
	plot_sizemat <- rbind(summsize1[,3], summsize1[,5])
	colnames(plot_sizemat) <- summsize1[,1]
	
	plot_agemat <- rbind(summage1[,3], summage1[,5])
	colnames(plot_agemat) <- summage1[,1]

	counter = counter+1
	plot(1:length(median2),type="n",xlim=c(0,length(median2)),ylim=0:1,xlab="Age",ylab=expression(Belief~italic(P)),cex.axis= plot.cex,cex.lab= plot.cex,col= col.bad.shade)
	polygon(c(1:length(median1),rev(1:length(median1))), c(iqrU1, rev(iqrL1)), col=col.good.shade, border= col.good.shade)
	polygon(c(1:length(median2),rev(1:length(median2))), c(iqrU2, rev(iqrL2)), col= col.bad.shade,border= col.bad.shade)
	lines(1:length(median1),median1,lwd=2,col=col.good.line)
	lines(1:length(median2),median2,lwd=2,lty=2,col=col.bad.line)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=0.8, adj=-0.25,cex=1.5)

	counter = counter+1		
	bp <- barplot(plot_sizemat,col=c(col.good.shade,col.bad.shade),border= c(col.good.line,col.bad.line), beside=T,xlab=expression(Size~italic(S)~at~maturity),ylab="Frequency",cex.axis=plot.cex,cex.lab=plot.cex,cex=plot.cex,ylim=ylim.freqSize,names.arg=rep(" ",length(plot_sizemat[1,])))
	par(new=T)	
	plot(0,0,type="n",axes=F,xlim=c(30,48),ylab="",xlab="")
	axis(1,cex.axis=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=0.8, adj=-0.25,cex=1.5)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=0.8, adj=-0.25,cex=1.5)
	box()
	
	counter = counter+1		
	bp <- barplot(plot_agemat,col=c(col.good.shade,col.bad.shade),border= c(col.good.line,col.bad.line),beside=T,xlab=expression(Age~at~maturity),ylab="Frequency",cex.axis= plot.cex,cex.lab=plot.cex,cex=plot.cex,ylim=ylim.freqAge,names.arg=rep(" ",length(plot_agemat[1,])))
	par(new=T)	
	plot(0,0,type="n",axes=F,xlim=c(60,141),ylab="",xlab="")
	axis(1,cex.axis=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=0.8, adj=-0.25,cex=1.5)
	box()
}


dev.off()

### PLOT EFFECT OF PRIORS AND EXPERIENCES

ylim.age.bad <- c(70,110)
ylim.age.good <- c(70,110)

ylim.size.bad <- c(25,40)
ylim.size.good <- c(40,55)

ylim.fit.bad <- c(0,0.1)
ylim.fit.good <- c(0.05,0.15)

use.adj=-0.35

# PLOT MATERNAL PRIORS
# get baselines sims: 

col.good = c("dodgerblue2","dodgerblue3","dodgerblue4") #c("grey50","grey40","grey30") #c("lightskyblue","grey60","royalblue")
col.bad = c("orange3","orange2","orange1") #c("grey70","grey80","grey90") #c("sienna3","grey60","lightsalmon")
plot.cex=2
# PLOT ALL
	
# data
for(parno in 3)
{
	plot_simno <- meandat$simno[meandat$parno==parno & meandat$age==0 & meandat$pfood==0]
	
	age1dat <- rbind(meandat$age1[meandat$simno %in% plot_simno],meandat$sdage1[meandat$simno %in% plot_simno])
	age2dat <- rbind(meandat$age2[meandat$simno %in% plot_simno],meandat$sdage2[meandat$simno %in% plot_simno])
	
	size1dat <- rbind(meandat$size1[meandat$simno %in% plot_simno],meandat$sdsize1[meandat$simno %in% plot_simno])
	size2dat <- rbind(meandat$size2[meandat$simno %in% plot_simno],meandat$sdsize2[meandat$simno %in% plot_simno])

	optact <- optactList[[parno]]
	opt.size.mat <- c(min(which(is.na(optact[,1]))),min(which(is.na(optact[,501]))))


	plot1dat <- rbind(meandat$fit1[meandat$simno %in% plot_simno],meandat$sdfit1[meandat$simno %in% plot_simno])
	plot2dat <- rbind(meandat$fit2[meandat$simno %in% plot_simno],meandat$sdfit2[meandat$simno %in% plot_simno])
}	



pdf(file="FigA3.pdf", width=7, height=4.5, family="Times", pointsize=9, useDingbat=FALSE)

par(mfrow=c(2,3),mai=c(0.6,0.6,0.2,0.2),oma=c(1,1,1.5,1.5))
	
counter=0
	
for(parno in 2)
{
	
	y.lim=ylim.age.bad
	
	counter=counter+1
	bpage2 <- barplot(age2dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Age at maturity",cex.lab=plot.cex,cex.names= plot.cex, col=col.bad,cex.axis=plot.cex)
	abline(h=0)
	arrows(bpage2, age2dat[1,]+ age2dat[2,]-y.lim[1], bpage2, age2dat[1,]-age2dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	axis(2,at=seq(0,y.lim[2]-y.lim[1],10),labels=seq(y.lim[1],y.lim[2],10),cex.axis=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1.5, adj=-0.25,cex=plot.cex*0.8)
	
	
	y.lim=ylim.size.bad

	counter=counter+1	
	bpsize2 <- barplot(size2dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Size at maturity",cex.lab= plot.cex,cex.names= plot.cex, col=col.bad,cex.axis=plot.cex)
	abline(h=opt.size.mat[1]-y.lim[1],lty=2)
	abline(h=0)
	arrows(bpsize2, size2dat[1,]+ size2dat[2,]-y.lim[1], bpsize2, size2dat[1,]-size2dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis= plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)

	y.lim=ylim.fit.bad

	counter=counter+1	
	bpplot2 <- barplot(plot2dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Fitness",cex.lab= plot.cex,cex.names= plot.cex, col=col.bad,cex.axis=plot.cex)
	arrows(bpplot2, plot2dat[1,]+ plot2dat[2,]-y.lim[1], bpplot2, plot2dat[1,]-plot2dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	axis(2,at=seq(0,y.lim[2]-y.lim[1],0.05),labels=seq(y.lim[1],y.lim[2],0.05),cex.axis= plot.cex)
	abline(h=0)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)
	

	y.lim=ylim.age.good
	
	counter=counter+1
	bpage1 <- barplot(age1dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Age at maturity",cex.lab= plot.cex,cex.names= plot.cex, col=col.good,cex.axis=plot.cex)
	abline(h=0)
	arrows(bpage1, age1dat[1,]+ age1dat[2,]-y.lim[1], bpage1, age1dat[1,]-age1dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	axis(2,at=seq(0,y.lim[2]-y.lim[1],10),labels=seq(y.lim[1],y.lim[2],10),cex.axis= plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)
	
	y.lim=ylim.size.good
	
	counter=counter+1	
	bpsize1 <- barplot(size1dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Size at maturity",cex.lab= plot.cex,cex.names= plot.cex, col=col.good,cex.axis=plot.cex)
	abline(h=opt.size.mat[2]-y.lim[1],lty=2)
	abline(h=0)
	arrows(bpsize1, size1dat[1,]+ size1dat[2,]-y.lim[1], bpsize1, size1dat[1,]-size1dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis= plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)

	y.lim=ylim.fit.good

	counter=counter+1	
	bpplot1 <- barplot(plot1dat[1,]-y.lim[1],names.arg=meandat$prior[meandat$simno %in% plot_simno],ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Fitness",cex.lab= plot.cex,cex.names= plot.cex, col=col.good,cex.axis=plot.cex)
	arrows(bpplot1, plot1dat[1,]+ plot1dat[2,]-y.lim[1], bpplot1, plot1dat[1,]-plot1dat[2,]-y.lim[1],
	       length = 0.05, # width of the arrowhead
	       angle = 90, # angle of the arrowhead
	       code = 3 # arrowhead in both ends
	       )
	abline(h=0)
	axis(2,at=seq(0,y.lim[2]-y.lim[1],0.05),labels=seq(y.lim[1],y.lim[2],0.05),cex.axis= plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line= 1.5, adj=-0.25,cex=plot.cex*0.8)
}		
mtext(expression(Starting~belief~italic(P)),side=1,outer=T,adj=0.5,line=-0.8,cex=plot.cex*0.8)

dev.off()




### PLOT FORAGING EXPERIENCES

# Plot poor/rich conditions = PLOT AGE, SIZE, FITNESS
plot.col=c("grey50","orange1","darkblue")
use.lwd=1.5
use.lty=9

# YLIM EXPERIENCES
ylim.age.bad <- c(80,100)
ylim.age.good <- c(90,110)

ylim.size.bad <- c(25,40)
ylim.size.good <- c(40,55)

ylim.fit.bad <- c(0.05,0.15)
ylim.fit.good <- c(0.05,0.15)

# data
parno = 3

plot_simno <- c(baseline_simno[parno],meandat$simno[meandat$parno==parno & meandat$pfood!=0])
	
plotage1 <- rbind(meandat$age1[meandat$simno %in% plot_simno],meandat$sdage1[meandat$simno %in% plot_simno])
plotage2 <- rbind(meandat$age2[meandat$simno %in% plot_simno],meandat$sdage2[meandat$simno %in% plot_simno])

plotsize1 <- rbind(meandat$size1[meandat$simno %in% plot_simno],meandat$sdsize1[meandat$simno %in% plot_simno])
plotsize2 <- rbind(meandat$size2[meandat$simno %in% plot_simno],meandat$sdsize2[meandat$simno %in% plot_simno])

plotfit1 <- rbind(meandat$fit1[meandat$simno %in% plot_simno],meandat$sdfit1[meandat$simno %in% plot_simno])
plotfit2 <- rbind(meandat$fit2[meandat$simno %in% plot_simno],meandat$sdfit2[meandat$simno %in% plot_simno])

	
pdf(file="FigA4.pdf", width=8, height=5, family="Times", pointsize=9, useDingbat=FALSE)

par(mfrow=c(2,3),mai=c(0.6,0.6,0.2,0.2),oma=c(1,1,1.5,1.5))
plot.cex=2
counter=0
	
## LOW-RESOURCE ENVIRONMENT

# AGE

y.lim=ylim.age.bad
plot2dat <- plotage2

counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Age at maturity",cex.lab=plot.cex)

# baseline
points(1, plot2dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot2dat[1,1]+ plot2dat[2,1]-y.lim[1],1, plot2dat[1,1]-plot2dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot2dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot2dat[1,2:7]-y.lim[1],type="o",pch=1,lty=55,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot2dat[1,2:7]+ plot2dat[2,2:7]-y.lim[1],2:7, plot2dat[1,2:7]-plot2dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=1,col=plot.col[2])	
# rich conditions
points((2:7)+0.1, plot2dat[1,8:13]-y.lim[1],type="o",pch=16,lty=1,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot2dat[1,8:13]+ plot2dat[2,8:13]-y.lim[1],(2:7)+0.1, plot2dat[1,8:13]-plot2dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[2])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1.5, adj=-0.25,cex=plot.cex*0.8)



# SIZE

y.lim=ylim.size.bad
plot2dat <- plotsize2

counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Size at maturity",cex.lab=plot.cex)

# baseline
points(1, plot2dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot2dat[1,1]+ plot2dat[2,1]-y.lim[1],1, plot2dat[1,1]-plot2dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot2dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot2dat[1,2:7]-y.lim[1],type="o",pch=1,lty=5,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot2dat[1,2:7]+ plot2dat[2,2:7]-y.lim[1],2:7, plot2dat[1,2:7]-plot2dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=5,col=plot.col[2])	
# rich conditions
points((2:7)+0.1, plot2dat[1,8:13]-y.lim[1],type="o",pch=16,lty=1,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot2dat[1,8:13]+ plot2dat[2,8:13]-y.lim[1],(2:7)+0.1, plot2dat[1,8:13]-plot2dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[2])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1.5, adj=-0.25,cex=plot.cex*0.8)



# FITNESS

y.lim=ylim.fit.bad
plot2dat <- plotfit2

## LOW-RESOURCE ENVIRONMENT
counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Fitness",cex.lab=plot.cex)

# baseline
points(1, plot2dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot2dat[1,1]+ plot2dat[2,1]-y.lim[1],1, plot2dat[1,1]-plot2dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot2dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot2dat[1,2:7]-y.lim[1],type="o",pch=1,lty=5,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot2dat[1,2:7]+ plot2dat[2,2:7]-y.lim[1],2:7, plot2dat[1,2:7]-plot2dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=5,col=plot.col[2])	
# rich conditions
points((2:7)+0.1, plot2dat[1,8:13]-y.lim[1],type="o",pch=16,lty=1,col=plot.col[2],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot2dat[1,8:13]+ plot2dat[2,8:13]-y.lim[1],(2:7)+0.1, plot2dat[1,8:13]-plot2dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[2])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],0.05),labels=seq(y.lim[1],y.lim[2],0.05),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1.5, adj=-0.25,cex=plot.cex*0.8)



## HIGH-RESOURCE ENVIRONMENT

# AGE

y.lim=ylim.age.good
plot1dat <- plotage1

counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Age at maturity",cex.lab=plot.cex)

# baseline
points(1, plot1dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot1dat[1,1]+ plot1dat[2,1]-y.lim[1],1, plot1dat[1,1]-plot1dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot1dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot1dat[1,2:7]-y.lim[1],type="o",pch=1,lty=5,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot1dat[1,2:7]+ plot1dat[2,2:7]-y.lim[1],2:7, plot1dat[1,2:7]-plot1dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=1,col=plot.col[3])	
# rich conditions
points((2:7)+0.1, plot1dat[1,8:13]-y.lim[1],type="o",pch=16,lty=1,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot1dat[1,8:13]+ plot1dat[2,8:13]-y.lim[1],(2:7)+0.1, plot1dat[1,8:13]-plot1dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[3])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1.5, adj=-0.25,cex=plot.cex*0.8)



# SIZE

y.lim=ylim.size.good
plot1dat <- plotsize1

counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Size at maturity",cex.lab=plot.cex)

# baseline
points(1, plot1dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot1dat[1,1]+ plot1dat[2,1]-y.lim[1],1, plot1dat[1,1]-plot1dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot1dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot1dat[1,2:7]-y.lim[1],type="o",pch=1,lty=5,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot1dat[1,2:7]+ plot1dat[2,2:7]-y.lim[1],2:7, plot1dat[1,2:7]-plot1dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=5,col=plot.col[3])	
# rich conditions
points((2:7)+0.1, plot1dat[1,8:13]-y.lim[1],type="o",pch=16,lty=1,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot1dat[1,8:13]+ plot1dat[2,8:13]-y.lim[1],(2:7)+0.1, plot1dat[1,8:13]-plot1dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[3])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],5),labels=seq(y.lim[1],y.lim[2],5),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1.5, adj=-0.25,cex=plot.cex*0.8)



# FITNESS

y.lim=ylim.fit.good
plot1dat <- plotfit1

counter=counter+1
plot(1:7,1:7,type="n",ylim=y.lim-y.lim[1],axes=F,xlab="",ylab="Fitness",cex.lab=plot.cex)

# baseline
points(1, plot1dat[1,1]-y.lim[1],pch=17,cex=plot.cex,col=plot.col[1])
arrows(1, plot1dat[1,1]+ plot1dat[2,1]-y.lim[1],1, plot1dat[1,1]-plot1dat[2,1]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[1])
abline(h= plot1dat[1,1]-y.lim[1],col=plot.col[1],lty=5)

# poor conditions
points(2:7, plot1dat[1,2:7]-y.lim[1],type="o",pch=1,lty=5,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows(2:7, plot1dat[1,2:7]+ plot1dat[2,2:7]-y.lim[1],2:7, plot1dat[1,2:7]-plot1dat[2,2:7]-y.lim[1],length=0.0,angle=90,code=3,lty=5,col=plot.col[3])	
# rich conditions
points((2:7)+0.1, plot1dat[1,8:13]-y.lim[1],type="o",pch=16,lty=1,col=plot.col[3],cex=plot.cex,lwd=use.lwd)
arrows((2:7)+0.1, plot1dat[1,8:13]+ plot1dat[2,8:13]-y.lim[1],(2:7)+0.1, plot1dat[1,8:13]-plot1dat[2,8:13]-y.lim[1],length=0.0,angle=90,code=3,col=plot.col[3])	

axis(2,at=seq(0,y.lim[2]-y.lim[1],0.05),labels=seq(y.lim[1],y.lim[2],0.05),cex.axis=plot.cex)
axis(1,at=1:7,labels=c("B",seq(0,50,10)),cex.axis=plot.cex)
	mtext(paste("(",letters[counter],")",sep=""),side=3, line=1.5, adj=-0.25,cex=plot.cex*0.8)


# add outer margin text
mtext("Age at onset of 10-step manipulated foraging experiences",side=1,outer=T,adj=0.5,line=-1,cex=plot.cex*0.7)


dev.off()





### FOR APPENDIX C: example trajectories

########
# do example trajectories from each environment in baseline and simulation experiments (randomly select five for each simulation) 
########

try.zlim=c(0.6,1)
bg.col = "grey90" # colour for individuals deciding to mature
try.col = rev(heat.colors(501)) # gray.colors(20, start=0.3, end=0.9)
plot.cex=1.1


pdf(file=paste("FigS4.SimExamplesBaseline.pdf",sep=""), width=6.7, height=3.7, family="Times", pointsize=9)
par(mfrow=c(2,5),mai=c(0.3,0.4,0.1,0.1),oma=c(1,1.5,2.2,0.5))	
counter=0			

for(parno in 2:3)
{
	optact=optactList[[parno]]
	info=infoList[[parno]]

	plotdat <- get(paste("simdat","ind",parno,0,0.5,"0.0",0,sep="_"))

	for(i in sample(1:100,5))
	{
		counter=counter+1
		image(matrix(0,80,51),axes=F,cex.main=1.5,xlim=c(0,1),ylim=c(0,1),zlim=c(0,0),col=bg.col)
		par(new=TRUE)
		image(t(optact),axes=F,cex.main=1.5,xlim=c(0,1),ylim=c(0,1),zlim=try.zlim,col=try.col,xlab="",ylab="",cex.lab=1)
		box()
		axis(1,seq(0,1,length.out=9),cex=1)  #,lab=letters[1:length(axTicks(1))])
		axis(2, at=seq(0,1,length.out=9), labels=seq(0,80,length.out=9),cex.lab=1)
				
		plotdf0 <- plotdat[plotdat$simno==1 & plotdat$it==i & plotdat$env==0,]
		plotdf1 <- plotdat[plotdat$simno==1 & plotdat$it==i & plotdat$env==1,]
		
		points((plotdf0$belief),(plotdf0$mass/80 - (1/80)),
			col="black",lwd=3,type="o",lty=1,pch=16,cex=0.8)  # or if want different lty, lty=j+2
			
		points((plotdf1$belief),(plotdf1$mass/80 - (1/80)),
			col="grey",lwd=3,type="o",lty=1,pch=16,cex=0.8)  # or if want different lty, lty=j+2
										
		mtext(paste("(",letters[counter],")",sep=""),side=3, line=0.8, adj=-0.25,cex=plot.cex)

	}
}		

mtext(expression(Size~italic(S)),side=2,line=0,adj=0.5,outer=T)
mtext(expression(Belief~italic(P)),side=1,line=0,adj=0.5,outer=T)

dev.off()

	
			

# then starting priors

for(parno in checkpars)
{
	pdf(file=paste("FigS5.SimExamplesPriors",parno,".pdf",sep=""), width=6.7, height=3.7, family="Times", pointsize=9)
	par(mfrow=c(2,5),mai=c(0.3,0.4,0.1,0.1),oma=c(1,1.5,2.2,0.5))	
	counter=0			

	optact=optactList[[parno]]
	info=infoList[[parno]]

	for(prior in c(0.1,0.9))
	{
		plotdat <- get(paste("simdat","ind",parno,0,prior,"0.0",0,sep="_"))
		
		for(i in sample(1:100,5))
		{
			counter=counter+1
			image(matrix(0,80,51),axes=F,cex.main=1.5,xlim=c(0,1),ylim=c(0,1),zlim=c(0,0),col=bg.col)
			par(new=TRUE)
			image(t(optact),axes=F,cex.main=1.5,xlim=c(0,1),ylim=c(0,1),zlim=try.zlim,col=try.col,xlab="",ylab="",cex.lab=1)
			box()
			axis(1,seq(0,1,length.out=9),cex=1)  #,lab=letters[1:length(axTicks(1))])
			axis(2, at=seq(0,1,length.out=9), labels=seq(0,80,length.out=9),cex.lab=1)
					
			plotdf0 <- plotdat[plotdat$simno==1 & plotdat$it==i & plotdat$env==0,]
			plotdf1 <- plotdat[plotdat$simno==1 & plotdat$it==i & plotdat$env==1,]
				
			points((plotdf0$belief),(plotdf0$mass/80 - (1/80)),
					col="black",lwd=3,type="o",lty=1,pch=16,cex=0.8)  # or if want different lty, lty=j+2
					
			points((plotdf1$belief),(plotdf1$mass/80 - (1/80)),
					col="grey",lwd=3,type="o",lty=1,pch=16,cex=0.8)  # or if want different lty, lty=j+2
					
			mtext(paste("(",letters[counter],")",sep=""),side=3, line=0.8, adj=-0.25,cex=plot.cex)
	
		}

	}	

	mtext(expression(Size~italic(S)),side=2,line=0,adj=0.5,outer=T)
	mtext(expression(Belief~italic(P)),side=1,line=0,adj=0.5,outer=T)

	dev.off()
}

# then experiences
for(parno in 2:3)
{
	for(pfood in c(0.1,0.9)) 
	
	{
		if(pfood==0.1) pdf(file=paste("FigS6.SimExamplesBadExperiences",parno,".pdf",sep=""), width=6.7, height=3.7, family="Times", pointsize=9)
		if(pfood==0.9) pdf(file=paste("FigS6.SimExamplesGoodExperiences",parno,".pdf",sep=""), width=6.7, height=3.7, family="Times", pointsize=9)

		par(mfrow=c(2,5),mai=c(0.3,0.4,0.1,0.1),oma=c(1,1.5,2.2,0.5))	
		counter=0			
	
		optact=optactList[[parno]]
		info=infoList[[parno]]
	
		for(age in c(0,50))
		{
			plotdat <- get(paste("simdat","ind",parno,0,0.5,pfood,age,sep="_"))
		
			for(i in sample(1:100,5))
			{
				counter=counter+1
				image(matrix(0,80,51),axes=F,cex.main=1.5,xlim=c(0,1),ylim=c(0,1),zlim=c(0,0),col=bg.col)
				par(new=TRUE)
				image(t(optact),axes=F,cex.main=1.5,xlim=c(0,1),ylim=c(0,1),zlim=try.zlim,col=try.col,xlab="",ylab="",cex.lab=1)
				box()
				axis(1,seq(0,1,length.out=9),cex=1)  #,lab=letters[1:length(axTicks(1))])
				axis(2, at=seq(0,1,length.out=9), labels=seq(0,80,length.out=9),cex.lab=1)
						
				plotdf0 <- plotdat[plotdat$simno==1 & plotdat$it==i & plotdat$env==0,]
				plotdf1 <- plotdat[plotdat$simno==1 & plotdat$it==i & plotdat$env==1,]
					
				points((plotdf0$belief),(plotdf0$mass/80 - (1/80)),
						col="black",lwd=3,type="o",lty=1,pch=16,cex=0.8)  # or if want different lty, lty=j+2
						
				points((plotdf1$belief),(plotdf1$mass/80 - (1/80)),
						col="grey",lwd=3,type="o",lty=1,pch=16,cex=0.8)  # or if want different lty, lty=j+2
						
				points((plotdf0$belief[plotdf0$age %in% age:(age+10)]),(plotdf0$mass[plotdf0$age %in% age:(age+10)]/80 - (1/80)),
						col="blue",lwd=3,type="o",lty=1,pch=16,cex=0.8)  # or if want different lty, lty=j+2
						
				points((plotdf1$belief[plotdf1$age %in% age:(age+10)]),(plotdf1$mass[plotdf1$age %in% age:(age+10)]/80 - (1/80)),
						col="blue",lwd=3,type="o",lty=1,pch=16,cex=0.8)  # or if want different lty, lty=j+2
					
				mtext(paste("(",letters[counter],")",sep=""),side=3, line=0.8, adj=-0.25,cex=plot.cex)	
			}	

		}
		

	mtext(expression(Size~italic(S)),side=2,line=0,adj=0.5,outer=T)
	mtext(expression(Belief~italic(P)),side=1,line=0,adj=0.5,outer=T)

	dev.off()
	}
}













