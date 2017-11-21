# Write code that will accomodate 2 or 3 species and all range shift scenarios 
#
##
################################
#
# Load these libraries
################################

library(MASS)
library(fields)

################################
#
# Function definitions
################################
#R equivalent for matlab's meshgrid
meshgrid=function(a,b) {
  list(
       x=outer(b*0,a,FUN="+"),
       y=outer(b,a*0,FUN="+")
       )
} 


#Make a species' intrinsic range using a Gaussian distribution. 
#Variables correspond to: 
#	Fr	Max height of distribution (max reproduction)
#	pks	Vector with: 
#		 1. Initial spatial mean (position)
#		 2. pk.ends: Ending spatial mean
#		 3. by: Rate at which mean moves 
#		If not a vector, it is assumed stationary (no change in position) 
#
#	Ds	Vector with:
#		 1. Spatial variance (width)
#		 2. Ds.end: ending spatial variance
#		 3. Rate at which variance changes
#		If not a vector, it is assumed stationary (no change in position) 
#
#	stc	spatial array with space and time coordinates
#	burns	time period where distribution is stationary
#	ngens	time period of change for distribution

make.range=function(Fr, pks, Ds, stc, ngens, burns, fast.by = FALSE, dnwidth=NULL) { 


np=dim(stc)[2]
nt=dim(stc)[1]
Fr.tmp = matrix(0, nt, np)

#Test if the peak and spatial width parameters are a single value
if(length(pks)>1){pk1 = pks[1]} else { 
			pk1 = pks} 

if(length(Ds)>1){Ds1 = Ds[1]} else { 
			Ds1 = Ds} 

#Burn phase
X1=stc[1:burns,(1:(np)),2]
G1D=exp(-((X1-pk1)^2)/(2*Ds1))
G1D[G1D<1e-7]=0
G1D=G1D/(max(G1D))*Fr
Fr.tmp[1:burns,]=G1D

#Test if it is the fast or normal scenario
# Note: fast scenario inherently assumes changing peak position
# BUT NO change in width. 

if(fast.by==FALSE) {
	#Moving phase
	#Test if mean and variance change and create appropriate variables if
	#they do
	if(length(pks)>1){ pk.end = pks[2]
				pk.by = pks[3]
				peak.stc=seq(pk1, pk.end, by=pk.by)[1:ngens]} else {
				peak.stc=c(matrix(pk1, ngens,1))}
	 

	if(length(Ds)>1){  Ds.end = Ds[2]
				Ds.by = Ds[3]
				Ds.stc=seq(Ds1, Ds.end, by=Ds.by)[1:ngens]} else {
				Ds.stc=c(matrix(Ds1, ngens,1))}
	  

	X2=stc[(1+burns):nt,(1:(np)),2]
	G1D=exp(-((X2-peak.stc)^2)/(2*Ds.stc))
	G1D[G1D<1e-7]=0
	G1D=G1D/(max(G1D))*Fr
	Fr.tmp[(1+burns):nt,]=G1D
} else { 

	pk.end = pks[2]
	pk.by = pks[3]
	peak.stc1=seq(pk1, pk.end, by=pk.by)[1:dnwidth]
	peak.stc2=matrix(pk.end, (ngens-dnwidth),1)
	peak.stc=c(peak.stc1,peak.stc2)

	Ds.stc=c(matrix(Ds1, ngens,1))

	X2=stc[(1+burns):nt,(1:(np)),2]
	G1D=exp(-((X2-peak.stc)^2)/(2*Ds.stc))
	G1D[G1D<1e-7]=0
	G1D=G1D/(max(G1D))*Fr
	Fr.tmp[(1+burns):nt,]=G1D

}


return(Fr.tmp)

}

#These two functions can be used to find upper and lower edges
#of the range distributions. Needs the peak position, variance
#of the distribution, and a tolerance specifying number
#of standard deviations to include. 
get.upper = function (np, burns, ngens, pks, Ds, D.tol= 3) {
	
	ngenst=burns+ngens
	#Test if the peak and spatial width parameters are a single value
	#(i.e. stationary or moving)
	if(length(pks)>1){
		pk1 = pks[1]
		pk.end = pks[2]
		pk.by = pks[3]
		peak.stc.burn=c(matrix(pk1, ngens,1))
		peak.stc.gens=seq(pk1, pk.end, by=pk.by)[1:burns]
		peak.stc = c(peak.stc.burn, peak.stc.gens)} else {
		pk1 = pks
		peak.stc=c(matrix(pk1, ngenst,1))
	}


	if(length(Ds)>1){
		Ds1 = Ds[1]
		Ds.end = Ds[2]
		Ds.by = Ds[3]
		Ds.stc.burn=c(matrix(Ds1, ngens,1))
		Ds.stc.gens=seq(Ds1, Ds.end, by=Ds.by)[1:ngens]
		Ds.stc= c(Ds.stc.burn, Ds.stc.gens)
		} else { 
				
		Ds1 = Ds
		Ds.stc=c(matrix(Ds1, ngens,1))
	} 
	
		sds.up=round(np/2)+peak.stc+sqrt(Ds.stc)*D.tol
		return(sds.up)

}

get.lower = function (np, burns, ngens, pks, Ds, D.tol= 3) {
	
	ngenst=burns+ngens
	#Test if the peak and spatial width parameters are a single value
	#(i.e. stationary or moving)
	if(length(pks)>1){
		pk1 = pks[1]
		pk.end = pks[2]
		pk.by = pks[3]
		peak.stc.burn=c(matrix(pk1, ngens,1))
		peak.stc.gens=seq(pk1, pk.end, by=pk.by)[1:burns]
		peak.stc = c(peak.stc.burn, peak.stc.gens)} else {
		pk1 = pks
		peak.stc=c(matrix(pk1, ngenst,1))
	}


	if(length(Ds)>1){
		Ds1 = Ds[1]
		Ds.end = Ds[2]
		Ds.by = Ds[3]
		Ds.stc.burn=c(matrix(Ds1, ngens,1))
		Ds.stc.gens=seq(Ds1, Ds.end, by=Ds.by)[1:ngens]
		Ds.stc= c(Ds.stc.burn, Ds.stc.gens)

		} else { 
				
		Ds1 = Ds
		Ds.stc=c(matrix(Ds1, ngens,1))
	} 
	
		sds.low=round(np/2)+peak.stc-sqrt(Ds.stc)*D.tol
		return(sds.low)

}

#This function is to get the fitness-density covariance for a single
#1D array of sites. It is based on Snyder 2008, Snyder and Chesson 2004.
#It takes as variables: 
#	Fr.spp 	The intrinsic range of invader
#	cc.spp 	The competition from the resident community
#	l1.spp 	The spatial average of intrinsic reproduction
#	D.spp 	The denominator of the derivative in the Taylor expansion
#	gr1.spp The invasion growth rate WITHOUT the cov(f,d)
#	sr.spp 	The survival of the invader. 
get.fd.cov = function (Fr.spp, cc.spp, l1.spp, D.spp, gr1.spp, sr.spp, a_rr){ 

	rmFr1=mean(Fr.spp) 
	rmcc2=mean(cc.spp) 
	ei1=l1.spp/D.spp*(Fr.spp/rmFr1)
	uijmuj=l1.spp*0.8/D.spp^2*(cc.spp/rmcc2)
	uijmuj[uijmuj<0]=0
	e1=(ei1-uijmuj)
	s_zeta=spectrum(e1,spans=3,taper=0, plot=F, na.action=na.omit)
	ls1=length(s_zeta$spec)
	#Rebuild dispersal kernel to accomodate width
	xx_pad=matrix(seq(-ls1,ls1))
	kd.tmp=a_rr/2*exp(-a_rr*abs(xx_pad))
	kd.tmp=kd.tmp/(sum(kd.tmp))
	fdc2=fft(fft(kd.tmp[(ls1+2):(ls1+1+ls1)])*fft(s_zeta$spec[1:ls1],inverse=T)/(1-fft(kd.tmp[(ls1+2):(ls1+1+ls1)])), inverse=T)/(ls1)
	#Then multiply by average igr
	agr1=(gr1.spp-sr.spp)

	return(Re(fdc2[1])*agr1)


}

#This function is to get the Time Derivative of the fitness-density 
#covariance for a single 1D array of sites. Math is based on same
#references as above. Variables are the same, but Fr.spp and cc.spp
#must now include 2 entries in the time direction 
#(i.e. 2 columns instead of 1)

get.fdcov.diff = function (Fr.spp, cc.spp, l1.spp, D.spp, gr1.spp, sr.spp,a_rr){ 

	#Stationary terms
	rmFr1=mean(Fr.spp) 
	rmcc2=mean(cc.spp) 
	ei1=l1.spp/D.spp*(Fr.spp/rmFr1)
	uijmuj=l1.spp/D.spp^2*(cc.spp/rmcc2)
	uijmuj[uijmuj<0]=0

	#Time derivatives
	dtmu=(cc.spp[2,]/mean(cc.spp[2,])-1)-(cc.spp[1,]/mean(cc.spp[1,])-1)
	dte=(Fr.spp[2,]/mean(Fr.spp[2,])-1)-(Fr.spp[1,]/mean(Fr.spp[1,])-1)

	#Spectrum, coherency
	sz=spectrum(cbind((ei1-uijmuj),(dte-dtmu)),spans=3,taper=0, plot=F, na.action=na.omit)
	ls1=length(sz$spec[,1]) 

	#Rebuild dispersal kernel to accomodate width
	xx_pad=matrix(seq(-ls1,ls1))
	kd.tmp=a_rr/2*exp(-a_rr*abs(xx_pad))
	kd.tmp=kd.tmp/(sum(kd.tmp))
	fdcdt=fft(fft(kd.tmp[(ls1+2):(ls1+1+ls1)])*fft(s_zeta$spec[1:ls1],inverse=T)/(1-fft(kd.tmp[(ls1+2):(ls1+1+ls1)])), inverse=T)/(ls1)

	return(Re(fdcdt[1]))


}

#######
#This block of functions is specifically for the FAST scenario
#It includes get.fast.igr, get.spread.rate, and 

#Get the IGR for a single time step
get.fast.igr = function (Fr.inv,cc,sr,a_rr ){

	l1.fast=mean(Fr.inv)
	cr_mean.fast=mean(cc)
	D.fast= 1+cr_mean.fast[1]
	
	#Non-linear competitive variance and the spatial storage effect
	#These terms represent the perturbation terms from Snyder 2008 
	var_mu_Us.fast = var(cc/mean(cc)-1)
	cov_e_mu_Us.fast = cov(Fr.inv/mean(Fr.inv)-1, cc/mean(cc)-1)
	
	#Components of the IGR
	Elam1.fast=l1.fast*(1/D.fast+(cr_mean.fast*var_mu_Us.fast)/D.fast^3-(cr_mean.fast*cov_e_mu_Us.fast)/D.fast^2)+sr[s]-1
	Elam2.fast=Elam1.fast^2
	gr1.n.fast = exp(Elam1.fast-0.5*Elam2.fast)
	
	#The fitness-density covariance 
	tryCatch( {cov_lam_vc.fast=get.fd.cov(Fr.inv,cc,l1.fast,D.fast,gr1.n.fast, sr[s],a_rr[s])} , error=function(e){} )
	
	#The full IGR
	gr1.fast = gr1.n.fast+cov_lam_vc.fast

	return(gr1.fast)

	}

#Calculate the population spread rate from an IGR. This is based on
#the math in 

get.spread.rate = function(gr1.fast,a_rr,sr) {
	u=seq(0,.1,0.0001)
	#For Lalpacian kernel 
	r1=gr1.fast*1/(1-(1/a_rr[1])^2*u^2) 
	r2=sr[1]*1
	cs=1/u*log(r1+r2)
	cs_all=min(cs,na.rm=T)/10
	
	return(cs_all)
	}

#Calculate a new fitness distribution for a single time step. 
get.fast.peak = function(peak.new,Fr,Ds,cs_all,stc,np){
	peak.new=peak.new+cs_all 
	max.new=Fr[1]

	X2.new=stc[n,(1:(np)),2]

	G1D.new=exp(-((X2.new-peak.new)^2)/(2*Ds[1]));
	G1D.new[G1D.new<1e-7]=0
	G1D.new=G1D.new/(max(G1D.new))*max.new

	return(G1D.new)

	}

##########################
#
#
###########################

#For naming files
f.name1=c("igr_3spp_s1000_t10000b")

########################################################
#Tunable lattice parameters (space and time)
#
########################################################
nspp=3 #number of species
ns=1000 #Lattice width/height -- This should be even
ngens=5000 #Number of generations 
burns=5000 #Let the resident establish
burns1=5000 #Initiate climate change

#Make these automatically:
xx=matrix(seq((-ns/2),(ns/2)))
np=length(xx)
ngenst=burns1+ngens

########################################################
#Tunable Species parameters
#
########################################################
#Calculate spatiotemporal Fr (fundamental niche)
###Maximum reproduction rates, corresponding to spatio-temporal ideal conditions

Fr=matrix(c(2,2,2))  

#In "fast" scenario: 
#Fr = matrix (c(4,2))

###############
#Spatial distriubtions
###############
# 3 species with different pariwise scenarios
#	1. With only peaks moving (peaks != peaks.end), corresponds to rang shift
#	2. With only Ds moving (Ds != Ds.end), corresponds to range expansion
#	3. With both 1 and 2, corresponds to range shift and expansion
#	4. With "fast" bys uncommented, corresponds to slow expansion. 
###############
#For mid to upper invasion (invader is in back) 
#	# 1. peaks
#	peaks=c(-50,250)
#	peaks.end=c(250,250)
#
#	# 2. Ds
#	Ds = c(np^2/1600, np^2/1600)
#	#Ds.end = c(np^2/25, np^2/1600)
#	Ds.end = c(np^2/1600, np^2/1600)
#
#	# 3. bys
#	fast.bys = FALSE
#	peaks.by=c( (abs(peaks[1]-peaks.end[1])/ngens ), (abs(peaks[2]-peaks.end[2])/ngens) )
#	Ds.by = c( (abs(Ds[1]-Ds.end[1])/ngens), (abs(Ds[2]-Ds.end[2])/ngens) )
#
#	# 4. fast bys
#	#fast.bys = TRUE
#	#deltaE=5
#	#dnwidth=floor((np/2)/deltaE)
#	#Ds.by=c((abs(pk1.end-pk2)/dnwidth),(abs(Ds[2]-Ds.end[2])/ngens) ) 

#For all 3 species with full lower/mid and lower/upper interaction -- 
#same starting place in space for lower 2. 
	# 1. peaks
	peaks=c(-50, -50, 150)
	peaks.end=c(-50, -50, 150); f.name2=c("_notrans")
	#peaks.end=c(150, -50, 150); f.name2=c("_trans")


	# 2. Ds
	Ds = c(np^2/1600, np^2/1600,np^2/1600)
	Ds.end = c(np^2/50, np^2/1600,np^2/1600); f.name3=c("_var")
	#Ds.end = c(np^2/1600, np^2/1600,np^2/1600); f.name3=c("_novar")

	# 3. bys
	fast.bys = FALSE; f.name4=c("_norm")
	peaks.by=c( (abs(peaks[1]-peaks.end[1])/ngens ), (abs(peaks[2]-peaks.end[2])/ngens),(abs(peaks[3]-peaks.end[3])/ngens) )
	Ds.by = c( (abs(Ds[1]-Ds.end[1])/ngens), (abs(Ds[2]-Ds.end[2])/ngens),(abs(Ds[3]-Ds.end[3])/ngens) )

	# 4. fast bys
	#fast.bys = TRUE; f.name4=c("_fast")
	#deltaE=5
	#dnwidth=floor((np/2)/deltaE)
	#peaks.by=c((abs(peaks.end[1]-peaks[1])/dnwidth),(abs(Ds[2]-Ds.end[2])/ngens),(abs(Ds[3]-Ds.end[3])/ngens) )
	#Ds.by = c( (abs(Ds[1]-Ds.end[1])/ngens), (abs(Ds[2]-Ds.end[2])/ngens),(abs(Ds[3]-Ds.end[3])/ngens) )


###Survival
sr=matrix(c(0.9,0.9,0.9)) 

###Competition coeffcients 
alphas=matrix( c(1,0.8,0.8, 0.8,1,0.8,0.8,0.8,1),3,3)

###Competition distance
#b_rr=1/(100*np) #Essentially global
b_rr=c(1,1,1)

###Dispersal distance
#a_rr=1/(100*np) #essentially global
a_rr=c(1,1,1)

###########################################################
# 
# Internal variables: 1D 
###########################################################


#Make the full array of space-time coordinates (for 1D space)
#meshgrid will produce an R data frame
stc.temp=meshgrid(1:ngenst,xx)
#Convert the coordinates into an array for easier access
stc=array(c(matrix(stc.temp$x,ngenst,np,byrow=T),matrix(stc.temp$y,ngenst,np,byrow=T)),dim=c(ngenst,np,2)) 
#Note: in 1d, things are a bit odd. stc.temp$x ends up as time, and $y as space. So in stc, stc[,,1] is time
# and stc[,,2] is space. stc[,,1] is the same as the column vector 1:ngenst repeated over np columns. Then
# stc[,,2] is space as a row vector (xx) repeated over ngenst rows. 

####Dispersal kernels and their Fourier transforms 
kd=matrix(0,np,nspp)
fkd=matrix(0,np,nspp)
for( s in 1:nspp){ 
kd[,s] = a_rr[s]/2*exp(-a_rr[s]*abs(xx))
kd[,s]=kd[,s]/(sum(kd[,s]))
fkd[,s]=fft(kd[,s])#/(np+1)
}

####Competition kernels and their Fourier transforms 
kc=matrix(0,np,nspp)
fkc=matrix(0,np,nspp)
for( s in 1:nspp){ 
kc[,s] = b_rr[s]/2*exp(-b_rr[s]*abs(xx))
kc[,s]=kc[,s]/(sum(kc[,s]))
fkc[,s]=fft(kc[,s])#/(np+1)
}

####Intrinsic ranges
Frs=array(c(matrix(0.00,ngenst,np),matrix(0.00,ngenst,np)),dim=c(ngenst,np,nspp)) 
for( s in 1:nspp) { 

if ( sum(peaks[s]-peaks.end[s]) != 0 ) { 
pks = c( peaks[s], peaks.end[s], peaks.by[s]) } else {
pks = peaks[s]} 

if ( sum(Ds[s]-Ds.end[s]) != 0) { 
Ds.tmp = c( Ds[s], Ds.end[s], Ds.by[s]) } else {
Ds.tmp = Ds[s]} 

if (fast.bys==FALSE) {
Frs[,,s] = make.range(Fr[s], pks, Ds.tmp, stc, ngens, burns) } else {
	#Do this for fast moving species, here species 1  
	if( s == 1){
		Frs[,,s] = make.range(Fr[s], pks, Ds.tmp, stc, ngens, burns, fast.bys,dnwidth)
	} else {
		#But the normal way for other species:
		Frs[,,s] = make.range(Fr[s], pks, Ds.tmp, stc, ngens, burns)
	}



	}
}

#############
#Map out upper and lower boundaries, set range for IGR measurements
#############
up1=matrix(0,ngenst,nspp)
low1=matrix(0,ngenst,nspp)

for(sa in 1:nspp){
if ( sum(peaks[sa]-peaks.end[sa]) != 0 ) { 
pks = c( peaks[sa], peaks.end[sa], peaks.by[sa]) } else {
pks = peaks[sa]} 

if ( sum(Ds[sa]-Ds.end[sa]) != 0) { 
Ds.tmp = c( Ds[sa], Ds.end[sa], Ds.by[sa]) } else {
Ds.tmp = Ds[sa]} 

up1[,sa ]=get.upper(np, burns,ngens, pks, Ds.tmp)
low1[,sa ]=get.lower(np, burns,ngens, pks, Ds.tmp)

}

low.tmp=min(round(low1),na.rm=T)
up.tmp=max(round(up1),na.rm=T)
#if (low.tmp<0){low.tmp=0}
#if (up.tmp>np){up.tmp=np}
low.tmp=0
up.tmp=np
low=matrix(low.tmp,ngenst,1)
up=matrix(up.tmp,ngenst,1)
wdth=round(up.tmp)-round(low.tmp)
sd.lim=min(which(up>np))
if(is.infinite(sd.lim)){sd.lim=ngenst}

###################################################################

###################################################################
# Key variables for output: the invasion growth rates, their components
# and the time derivatives of igr. 

#Components of IGR
l1=matrix(0,ngenst,nspp)
y=matrix(0,ngenst,nspp)
cr_mean=matrix(0,ngenst,nspp)
D = matrix(0,ngenst,nspp)
var_mu_Us=matrix(0,ngenst,nspp) 
cov_e_mu_Us=matrix(0,ngenst,nspp)
cov_lam_vc=matrix(0,ngenst,nspp)
Elam1=matrix(0,ngenst,nspp) 
Elam2=matrix(0,ngenst,nspp)
gr1.n=matrix(0,ngenst,nspp)
gr1=matrix(0,ngenst,nspp)

#Time derivatives
dt_l1=matrix(0,ngenst,nspp)
dt_y=matrix(0,ngenst,nspp)
dt_cr_mean=matrix(0,ngenst,nspp)
dt_D = matrix(0,ngenst,nspp)
cov_dtmu_mu_Us=matrix(0,ngenst,nspp)
cov_e_dtmu_Us=matrix(0,ngenst,nspp)
cov_dte_mu_Us=matrix(0,ngenst,nspp)
cov_lam_vcdt1=matrix(0,ngenst,nspp)
dtElam1=matrix(0,ngenst,nspp)
dtElam2=matrix(0,ngenst,nspp)
dtgr1.n=matrix(0,ngenst,nspp)
dtgr1=matrix(0,ngenst,nspp)

###################################################################

###################################################################
#These are for each of the pairwise scenarios. Same variables, but for 
#1 vs 2, 1 vs 3, etc. 

#Components of IGR
l1.p=matrix(0,ngenst,nspp*2)
y.p=matrix(0,ngenst,nspp*2)
cr_mean.p=matrix(0,ngenst,nspp*2)
D.p = matrix(0,ngenst,nspp*2)
var_mu_Us.p=matrix(0,ngenst,nspp*2) 
cov_e_mu_Us.p=matrix(0,ngenst,nspp*2)
cov_lam_vc.p=matrix(0,ngenst,nspp*2)
Elam1.p=matrix(0,ngenst,nspp*2) 
Elam2.p=matrix(0,ngenst,nspp*2)
gr1.n.p=matrix(0,ngenst,nspp*2)
gr1.p=matrix(0,ngenst,nspp*2)

#Time derivatives
dt_l1.p=matrix(0,ngenst,nspp*2)
dt_y.p=matrix(0,ngenst,nspp*2)
dt_cr_mean.p=matrix(0,ngenst,nspp*2)
dt_D.p = matrix(0,ngenst,nspp*2)
cov_dtmu_mu_Us.p=matrix(0,ngenst,nspp*2)
cov_e_dtmu_Us.p=matrix(0,ngenst,nspp*2)
cov_dte_mu_Us.p=matrix(0,ngenst,nspp*2)
cov_lam_vcdt1.p=matrix(0,ngenst,nspp*2)
dtElam1.p=matrix(0,ngenst,nspp*2)
dtElam2.p=matrix(0,ngenst,nspp*2)
dtgr1.n.p=matrix(0,ngenst,nspp*2)
dtgr1.p=matrix(0,ngenst,nspp*2)

#Outermost loop: treat each species as invader

#Keep track of total interspecific competition of species when it is invader
cc=array(c(matrix(0,ngenst,np),matrix(0,ngenst,np)),dim=c(ngenst,np,nspp)) 

for ( s in 1:nspp) { 

s.index= 1:nspp
s.index = s.index[-s]

#The population matrix. Columns represent space, rows are time. 
nr=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 
#Initialize residents
for ( sa in s.index) { nr[1,,sa] = matrix(0.01,1,np)}

#Equilibrium of residents
nr2_eq_all = array(c(matrix(0,ngenst,np),matrix(0,ngenst,np)),dim=c(ngenst,np,(nspp))) 

###############################################
#Population dynamics 
#
###############################################
#
#
# Burn first

for (n in 1:(burns)) {
#Seedlings germinate and produce seeds at rate Fr, weighted by competition
fp = nr[n,,] 
for (sa in 1:nspp){ fp[,sa] = fft(fp[,sa]) } 

#Competition from each pairwise combination of species
cr=matrix(0,np,nspp^2)
cr_tot = matrix(0,np,nspp)
for (sa in 1:nspp) { 
	for (sb in 1:nspp) {
			cr[,(sa*(sa-1)+sb)] = Re(fft((alphas[sa,sb]*fp[,sb]*fkc[,sb]),inverse=T)/(np+1))
			cr[,(sa*(sa-1)+sb)] = c(cr[(ceiling(np/2):(np)),(sa*(sa-1)+sb)], cr[(1:floor(np/2)),(sa*(sa-1)+sb)])
			cr_tot[,sa] = cr_tot[,sa] + cr[ ,(sa*(sa-1)+sb)]	
	}
}

#Store this value of competition for the invader to calculate invasion growth rates
cc[n,,s] = cr_tot[,s]

#competition-weighted seed production
lam = matrix(0,np,nspp)

#Seeds disperse
fd = matrix(0,np,nspp)
nr_disp = matrix(0,np,nspp)

for(sa in 1:nspp) { 
	#competition-weighted seed production
	lam[,sa]=(Frs[n,,sa])/(1+cr_tot[,sa])
	lam[,sa]= lam[,sa]*(nr[n,,sa]>1e-34)
	#If both species have zero germination, then NAs are produced: 
	lam[,sa][is.na(lam[,sa]) ] = 0

	#Seeds disperse
	fd[,sa]=fft(lam[,sa]*nr[n,,sa])
	nr_disp[,sa] = Re(fft( (fd[,sa]*fkd[,sa]), inverse=T)/(np+1))
	nr_disp[,sa] = c(nr_disp[(ceiling(np/2):(np)),sa],nr_disp[(1:floor(np/2)),sa])

	#Final population step
	nr[(n+1),,sa] = nr_disp[,sa] + nr[n,,sa]*sr[sa]

}

#This is the solved version of each resident's equilibrium density when it is alone. 
for( t in 1:(length(s.index)) ) {
sa=s.index[t]
nr2_eq_all[n,,sa]=Re(fft(fkc[,sa]*fkd[,sa]*fft(Frs[n,,sa])/(1-sr[sa])-1,inverse=T)/(np+0.2*np))
}

}

########################################
#Invasion -- Do not do this if the goal is to calculate IGRs! 
########################################
#nr[(burns), ,s]= 1e-6
########################################
#If "FAST" scenario is chosen
if (fast.bys==TRUE){
		#Remake the range of fast species
		if ( sum(peaks[1]-peaks.end[1]) != 0 ) { 
		pks = c( peaks[1], peaks.end[1], peaks.by[1]) } else {
		pks = peaks[1]} 

		if ( sum(Ds[1]-Ds.end[1]) != 0) { 
		Ds.tmp = c( Ds[1], Ds.end[1], Ds.by[1]) } else {
		Ds.tmp = Ds[1]} 

		Frs[,,1] = make.range(Fr[1], pks, Ds.tmp, stc, ngens, burns, fast.bys,dnwidth)
		#Declare new variables for range 
		peak.new=peaks[1]
		Fr.inv=Frs[,,1]		
}


for (n in burns:(ngenst-1)) {
#Seedlings germinate and produce seeds at rate Fr, weighted by competition
fp = nr[n,,] 
for (sa in 1:nspp){ fp[,sa] = fft(fp[,sa]) } 

#Competition from each pairwise combination of species
cr=matrix(0,np,nspp^2)
cr_tot = matrix(0,np,nspp)
for (sa in 1:nspp) { 
	for (sb in 1:nspp) {
			cr[,(sa*(sa-1)+sb)] = Re(fft((alphas[sa,sb]*fp[,sb]*fkc[,sb]),inverse=T)/(np+1))
			cr[,(sa*(sa-1)+sb)] = c(cr[(ceiling(np/2):(np)),(sa*(sa-1)+sb)], cr[(1:floor(np/2)),(sa*(sa-1)+sb)])
			cr_tot[,sa] = cr_tot[,sa] + cr[ ,(sa*(sa-1)+sb)]	
	}
}

#Store this value of competition for the invader to calculate invasion growth rates
cc[n,,s] = cr_tot[,s]

#competition-weighted seed production
lam = matrix(0,np,nspp)

#Seeds disperse
fd = matrix(0,np,nspp)
nr_disp = matrix(0,np,nspp)


for(sa in 1:nspp) { 
	#competition-weighted seed production
	lam[,sa]=(Frs[n,,sa])/(1+cr_tot[,sa])
	lam[,sa]= lam[,sa]*(nr[n,,sa]>1e-34)
	#If both species have zero germination, then NAs are produced: 
	lam[,sa][is.na(lam[,sa]) ] = 0

	#Seeds disperse
	fd[,sa]=fft(lam[,sa]*nr[n,,sa])
	nr_disp[,sa] = Re(fft( (fd[,sa]*fkd[,sa]), inverse=T)/(np+1))
	nr_disp[,sa] = c(nr_disp[(ceiling(np/2):(np)),sa],nr_disp[(1:floor(np/2)),sa])

	#Final population step
	nr[(n+1),,sa] = nr_disp[,sa] + nr[n,,sa]*sr[sa]

}

#This is the solved version of each resident's equilibrium density when it is alone. 
for( t in 1:(length(s.index)) ) {
sa=s.index[t]
nr2_eq_all[n,,sa]=Re(fft(fkc[,sa]*fkd[,sa]*fft(Frs[n,,sa])/(1-sr[sa])-1,inverse=T)/(np+0.2*np))
}

#If the "FAST" scenario is chosen, need to do a special set of 
#calculations of the invader.
if(fast.bys==TRUE){   

	#Calculate the IGR for species 1 for only this time step
	gr1.fast=get.fast.igr(Fr.inv[n,low[n]:up[n]],cc[n,low[n]:up[n],1],sr,a_rr )
	#Now use the IGR to calculate the spread rate
	cs_all = get.spread.rate(gr1.fast,a_rr,sr)	

	# Make the new intrinsic fitness distribution for the next timestep
	# First, make a "pretend" fitness distribution based on how far 
	# the population of the invader can spread in a single time step
	Fr.inv[n+1,] = get.fast.peak(peak.new,Fr,Ds,cs_all,stc,np)
	# Then filter the actual intrinsic range according to the amount 
	# of overlap
	Frs[n+1,,1] = Frs[n+1,,1]*as.numeric(Fr.inv[n+1,]>1e-2) 

	}


}


###############################################
# Invasion growth rates -- using simulation data
#
###############################################

#These are all of the necessary terms to calculate the standard invasion growth rates

#Calculate standard spatial terms
for (t in 1:sd.lim) { 
	
	l1[t,s]=mean(Frs[t,low[t]:up[t],s])
	y[t,s]=mean(nr[t,low[t]:up[t],s])
	cr_mean[t,s]=mean(cc[t,low[t]:up[t],s])
	D[t,s]= 1+cr_mean[t,s]
	
	#Non-linear competitive variance and the spatial storage effect
	#These terms represent the perturbation terms from Snyder 2008 
	var_mu_Us[t,s] = var(cc[t,low[t]:up[t],s]/mean(cc[t,low[t]:up[t],s])-1)
	cov_e_mu_Us[t,s] = cov(Frs[t,low[t]:up[t],s]/mean(Frs[t,low[t]:up[t],s])-1, cc[t,low[t]:up[t],s]/mean(cc[t,low[t]:up[t],s])-1)
	
	#Components of the IGR
	Elam1[t,s]=l1[t,s]*(1/D[t,s]+(cr_mean[t,s]*var_mu_Us[t,s])/D[t,s]^3-(cr_mean[t,s]*cov_e_mu_Us[t,s])/D[t,s]^2)+sr[s]-1
	Elam2.p[t,s]=(l1.p[t,s]*(1/D.p[t,s]) +sr[s]-1)^2+2*(l1.p[t,s]*(1/D.p[t,s]) +sr[s]-1)*
			(var_mu_Us.p[t,s]-cov_e_mu_Us.p[t,s])/D.p[t,s]^4
	gr1.n[t,s] = exp(Elam1[t,s]-0.5*Elam2[t,s])
	
	#The fitness-density covariance 
	tryCatch( {cov_lam_vc[t,s]=get.fd.cov(Frs[t,low[t]:up[t],s],cc[t,low[t]:up[t],s],l1[t,s],D[t,s],gr1.n[t,s], sr[s],a_rr[s])} , error=function(e){} )
	
	#The full IGR
	gr1[t,s] = gr1.n[t,s]+cov_lam_vc[t,s]
} #End spatial terms


#Calculate time derivatives

for (t in 2:sd.lim){ 

	dt_l1[t,s]=l1[t,s]-l1[t-1,s]
	dt_y[t,s] = y[t,s]-y[t-1,s]
	dt_cr_mean[t,s]=cr_mean[t,s] - cr_mean[t-1,s]
	dt_D[t,s] = D[t,s] - D[t-1,s]

	#Time derivatives of competitive variance and SSE
	cc.diff=(cc[t,low[t]:up[t],s]/mean(cc[t,low[t]:up[t],s])-1)-(cc[(t-1),low[(t-1)]:up[(t-1)],s]/mean(cc[(t-1),low[(t-1)]:up[(t-1)],s])-1)
	Frs.diff=(Frs[t,low[t]:up[t],s]/mean(Frs[t,low[t]:up[t],s])-1)-(Frs[(t-1),low[(t-1)]:up[(t-1)],s]/mean(Frs[(t-1),low[(t-1)]:up[(t-1)],s])-1)
	cov_dtmu_mu_Us[t,s]= cov(cc[t,low[t]:up[t],s]/mean(cc[t,low[t]:up[t],s])-1,cc.diff)
	cov_e_dtmu_Us[t,s]=cov( (Frs[t,low[t]:up[t],s]/mean(Frs[t,low[t]:up[t],s])-1), cc.diff )
	cov_dte_mu_Us[t,s]=cov( Frs.diff, (cc[t,low[t]:up[t],s]/mean(cc[t,low[t]:up[t],s])-1) )

	dtElam1[t,s]=dt_l1[t,s]*(1/D[t,s]+(cr_mean[t,s]^2*var_mu_Us[t,s])/D[t,s]^3-(cr_mean[t,s]*cov_e_mu_Us[t,s])/D[t,s]^2)+(l1[t,s]*(2*cr_mean[t,s]^2*cov_dtmu_mu_Us[t,s])/D[t,s]^3-l1[t,s]*cr_mean[t,s]*(cov_e_dtmu_Us[t,s]+cov_dte_mu_Us[t,s])/D[t,s]^2)
	dtElam2[t,s]=2*Elam1[t,s]*dtElam1[t,s]
	dtgr1.n[t,s]=dtElam1[t,s]-0.5*dtElam2[t,s]

	#Time derivative of fitness-density covariance 
	#cov_lam_vcdt1=2*get.fdcov.diff(Frs[((t-1):t),low[t]:up[t],s],cc[((t-1):t),low[t]:up[t],s], l1[t],D[t],gr1.n[t], sr[s],a_rr[s] )
	#cov_lam_vcdt2=get.fdcov.diff( )



}#End time derivatives


###############################################
# Invasion growth rates for pairwise cases!  -- using simulation data
#
###############################################

#Outer loop through residents 

for (rs in 1:(nspp-1)){
#These are all of the necessary terms to calculate the standard invasion growth rates

enter1 = 2*(s-1)+rs
rez.no = s.index[rs]

#Calculate standard spatial terms
for (t in 1:sd.lim) { 
	
	cc.p=Re(fft((alphas[s,rez.no]*fft(nr2_eq_all[t,,rez.no])*fkc[,1]),inverse=T)/(np+1))
	cc.p=c(cc.p[ceiling(np/2):(np)],cc.p[1:floor(np/2)])
	cc.p=cc.p[low[t]:up[t]]
	l1.p[t,enter1]=mean(Frs[t,low[t]:up[t],s])
	y.p[t,enter1]=mean(nr2_eq_all[t,,rez.no])
	cr_mean.p[t,enter1]=mean(cc.p)
	D.p[t,enter1]= 1+cr_mean.p[t,enter1]
	
	#Non-linear competitive variance and the spatial storage effect
	#These terms represent the perturbation terms from Snyder 2008 
	var_mu_Us.p[t,enter1] = var(cc.p/mean(cc.p)-1)
	cov_e_mu_Us.p[t,enter1] = cov(Frs[t,low[t]:up[t],s]/mean(Frs[t,low[t]:up[t],s])-1, cc.p/mean(cc.p)-1)
	
	#Components of the IGR
	Elam1.p[t,enter1]=l1.p[t,enter1]*(1/D.p[t,enter1]+(cr_mean.p[t,enter1]*var_mu_Us.p[t,enter1])/D.p[t,enter1]^3-(cr_mean.p[t,enter1]*cov_e_mu_Us.p[t,enter1])/D.p[t,enter1]^2)+sr[s]-1
	Elam2.p[t,enter1]=(l1.p[t,enter1]*(1/D.p[t,enter1]) +sr[s]-1)^2+2*(l1.p[t,enter1]*(1/D.p[t,enter1]) +sr[s]-1)*
			(var_mu_Us.p[t,enter1]-cov_e_mu_Us.p[t,enter1])/D.p[t,enter1]^4
	gr1.n.p[t,enter1] = exp(Elam1.p[t,enter1]-0.5*Elam2.p[t,enter1])
	
	#The fitness-density covariance 
	tryCatch( {cov_lam_vc.p[t,enter1]=get.fd.cov(Frs[t,low[t]:up[t],s], cc.p, l1.p[t,enter1], D.p[t,enter1], gr1.n.p[t,enter1], sr[s], a_rr[s])} , error=function(e){} )
	
	#The full IGR
	gr1.p[t,enter1] = gr1.n.p[t,enter1]+cov_lam_vc.p[t,enter1]
} #End spatial terms


#Calculate time derivatives

for (t in 2:sd.lim){ 
	
	cc.p=Re(fft((alphas[s,rez.no]*fft(nr2_eq_all[t,,rez.no])*fkc[,1]),inverse=T)/(np+1))
	cc.p=c(cc.p[ceiling(np/2):(np)],cc.p[1:floor(np/2)])
	cc.p=cc.p[low[t]:up[t]]
	cc.pt1=Re(fft((alphas[s,rez.no]*fft(nr2_eq_all[(t-1),,rez.no])*fkc[,1]),inverse=T)/(np+1))
	cc.pt1=c(cc.p[ceiling(np/2):(np)],cc.p[1:floor(np/2)])
	cc.pt1=cc.p[low[(t-1)]:up[(t-1)]]
	dt_l1.p[t,enter1]=l1.p[t,enter1]-l1.p[t-1,enter1]
	dt_y.p[t,enter1] = y.p[t,enter1]-y.p[t-1,enter1]
	dt_cr_mean.p[t,enter1]=cr_mean.p[t,enter1] - cr_mean.p[t-1,enter1]
	dt_D.p[t,enter1] = D.p[t,enter1] - D.p[t-1,enter1]

	#Time derivatives of competitive variance and SSE
	cc.diff=(cc.p/mean(cc.p)-1)-(cc.pt1/mean(cc.pt1)-1)
	Frs.diff=(Frs[t,low[t]:up[t],s]/mean(Frs[t,low[t]:up[t],s])-1)-(Frs[(t-1),low[(t-1)]:up[(t-1)],s]/mean(Frs[(t-1),low[(t-1)]:up[(t-1)],s])-1)
	cov_dtmu_mu_Us.p[t,enter1]= cov(cc.p/mean(cc.p)-1,cc.diff)
	cov_e_dtmu_Us.p[t,enter1]=cov( (Frs[t,low[t]:up[t],s]/mean(Frs[t,low[t]:up[t],s])-1), cc.diff )
	cov_dte_mu_Us.p[t,enter1]=cov( Frs.diff, (cc[t,low[t]:up[t],s]/mean(cc[t,low[t]:up[t],s])-1) )

	dtElam1.p[t,enter1]=dt_l1.p[t,s]*(1/D.p[t,enter1]+(cr_mean.p[t,enter1]^2*var_mu_Us.p[t,enter1])/D.p[t,enter1]^3-(cr_mean.p[t,enter1]*cov_e_mu_Us.p[t,enter1])/D.p[t,enter1]^2)+(l1.p[t,enter1]*(2*cr_mean.p[t,enter1]^2*cov_dtmu_mu_Us.p[t,enter1])/D.p[t,enter1]^3-l1.p[t,enter1]*cr_mean.p[t,enter1]*(cov_e_dtmu_Us.p[t,enter1]+cov_dte_mu_Us.p[t,enter1])/D.p[t,enter1]^2)
	dtElam2.p[t,enter1]=2*Elam1.p[t,enter1]*dtElam1.p[t,enter1]
	dtgr1.n.p[t,enter1]=dtElam1.p[t,enter1]-0.5*dtElam2.p[t,enter1]

	#Time derivative of fitness-density covariance 
	#cov_lam_vcdt1=2*get.fdcov.diff(Frs[((t-1):t),low[t]:up[t],s],cc[((t-1):t),low[t]:up[t],s], l1[t],D[t],gr1.n[t], sr[s],a_rr[s] )
	#cov_lam_vcdt2=get.fdcov.diff( )



}#End time derivatives

}#End outer species loop through residents

} #End outer loop through species

###########################################################
#Data saving
###########################################################

#Save all of the meaningful invasion growth rate variables
#in a file with an informative but awkwardly long file name. 
#This includes information about the temporal and spatial extent,
#and the details of the scenario type. 

file.name = (paste(f.name1,f.name2,f.name3,f.name4, ".var",sep=""))
save(file=file.name, "l1", "D", "var_mu_Us","cov_e_mu_Us",
"cov_lam_vc", "Elam1", "gr1.n", "gr1",
"dt_l1","dt_cr_mean", "dt_D", "cov_dtmu_mu_Us","cov_e_dtmu_Us","cov_dte_mu_Us",
"cov_lam_vcdt1", "dtElam1","dtElam2","dtgr1.n", "dtgr1", "l1.p",
"y.p", "cr_mean.p", "D.p", "var_mu_Us.p", "cov_e_mu_Us.p","cov_lam_vc.p", 
"Elam1.p","Elam2.p", "gr1.n.p", "gr1.p", "dt_l1.p", "dt_y.p" ,"dt_cr_mean.p",
"dt_D.p", "cov_dtmu_mu_Us.p","cov_e_dtmu_Us.p", "cov_dte_mu_Us.p", "cov_lam_vcdt1.p",
"dtElam1.p", "dtElam2.p", "dtgr1.n.p", "dtgr1.p"
)
