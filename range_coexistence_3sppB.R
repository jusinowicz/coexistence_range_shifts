# Write code that will accomodate 2 or 3 species and all range shift scenarios 
#
##
#=============================================================================
# Load these libraries
#=============================================================================
library(MASS)
library(fields)
source("./range_coexistence_functionsB.R")

#=============================================================================
#For naming files
#=============================================================================
f.name1=c("igr_3spp_s1000_t10000B")

#=============================================================================
#Tunable lattice parameters (space and time)
#=============================================================================
nspp=3 #number of species
ns=1000 #Lattice width/height -- This should be even
ngens=4 #Number of generations for environmental change
iconfig=1 #NUmber of interattions of initial configuration of species

#Make these automatically:
xx=matrix(seq((-ns/2),(ns/2)))
np=length(xx)
ngenst=iconfig+ngens

#=============================================================================
#Tunable Species parameters
#=============================================================================
#Calculate spatiotemporal Fr (fundamental niche)
###Maximum reproduction rates, corresponding to spatio-temporal ideal conditions

Fr=matrix(c(2,2,2))  

#In "fast" scenario: 
#Fr = matrix (c(4,2))

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

	# # 4. fast bys
	# fast.bys = TRUE; f.name4=c("_fast")
	# deltaE=5
	# dnwidth=floor((np/2)/deltaE)
	# peaks.by=c((abs(peaks.end[1]-peaks[1])/dnwidth),(abs(Ds[2]-Ds.end[2])/ngens),(abs(Ds[3]-Ds.end[3])/ngens) )
	# Ds.by = c( (abs(Ds[1]-Ds.end[1])/ngens), (abs(Ds[2]-Ds.end[2])/ngens),(abs(Ds[3]-Ds.end[3])/ngens) )


#=============================================================================
# Internal variables: 1D 
#=============================================================================


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
Frs[,,s] = make.range(Fr[s], pks, Ds.tmp, stc, ngens, iconfig) } else {
	#Do this for fast moving species, here species 1  
	if( s == 1){
		Frs[,,s] = make.range(Fr[s], pks, Ds.tmp, stc, ngens, iconfig, fast.bys,dnwidth)
	} else {
		#But the normal way for other species:
		Frs[,,s] = make.range(Fr[s], pks, Ds.tmp, stc, ngens, iconfig)
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

up1[,sa ]=get.upper(np, iconfig,ngens, pks, Ds.tmp)
low1[,sa ]=get.lower(np, iconfig,ngens, pks, Ds.tmp)

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

#=============================================================================
# Key variables for output: the invasion growth rates, their components
#=============================================================================
#Simulation IGR
ldg.sim = matrix(0,ngenst, nspp)
ldgIMP.sim = matrix(0,ngenst, nspp)
covIMP.sim = matrix(0,ngenst, nspp)

#Components of IGR
l1=matrix(0,ngenst,nspp)
D = matrix(0,ngenst,nspp)
var_mu_Us=matrix(0,ngenst,nspp) 
cov_e_mu_Us=matrix(0,ngenst,nspp)
cov_lam_vc=matrix(0,ngenst,nspp)
cov_lam_vc2=matrix(0,ngenst,nspp)
Elam1=matrix(0,ngenst,nspp) 
Elam2=matrix(0,ngenst,nspp)
gr1.n=matrix(0,ngenst,nspp)
gr1=matrix(0,ngenst,nspp)

#Components kept as lists
y.full=NULL
w.eq.full=NULL

###################################################################
#These are for each of the pairwise scenarios. Same variables, but for 
#1 vs 2, 1 vs 3, etc. 

#Components of IGR
l1.p=matrix(0,ngenst,nspp*(nspp-1))
y.p=matrix(0,ngenst,nspp*(nspp-1))
w.eq.p=matrix(0,ngenst,nspp*(nspp-1))
D.p = matrix(0,ngenst,nspp*(nspp-1))
var_mu_Us.p=matrix(0,ngenst,nspp*(nspp-1)) 
cov_e_mu_Us.p=matrix(0,ngenst,nspp*(nspp-1))
cov_lam_vc.p=matrix(0,ngenst,nspp*(nspp-1))
cov_lam_vc.p2=matrix(0,ngenst,nspp*(nspp-1))
Elam1.p=matrix(0,ngenst,nspp*(nspp-1)) 
Elam2.p=matrix(0,ngenst,nspp*(nspp-1))
gr1.n.p=matrix(0,ngenst,nspp*(nspp-1))
gr1.p=matrix(0,ngenst,nspp*(nspp-1))

#=============================================================================
#Outermost loop: treat each species as invader
#=============================================================================

#Keep track of total interspecific competition of species when it is invader
cc=array(c(matrix(0,ngenst,np),matrix(0,ngenst,np)),dim=c(ngenst,np,nspp)) 

for ( s in 1:nspp) { 

	s.index= 1:nspp
	s.index = s.index[-s]

	#Equilibrium populations of residents. Columns represent space, rows are time. 
	nr=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 
	#Theoretical pairwise equilibrium of residents
	nr2_eq_all = array(c(matrix(0,ngenst,np),matrix(0,ngenst,np)),dim=c(ngenst,np,(nspp))) 

	#Initialize residents
	for ( sa in s.index) { nr[1,,sa] = matrix(0.01,1,np)}

	
	#=============================================================================
	#Invader and resident stationary distributions
	#=============================================================================


	for( t in 1: ngenst){ 

		#If the "FAST" scenario is chosen, need to do a special set of 
		#calculations of the invader.
		if(fast.bys==TRUE ){   
			#Declare new variables for range 
	

			if (t==(iconfig)){
				#Remake the range of fast species
				if ( sum(peaks[1]-peaks.end[1]) != 0 ) { 
				pks = c( peaks[1], peaks.end[1], peaks.by[1]) } else {
				pks = peaks[1]} 

				if ( sum(Ds[1]-Ds.end[1]) != 0) { 
				Ds.tmp = c( Ds[1], Ds.end[1], Ds.by[1]) } else {
				Ds.tmp = Ds[1]} 
				
				Frs[,,1] = make.range(Fr[1], pks, Ds.tmp, stc, ngens, iconfig, fast.bys,dnwidth)
				#Declare new variables for range 
				peak.new=peaks[1]
				Fr.inv=Frs[,,1]		
			} 

			if (t >(iconfig)){
				#Calculate the IGR for species 1 for only this time step
				gr1.fast=get.fast.igr(Frs[(t-1),,],nr[(t-1),,], sr, alphas, kd,kc,n.inv=1 )
				#Now use the IGR to calculate the spread rate
				cs_all = get.spread.rate(gr1.fast,a_rr,sr)	

				# Make the new intrinsic fitness distribution for the next timestep
				# First, make a "pretend" fitness distribution based on how far 
				# the population of the invader can spread in a single time step
				Fr.inv[t,] = get.fast.peak(peak.new,Fr,Ds,cs_all,stc[(t-1),,],np)
				# Then filter the actual intrinsic range according to the amount 
				# of overlap
				Frs[t,,1] = Frs[t,,1]*as.numeric(Fr.inv[t,]>1e-2) 
			}
		}
		
		
		#Get the equilibrium of the resident community when invader is absent
		nr[t,,(s.index)] = get.res.eq(Frs[t,,],s.index,sr,alphas, kd,kc,fast=TRUE )
		
		#This is the solved version of each resident's equilibrium density when it is alone. 
		for( sa in 1:(length(s.index)) ) {
			sb=s.index[sa]
			nr2_eq_all[t,,sb] = get.rsd.pair(Frs[t,,sb],sr[sb], alphas[s,sb],kd[,sb],kc[,sb])
			}


		#Get the low-density equilibrium density of the invader against the resident community
		nr[t,,s] = get.inv.ldeq(Frs[t,,], nr[t,,], s, sr,alphas, kd,kc,fast=TRUE )
	

	print(t)
	} 

	#=============================================================================
	# Low-density growth rates -- using simulation data, all 3 spp
	#=============================================================================
	#Simple numerical LDG: use the invader low-density equilibrium density, 
	#then let it invade against the community: 
	
	for (t in 1:ngenst) { 
		inv.one =pop_lg(Frs[t,,],nr[t,,], sr, alphas, kd,kc )[,s]
		ldg.sim[t,s]=mean(inv.one)/mean(nr[t,,s])

	}

	#Calculate the spatially implicit portion and the fitness-density covariance
	#separately, still using simulations. This just uses a "global" dispersal kernel
	####Dispersal kernels and their Fourier transforms 
	inv.id =1e-5

	kdg=matrix(0,np,nspp)
	fkdg=matrix(0,np,nspp)
	for( sa in 1:nspp){ 
	kdg[,sa] = 1/np
	fkdg[,sa]=fft(kdg[,sa])#/(np+1)
	}
	
	nr2=nr
	#Make the invader spatially homogenous, at low density
	nr2[,,s] = matrix(inv.id,ngenst,np)


	for (t in 1:ngenst) { 
		inv.oneIMP =pop_lg(Frs[t,,], nr2[t,,], sr, alphas, kdg,kc )[,s]
		ldgIMP.sim[t,s]=mean(inv.oneIMP)/mean(nr2[t,,s])

	}

	covIMP.sim[,s] = ldg.sim[,s]-ldgIMP.sim[,s]

	#Calculate the LDG for species vs. community by component

	y = matrix(0,sd.lim, nspp-1)
	w.eq = matrix(0,sd.lim, nspp-1)

	#Calculate standard spatial terms
	for (t in 1:ngenst) { 
		
		sp.ids =c(s, s.index) #IDs of resident, and invaders (this way uses all species)
		nrs= nr[t,,s.index] #Residents 
		
		#isd.tmp = get.isd(Frs[t,,s], nrs,s, alphas, sr[s], kc, kd[,s])

		y[t,] = colMeans(nrs,na.rm=T) #Means of resident densities 
		w.eq[t,] = alphas[s,s.index]*y[t,] #Weights for the LDG

		l1[t,s]=mean(Frs[t,,s])
		D[t,s]= 1+sum(w.eq[t,])

		#Calculate the standardized competition from residents
		muj=nrs/(matrix(y[t,],np,( nspp-1),byrow=T))-1
		uijmuj = matrix(0,np,nspp-1)
		for(a in 1:(nspp-1)) {
				ns = s.index[a]
				uijmuj[,a] = convolve(muj[,a],kc[,ns])
				uijmuj[,a]= c(uijmuj[ceiling(np/2):(np),a], uijmuj[1:floor(np/2),a] )
				uijmuj[,a] = w.eq[t,a]*(uijmuj[,a]+abs(min(uijmuj[,a])))

		}

		uijmuj[uijmuj<0]=0
		#uijmuj = uijmuj*(matrix(c(w.eq[t,]/D[t,s]),np,nspp-1,byrow=T))

		#Total (standardized) competition experienced by invader used to calculate invasion growth rates
		cc[t,,s] = apply(uijmuj,1,sum)/D[t,s]
		
		#Non-linear competitive variance and the spatial storage effect
		#These terms represent the perturbation terms from Snyder 2008 
		var_mu_Us[t,s] = var(cc[t,,s])
		cov_e_mu_Us[t,s] = cov(Frs[t,,s]/mean(Frs[t,,s])-1, cc[t,,s])

		Elam1[t,s]=l1[t,s]/D[t,s]*(1+(var_mu_Us[t,s])-
				(cov_e_mu_Us[t,s]))+sr[s]-1

		Elam2[t,s]=(l1[t,s]*(1/D[t,s]) +sr[s]-1)^2+2*(l1[t,s]*(1/D[t,s]) +sr[s]-1)*
			(var_mu_Us[t,s]-cov_e_mu_Us[t,s])/D[t,s]^4

		
		# Elam1[t,s]=l1[t,s]*(1/D[t,s]+(var_mu_Us[t,s])/D[t,s]^3-
		# 		(cov_e_mu_Us[t,s])/D[t,s]^2)+sr[s]-1

		# Elam2[t,s]=(l1[t,s]*(1/D[t,s]) +sr[s]-1)^2+2*(l1[t,s]*(1/D[t,s]) +sr[s]-1)*
		# 	(var_mu_Us[t,s]-cov_e_mu_Us[t,s])/D[t,s]^4

		gr1.n[t,s] = exp(Elam1[t,s]-0.5*Elam2[t,s])
		
		#The fitness-density covariance 
		tryCatch( {cov_lam_vc[t,s]=get.fd.cov(Frs[t,,s],sp.ids, 
			nr[t,,],sr[s], alphas[s,], a_rr[s], kc)} , error=function(e){} )
								
		cov_lam_vc2[t,s]=get.fd.cov2(Frs[t,,s],sp.ids, 
			nr[t,,],sr[s], alphas, kd, kc)		
		#The full LDG
		#gr1[t,s] = gr1.n[t,s]+cov_lam_vc[t,s]
		gr1[t,s] = gr1.n[t,s]+covIMP.sim[t,s]

	} #End spatial terms

	plot(gr1.n[,s])
	points(ldgIMP.sim[,s],col="red")

	plot(gr1[,s])
	points(ldg.sim[,s],col="red")

	y.full=c(y.full, list(y) )
	w.eq.full=c(w.eq.full,list(w.eq))
	#=============================================================================
	# Low-density growth rates for pairwise cases!  -- using theoretical 2spp equilibrium
	#=============================================================================


	#Outer loop through residents 

	for (rs in 1:(nspp-1)){
		#These are all of the necessary terms to calculate the standard invasion growth rates

		enter1 = (nspp-1)*(s-1)+rs
		rez.no = s.index[rs]

		#Calculate standard spatial terms
		for (t in 1:sd.lim) { 

			l1.p[t,enter1]=mean(Frs[t,low[t]:up[t],s])
			y.p[t,enter1]=mean(nr2_eq_all[t,,rez.no])#Means of resident densities 
			w.eq.p[t,enter1]=alphas[s,rez.no]*y.p[t,enter1]#Weights for the LDG
			D.p[t,enter1]= 1+sum(w.eq.p[t,enter1])

			#Calculate the standardized competition from residents
			nrs.p=Re(fft((alphas[s,rez.no]*fft(nr2_eq_all[t,,rez.no])*fkc[,rez.no]),inverse=T)/(np+1))
			nrs.p=c(nrs.p[ceiling(np/2):(np)],nrs.p[1:floor(np/2)])
		
			#Calculate the standardized competition from residents
			muj.p=nrs.p/y.p[t,enter1]-1	
			uijmuj.p = convolve(muj.p,kc[,rez.no])
			uijmuj.p= c(uijmuj.p[ceiling(np/2):(np)], uijmuj.p[1:floor(np/2)] )
			uijmuj.p = w.eq.p[t,enter1]*(uijmuj.p+abs(min(uijmuj.p)))	
			uijmuj.p[uijmuj.p<0]=0

			cc.p=uijmuj.p

			#Non-linear competitive variance and the spatial storage effect
			#These terms represent the perturbation terms from Snyder 2008 
			var_mu_Us.p[t,enter1] = var(cc.p)
			cov_e_mu_Us.p[t,enter1] = cov(Frs[t,low[t]:up[t],s]/mean(Frs[t,low[t]:up[t],s])-1, cc.p)
			
			#First-order approximation
			Elam1.p[t,enter1]=l1.p[t,enter1]/(D.p[t,enter1])*( 1+1/(D.p[t,enter1])^2*var_mu_Us.p[t,enter1]-
					(1/D.p[t,enter1])*cov_e_mu_Us.p[t,enter1] )+sr[s]-1
			#Second-order approximation
			Elam2.p[t,enter1]=(l1.p[t,enter1]*(1/(D.p[t,enter1])) +
				sr[s]-1)^2+2*(l1.p[t,enter1]*(1/(D.p[t,enter1])) +
				sr[s]-1)*(var_mu_Us.p[t,enter1]-cov_e_mu_Us.p[t,enter1])/(D.p[t,enter1])^4
					
			gr1.n.p[t,enter1] = exp(Elam1.p[t,enter1]-0.5*Elam2.p[t,enter1])
			
			#The fitness-density covariance 
			#Method 1: 
			tryCatch( {cov_lam_vc.p[t,enter1]=get.fd.cov.p(Frs[t,,s], nr2_eq_all[t,,rez.no], 
				sr[s], alphas[s,rez.no], a_rr[s], kc[,s])} , error=function(e){} )
			#Method 2: 
			#calculate the invader stationary distribution: 
			cov_lam_vc.p2[t,enter1]=get.fd.cov.p2(Frs[t,low[t]:up[t],s],
				nr2_eq_all[t,,rez.no], sr[s], alphas[s,rez.no], kd[,s], kc[,s])
			

			#The full LDG
			gr1.p[t,enter1] = gr1.n.p[t,enter1]+cov_lam_vc.p[t,enter1]

		} #End spatial terms


	}#End outer species loop through residents
 print(s)

} #End main loop through species

#=============================================================================
#Data saving
#=============================================================================


#Save all of the meaningful invasion growth rate variables
#in a file with an informative but awkwardly long file name. 
#This includes information about the temporal and spatial extent,
#and the details of the scenario type. 

file.name = (paste(f.name1,f.name2,f.name3,f.name4, ".var",sep=""))
save(file=file.name, "l1", "D", "var_mu_Us","cov_e_mu_Us",
"cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", "w.eq.full", "l1.p",
"y.p", "w.eq.p", "D.p", "var_mu_Us.p", "cov_e_mu_Us.p","cov_lam_vc.p", "cov_lam_vc.p2",
"Elam1.p","Elam2.p", "gr1.n.p", "gr1.p","ldg.sim","ldgIMP.sim","covIMP.sim")
