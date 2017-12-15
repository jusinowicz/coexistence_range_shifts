###########################################################
#Figure making
###########################################################
#
#==============================================================================
#	Figure 1: Just make some ranges for N species and their 
#	expansions for a picture. Should encompass 3 main scenarios: 
#	A. Range translocation, 
#	B. Range expansion
#	C. Translocation and expansion
#==============================================================================

source("./range_coexistence_functionsB.R")

nspp=8 #Number of species
ns=800 #Lattice width/height -- This should be even
ngens=2 #Number of generations for environmental change
iconfig=1 #NUmber of interattions of initial configuration of species

#Make these automatically:
xx=matrix(seq((-ns/2),(ns/2)))
np=length(xx)
ngenst=iconfig+ngens

#Make the full array of space-time coordinates (for 1D space)
#meshgrid will produce an R data frame
stc.temp=meshgrid(1:ngenst,xx)
#Convert the coordinates into an array for easier access
stc=array(c(matrix(stc.temp$x,ngenst,np,byrow=T),matrix(stc.temp$y,ngenst,np,byrow=T)),dim=c(ngenst,np,2)) 


###############
#Distribution of intrinsic ranges
###############
# N species with random starting and ending mean positions, and range shapes 
# (both variance and distribution types)

Fr=matrix(1,nspp,1) #Peak height

# This way divides species into an upper and lower group, then lets members of 
# lower group move into upper group. 
	# 1. peaks (spatial mean)
	peaks=c(-200, -100, -90, -50,  50, 100, 150, 200)-100
	peaks.end=c(-160, -100, 200, -50, 50, 100, 150, 220)+100
	# 2. Ds (spatial variance)
	Ds = matrix( (ns^2/1000),nspp,1)
	Ds.end = c((ns^2/4000),(ns^2/1000),(ns^2/200), (ns^2/200), (ns^2/1000),(ns^2/4000) ,(ns^2/500), (ns^2/4000))

	# 3. bys
	fast.bys = FALSE
	peaks.by=(peaks.end-peaks)/ngens 
	Ds.by = (Ds.end-Ds)/ngens


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


color.use = c("gray30","red", "darkgreen", "gray30","blue", "gray50","gray30","gray60" )

fig.name = paste("sim_rangeshifts_8sppB.pdf",sep="")
pdf(file=fig.name, family='Helvetica', pointsize=16)

par(mfrow=c(2,1))
plot(Frs[1,,1],t="l",xlab="Location", ylab="Fitness", cex.lab=1.7,cex.axis=1.2,ylim=c(0,3),
		col=color.use[1],xaxs="i",xaxt='n',yaxs="i",yaxt='n',bty="n",lwd = 3)
abline(h=0)
for(j in 2:nspp) {
	lines( Frs[1,,j],col = color.use[j],lwd=3)
}

plot(Frs[3,,1],t="l",xlab="Location", ylab="Fitness", cex.lab=1.7,cex.axis=1.2,ylim=c(0,3),
		col=color.use[1],xaxs="i",xaxt='n',yaxs="i",yaxt='n',bty="n",lwd=3)
abline(h=0)
for(j in 2:nspp) {
	lines( Frs[3,,j],col = color.use[j],lwd=3)
}

dev.off()

#==============================================================================
#Figure 3
#
#
#==============================================================================

#This requires user input!

variable.list=list("l1", "D", "var_mu_Us","cov_e_mu_Us",
"cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", "w.eq.full", "l1.p",
"y.p", "w.eq.p", "D.p", "var_mu_Us.p", "cov_e_mu_Us.p","cov_lam_vc.p", "cov_lam_vc.p2",
"Elam1.p","Elam2.p", "gr1.n.p", "gr1.p","ldg.sim","ldgIMP.sim","covIMP.sim")


#List of scenario file names 
file.name.list= list("igr_3spp_s1000_t10000B_trans_novar_norm.var", "igr_3spp_s1000_t10000B_notrans_var_norm.var",
	"igr_3spp_s1000_t10000B_trans_var_norm.var", "igr_3spp_s1000_t10000B_trans_novar_fast.var" )

#In order to deal with each of the scenarios simultaneously, it is necessary
#to systematically rename variables based on their scenario. 
#Assosciated prefixes for renaming
prefix.list= list("tr","var","both","fast")

#Rename the files from each scenario
var.length=length(variable.list)

for (g in 1:4){
	load(file.name.list[[g]])
	for( f in 1:(var.length)) {      
		new.name=paste(prefix.list[[g]],variable.list[[f]],sep="")	
		assign(new.name,eval(as.name(variable.list[[f]])))

	}	
}


#Set this up 

################################################################
#Figure 3A: Non-Arch (linear) figures showing how the igr of each 
# scenario moves according to the benefit of either range width 
# or range overlap.
# This version for the 2SPP picture. 
# Components STANDARDIZED BY MEAN coefficients
# Figures are a SINGLE PLOT
################################################################
#Set these variables manually (hopefully you kept good notes)
iconfig=1
ngens=100
ngenst=iconfig+ngens
nspp=3
scale=10
#tol2=c(100,100,100,50)
tol2=c(100,100,100,70)

#Skipping columns
uu=c(1,2,3,5)
#Different limits
xlims = matrix(c(-0.0, 1.5, -0.0,1.5,-0.0,1.5, -0.0,1.5),4,2,byrow=T)
x.axx = matrix(c(0.2,2.4,0.4,0.2,2.4,0.4,0.2,2.4,0.4, 0.2,1.4,0.2 ),4,3,byrow=T)

#fig.name = paste("all4pairs_linear_widthOverlap_2Dgr1.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#Set this up as a panel of 4, or see below for individual files. 
par(mfrow=c(1,4),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
#This is a constant to scale the range overlap
min_overlap1=0
#Find the minimum value across the 3 normal scenarios
for( g in 1:3){ 
	for(u in 1:6){

		#Identify the variable names by scenario, then assign them to the right
		#variable for plotting
		cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[18]] ,sep="")
		cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[20]] ,sep="")
		l1.name=paste(prefix.list[[g]],variable.list[[13]] ,sep="")
		D.name = paste(prefix.list[[g]],variable.list[[16]] ,sep="")

		assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
		assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
		assign("l1",eval(as.name(l1.name)) )
		assign("D",eval(as.name(D.name)) )

		#Create a mean standardization coefficient to simplify presentation. 
		C1 =  median((l1[1,u]/(D[1,u])^3), na.rm=T)
		C2 =  median(l1[1,u]/(D[1,u])^2, na.rm=T )
		
		min_overlap1 = min(c (min_overlap1,(-C2*cov_e_mu_Us[tol2[u],u]+cov_lam_vc[tol2[u],u]),(-C2*cov_e_mu_Us[1,u]+cov_lam_vc[1,u])),na.rm=T)
	}
}

#Offset points in the y direction for readibility 	
j1 = jitter((1:4)*0, 0.5)

#Now make the plots: 
for( g in 1:4){ 

	color.use=list("darkgreen","darkgreen","red", "red","blue","blue")
	up.low=list("upper","lower")

	#Or just switch between devices for plotting
	#dev.set(g)
	#par(mfrow=c(1,4), mar=c(5,6,6,2))

	ylim=c(0,0.3)
	xlim=xlims[g,]
	#plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Range width", xlab="Range overlap",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3 )
	if (g ==1){
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Benefit of Competitor Aggregation", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n' )
	axis(side=1, labels=T, at=seq(x.axx[g,1], x.axx[g,2],x.axx[g,3]))
	axis(side=2, labels=T, at=seq(0.1,0.3,0.1))

	} else { 
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n')
	axis(side=1, labels=T, at=seq(x.axx[g,1], x.axx[g,2],x.axx[g,3]))

		} 

	for (ug in 1:4){ 
			u=uu[ug]
			
			#Identify the variable names by scenario, then assign them to the right
			#variable for plotting
			gr.name=paste(prefix.list[[g]],variable.list[[24]] ,sep="")
			var_mu_Us.name=paste(prefix.list[[g]],variable.list[[17]] ,sep="")
			cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[18]] ,sep="")
			cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[20]] ,sep="")
			#For scaling
			l1.name=paste(prefix.list[[g]],variable.list[[13]] ,sep="")
			w.eq.name = paste(prefix.list[[g]],variable.list[[15]] ,sep="")
			D.name = paste(prefix.list[[g]],variable.list[[16]] ,sep="")

			assign("gr",eval(as.name(gr.name)) )
			assign("var_mu_Us",eval(as.name(var_mu_Us.name)) )
			assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
			assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
			assign("l1",eval(as.name(l1.name)) )
			assign("D",eval(as.name(D.name)) )
			assign("cr_mean",D-1 )

			#This line represents persistence in the absence of spatial mechanisms
			#mlin.non1 =1-(l1[iconfig,u]/D[iconfig,u]+0.9)
			mlin.non1 =0.1
			lin.non1= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non1
			lines((seq(xlim[1],xlim[2],by=0.1)),lin.non1,col=color.use[[u]],lty = 2,lwd=0.5)

			#Create a mean standardization coefficient to simplify presentation. 
			#C1 =  median((l1[,u]*cr_mean[,u]^2/(D[,u])^3), na.rm=T)
			#C2 =  median(l1[,u]*cr_mean[,u]/(D[,u])^2, na.rm=T )

			C1 = (l1[1,u]/(D[1,u])^3)
			C2 = l1[1,u]/(D[1,u])^2


			#Spatial mechanisms add this much
			width_y1 = C1*var_mu_Us[iconfig,u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap1 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x1 = -min_overlap1+((-C2*cov_e_mu_Us[iconfig,u])+cov_lam_vc[iconfig,u])
			#Point where lines intercept: x = (b1-b2)/2
			#non.x1 = (mlin.non1 - (width_y1 - overlap_x1) )/2
			#non.y1 = -non.x1+ mlin.non1 
			#Intersection through origin instead:
			non.x1 = -mlin.non1/(-1-(width_y1/overlap_x1) )
			non.y1 = -non.x1+ mlin.non1 
			#points(non.x1,non.y1,col=color.use[[u]], pch=16)
			points(overlap_x1,width_y1+j1[ug],col=color.use[[u]] )
			#segments(non.x1,non.y1,overlap_x1,width_y1,col=color.use[[u]],lty=2)			


			#Now the later point
			#mlin.non2 =1-(l1[(ngenst-tol2[g]),u]/D[(ngenst-tol2[g]),u]+0.9)
			mlin.non2 =0.1
			lin.non2= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non2
			lines((seq(xlim[1],xlim[2],by=0.1)),y=lin.non2,col=color.use[[u]],lty=2,lwd=0.5  )

			#Spatial mechanisms add this much
			width_y2 = C1*var_mu_Us[(tol2[g]),u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap2 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x2 = -min_overlap1+((-C2*cov_e_mu_Us[(tol2[g]),u])+cov_lam_vc[(tol2[g]),u])
			non.x2 = -mlin.non2/(-1-(width_y2/overlap_x2) )
			non.y2 = -non.x2+ mlin.non2 
			#points(non.x2,non.y2 )
			points(overlap_x2,width_y2+j1[ug],col=color.use[[u]])
			#segments(non.x2,non.y2,overlap_x2,width_y2,col=color.use[[u]],lty=2)			

			#segments(non.x1,non.y1,new.x2,new.y2,col=color.use[[u]],lty=2)
			arrows(overlap_x1,width_y1+j1[ug],overlap_x2,width_y2+j1[ug],col=color.use[[u]],length=0.1,lwd=2)
		} 
}

dev.off()


################################################################
#Figure 3B: Non-Arch (linear) figures showing how the igr of each 
# scenario moves according to the benefit of either range width 
# or range overlap.
# This version for the full 3spp picture. 
# Components STANDARDIZED BY MEAN coefficients
# Figures are a SINGLE PLOT
################################################################
#Set these variables manually (hopefully you kept good notes)
iconfig=1
ngens=100
ngenst=iconfig+ngens
nspp=3
scale=10
#tol2=c(100,100,100,100)
tol2=c(100,100,100,70)

xlims = matrix(c(-0.0, 2.4, -0.0,2.4,-0.0,2.4, -0.0,2.4),4,2,byrow=T)
x.axx = matrix(c(0.2,2.4,0.4,0.2,2.4,0.4,0.2,2.4,0.4, 0.2,2.4,0.4 ),4,3,byrow=T)


#fig.name = paste("all4_linear_widthOverlap_2Dgr1.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#Set this up as a panel of 4, or see below for individual files. 
par(mfrow=c(1,4),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
#This is a constant to scale the range overlap
min_overlap1=0
#Find the minimum value across the 3 normal scenarios
for( g in 1:3){ 
	for(u in 1:3){
		#Identify the variable names by scenario, then assign them to the right
		#variable for plotting
		cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[4]] ,sep="")
		cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[27]] ,sep="")
		l1.name=paste(prefix.list[[g]],variable.list[[1]] ,sep="")
		D.name = paste(prefix.list[[g]],variable.list[[2]] ,sep="")

		assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
		assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
		assign("l1",eval(as.name(l1.name)) )
		assign("D",eval(as.name(D.name)) )

		#Create a mean standardization coefficient to simplify presentation. 
		C1 =  median((l1[1,u]/(D[1,u])^3), na.rm=T)
		C2 =  median(l1[1,u]/(D[1,u])^2, na.rm=T )

		min_overlap1 = min(c (min_overlap1,(-C2*cov_e_mu_Us[tol2[u],u]+cov_lam_vc[tol2[u],u]),(-C2*cov_e_mu_Us[1,u]+cov_lam_vc[1,u])),na.rm=T)
	}
}

#Now make the plots: 
for( g in 1:4){ 

	color.use=list("darkgreen","red","blue")
	up.low=list("upper","lower")

	#Or just switch between devices for plotting
	#dev.set(g)
	#par(mfrow=c(1,4), mar=c(5,6,6,2))

	ylim=c(0,0.6)
	#xlim=c(0,0.8)
	xlim=xlims[g,]
	#plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Range width", xlab="Range overlap",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3 )
	if (g ==1){
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Benefit of Competitor Aggregation", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n' )
	axis(side=1, labels=T, at=seq(x.axx[g,1], x.axx[g,2],x.axx[g,3]))
	axis(side=2, labels=T, at=seq(0.1,0.6,0.2))

	} else { 
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n')
	axis(side=1, labels=T, at=seq(x.axx[g,1], x.axx[g,2],x.axx[g,3]))

		} 

	for (u in 1:3){ 
			
			

			#Identify the variable names by scenario, then assign them to the right
			#variable for plotting
			gr.name=paste(prefix.list[[g]],variable.list[[10]] ,sep="")
			var_mu_Us.name=paste(prefix.list[[g]],variable.list[[3]] ,sep="")
			cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[4]] ,sep="")
			cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[27]] ,sep="")
			#For scaling
			l1.name=paste(prefix.list[[g]],variable.list[[1]] ,sep="")
			w.eq.name = paste(prefix.list[[g]],variable.list[[12]] ,sep="")
			D.name = paste(prefix.list[[g]],variable.list[[2]] ,sep="")

			assign("gr",eval(as.name(gr.name)) )
			assign("var_mu_Us",eval(as.name(var_mu_Us.name)) )
			assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
			assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
			assign("l1",eval(as.name(l1.name)) )
			assign("D",eval(as.name(D.name)) )
			assign("w.eq",eval(as.name(w.eq.name))  )

			#This line represents persistence in the absence of spatial mechanisms
			#mlin.non1 =1-(l1[iconfig,u]/D[iconfig,u]+0.9)
			mlin.non1 =0.1
			lin.non1= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non1
			lines((seq(xlim[1],xlim[2],by=0.1)),lin.non1,col=color.use[[u]],lty=2,lwd=0.5  )

			#Create a mean standardization coefficient to simplify presentation. 
			#C1 =  median((l1[,u]/(D[,u])^3), na.rm=T)
			#C2 =  median(l1[,u]/(D[,u])^2, na.rm=T )

			C1 = (l1[1,u]/(D[1,u])^3)
			C2 = l1[1,u]/(D[1,u])^2

			#Spatial mechanisms add this much
			width_y1 = C1*var_mu_Us[iconfig,u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap1 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x1 = -min_overlap1+((-C2*cov_e_mu_Us[iconfig,u])+cov_lam_vc[iconfig,u])
			#Point where lines intercept: x = (b1-b2)/2
			#non.x1 = (mlin.non1 - (width_y1 - overlap_x1) )/2
			#non.y1 = -non.x1+ mlin.non1 
			#Intersection through origin instead:
			non.x1 = -mlin.non1/(-1-(width_y1/overlap_x1) )
			non.y1 = -non.x1+ mlin.non1 
			#points(non.x1,non.y1,col=color.use[[u]], pch=16)
			points(overlap_x1,width_y1,col=color.use[[u]] )
			#segments(non.x1,non.y1,overlap_x1,width_y1,col=color.use[[u]],lty=2)			


			#Now the later point
			#mlin.non2 =1-(l1[(ngenst-tol2[u]),u]/D[(ngenst-tol2[u]),u]+0.9)
			mlin.non2 = 0.1
			lin.non2= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non2
			lines((seq(xlim[1],xlim[2],by=0.1)),lin.non2,col=color.use[[u]] ,lty=2,lwd=0.5  )

			#Spatial mechanisms add this much
			width_y2 = C1*var_mu_Us[(tol2[g]),u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap2 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x2 = -min_overlap1+((-C2*cov_e_mu_Us[(tol2[g]),u])+cov_lam_vc[(tol2[g]),u])
			non.x2 = -mlin.non2/(-1-(width_y2/overlap_x2) )
			non.y2 = -non.x2+ mlin.non2 
			#points(non.x2,non.y2 )
			points(overlap_x2,width_y2,col=color.use[[u]])
			#segments(non.x2,non.y2,overlap_x2,width_y2,col=color.use[[u]],lty=2)			

			#segments(non.x1,non.y1,new.x2,new.y2,col=color.use[[u]],lty=2)
			arrows(overlap_x1,width_y1,overlap_x2,width_y2,col=color.use[[u]],length=0.1,lwd=2)
		} 
}
dev.off()

###############################################################################
#Use this to create plots of the intrinsic ranges and their shifts
#for each scenario used in this figure
#Go through each scenario manually by commenting/uncommenting below
###############################################################################
fig.name = paste("all4_ranges_2Dgr1.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, family='Helvetica', pointsize=16)
###############################################################################
par(mfrow=c(2,4),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
###############################################################################
color.use=list("darkgreen","red","blue")

nspp=3 #number of species
ns=1000 #Lattice width/height -- This should be even
ngens=100 #Number of generations for environmental change
iconfig=1 #NUmber of interattions of initial configuration of species
ngenst=iconfig+ngens
scale=10

#tol2=c(100,100,100,100)
tol2=c(100,100,100,70)
Fr=matrix(c(2,2,2))  

#Make these automatically:
xx=matrix(seq((-ns/2),(ns/2)))
np=length(xx)
ngenst=iconfig+ngens
#Make the full array of space-time coordinates (for 1D space)
#meshgrid will produce an R data frame
stc.temp=meshgrid(1:ngenst,xx)
#Convert the coordinates into an array for easier access
stc=array(c(matrix(stc.temp$x,ngenst,np,byrow=T),matrix(stc.temp$y,ngenst,np,byrow=T)),dim=c(ngenst,np,2)) 
#Note: in 1d, things are a bit odd. stc.temp$x ends up as time, and $y as space. So in stc, stc[,,1] is time
# and stc[,,2] is space. stc[,,1] is the same as the column vector 1:ngenst repeated over np columns. Then
# stc[,,2] is space as a row vector (xx) repeated over ngenst rows. 

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
	#peaks.end=c(-50, -50, 150); f.name2=c("_notrans")
	peaks.end=c(150, -50, 150); f.name2=c("_trans")

	# 2. Ds
	Ds = c(np^2/1600, np^2/1600,np^2/1600)
	#Ds.end = c(np^2/50, np^2/1600,np^2/1600); f.name3=c("_var")
	Ds.end = c(np^2/1600, np^2/1600,np^2/1600); f.name3=c("_novar")

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

j1 = jitter((1:4), 1)
plot( (seq( (-ns/2),(ns/2)))[Frs[1,,1]>1e-2]+j1[1]*15, (Frs[1,,1][Frs[1,,1]>1e-2]), xlim=c((-ns/2),(ns/2)), 
	type="l", ylab="", xlab="", cex.lab=1.7,cex.axis=1.2, col=color.use[[1]],
	xaxs="i",xaxt='n',yaxs="i",yaxt='n',bty="n",lwd = 2 )
abline(h=0)
for(j in 2:nspp) {
	lines((seq( (-ns/2),(ns/2)))[Frs[1,,j]>1e-2], (Frs[1,,j][Frs[1,,j]>1e-2]),col = color.use[[j]],lwd=2)
}
lines((seq( (-ns/2),(ns/2)))[Frs[tol2[j],,1]>1e-2]-j1[1]*15, (Frs[tol2[j],,1][Frs[tol2[j],,1]>1e-2]),col = color.use[[1]],lwd=2,lty=2) 



dev.off()



#==============================================================================
# Figure 4 (?) A figure showing shifts in a multi-(N) species community
#
#==============================================================================
#load("igr_10spp_s1000_t10000c_trans_var_norm.var")
variable.list=list("frs_all", "l1", "D", "var_mu_Us","cov_e_mu_Us",
"cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", "w.eq.full", "l1.p",
"y.p", "w.eq.p", "D.p", "var_mu_Us.p", "cov_e_mu_Us.p","cov_lam_vc.p", "cov_lam_vc.p2",
"Elam1.p","Elam2.p", "gr1.n.p", "gr1.p","ldg.sim","ldgIMP.sim","covIMP.sim")

#List of scenario file names 
file.name.list= list("igr_10spp_s1000_t10000c_trans_var_norm.var")

iconfig=1
ngens=5
ngenst=iconfig+ngens
nspp=10
scale=1
tol2=3

color.use = c("gray30","red", "darkgreen", "gray50","blue", "brown","aquamarine3","darkorange","blueviolet","chocolate" )


fig.name = paste("allN_linear_widthOverlap_2Dgr1F.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#Set this up as a panel of 4, or see below for individual files. 
par(mfrow=c(1,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
#This is a constant to scale the range overlap
min_overlap1=0
#Find the minimum value across the 3 normal scenarios
 
	for(u in 1:nspp){
		#Identify the variable names by scenario, then assign them to the right
		#variable for plotting
		cov_e_mu_Us.name=paste(variable.list[[5]] ,sep="")
		cov_lam_vc.name=paste(variable.list[[28]] ,sep="")
		#For scaling
		l1.name=paste(variable.list[[2]] ,sep="")
		w.eq.name = paste(variable.list[[13]] ,sep="")
		D.name = paste(variable.list[[3]] ,sep="")

		assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
		assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
		assign("l1",eval(as.name(l1.name)) )
		assign("D",eval(as.name(D.name)) )

		#Create a mean standardization coefficient to simplify presentation. 
		C1 =  median((l1[1,u]/(D[1,u])^3), na.rm=T)
		C2 =  median(l1[1,u]/(D[1,u])^2, na.rm=T )

		min_overlap1 = min(c (min_overlap1,(-C2*cov_e_mu_Us[tol2[u],u]+cov_lam_vc[tol2[u],u]),(-C2*cov_e_mu_Us[1,u]+cov_lam_vc[1,u])),na.rm=T)
	}

xlims = matrix(c(-0.0, 3.6),byrow=T)
x.axx = matrix(c(0.2,3.6,0.4),byrow=T)
ylim=c(0,0.1)
#xlim=c(0,0.8)
xlim=xlims

plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Benefit of Competitor Aggregation", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n' )
	axis(side=1, labels=T, at=seq(x.axx[1], x.axx[2],x.axx[3]))
	axis(side=2, labels=T, at=seq(ylim[2]/2,ylim[2]/2,ylim[2]/2))

j1 = jitter((1:nspp)*0, 0.01)

for (u in 1:nspp){ 
			
			

			#Identify the variable names by scenario, then assign them to the right
			#variable for plotting
			gr.name=paste(variable.list[[11]] ,sep="")
			var_mu_Us.name=paste(variable.list[[4]] ,sep="")
			cov_e_mu_Us.name=paste(variable.list[[5]] ,sep="")
			cov_lam_vc.name=paste(variable.list[[28]] ,sep="")
			#For scaling
			l1.name=paste(variable.list[[2]] ,sep="")
			w.eq.name = paste(variable.list[[13]] ,sep="")
			D.name = paste(variable.list[[3]] ,sep="")

			assign("gr",eval(as.name(gr.name)) )
			assign("var_mu_Us",eval(as.name(var_mu_Us.name)) )
			assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
			assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
			assign("l1",eval(as.name(l1.name)) )
			assign("D",eval(as.name(D.name)) )
			assign("w.eq",eval(as.name(w.eq.name))  )

			#This line represents persistence in the absence of spatial mechanisms
			#mlin.non1 =1-(l1[iconfig,u]/D[iconfig,u]+0.9)
			mlin.non1 =0.1
			lin.non1= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non1
			lines((seq(xlim[1],xlim[2],by=0.1)),lin.non1,col=color.use[[u]],lty=2,lwd=0.5  )

			#Create a mean standardization coefficient to simplify presentation. 
			#C1 =  median((l1[,u]/(D[,u])^3), na.rm=T)
			#C2 =  median(l1[,u]/(D[,u])^2, na.rm=T )

			C1 = (l1[1,u]/(D[1,u])^3)
			C2 = l1[1,u]/(D[1,u])^2

			#Spatial mechanisms add this much
			width_y1 = C1*var_mu_Us[iconfig,u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap1 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x1 = -min_overlap1+((-C2*cov_e_mu_Us[iconfig,u])+cov_lam_vc[iconfig,u])
			#Point where lines intercept: x = (b1-b2)/2
			#non.x1 = (mlin.non1 - (width_y1 - overlap_x1) )/2
			#non.y1 = -non.x1+ mlin.non1 
			#Intersection through origin instead:
			non.x1 = -mlin.non1/(-1-(width_y1/overlap_x1) )
			non.y1 = -non.x1+ mlin.non1 
			#points(non.x1,non.y1,col=color.use[[u]], pch=16)
			points(overlap_x1,width_y1+j1[u],col=color.use[[u]] )
			#segments(non.x1,non.y1,overlap_x1,width_y1,col=color.use[[u]],lty=2)			


			#Now the later point
			#mlin.non2 =1-(l1[(ngenst-tol2[u]),u]/D[(ngenst-tol2[u]),u]+0.9)
			mlin.non2 = 0.1
			lin.non2= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non2
			lines((seq(xlim[1],xlim[2],by=0.1)),lin.non2,col=color.use[[u]] ,lty=2,lwd=0.5  )

			#Spatial mechanisms add this much
			width_y2 = C1*var_mu_Us[(tol2),u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap2 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x2 = -min_overlap1+((-C2*cov_e_mu_Us[(tol2),u])+cov_lam_vc[(tol2),u])
			non.x2 = -mlin.non2/(-1-(width_y2/overlap_x2) )
			non.y2 = -non.x2+ mlin.non2 
			#points(non.x2,non.y2 )
			points(overlap_x2,width_y2+j1[u],col=color.use[[u]])
			#segments(non.x2,non.y2,overlap_x2,width_y2,col=color.use[[u]],lty=2)			

			#segments(non.x1,non.y1,new.x2,new.y2,col=color.use[[u]],lty=2)
			arrows(overlap_x1,width_y1+j1[u],overlap_x2,width_y2+j1[u],col=color.use[[u]],length=0.1,lwd=2)
		} 

dev.off


fig.name = paste("allN_ranges_2Dgr1F.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, family='Helvetica', pointsize=16)
###############################################################################
par(mfrow=c(2,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))

j1 = jitter((1:nspp), 1)

ns = np
plot( (seq( (-ns/2),(ns/2)))[Frs[1,,1]>1e-2], (Frs[1,,1][Frs[1,,1]>1e-2]), xlim=c((-ns/2),(ns/2)), 
	ylim=c(0,6.5), type="l", ylab="", xlab="", cex.lab=1.7,cex.axis=1.2, col=color.use[[1]],
	xaxs="i",xaxt='n',yaxs="i",yaxt='n',bty="n",lwd = 2 )
abline(h=0)
for(j in 1:nspp) {
	lines((seq( (-ns/2),(ns/2)))[Frs[1,,j]>1e-2], (Frs[1,,j][Frs[1,,j]>1e-2]),col = color.use[[j]],lwd=2)
	lines((seq( (-ns/2),(ns/2)))[Frs[tol2,,j]>1e-2], (Frs[tol2,,j][Frs[tol2,,j]>1e-2]),col = color.use[[j]],lwd=2,lty=2) 

}

dev.off()
#==============================================================================
# Not figures, but useful for plotting
#
#==============================================================================

################################################################
#Standard IGR and components versus Environmental change. 
#For full 3spp community
#
################################################################
#Set these variables manually (hopefully you kept good notes)
iconfig=1
ngens=100
ngenst=iconfig+ngens
nspp=3
scale=10


#In four seperate files/plots: 

for( g in 1:3){ 

	#fig.name = paste(prefix.list[[g]],up.low[[u]],"_3sppscaled_ricomponents.jpg",sep="")
	#jpeg( file=fig.name,height=4.5, width=4.75, units="in", family='Helvetica', quality=100, res=300, pointsize=12)

	color.use=list("darkgreen","red","blue")
	up.low=list("upper","lower")

	u=1
	plot( gr[,u],t="n",xlab="Generations", ylab="Invasion growth rate", cex.lab=1.7,cex.axis=1.2,ylim=c(-0.5,2),col="black",xaxs="i",xaxt='n',bty="L")
	#axis(1, at=c(1000,2000,3000,4000,5000),label=c(10,20,30,40,50),cex.axis=1.2)
	#axis(1, at=c(8000,8500,9000,9500,10000),label=c(10,20,30,40,50),cex.axis=1.2)


	for (u in 1:3){ 
	
		#Identify the variable names by scenario, then assign them to the right
		#variable for plotting
		gr.name=paste(prefix.list[[g]],variable.list[[10]] ,sep="")
		var_mu_Us.name=paste(prefix.list[[g]],variable.list[[3]] ,sep="")
		cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[4]] ,sep="")
		cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[27]] ,sep="")
		#For scaling
		l1.name=paste(prefix.list[[g]],variable.list[[1]] ,sep="")
		w.eq.name = paste(prefix.list[[g]],variable.list[[12]] ,sep="")
		D.name = paste(prefix.list[[g]],variable.list[[2]] ,sep="")

		assign("gr",eval(as.name(gr.name)) )
		assign("var_mu_Us",eval(as.name(var_mu_Us.name)) )
		assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
		assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
		assign("l1",eval(as.name(l1.name)) )
		assign("D",eval(as.name(D.name)) )
		assign("w.eq",eval(as.name(w.eq.name))  )


		#Scaling of points
		pnts1=seq(1,ngenst, by = scale)
		#pnts2=seq(1,ngenst, by = scale)
		#pnts1=seq(burns,ngenst, by = 400)
		#pnts2=seq(1,ngenst, by = 10)

		#Species 1 as invader against residents
		#lines(gr[burns:(ngenst-10),u],col=color.use[[u]])
		points(pnts2, gr[pnts1,u],col=color.use[[u]])
		lines(-l1[,u]/(D[,u])^2*cov_e_mu_Us[,u], lty=2,col=color.use[[u]])
		lines(cov_lam_vc[,u], lty=2,col=color.use[[u]])
		points(pnts2,-l1[pnts1,u]/(D[pnts1,u])^2*cov_e_mu_Us[pnts1,u], pch=17,col=color.use[[u]])
		points(pnts2,cov_lam_vc[pnts1,u], pch=19,col=color.use[[u]])
		lines(l1[,u]/(D[,u])^3*var_mu_Us[,u], lty=2,col=color.use[[u]])

	
		
	}

	#dev.off()
}


################################################################
#Standard IGR and components versus Environmental change. 
#	Only for 1 vs. 2 and 1 vs.3 
#	Scaled
################################################################
#Set these variables manually (hopefully you kept good notes)
iconfig=1
ngens=100
ngenst=iconfig+ngens
nspp=3
scale=1

#In four seperate files/plots: 

for( g in 1:4){ 
	color.use=list("red","blue")
	up.low=list("upper","lower")
	for (u in 1:2){ 
		#Dynamic file naming
		#fig.name = paste(prefix.list[[g]],up.low[[u]],"_scaled_ricomponents.jpg",sep="")
		#jpeg( file=fig.name,height=4.5, width=4.75, units="in", family='Helvetica', quality=100, res=300, pointsize=12)

		#Or just switch between devices for plotting
		#dev.set(g)
		#par(mar=c(5,6,6,2))

		#Identify the variable names by scenario, then assign them to the right
		#variable for plotting
		gr.name=paste(prefix.list[[g]],variable.list[[10]] ,sep="")
		var_mu_Us.name=paste(prefix.list[[g]],variable.list[[3]] ,sep="")
		cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[4]] ,sep="")
		cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[27]] ,sep="")
		#For scaling
		l1.name=paste(prefix.list[[g]],variable.list[[1]] ,sep="")
		w.eq.name = paste(prefix.list[[g]],variable.list[[12]] ,sep="")
		D.name = paste(prefix.list[[g]],variable.list[[2]] ,sep="")

		assign("gr",eval(as.name(gr.name)) )
		assign("var_mu_Us",eval(as.name(var_mu_Us.name)) )
		assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
		assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
		assign("l1",eval(as.name(l1.name)) )
		assign("cr_mean",eval(as.name(cr_mean.name)) )
		assign("D",eval(as.name(D.name)) )

		#Scaling of points
		pnts1=seq(iconfig,ngenst, by = scale)
		pnts2=seq(1,(1+iconfig), by = scale)

		#Species 1 as invader against residents
		plot(gr[iconfig:(ngenst-10),u],t="l",xlab="Generations", ylab="Invasion growth rate", cex.lab=1.7,cex.axis=1.2,ylim=c(-0.5,2),col="black",xaxs="i",xaxt='n')
		axis(1, at=c(1000,2000,3000,4000,5000),label=c(10,20,30,40,50),cex.axis=1.2)
		lines(-l1[iconfig:(ngenst-10),u]*cr_mean[iconfig:(ngenst-10),u]/(D[iconfig:(ngenst-10),u])^2*cov_e_mu_Us[iconfig:(ngenst-10),u], lty=2,col="black")
		lines(cov_lam_vc[iconfig:(ngenst-10),u], lty=2,col="black")
		points(pnts2,-l1[pnts1,u]*cr_mean[pnts1,u]/(D[pnts1,u])^2*cov_e_mu_Us[pnts1,u], pch=17,col="black")
		points(pnts2,cov_lam_vc[pnts1,u], pch=19,col="black")
		lines(l1[iconfig:(ngenst-10),u]*cr_mean[iconfig:(ngenst-10),u]^2/(D[iconfig:(ngenst-10),u])^3*var_mu_Us[iconfig:(ngenst-10),u], lty=2,col="black")

		#Now each resident species as invader against species 1
		lines(gr[ iconfig:(ngenst-10),(2*u+1)],col=color.use[[u]])
		lines(-l1[iconfig:(ngenst-10),(2*u+1)]*cr_mean[iconfig:(ngenst-10),(2*u+1)]/(D[iconfig:(ngenst-10),(2*u+1)])^2*cov_e_mu_Us[ iconfig:(ngenst-10),(2*u+1)], lty=2,col=color.use[[u]])
		points(pnts2,cov_lam_vc[pnts1,(2*u+1)], pch=19,col=color.use[[u]])
		points(pnts2,-l1[pnts1,(2*u+1)]*cr_mean[pnts1,(2*u+1)]/(D[pnts1,(2*u+1)])^2*cov_e_mu_Us[pnts1,(2*u+1)], pch=17,col=color.use[[u]])
		lines(cov_lam_vc[ iconfig:(ngenst-10),(2*u+1)], lty=2,col=color.use[[u]])
		lines(l1[iconfig:(ngenst-10),(2*u+1)]*cr_mean[iconfig:(ngenst-10),(2*u+1)]^2/(D[iconfig:(ngenst-10),(2*u+1)])^3*var_mu_Us[ iconfig:(ngenst-10),(2*u+1)], lty=2,col=color.use[[u]])

		dev.off()
	}
}