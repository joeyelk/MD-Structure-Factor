import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator
import math
import tqdm
from tqdm import trange
import datetime
import time


mainlabel=""

units="$\AA^{-1}$"				#units of q vector

xunits=units
yunits=units
zunits=units

contours=200					#number of contours for contour plots
DPI=300							#resolution in DPI of images
format=".png"					#format to save images as
text="Structure Factor"			#label


#SHOW_PLOTS=False  #show plots in addition to saving.  This is useful for cropping or debugging

PLOT_EWALDS=True 	#enable ewald-corrected SF plots
savelog=True		#save log(SF)
savelin=True		#save SF
NBINSRAD=0			#Number of bins in Radial plot

FP_THRESHOLD=1.0E-12	#used to check for floating point zero value-anything less than this is considered zero (used in determining whether to use fast interpolation for ortho grid)


normplot=1 			#normalize plots

title_fontsize=9

path=""				


def make_flat_plot(D,xr,yr,zr):
	if len(xr)!=1 and len(yr)!=1 and len(zr)!=1:
		print "error in make_flat_plot!  one of these lengths must be 1"
		exit()

	for ix in xr:
		for iy in yr:
			for iz in zr:
				r=D[ix,iy,iz,:4]
				print r
				pts.append((r))



def pl(title,obj):		
	delim="="*20
	print delim,title,delim
	print obj
	print
	#print delim,"end ",title,delim


	
def pli(obj,title=""):
	pl(title,obj)
	buf=raw_input("enter q to quit, anything else to continue")
	if buf=='q':
		exit()
def ple(title,obj):
	pl(title,obj)
	exit()

		

def csplot_wlog(X,Y,Z,contours,lab,xlab,ylab,**kwargs):

	# kwargs['vmax']=np.amax(Z)
	# kwargs['vmin']=np.amin(Z)

	csplot(X,Y,Z,contours,lab,xlab,ylab,**kwargs)
	
	# kwargs['vmax']=np.amax(np.log10(Z))
	# kwargs['vmin']=np.amin(np.log10(Z))
	
	# csplot(X,Y,np.log(Z),contours,"log_"+lab,xlab,ylab,**kwargs)
	csplot(X,Y,np.log(Z),contours,"log_"+lab,xlab,ylab,**kwargs)


def csplot(X,Y,Z,contours,lab,xlab,ylab,**kwargs):
	
	title=lab+" S("+xlab+","+ylab+")"
	fname=lab+"_"+xlab+"_"+ylab
	fig, ax = plt.subplots()
	plt.suptitle(title)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	
	if normplot==1:
		cax=plt.contourf(X,Y,Z/np.amax(Z),contours,**kwargs)
	else:
		cax=plt.contourf(X,Y,Z,contours,**kwargs)
	

	#ax.set_aspect((np.amax(Y)-np.amin(Y))/(np.amax(X)-np.amin(X)))
	# ax.set_aspect('auto')
	cbar = fig.colorbar(cax)
	
	plt.savefig(path+fname+format,dpi=DPI)
	#print
	#print "saving ",path+fname+format
	# 
	plt.clf()
	
	#ax.pcolormesh(XMG[:-1,:], YMG[:-1,:], Hrz.T,vmin=0.0,vmax=np.amax(Hrz))
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
	
	#plt.savefig(path+"rzplot"+format,dpi=DPI)
	#plt.clf()
	

#plot slice through structure factor
def sfplot(data,**kwargs):



	if not os.path.exists(path):
		os.makedirs(path)
	cspos=0.0
	la=[]
	
	lb=0
	
	an=['x','y','z']  #axes names
	for i in xrange(data.shape[2]-1):
		if np.unique(data[...,i]).size>1:
			la.append(i)
		else:
			lb=i
			cspos=data[0,0,i]
			

	title=mainlabel+"\n"+text+"\n"+an[lb]+"="+str(round(cspos,2)) + zunits
	ltitle=mainlabel+"\n"+"log "+text+"\n"+an[lb]+"="+str(round(cspos,2)) + zunits
	
	xlab=an[la[0]]
	ylab=an[la[1]]
	
	filename=path+an[lb]+"="+str(round(cspos,2))

	
	xlab+="("+xunits+")"
	ylab+="("+yunits+")"
	
	if savelog:
		plt.suptitle(ltitle,fontsize=title_fontsize)
		plt.xlabel(xlab)
		plt.ylabel(ylab)
		plt.contourf(data[...,la[0]],data[...,la[1]],np.log(data[...,3]),contours,**kwargs)
		
		plt.savefig(filename+"_log"+format,dpi=DPI)
		plt.clf()
	
	if savelin:
		plt.suptitle(title,fontsize=title_fontsize)
		plt.xlabel(xlab)
		plt.ylabel(ylab)
		plt.contourf(data[...,la[0]],data[...,la[1]],data[...,3],contours,**kwargs)
		plt.savefig(filename+format,dpi=DPI)
		plt.clf()

def radial_integrate(D,Nbins,outputname):



	#X =D[:,0,0,0]
	#Y =D[0,:,0,1]
	#Z =D[0,0,:,2]
	
	
	
	SF=D[:,:,:,3]
	#SF[SF.shape[0]/2,SF.shape[1]/2,SF.shape[2]/2]=0.0
	
	#for i in tqdm.tqdm(xrange(100)):
		#SF[np.unravel_index(np.argmax(SF),(SF.shape))]=0.0
		
	#print SF.shape[0]/2,SF.shape[1]/2,SF.shape[2]/2
	#SF[127,126,127]=0.0
	# 
	
	
	print D[:,0,0,0]
	print "done!"
	# print "done!"
	# exit()
	
	R=(D[:,:,:,0]**2).astype(np.float16)+(D[:,:,:,1]**2).astype(np.float16)+(D[:,:,:,2]**2).astype(np.float16)
	H,E=np.histogram(R,bins=Nbins,weights=SF)
	Hc,E=np.histogram(R,bins=Nbins)
	Hc=np.where(Hc!=0,Hc,1.0)
	H/=Hc
	H[:1]=0.0
	H/=np.amax(H)
	
	#print np.argmax(H)
	
	
	plt.plot(E[:-1],H)
	plt.xlim(0,5)
	plt.savefig(outputname,dpi=DPI)
	
	
	#print H
	#print type(H)
	#print H.shape
	#print H.shape
	#print R
	#print R.shape
	
		
	# ES = RegularGridInterpolator((X, Y, Z), SF)	
	
	# pts=[]
	# for x in np.linspace(np.amin(X),np.amax(X),X.size*Ni):
		# for y in np.linspace(np.amin(Y),np.amax(Y),Y.size*Ni):
			# for z in np.linspace(np.amin(Z),np.amax(Z),Z.size*Ni):
				# pts.append((x,y,z))
			
	# A=ES(pts)
	
	



def Plot_Ewald_Sphere_Correction_old(D,wavelength_angstroms):  #pass full 3d data,SF,wavelength in angstroms
	
	X =D[:,0,0,0]
	Y =D[0,:,0,1]
	Z =D[0,0,:,2]
	SF=D[:,:,:,3]
	
	K_ES=2.0*math.pi/wavelength_angstroms  #calculate k for incident xrays in inverse angstroms
	
	ES = RegularGridInterpolator((X, Y, Z), SF)		
	
	pts=[]
	for ix in xrange(D.shape[0]):
		xsq=X[ix]**2.0
		for iy in xrange(D.shape[1]):
			R=np.sqrt(xsq+Y[iy]**2.0)
			theta=np.arctan(R/K_ES)
			xnew=X[ix]*np.cos(theta)
			ynew=Y[iy]*np.cos(theta)
			znew=K_ES*(1.0-np.cos(theta))
			pts.append((X[ix],Y[iy],xnew,ynew,znew))

		
	
	pts=np.asarray(pts)
	
	EWD=ES(pts[:,2:])
	EWD=EWD.reshape(D.shape[0],D.shape[1])
	plt.contourf(D[:,:,0,0],D[:,:,0,1],EWD,200,interpolation=interp)
	
	plt.savefig("EWxy.png",dpi=300)
	plt.clf()
	
	plt.contourf(D[:,:,0,0],D[:,:,0,1],np.log(EWD),200,interpolation=interp)
	
	plt.savefig("EWxylog.png",dpi=300)
	plt.clf()
	
def Plot_Ewald_Sphere_Correction(D,wavelength_angstroms,ucell=[],**kwargs):  #pass full 3d data,SF,wavelength in angstroms

	
	
	print D.shape

	if not os.path.exists(path):
		os.makedirs(path)
	
	X =D[:,0,0,0]
	Y =D[0,:,0,1]
	Z =D[0,0,:,2]
	SF=D[:,:,:,3]
	
	K_ES=2.0*math.pi/wavelength_angstroms  #calculate k for incident xrays in inverse angstroms

	
	ES = RegularGridInterpolator((X, Y, Z), SF,bounds_error=False)		
	
	xypts=[]
	for ix in xrange(D.shape[0]):
		xsq=X[ix]**2.0
		for iy in xrange(D.shape[1]):
			theta=np.arctan(np.sqrt(xsq+Y[iy]**2.0)/K_ES)
			xypts.append((X[ix]*np.cos(theta),Y[iy]*np.cos(theta),K_ES*(1.0-np.cos(theta))))
			
	xzpts=[]
	for ix in xrange(D.shape[0]):
		xsq=X[ix]**2.0
		for iz in xrange(D.shape[2]):
			theta=np.arctan(np.sqrt(xsq+Z[iz]**2.0)/K_ES)
			xzpts.append((X[ix]*np.cos(theta),K_ES*(1.0-np.cos(theta)),Z[iz]*np.cos(theta)))
	
	yzpts=[]
	for iy in xrange(D.shape[1]):
		ysq=Y[iy]**2.0
		for iz in xrange(D.shape[2]):
			theta=np.arctan(np.sqrt(ysq+Z[iz]**2.0)/K_ES)
			yzpts.append((K_ES*(1.0-np.cos(theta)),Y[iy]*np.cos(theta),Z[iz]*np.cos(theta)))
	
	xypts=np.asarray(xypts)
	xzpts=np.asarray(xzpts)
	yzpts=np.asarray(yzpts)
	
	EWDxy=ES(xypts)
	EWDxz=ES(xzpts)
	EWDyz=ES(yzpts)
	
	EWDxy=EWDxy.reshape(D.shape[0],D.shape[1])
	EWDxz=EWDxz.reshape(D.shape[0],D.shape[2])
	EWDyz=EWDyz.reshape(D.shape[1],D.shape[2])
	
	title="Ewald Corrected Structure Factor \n $\lambda=$"+str(wavelength_angstroms)+" $\AA$   $k_{ew}=$"+str(round(K_ES,2))+" $\AA^{-1}$"
	ltitle='log ' + title
	
	xlab='x ('+units + ")"
	ylab='y ('+units + ")"
	zlab='z ('+units + ")"
	
	fname="Ewald_"	
	
	plt.suptitle(title)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.contourf(D[:,:,0,0],D[:,:,0,1],EWDxy,contours,**kwargs)
	plt.savefig(path+fname+"xy"+format,dpi=DPI)
	plt.clf()
	
	plt.suptitle(ltitle)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.contourf(D[:,:,0,0],D[:,:,0,1],np.log(EWDxy),contours,**kwargs)
	plt.savefig(path+fname+"xylog"+format,dpi=DPI)
	plt.clf()
	
	plt.suptitle(title)
	plt.xlabel(xlab)
	plt.ylabel(zlab)
	plt.contourf(D[:,0,:,0],D[:,0,:,2],EWDxz,contours,**kwargs)
	plt.savefig(path+fname+"xz"+format,dpi=DPI)
	plt.clf()
	
	plt.suptitle(ltitle)
	plt.xlabel(xlab)
	plt.ylabel(zlab)
	plt.contourf(D[:,0,:,0],D[:,0,:,2],np.log(EWDxz),contours,**kwargs)
	plt.savefig(path+fname+"xzlog"+format,dpi=DPI)
	plt.clf()
	
	plt.suptitle(title)
	plt.xlabel(ylab)
	plt.ylabel(zlab)
	plt.contourf(D[0,:,:,1],D[0,:,:,2],EWDyz,contours,**kwargs)
	plt.savefig(path+fname+"yz"+format,dpi=DPI)
	plt.clf()
	
	plt.suptitle(ltitle)
	plt.xlabel(ylab)
	plt.ylabel(zlab)
	plt.contourf(D[0,:,:,1],D[0,:,:,2],np.log(EWDyz),contours,**kwargs)
	plt.savefig(path+fname+"yzlog"+format,dpi=DPI)
	plt.clf()



def PLOT_RAD_NEW(D,wavelength_angstroms,ucell,**kwargs):
	if not os.path.exists(path):
		os.makedirs(path)
		
	X=D[:,0,0,0]
	Y=D[0,:,0,1]
	Z=D[0,0,:,2]
	SF=D[...,3]
	
	# print ucell
	# exit()
		
	ES = RegularGridInterpolator((X, Y, Z), SF,bounds_error=False)
	
	THETA_BINS_PER_INV_ANG=20.
	MIN_THETA_BINS=10		#minimum allowed bins
	RBINS=400		
	NLEVELS=200		#number of levels for contour plots
	
	# MIN=0.5
	# MAX=1.5
	
	ZBINS=Z.shape[0]	#400
	
	XR=(X[-1]-X[0])*ucell[0][0]
	YR=(Y[-1]-Y[0])*ucell[1][1]
	
	Rmax=min(XR,YR)/2.0
	Rmax*=0.95
	
	rarr,rspace=np.linspace(0.0,Rmax,RBINS,retstep=True)
	zar=np.linspace(Z[0],Z[-1],ZBINS)
	
	oa=np.zeros((rarr.shape[0],zar.shape[0]))			
	circ=2.*np.pi*rarr											#circumference
	
	a1=ucell[0]
	a2=ucell[1]
	a3=ucell[2]
		
	b1=(np.cross(a2,a3))/(np.dot(a1,np.cross(a2,a3)))#
	b2=(np.cross(a3,a1))/(np.dot(a2,np.cross(a3,a1)))#*2.0*math.pi
	b3=(np.cross(a1,a2))/(np.dot(a3,np.cross(a1,a2)))#*2.0*math.pi 
	
	
	for ir in trange(rarr.shape[0]):
		
		NTHETABINS=max(int(THETA_BINS_PER_INV_ANG*circ[ir]),MIN_THETA_BINS)	#calculate number of bins at this r
		thetas=np.linspace(0.0,np.pi*2.0,NTHETABINS,endpoint=False)		#generate theta array
		
		t,r,z=np.meshgrid(thetas,rarr[ir],zar)							#generate grid of cylindrical points
		
		xar=r*np.cos(t)													#set up x,y coords
		yar=r*np.sin(t)	
		
		pts=np.vstack((xar.ravel(),yar.ravel(),z.ravel())).T		#reshape for interpolation
		
		# print pts.shape
		# exit()
		#MCpts=tm2(pts,ucell)
		MCpts=mc_inv(pts,ucell)
		# MCpts=tm3(pts,b1,b2,b3)								#transform to monoclinic cell
		# MCpts=pts
		
		oa[ir,:]=np.average(ES(MCpts).reshape(r.shape),axis=1)	#store average values in final array
		# print np.amin(oa)
		# print np.amax(oa)
		# exit()
		
	
	# oa=np.where(oa.isnan
	# oa=np.nan_to_num(oa)
	mn=np.nanmin(oa)
	print oa
	print np.nanmin(oa)
	print np.amax(oa)
	oa=np.where(np.isnan(oa),mn,oa)
	
	print np.amin(oa)
	print np.amax(oa)
	
	# exit()
	
	rad_avg=np.average(oa)	
	
	oa/=rad_avg		#normalize
	
	
	
	#set up data for contourf plot
	
	final=np.append(oa[::-1,:],oa[1:],axis=0)		#SF
	rfin=np.append(-rarr[::-1],rarr[1:])			#R
	zfin=np.append(z[:,0,:],z[1:,0,:],axis=0)		#Z
	
	unitlab='($\AA^{-1}$)'				#Angstroms
	
	logfinal=np.log(final)
	
	MIN=np.amin(final)
	MAX=np.amax(final)
	
	lvls=np.linspace(MIN,MAX,NLEVELS)			#contour levels
	
	
	plt.contourf(rfin,zfin[0],final.T,levels=lvls,cmap='jet')
	plt.colorbar()
	
	plt.title('S(r,z)')
	plt.xlabel('r '+unitlab)
	plt.ylabel('z '+unitlab)
	
	plt.savefig('new_rzplot.png')
	plt.clf()
	
	
	
	cs = plt.contourf(rfin, zfin[0], final.T, levels=lvls, cmap='jet', extend='both')
	cs.cmap.set_under('k')
	cs.set_clim(MIN, 0.1*MAX)
	plt.title('S(r,z)')
	plt.xlabel('r '+unitlab)
	plt.ylabel('z '+unitlab)
	plt.colorbar()
	plt.savefig('cs.png')
	plt.clf()
	
	
	
	
	plt.contourf(rfin,zfin[0],final.T,levels=lvls,cmap='jet')
	plt.colorbar()
	
	plt.title('S(r,z)')
	plt.xlabel('r '+unitlab)
	plt.ylabel('z '+unitlab)
	plt.savefig('new_rzplot2.png')
	plt.clf()
	
	
	
	
	lglvls=np.linspace(np.amin(logfinal),np.amax(logfinal),NLEVELS)
	
	plt.contourf(rfin,zfin[0],logfinal.T,levels=lglvls,cmap='jet')	
	plt.colorbar()
	
	plt.title('ln(S(r,z))')
	plt.xlabel('r '+unitlab)
	plt.ylabel('z '+unitlab)
	plt.savefig('new_log_rzplot.png')
	plt.clf()
	
	
	
	x2=np.linspace(-Rmax,Rmax,RBINS*2-1)
	z2=np.linspace(Z[0],Z[-1],RBINS)
	
	xg2,yg2,zg2=np.meshgrid(x2,np.asarray(0),z2)
	pts=np.vstack((xg2.ravel(),yg2.ravel(),zg2.ravel())).T
	out2=ES(pts).reshape(xg2.shape[1],xg2.shape[2])
	
	o2n=out2[:,:]/rad_avg
	
	plt.contourf(xg2[0,:,:],zg2[0,:,:],o2n,levels=lvls,cmap='jet')
	
	
	plt.xlabel('x '+unitlab)
	plt.ylabel('z '+unitlab)
	plt.title('S(x,z)|$_{y=0}$')
	
	plt.colorbar()
	plt.savefig('new_xzplot.png')
	plt.clf()
	
	
	x2=np.linspace(-Rmax,Rmax,RBINS*2-1)
	y2=np.linspace(-Rmax,Rmax,RBINS*2-1)
	
	xg2,yg2,zg2=np.meshgrid(x2,y2,np.asarray(0))
	pts=np.vstack((xg2.ravel(),yg2.ravel(),zg2.ravel())).T
	out2=ES(pts).reshape(xg2.shape[0],xg2.shape[1])
	
	o2n=out2[:,:]/np.average(out2)
	
	lvlsxy=np.linspace(np.amin(o2n),np.amax(o2n),NLEVELS)			#contour levels
	
	plt.contourf(xg2[:,:,0],yg2[:,:,0],o2n,levels=lvlsxy,cmap='jet')
	
	
	plt.xlabel('x '+unitlab)
	plt.ylabel('y '+unitlab)
	plt.title('S(x,y)')#|$_{y=0}$')
	
	plt.colorbar()
	plt.savefig('new_xyplot.png')
	plt.clf()
	
	
	
	if False:
		dif=o2n-final
		lvls2=np.linspace(-0.4,0.4,100)
		
		plt.contourf(xg2[0,:,:],zg2[0,:,:],dif,levels=lvls2,cmap='seismic')
		plt.xlabel('x,r '+unitlab)
		plt.ylabel('z '+unitlab)
		plt.title('S(r,z)-S(x,z)|$_{y=0}$')
		
		
		plt.colorbar()
		plt.savefig('difference.png')

	
def PLOT_EW_NEW_OLD(D,wavelength_angstroms,ucell,**kwargs):
	if not os.path.exists(path):
		os.makedirs(path)
		
	X=D[:,0,0,0]
		
	Y=D[0,:,0,1]
	Z=D[0,0,:,2]
	SF=D[...,3]
	
	
	
	
	# print X
	# exit()
	
	ES = RegularGridInterpolator((X, Y, Z), SF,bounds_error=False)
	
	
	THETABINS=400
	RBINS=400
	ZBINS=400
	
	XR=(X[-1]-X[0])*ucell[0][0]
	YR=(Y[-1]-Y[0])*ucell[1][1]
	
	Rmax=min(XR,YR)/2.0
	Rmax*=0.95
	print Rmax
	
	thetas=np.linspace(0.0,np.pi*2.0,THETABINS,endpoint=False)
	rarr,rspace=np.linspace(0.0,Rmax,RBINS,retstep=True)
	# rarr+=rspace
	zar=np.linspace(Z[0],Z[-1],ZBINS)
	
	t,r,z=np.meshgrid(thetas,rarr,zar)
	
	
	# print t[:,0,:]
	# exit()
	# print thetas
	# exit()
	xar=r*np.cos(t)
	yar=r*np.sin(t)	
	
	# pts=np.asarray([xar,yar,zar])
	
	print xar.shape
	print yar.shape
	print z.shape
	
	pts=np.vstack((xar.ravel(),yar.ravel(),z.ravel())).T
	
	MCpts=to_monoclinic(pts,ucell)
	
	print pts.shape
	# print "interp"
	# exit()
	NSPLIT=4
	if NSPLIT>1:
		out=np.asarray([])
		LAR=np.array_split(MCpts,NSPLIT)
		print "processing split array"
		for i in tqdm.tqdm(LAR):
			out=np.append(out,ES(i))
	else:
		out=ES(MCpts)
	# print MCpts.shape
	# exit()
	print out.shape
	
	
	o=out.reshape(r.shape)
	print np.amax(o)
	
	oa=np.average(o,axis=1)	#average over thetas
	
	print np.amax(oa)
	# exit()
	print oa.shape
	# print rfin.shape
	
	# rad_max=np.amax(oa)
	rad_max=np.average(oa)
	
	oa/=rad_max
	final=np.append(oa[::-1,:],oa[1:],axis=0)
	rfin=np.append(-r[::-1,0,:],r[1:,0,:],axis=0)
	zfin=np.append(z[:,0,:],z[1:,0,:],axis=0)
	
	# print rfin[:,0]
	# exit()
	
	print rfin.shape
	print zfin.shape
	
	print np.amax(final)
	print np.amin(final)
	
	unitlab='($\AA^{-1}$)'
	
	NLEVELS=200
	lvls=np.linspace(0.5,1.5,NLEVELS)
	plt.contourf(rfin,zfin,final,levels=lvls,cmap='seismic')
	plt.colorbar()
	
	plt.title('S(r,z)')
	plt.xlabel('r '+unitlab)
	plt.ylabel('z '+unitlab)
	plt.savefig('new_rzplot.png')
	plt.clf()
	
	x2=np.linspace(-Rmax,Rmax,RBINS*2-1)
	z2=np.linspace(Z[0],Z[-1],RBINS)
	
	
	
	xg2,yg2,zg2=np.meshgrid(x2,np.asarray(0),z2)
	pts=np.vstack((xg2.ravel(),yg2.ravel(),zg2.ravel())).T
	out2=ES(pts).reshape(xg2.shape[1],xg2.shape[2])
	
	print out2
	print out2.shape
	print xg2.shape
	
	o2n=out2[:,:]/rad_max
	
	plt.contourf(xg2[0,:,:],zg2[0,:,:],o2n,levels=lvls,cmap='seismic')
	
	
	plt.xlabel('x '+unitlab)
	plt.ylabel('z '+unitlab)
	plt.title('S(x,z)|$_{y=0}$')
	# plt.contourf(rfin,zfin,final,levels=lvls)
	plt.colorbar()
	plt.savefig('new_xzplot.png')
	plt.clf()
	# print o2n.shape
	# print final.shape
	# exit()
	
	dif=o2n-final
	lvls2=np.linspace(-0.4,0.4,100)
	
	plt.contourf(xg2[0,:,:],zg2[0,:,:],dif,levels=lvls2,cmap='seismic')
	plt.xlabel('x,r '+unitlab)
	plt.ylabel('z '+unitlab)
	plt.title('S(r,z)-S(x,z)|$_{y=0}$')
	
	
	plt.colorbar()
	plt.savefig('difference.png')
	
	
	
	
	# plt.contourf(r[:,0,:],z[:,0,:],oa,100)
	# plt.contourf(-r[1:,0,:],z[1:,0,:],oa[1:,...],100)
	
	
	# plt.show()
		
	# for ir in np.linspace(0.0,10.0,50):
		# print ir
	exit()
	
	
	
	
	# buf=np.asarray([[1.,0,0],[0,1.,0],[0,0,1.],[1.,1.123,1.567]])
	# a=to_monoclinic(buf,ucell)
	# b=to_cartesian(a,ucell)
	
	
	
	print ucell[0]
	print ucell[1]
	print ucell[2]
	exit()
	
def mc_inv(D,ucell):
		
	a1=ucell[0]
	a2=ucell[1]
	a3=ucell[2]
	

	
	b1=(np.cross(a2,a3))/(np.dot(a1,np.cross(a2,a3)))#
	b2=(np.cross(a3,a1))/(np.dot(a2,np.cross(a3,a1)))#*2.0*math.pi
	b3=(np.cross(a1,a2))/(np.dot(a3,np.cross(a1,a2)))#*2.0*math.pi 
	
	
	b_inv=np.linalg.inv(np.vstack((b1,b2,b3)))
	Dnew=np.zeros_like(D)
	
	X=D[...,0]
	Y=D[...,1]
	Z=D[...,2]
	
	for ix in xrange(D.shape[0]):			
		Dnew[ix,0:3]+=X[ix]*b_inv[0]  #(X[ix]-X[X.shape[0]/2])*b1
		
	for iy in xrange(D.shape[0]):			
		Dnew[iy,0:3]+=Y[iy]*b_inv[1]  #(Y[iy]-Y[Y.shape[0]/2])*b2
		
	for iz in xrange(D.shape[0]):			
		Dnew[iz,0:3]+=Z[iz]*b_inv[2]  #(Z[iz]-Z[Z.shape[0]/2])*b3
	return Dnew
	
def tm3(D,b1,b2,b3):
	
	Dnew=np.zeros_like(D)
	
	X=D[...,0]
	Y=D[...,1]
	Z=D[...,2]
	
	for ix in xrange(D.shape[0]):			
		Dnew[ix,0:3]+=X[ix]*b1  #(X[ix]-X[X.shape[0]/2])*b1
		
	for iy in xrange(D.shape[0]):			
		Dnew[iy,0:3]+=Y[iy]*b2  #(Y[iy]-Y[Y.shape[0]/2])*b2
		
	for iz in xrange(D.shape[0]):			
		Dnew[iz,0:3]+=Z[iz]*b3  #(Z[iz]-Z[Z.shape[0]/2])*b3
	return Dnew

def tm2(D,ucell):
	a1=ucell[0]
	a2=ucell[1]
	a3=ucell[2]
	

	
	b1=(np.cross(a2,a3))/(np.dot(a1,np.cross(a2,a3)))#
	b2=(np.cross(a3,a1))/(np.dot(a2,np.cross(a3,a1)))#*2.0*math.pi
	b3=(np.cross(a1,a2))/(np.dot(a3,np.cross(a1,a2)))#*2.0*math.pi 
	
	

	
	Dnew=np.zeros_like(D)
	
	X=D[...,0]
	Y=D[...,1]
	Z=D[...,2]
	
	for ix in xrange(D.shape[0]):			
		Dnew[ix,0:3]+=X[ix]*b1  #(X[ix]-X[X.shape[0]/2])*b1
		
	for iy in xrange(D.shape[0]):			
		Dnew[iy,0:3]+=Y[iy]*b2  #(Y[iy]-Y[Y.shape[0]/2])*b2
		
	for iz in xrange(D.shape[0]):			
		Dnew[iz,0:3]+=Z[iz]*b3  #(Z[iz]-Z[Z.shape[0]/2])*b3
	return Dnew
	
	
def to_monoclinic(coords,ucell):		#monoclinic for now
	out=coords.copy()
	out[...,1]/=ucell[1,1]
	out[...,0]-=out[...,1]*ucell[1,0]
	return out

	# T[...,1]=T[...,1]/np.sin(theta)
	# T[...,0]=T[...,0]-T[...,1]*np.cos(theta)
	
	
def to_cartesian(coords,ucell):
	out=np.zeros_like(coords)
	
	out=np.einsum('...i,ij->...j',coords,ucell)
	# out=np.einsum('...i,ij->...j',coords[...,0],ucell[0])
	return out
	
	
	
	
def Plot_Ewald_triclinic(D,wavelength_angstroms,ucell,**kwargs):  #pass full 3d data,SF,wavelength in angstroms
		
	PLOT_RAD_NEW(D,wavelength_angstroms,ucell,**kwargs)
	exit()
		
	if not os.path.exists(path):
		os.makedirs(path)
		
		
	# print path
	# exit()
		
	X =D[:,0,0,0].copy()
	Y =D[0,:,0,1].copy()
	Z =D[0,0,:,2].copy()
	
	# NBINSRAD=0
	NBINSZ=1*D[0,0,:,2].size
	ZBNS=np.linspace(Z[0],Z[-1],NBINSZ)
	
	# print NBINSRAD
	# exit()
	
	if NBINSRAD>0:
		XBNSRD=np.linspace(-NBINSRAD,NBINSRAD,num=NBINSRAD*2)
		XBNSRD=np.sqrt(np.abs(XBNSRD))*np.sign(XBNSRD)
		XBNSRD*=(X[-1]/XBNSRD[-1])
		
	else:
		XBNSRD=X
		# XBNSRD=np.linspace(X[0],X[-1],200)
		# ZBNS=np.linspace(Z[0],Z[-1],200)
		print 'setting XBNSRD=',X
		# exit()
	
	# print XBNSRD
	# exit()
	dx1=X[1+X.shape[0]/2]-X[X.shape[0]/2]
	
	
	# XBNSRD=np.sqrt(np.abs(X))*np.sign(X)#np.sqrt(X[X.shape[0]/2:])
	
	# XBNSRD=X*4
	
	# print XBNSRD
	# print "="*50
	# print X[X.shape[0]/2:]
	# exit()
	
	
	SF=D[:,:,:,3]

	a1=ucell[0]
	a2=ucell[1]
	a3=ucell[2]
	
	b1=(np.cross(a2,a3))/(np.dot(a1,np.cross(a2,a3)))#
	b2=(np.cross(a3,a1))/(np.dot(a2,np.cross(a3,a1)))#*2.0*math.pi
	b3=(np.cross(a1,a2))/(np.dot(a3,np.cross(a1,a2)))#*2.0*math.pi 
	
	Dnew=np.zeros_like(D)
	
	for ix in trange(D.shape[0]):			
		Dnew[ix,:,:,0:3]+=X[ix]*b1  #(X[ix]-X[X.shape[0]/2])*b1
		
	for iy in trange(D.shape[1]):			
		Dnew[:,iy,:,0:3]+=Y[iy]*b2  #(Y[iy]-Y[Y.shape[0]/2])*b2
		
	for iz in trange(D.shape[2]):			
		Dnew[:,:,iz,0:3]+=Z[iz]*b3  #(Z[iz]-Z[Z.shape[0]/2])*b3
	
	D[...,:3]=Dnew[...,:3]
	

	K_ES=2.0*math.pi/wavelength_angstroms  #calculate k for incident xrays in inverse angstroms
	
	#https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html#scipy.interpolate.RegularGridInterpolator
	#Notes
	#Contrary to LinearNDInterpolator and NearestNDInterpolator, RegularGridInterpolator class avoids expensive triangulation of the input data by taking advantage of the regular grid structure.
	#this is why this style of interpolation is so slow
	
	
	XGD=D[:,:,:,0]  #X spatial grid view
	YGD=D[:,:,:,1]
	ZGD=D[:,:,:,2]
	VGD=D[:,:,:,3]
	
	DC=D[:,:,:,0:3]
	
	DR=DC.reshape(DC.size/3,3)
		
	#check if fast interpolation can be used
	lbuf=True
	for i in xrange(3):
		for j in xrange(i+1,3):
			if ucell[i,j]>FP_THRESHOLD or ucell[j,i]>FP_THRESHOLD:	#floating point rounding
				lbuf=False
				
	print "Interpolating grid..."
	
	print ucell
	
	if ucell[0,0]==ucell[1,1] and ucell[0,0]==ucell[2,2] and lbuf:
	# if True:
		print "using fast interpolation for orthorhombic cell"
		ES = RegularGridInterpolator((X, Y, Z), SF,bounds_error=False)  #can be used in case of orthorhombic unit cell
	else:
		print "Interpolating non-orthorhombic cell"
		dtime=480.0*XGD.size/(98*98*99)  #empirical time estimate
		
		print "interpolation time estimate: ", round(dtime/60,1) , " minutes, finishing around ", (datetime.datetime.now()+datetime.timedelta(seconds=dtime)).strftime('%I:%M %p')
		#ES = griddata((XGD.reshape(XGD.size/3,3), YGD.reshape(YGD.size/3,3), ZGD.reshape(ZGD.size/3,3)), VGD.reshape(ZGD.size/3,3), (grid_x,grid_y,grid_z),method='linear')
		start=time.time()
		#ES = griddata(DC.reshape(DC.size/3,3), VGD.ravel(), (grid_x,grid_y,grid_z),method='linear')  #485s
		coords=zip(XGD.ravel(),YGD.ravel(),ZGD.ravel())
		
		if False:
			ES= LinearNDInterpolator(coords,VGD.ravel())
		else:
			ES= NearestNDInterpolator(coords,VGD.ravel())		#testing this function
		
		
		end=time.time()	
		print "interpolation finished, taking", end-start, "seconds"


	xyzpts=np.asarray([])
	print "setting up points for radial integration"
	#There's a better/faster way of doing this, but it doesn't take too much time as-is
	Scale=1
	
	if False:
		for ix in trange(D.shape[0]):
			#xsq=X[ix]**2.0
			for iy in xrange(D.shape[1]):
				#r=np.sqrt(Xsq+Y[iy]**2.0)
				for iz in xrange(D.shape[2]):
					#xyzpts.append((X[ix],Y[iy],Z[iz]))
					xyzpts.append((D[ix,iy,iz,0],D[ix,iy,iz,1],D[ix,iy,iz,2]))
					#rzpts.append((X[ix]*np.cos(theta),Y[iy]*np.cos(theta),K_ES*(1.0-np.cos(theta))))
	else:
		XPTS=np.linspace(D[0,0,0,0],D[-1,0,0,0],Scale*D.shape[0],dtype=np.float16)
		YPTS=np.linspace(D[0,0,0,1],D[0,-1,0,1],Scale*D.shape[1],dtype=np.float16)
		ZPTS=np.linspace(D[0,0,0,2],D[0,0,-1,2],Scale*D.shape[2],dtype=np.float16)
		print "mesh"
		xyzpts=np.meshgrid(XPTS,YPTS,ZPTS)
		print "stack"
		xyzpts=np.stack(xyzpts,-1).reshape(-1,3)
		print "done"
		# print type(xyzpts)
		# print xyzpts.shape
		# exit()
		# for ix in tqdm.tqdm(np.nditer(XPTS)):
			# lbuf=[]
			# for iy in np.nditer(YPTS):
				# for iz in np.nditer(ZPTS):
					# lbuf.append((ix,iy,iz))
			# xyzpts=np.append(xyzpts,np.asarray(lbuf))
	
		
		
		
		
		# xyzpts=np.vstack((XPTS,YPTS,ZPTS))	#make an array from these

	
		# for ix in trange(D.shape[0]*Scale):
			#xsq=X[ix]**2.0
			# for iy in xrange(D.shape[1]*Scale):
				#r=np.sqrt(Xsq+Y[iy]**2.0)
				# for iz in xrange(D.shape[2]*Scale):
					#xyzpts.append((X[ix],Y[iy],Z[iz]))
					# xyzpts.append((D[ix,iy,iz,0],D[ix,iy,iz,1],D[ix,iy,iz,2]))
		# Ninterp=3.0
		# for ix in trange(D.shape[0]*Ninterp):
			# xval=(ix/Ninterp)*D[]
			# for iy in xrange(D.shape[1]*Ninterp):
				# for iz in xrange(D.shape[2]*Ninterp):
					# xyzpts.append((D[ix,iy,iz,0],D[ix,iy,iz,1],D[ix,iy,iz,2]))




	# xyzpts=np.asarray(xyzpts)
	# print xyzpts.shape
	# exit()
	print xyzpts.shape
	
	NSP=20
	NSP=np.minimum(NSP,xyzpts.shape[0])		#split into at most 20 chunks before processing to limit memory usage
	
	xyzpieces=np.array_split(xyzpts,NSP)
	# print type(xyzpieces)
	# exit()
	EWDxyz=np.asarray([])
	print "interpolating"
	for i in tqdm.tqdm(xyzpieces):
		# print "A"
		buf=ES(i)
		# print 
		# print i.shape
		# print "B"
		EWDxyz=np.append(EWDxyz,buf,axis=0)
	#EWDxyz=ES(xyzpts)
	print "EWD done"
	
	# print type(ES)
	# exit()
	
	#print rzpts.shape
	
	#EWDxyz=EWDxyz.reshape(xyzpts.shape[0]/D.shape[2],D.shape[2])
	rpts=np.sqrt(xyzpts[:,0]**2.0+xyzpts[:,1]**2.0)
		
	
	Hcount,XEC,YEC=np.histogram2d(rpts,xyzpts[:,2],bins=(XBNSRD,ZBNS))
	
	Hval,XEV,YEV=np.histogram2d(rpts,xyzpts[:,2],weights=EWDxyz,normed=False,bins=(XBNSRD,ZBNS))
	
	# Hval,XEV,YEV=np.histogram2d(rpts,xyzpts[:,2],normed=False,bins=(XBNSRD,Z))

	switch1=True
	
	
	if switch1:
		Hcount=np.where(Hcount==0,1,Hcount)
		
	Hrz=Hval/Hcount
	
	# print Hrz[:,Hrz.shape[1]/2+4]
	# raw_input("Hrz")
	
	if not switch1:
		Hrz = np.ma.masked_invalid(Hrz)
	
	S1=np.sum(Hrz)
	S3=np.sum(Hrz[Hrz.shape[0]/2,:])
	print S3
	# exit()
	
	Condition1=False		#Need to figure this out-when should this be true?
	
	if Condition1:
		for ir in xrange(1,Hrz.shape[0]/2-1):
			Hrz[-ir+Hrz.shape[0]/2,:]=Hrz[ir+Hrz.shape[0]/2,:]   #this needs to be tested for both even and odd numbers of bins 
	else:
		for ir in xrange(1,Hrz.shape[0]/2-1):
			Hrz[-ir+2+Hrz.shape[0]/2,:]=Hrz[ir+Hrz.shape[0]/2,:]   #this needs to be tested for both even and odd numbers of bins 
	S2=np.sum(Hrz)
	# print S2/S1
	# print (S2+2.*S3)/S1
	
	# exit()
	XMG, YMG = np.meshgrid(XEV, YEV)#-eyev)
	#dx1=0.0#XMG[0,1]-XMG[0,0]
	
	print dx1
	print dx1/2.0
	# raw_input("inpt")
	
	print XMG[0]
	# for i in xrange(1,XMG.shape[1]):
		# print XMG[0,i]-XMG[0,i-1]
	
	# exit()
	# print XMG[:-1,0]
	# exit()
	
	plt.pcolormesh(XMG[:-1,:]-dx1/2.0, YMG[:-1,:], Hrz.T,vmin=0.0,vmax=np.amax(Hrz))
	plt.savefig(path+"rzplot_1"+format,dpi=DPI)
	plt.clf()
	# exit()
	
	
	
	mn=np.amin(Hrz[np.nonzero(Hrz)])
	Hbuf=np.where(Hrz>0.0,Hrz,mn)
	# print np.amin(Hbuf),mn
	# exit()
	Log_HRZ=np.log10(Hbuf)
		
	plt.pcolormesh(XMG[:-1,:]-dx1/2.0, YMG[:-1,:], Log_HRZ.T,vmin=np.amin(Log_HRZ),vmax=np.amax(Log_HRZ),cmap='nipy_spectral')
	
	plt.colorbar()
	plt.savefig(path+"_log_rzplot"+format,dpi=DPI)
	plt.clf()
	
	
	Nx=D.shape[0]
	Ny=D.shape[1]
	Nz=D.shape[2]
	
	
#==============flat and Ewald-corrected plots=================
	
	# def get_pts(D,i0,i1,i2cross):
		# for i in trange(D.shape[i0]):
			# for j in xrange(D.shape[i1]):
				# p1=D[ix,iy,i2cross,0]
	
	
	
	
	# make_flat_plot(D,[0.0],xrange(D.shape[0]),[])
	
	
	xypts=[]
	xyflat=[]
	print "xy"
	
	for ix in trange(D.shape[0]):
		for iy in xrange(D.shape[1]):
			xp=D[ix,iy,Nz/2,0]
			yp=D[ix,iy,Nz/2,1]
			
			
			# if np.abs(xp)<0.2 and np.abs(yp)<0.2:
				# print xp,yp
			
			# zp=D[ix,iy,Nz/2,2]
			# zp=0.0
			# print zp
			# exit()
			theta=np.arctan(np.sqrt(xp**2.0+yp**2.0)/K_ES)
			xypts.append((xp*np.cos(theta),yp*np.cos(theta),K_ES*(1.0-np.cos(theta))))
			xyflat.append((xp,yp,0.0))
			# xyflat.append((xp,yp,0.0))
			
			
	xzpts =[]
	xzflat=[]
	
	print "xz"
	
	for ix in trange(D.shape[0]):
		for iz in xrange(D.shape[2]):
			xp=D[ix,Ny/2,iz,0]
			zp=D[ix,Ny/2,iz,2]
			theta=np.arctan(np.sqrt(xp**2.0+yp**2.0)/K_ES)
			xzpts.append((xp*np.cos(theta),K_ES*(1.0-np.cos(theta)),zp*np.cos(theta)))
			xzflat.append((xp,0.0,zp))

	
	

					
	
	print "yz"
	yzpts=[]
	yzflat=[]
	
	print D.shape
	# exit()
	
	for iy in trange(D.shape[1]):
		# ysq=Y[iy]**2.0
		for iz in xrange(D.shape[2]):
			# yp=D[Nz/2,iy,iz,1]
			# zp=D[Nz/2,iy,iz,2]
			yp=D[Nx/2,iy,iz,1]		#fixed index
			zp=D[Nx/2,iy,iz,2]
			
			
			
			
			theta=np.arctan(np.sqrt(yp**2.0+zp**2.0)/K_ES)
			yzpts.append((K_ES*(1.0-np.cos(theta)),yp*np.cos(theta),zp*np.cos(theta)))
			
			yzflat.append((0.0,yp,zp))
	
	xypts=np.asarray(xypts)
	xzpts=np.asarray(xzpts)
	yzpts=np.asarray(yzpts)
	
	xyflat=np.asarray(xyflat)
	xzflat=np.asarray(xzflat)
	yzflat=np.asarray(yzflat)
		
	EWDxy=ES(xypts)
	EWDxz=ES(xzpts)
	EWDyz=ES(yzpts)
	
	EWDxyflat=ES(xyflat)
	EWDxzflat=ES(xzflat)
	EWDyzflat=ES(yzflat)
	
	
	
	EWDxy=EWDxy.reshape(D.shape[0],D.shape[1])
	EWDxz=EWDxz.reshape(D.shape[0],D.shape[2])
	EWDyz=EWDyz.reshape(D.shape[1],D.shape[2])
	
	EWDxyflat=EWDxyflat.reshape(D.shape[0],D.shape[1])
	EWDxzflat=EWDxzflat.reshape(D.shape[0],D.shape[2])
	EWDyzflat=EWDyzflat.reshape(D.shape[1],D.shape[2])
	
	
	
	title="Ewald Corrected Structure Factor \n $\lambda=$"+str(wavelength_angstroms)+" $\AA$   $k_{ew}=$"+str(round(K_ES,2))+" $\AA^{-1}$"
	ltitle='log ' + title
	
	xlab='x ('+units + ")"
	ylab='y ('+units + ")"
	zlab='z ('+units + ")"
	
	fname="Ewald_"	

	

	iz=0
	plt.suptitle(title)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.contourf(D[:,:,iz,0],D[:,:,iz,1],EWDxy,contours,**kwargs)
	# plt.contourf(XGT,YGT,EWDxy,contours,**kwargs)
	plt.savefig(path+fname+"xy"+str(iz)+format,dpi=DPI)
	plt.clf()

	

	
	lax=['x','y','z']
	
	ewlab="Ewald"
	flab="Flat"

	
	iax1=0
	iax2=1
	

	
	EWDxy=np.ma.masked_invalid(EWDxy)
	EWDxyflat=np.ma.masked_invalid(EWDxyflat)
	
	EWDxz=np.ma.masked_invalid(EWDxz)
	EWDxzflat=np.ma.masked_invalid(EWDxzflat)
	
	EWDyz=np.ma.masked_invalid(EWDyz)
	EWDyzflat=np.ma.masked_invalid(EWDyzflat)
	
	#kwargs={'cmap':'jet'}
	# vx=np.amax(D)
	# vn=np.amin(D)
	
	
	if PLOT_EWALDS:
		csplot_wlog(D[:,:,Nz/2+1,iax1],D[:,:,Nz/2+1,iax2],EWDxy,    contours,ewlab,lax[iax1],lax[iax2],**kwargs)
	csplot_wlog(D[:,:,Nz/2+1,iax1],D[:,:,Nz/2+1,iax2],EWDxyflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)
	
	iax1=0
	iax2=2
	if PLOT_EWALDS:
		csplot_wlog(D[:,Ny/2,:,iax1],D[:,Ny/2,:,iax2],EWDxz,    contours,ewlab,lax[iax1],lax[iax2],**kwargs)
	csplot_wlog(D[:,Ny/2,:,iax1],D[:,Ny/2,:,iax2],EWDxzflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)
	
	iax1=1
	iax2=2
	if PLOT_EWALDS:
		csplot_wlog(D[Nx/2,:,:,iax1],D[Nx/2,:,:,iax2],EWDyz,    contours,ewlab,lax[iax1],lax[iax2],**kwargs)
	csplot_wlog(D[Nx/2,:,:,iax1],D[Nx/2,:,:,iax2],EWDyzflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)

	
def Plot_Ewald_triclinic_old(D,wavelength_angstroms,ucell,**kwargs):  #pass full 3d data,SF,wavelength in angstroms
	
	# print D.shape

	# for i in xrange(4):
		# pl(i,np.amax(D[:,:,:,i]))
	
	# print "exiting"
	# exit()
	
	if not os.path.exists(path):
		os.makedirs(path)
		
	X =D[:,0,0,0]
	Y =D[0,:,0,1]
	Z =D[0,0,:,2]
	SF=D[:,:,:,3]
	#pli('Y0',Y)
	a1=ucell[0]
	a2=ucell[1]
	a3=ucell[2]
	
	b1=(np.cross(a2,a3))/(np.dot(a1,np.cross(a2,a3)))#*2.0*math.pi
	b2=(np.cross(a3,a1))/(np.dot(a2,np.cross(a3,a1)))#*2.0*math.pi
	b3=(np.cross(a1,a2))/(np.dot(a3,np.cross(a1,a2)))#*2.0*math.pi 
	
	Dnew=np.zeros_like(D)
	
	for ix in trange(D.shape[0]):			
		Dnew[ix,:,:,0:3]+=X[ix]*b1  #(X[ix]-X[X.shape[0]/2])*b1
		
	for iy in trange(D.shape[1]):			
		Dnew[:,iy,:,0:3]+=Y[iy]*b2  #(Y[iy]-Y[Y.shape[0]/2])*b2
		
	for iz in trange(D.shape[2]):			
		Dnew[:,:,iz,0:3]+=Z[iz]*b3  #(Z[iz]-Z[Z.shape[0]/2])*b3
	
	D=Dnew
	
	#pli('Y',Y)
	
	# X =Dnew[:,0,0,0]
	# Y =Dnew[0,:,0,1]
	# Z =Dnew[0,0,:,2]
	K_ES=2.0*math.pi/wavelength_angstroms  #calculate k for incident xrays in inverse angstroms
	
	ES = RegularGridInterpolator((X, Y, Z), SF,bounds_error=False)
	
	pli('warning! RegularGridInterpolator will produce incorrect results on a triclinic grid.  This issue is being worked on.')
	
	#angle averaging
	xyzpts=[]
	print "setting up points for radial integration"
	#There's a better/faster way of doing this, but it doesn't take too much time as-is
	if True:
		for ix in trange(D.shape[0]):
			#xsq=X[ix]**2.0
			for iy in xrange(D.shape[1]):
				#r=np.sqrt(Xsq+Y[iy]**2.0)
				for iz in xrange(D.shape[2]):
					#xyzpts.append((X[ix],Y[iy],Z[iz]))
					xyzpts.append((D[ix,iy,iz,0],D[ix,iy,iz,1],D[ix,iy,iz,2]))
		
					#rzpts.append((X[ix]*np.cos(theta),Y[iy]*np.cos(theta),K_ES*(1.0-np.cos(theta))))
	else:
		pass
		# Ninterp=3.0
		# for ix in trange(D.shape[0]*Ninterp):
			# xval=(ix/Ninterp)*D[]
			# for iy in xrange(D.shape[1]*Ninterp):
				# for iz in xrange(D.shape[2]*Ninterp):
					# xyzpts.append((D[ix,iy,iz,0],D[ix,iy,iz,1],D[ix,iy,iz,2]))




	xyzpts=np.asarray(xyzpts)
	EWDxyz=ES(xyzpts)
	
	#print rzpts.shape
	
	#EWDxyz=EWDxyz.reshape(xyzpts.shape[0]/D.shape[2],D.shape[2])
	rpts=np.sqrt(xyzpts[:,0]**2.0+xyzpts[:,1]**2.0)
	# for i in xrange(rpts.shape[0]):
	
		# print EWDxyz[i],rpts[i]
	#  
	
	
	# print rpts.shape
	#print rpts.reshape(rpts.shape[0]).shape
	
	#rpts=np.reshape[rpts.shape[0]]
	
	
	#print [xyzpts[i,0] ]
	#for i in xrange(xyzpts.shape[0]):
		#print xyzpts[i,0] 
	#a=raw_input("input")
	
	
	#print rpts
	#a=raw_input("input")
	#print xyzpts[:,2]
	#a=raw_input("input")
	Hcount,XEC,YEC=np.histogram2d(rpts,xyzpts[:,2],bins=(X,Z))
	
	# pl('hc',Hcount)
	# pl('hc shape',Hcount.shape)
	# pli('XEC',XEC)
	# pli('YEC',YEC)
	
	# for i in xrange(Hcount.shape[0]):
		# print i,np.amax(Hcount[i,:])
	# a=raw_input("input")
	# for i in xrange(Hcount.shape[1]):
		# print i,np.amax(Hcount[:,i])
	# 
	
	# print Hcount
	# print np.amax(Hcount)
	# 
	
	# for i in xrange(Hcount.shape[0]):
		# print i,np.amax(Hcount[i,:])
	# print Hcount.shape[0]
#		pli('hc '+str(i),Hcount[i,:])
	
	# plt.plot(Hcount[Hcount.shape[0]/2,:])
	# plt.show()
	# exit()
	
	
	
	Hval,XEV,YEV=np.histogram2d(rpts,xyzpts[:,2],weights=EWDxyz,normed=True,bins=(X,Z))
	# a=raw_input("hvalx")
	# for i in xrange(Hval.shape[0]):
		# print i,np.amax(Hval[i,:])
	# a=raw_input("hvalz")
	# for i in xrange(Hval.shape[1]):
		# print i,np.amax(Hval[:,i])
	
	Hrz=Hval/Hcount
	Hrz = np.ma.masked_invalid(Hrz)
	
	#for ir in xrange(Hrz.shape[0]):
	for ir in xrange(0,Hrz.shape[0]/2):
		#Hrz[:,-iz+Hrz.shape[1]/2]=Hrz[:,iz+Hrz.shape[1]/2]
		Hrz[-ir+Hrz.shape[0]/2,:]=Hrz[ir+1+Hrz.shape[0]/2,:]   #this needs to be tested for both even and odd numbers of bins 
		#Hrz[ir+Hrz.shape[1]/2]=1.0
	# print Hcount
	# a=raw_input("input")
	# print Hval
	# a=raw_input("input")
	# print Hrz
	#ax = plt.add_subplot(111, title='pcolormesh: actual edges',
	#	aspect='equal')
	
	
	# exev=(XEV[1]-XEV[0])/0.5
	# eyev=(YEV[1]-YEV[0])/0.5
	XMG, YMG = np.meshgrid(XEV, YEV)#-eyev)
	# XMG, YMG = np.meshgrid(XEV, YEV)
	# print XMG
	# print YMG
	#  
	#tst=np.where(YMG[:-1,:-1]>0.0,np.log(Hrz),0.0)
	#plt.pcolormesh(XMG[:-1], YMG[:-1], np.log(Hrz))
	
	
	#fig, ax = plt.subplots()
	# print XMG[:-1,:]
	# print XMG[:-1,:].shape
	# exit()
	plt.pcolormesh(XMG[:-1,:], YMG[:-1,:], Hrz.T,vmin=0.0,vmax=np.amax(Hrz))
	
	
	
	#plt.contourf(Hrz.T)
	#plt.savefig(path+"rzplt2")
	
	
	
	#cax=ax.pcolormesh(XMG[:-1,:], YMG[:-1,:], Hrz.T)#,vmin=0.0)#,vmax=np.amax(Hrz))
	#cax=ax.contourf(XMG[:-1,:-1], YMG[:-1,:-1], Hrz.T,**kwargs)
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
	#cbar = fig.colorbar(cax)
	plt.savefig(path+"rzplot"+format,dpi=DPI)
	
	plt.clf()
	
	plt.pcolormesh(XMG[:-1,:], YMG[:-1,:], np.log10(Hrz.T),vmin=np.amin(np.log10(Hrz)),vmax=np.amax(np.log10(Hrz)))
	
	#plt.contourf(Hrz.T)
	#plt.savefig(path+"rzplt2")
	# 
	
	
	#cax=ax.pcolormesh(XMG[:-1,:], YMG[:-1,:], Hrz.T)#,vmin=0.0)#,vmax=np.amax(Hrz))
	#cax=ax.contourf(XMG[:-1,:-1], YMG[:-1,:-1], Hrz.T,**kwargs)
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
	#cbar = fig.colorbar(cax)
	plt.savefig(path+"_log_rzplot"+format,dpi=DPI)
	
	plt.clf()
	
	
	# 
	
	# fig, ax = plt.subplots()
	# cax=ax.pcolormesh(XMG[:-1,:], YMG[:-1,:], np.log(Hrz.T))
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
	# cbar = fig.colorbar(cax)
	# plt.savefig(path+"rzplot_log"+format,dpi=DPI)
	# plt.clf()
	
	
	
	
	
	
	
	# print Hcount.shape
	# print "HC"
	#  
	
	
	# print EWDxyz.shape
	# 
	
	
	xypts=[]
	xyflat=[]
	for ix in xrange(D.shape[0]):
		xsq=X[ix]**2.0
		for iy in xrange(D.shape[1]):
			theta=np.arctan(np.sqrt(xsq+Y[iy]**2.0)/K_ES)
			xypts.append((X[ix]*np.cos(theta),Y[iy]*np.cos(theta),K_ES*(1.0-np.cos(theta))))
			xyflat.append((X[ix],Y[iy],0.0))
			
			
	xzpts =[]
	xzflat=[]
	for ix in xrange(D.shape[0]):
		xsq=X[ix]**2.0
		for iz in xrange(D.shape[2]):
			theta=np.arctan(np.sqrt(xsq+Z[iz]**2.0)/K_ES)
			xzpts.append((X[ix]*np.cos(theta),K_ES*(1.0-np.cos(theta)),Z[iz]*np.cos(theta)))
			xzflat.append((X[ix],0.0,Z[iz]))
	
	yzpts=[]
	yzflat=[]
	for iy in xrange(D.shape[1]):
		ysq=Y[iy]**2.0
		for iz in xrange(D.shape[2]):
			theta=np.arctan(np.sqrt(ysq+Z[iz]**2.0)/K_ES)
			yzpts.append((K_ES*(1.0-np.cos(theta)),Y[iy]*np.cos(theta),Z[iz]*np.cos(theta)))
			yzflat.append((0.0,Y[iy],Z[iz]))
	
	xypts=np.asarray(xypts)
	xzpts=np.asarray(xzpts)
	yzpts=np.asarray(yzpts)
	
	xyflat=np.asarray(xyflat)
	xzflat=np.asarray(xzflat)
	yzflat=np.asarray(yzflat)
	
	
	
	EWDxy=ES(xypts)
	EWDxz=ES(xzpts)
	EWDyz=ES(yzpts)
	
	EWDxyflat=ES(xyflat)
	EWDxzflat=ES(xzflat)
	EWDyzflat=ES(yzflat)
	
	
	
	EWDxy=EWDxy.reshape(D.shape[0],D.shape[1])
	EWDxz=EWDxz.reshape(D.shape[0],D.shape[2])
	EWDyz=EWDyz.reshape(D.shape[1],D.shape[2])
	
	EWDxyflat=EWDxyflat.reshape(D.shape[0],D.shape[1])
	EWDxzflat=EWDxzflat.reshape(D.shape[0],D.shape[2])
	EWDyzflat=EWDyzflat.reshape(D.shape[1],D.shape[2])
	
	
	
	title="Ewald Corrected Structure Factor \n $\lambda=$"+str(wavelength_angstroms)+" $\AA$   $k_{ew}=$"+str(round(K_ES,2))+" $\AA^{-1}$"
	ltitle='log ' + title
	
	xlab='x ('+units + ")"
	ylab='y ('+units + ")"
	zlab='z ('+units + ")"
	
	fname="Ewald_"	
	
	plt.suptitle(title)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	plt.contourf(D[:,:,0,0],D[:,:,0,1],EWDxy,contours,**kwargs)
	plt.savefig(path+fname+"xy"+format,dpi=DPI)
	plt.clf()

	
	Nx=D.shape[0]
	Ny=D.shape[1]
	Nz=D.shape[2]
	
	
	lax=['x','y','z']
	
	ewlab="Ewald"
	flab="Flat"
	
	iax1=0
	iax2=1
	
	
	if PLOT_EWALDS:
		csplot_wlog(D[:,:,Nz/2+1,iax1],D[:,:,Nz/2+1,iax2],EWDxy,    contours,ewlab,lax[iax1],lax[iax2],**kwargs)
	csplot_wlog(D[:,:,Nz/2+1,iax1],D[:,:,Nz/2+1,iax2],EWDxyflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)
	
	iax1=0
	iax2=2
	if PLOT_EWALDS:
		csplot_wlog(D[:,Ny/2,:,iax1],D[:,Ny/2,:,iax2],EWDxz,    contours,ewlab,lax[iax1],lax[iax2],**kwargs)
	csplot_wlog(D[:,Ny/2,:,iax1],D[:,Ny/2,:,iax2],EWDxzflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)
	
	iax1=1
	iax2=2
	if PLOT_EWALDS:
		csplot_wlog(D[Nx/2,:,:,iax1],D[Nx/2,:,:,iax2],EWDyz,    contours,ewlab,lax[iax1],lax[iax2],**kwargs)
	csplot_wlog(D[Nx/2,:,:,iax1],D[Nx/2,:,:,iax2],EWDyzflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)
	
	
	
	#csplot(D[:,:,0,0],D[:,0,:,2],EWDxy,contours,"Ewald","x","z",**kwargs)
	
	
	
	#  
	#csplot(X,Y,Z,contours,xlab,ylab,**kwargs):
	
	
	# plt.suptitle(ltitle)
	# plt.xlabel(xlab)
	# plt.ylabel(ylab)
	# plt.contourf(D[:,:,0,0],D[:,:,0,1],np.log(EWDxy),contours,**kwargs)
	# plt.savefig(path+fname+"xylog"+format,dpi=DPI)
	# plt.clf()
	
	# plt.suptitle(title)
	# plt.xlabel(xlab)
	# plt.ylabel(zlab)
	# plt.contourf(D[:,0,:,0],D[:,0,:,2],EWDxz,contours,**kwargs)
	# plt.savefig(path+fname+"xz"+format,dpi=DPI)
	# plt.clf()
	
	# plt.suptitle(ltitle)
	# plt.xlabel(xlab)
	# plt.ylabel(zlab)
	# plt.contourf(D[:,0,:,0],D[:,0,:,2],np.log(EWDxz),contours,**kwargs)
	# plt.savefig(path+fname+"xzlog"+format,dpi=DPI)
	# plt.clf()
	
	# plt.suptitle(title)
	# plt.xlabel(ylab)
	# plt.ylabel(zlab)
	# plt.contourf(D[0,:,:,1],D[0,:,:,2],EWDyz,contours,**kwargs)
	# plt.savefig(path+fname+"yz"+format,dpi=DPI)
	# plt.clf()
	
	# plt.suptitle(ltitle)
	# plt.xlabel(ylab)
	# plt.ylabel(zlab)
	# plt.contourf(D[0,:,:,1],D[0,:,:,2],np.log(EWDyz),contours,**kwargs)
	# plt.savefig(path+fname+"yzlog"+format,dpi=DPI)
	# plt.clf()
	

	
	# title="Flat Structure Factor \n $\lambda=$"+str(wavelength_angstroms)+" $\AA$   $k_{ew}=$"+str(round(K_ES,2))+" $\AA^{-1}$"
	# ltitle='log ' + title
	
	# xlab='x ('+units + ")"
	# ylab='y ('+units + ")"
	# zlab='z ('+units + ")"
	
	# fname="Flat_"	
	
	# plt.suptitle(title)
	# plt.xlabel(xlab)
	# plt.ylabel(ylab)
	# plt.contourf(D[:,:,0,0],D[:,:,0,1],EWDxyflat,contours,**kwargs)
	# plt.savefig(path+fname+"xy"+format,dpi=DPI)
	# plt.clf()
	
	# plt.suptitle(ltitle)
	# plt.xlabel(xlab)
	# plt.ylabel(ylab)
	# plt.contourf(D[:,:,0,0],D[:,:,0,1],np.log(EWDxyflat),contours,**kwargs)
	# plt.savefig(path+fname+"xylog"+format,dpi=DPI)
	# plt.clf()
	
	# plt.suptitle(title)
	# plt.xlabel(xlab)
	# plt.ylabel(zlab)
	# plt.contourf(D[:,0,:,0],D[:,0,:,2],EWDxzflat,contours,**kwargs)
	# plt.savefig(path+fname+"xz"+format,dpi=DPI)
	# plt.clf()
	
	# plt.suptitle(ltitle)
	# plt.xlabel(xlab)
	# plt.ylabel(zlab)
	# plt.contourf(D[:,0,:,0],D[:,0,:,2],np.log(EWDxzflat),contours,**kwargs)
	# plt.savefig(path+fname+"xzlog"+format,dpi=DPI)
	# plt.clf()
	
	# plt.suptitle(title)
	# plt.xlabel(ylab)
	# plt.ylabel(zlab)
	# plt.contourf(D[0,:,:,1],D[0,:,:,2],EWDyzflat,contours,**kwargs)
	# plt.savefig(path+fname+"yz"+format,dpi=DPI)
	# plt.clf()
	
	# plt.suptitle(ltitle)
	# plt.xlabel(ylab)
	# plt.ylabel(zlab)
	# plt.contourf(D[0,:,:,1],D[0,:,:,2],np.log(EWDyzflat),contours,**kwargs)
	# plt.savefig(path+fname+"yzlog"+format,dpi=DPI)
	# plt.clf()
	
	
	
