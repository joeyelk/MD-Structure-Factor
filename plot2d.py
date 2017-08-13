import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
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

normplot=1 			#normalize plots

title_fontsize=9

path=""				

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
	csplot(X,Y,Z,contours,lab,xlab,ylab,**kwargs)
	csplot(X,Y,np.log(Z),contours,"log_"+lab,xlab,ylab,**kwargs)


def csplot(X,Y,Z,contours,lab,xlab,ylab,**kwargs):
	
	title=lab+" S("+xlab+","+ylab+")"
	fname=lab+"_"+xlab+"_"+ylab
	fig, ax = plt.subplots()
	plt.suptitle(title)
	plt.xlabel(xlab)
	plt.ylabel(ylab)
	
	if normplot==1:
		cax=plt.contourf(X,Y,Z/np.amax(Z),contours,vmin=0.0,vmax=1.0,**kwargs)
	else:
		cax=plt.contourf(X,Y,Z,contours,**kwargs)
	

	#ax.set_aspect((np.amax(Y)-np.amin(Y))/(np.amax(X)-np.amin(X)))
	# ax.set_aspect('auto')
	#cbar = fig.colorbar(cax)
	
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

def Plot_Ewald_triclinic(D,wavelength_angstroms,ucell,**kwargs):  #pass full 3d data,SF,wavelength in angstroms
		
	if not os.path.exists(path):
		os.makedirs(path)
		
	X =D[:,0,0,0].copy()  
	Y =D[0,:,0,1].copy()
	Z =D[0,0,:,2].copy()
	
	
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
			if ucell[i,j]!=0 or ucell[j,i]!=0:
				lbuf=False
				
	print "Interpolating grid..."
	
	if ucell[0,0]==ucell[1,1] and ucell[0,0]==ucell[2,2] and lbuf:
	# if True:
		ES = RegularGridInterpolator((X, Y, Z), SF,bounds_error=False)  #can be used in case of orthorhombic unit cell
	else:
		
		dtime=480.0*XGD.size/(98*98*99)  #empirical time estimate
		
		print "interpolation time estimate: ", round(dtime/60,1) , " minutes, finishing around ", (datetime.datetime.now()+datetime.timedelta(seconds=dtime)).strftime('%I:%M %p')
		#ES = griddata((XGD.reshape(XGD.size/3,3), YGD.reshape(YGD.size/3,3), ZGD.reshape(ZGD.size/3,3)), VGD.reshape(ZGD.size/3,3), (grid_x,grid_y,grid_z),method='linear')
		start=time.time()
		#ES = griddata(DC.reshape(DC.size/3,3), VGD.ravel(), (grid_x,grid_y,grid_z),method='linear')  #485s
		coords=zip(XGD.ravel(),YGD.ravel(),ZGD.ravel())
		
		ES= LinearNDInterpolator(coords,VGD.ravel())
		end=time.time()	
		print "interpolation finished, taking", end-start, "seconds"


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
	
	Hcount,XEC,YEC=np.histogram2d(rpts,xyzpts[:,2],bins=(X,Z))
	
	Hval,XEV,YEV=np.histogram2d(rpts,xyzpts[:,2],weights=EWDxyz,normed=True,bins=(X,Z))

	switch1=True
	
	
	
	if switch1:
		Hcount=np.where(Hcount==0,1,Hcount)
		
	Hrz=Hval/Hcount
	
	if not switch1:
		Hrz = np.ma.masked_invalid(Hrz)
	
	for ir in xrange(0,Hrz.shape[0]/2-1):
	
		Hrz[-ir+Hrz.shape[0]/2,:]=Hrz[ir+2+Hrz.shape[0]/2,:]   #this needs to be tested for both even and odd numbers of bins 
	
	XMG, YMG = np.meshgrid(XEV, YEV)#-eyev)
	dx1=XMG[0,1]-XMG[0,0]
	
	print dx1
	print dx1/2.0
	# raw_input("inpt")
	
	plt.pcolormesh(XMG[:-1,:]-dx1/2.0, YMG[:-1,:], Hrz.T,vmin=0.0,vmax=np.amax(Hrz))
	plt.savefig(path+"rzplot"+format,dpi=DPI)
	plt.clf()
	
	plt.pcolormesh(XMG[:-1,:]-dx1/2.0, YMG[:-1,:], np.log10(Hrz.T),vmin=np.amin(np.log10(Hrz)),vmax=np.amax(np.log10(Hrz)))
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
	
	xypts=[]
	xyflat=[]
	print "xy"
	
	for ix in trange(D.shape[0]):
		for iy in xrange(D.shape[1]):
			xp=D[ix,iy,Nz/2,0]
			yp=D[ix,iy,Nz/2,1]
			theta=np.arctan(np.sqrt(xp**2.0+yp**2.0)/K_ES)
			xypts.append((xp*np.cos(theta),yp*np.cos(theta),K_ES*(1.0-np.cos(theta))))
			xyflat.append((xp,yp,0.0))
			
			
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
	for iy in trange(D.shape[1]):
		# ysq=Y[iy]**2.0
		for iz in xrange(D.shape[2]):
			yp=D[Nz/2,iy,iz,1]
			zp=D[Nz/2,iy,iz,2]
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
	plt.pcolormesh(XMG[:-1,:], YMG[:-1,:], Hrz.T,vmin=0.0,vmax=np.amax(Hrz))
	
	#plt.contourf(Hrz.T)
	#plt.savefig(path+"rzplt2")
	# 
	
	
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
	
	
	
