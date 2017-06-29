import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import RegularGridInterpolator
import math
import tqdm
from tqdm import trange
mainlabel=""

units="$\AA^{-1}$"

xunits=units
yunits=units
zunits=units

contours=200
DPI=300
format=".png"
text="Structure Factor"

savelog=True
savelin=True

title_fontsize=9

path=""


		

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
	#exit()
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
	# B=
	



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
	
	print D.shape

	if not os.path.exists(path):
		os.makedirs(path)
		
	X =D[:,0,0,0]
	Y =D[0,:,0,1]
	Z =D[0,0,:,2]
	SF=D[:,:,:,3]
	
	a1=ucell[0]
	a2=ucell[1]
	a3=ucell[2]
	
	b1=(np.cross(a2,a3))/(np.dot(a1,np.cross(a2,a3)))#*2.0*math.pi
	b2=(np.cross(a3,a1))/(np.dot(a2,np.cross(a3,a1)))#*2.0*math.pi
	b3=(np.cross(a1,a2))/(np.dot(a3,np.cross(a1,a2)))#*2.0*math.pi 
	

	Dnew=np.zeros_like(D)

	
	for ix in trange(D.shape[0]):			
		Dnew[ix,:,:,0:3]+=X[ix]*b1
	for iy in trange(D.shape[1]):			
		Dnew[:,iy,:,0:3]+=Y[iy]*b2
	for iz in trange(D.shape[2]):			
		Dnew[:,:,iz,0:3]+=Z[iz]*b3
	
	D=Dnew
	
	K_ES=2.0*math.pi/wavelength_angstroms  #calculate k for incident xrays in inverse angstroms
	
	ES = RegularGridInterpolator((X, Y, Z), SF,bounds_error=False)		
	
	
	#angle averaging
	xyzpts=[]
	print "xyz points"
	for ix in trange(D.shape[0]):
		#xsq=X[ix]**2.0
		for iy in xrange(D.shape[1]):
			#r=np.sqrt(Xsq+Y[iy]**2.0)
			for iz in xrange(D.shape[2]):
				xyzpts.append((X[ix],Y[iy],Z[iz]))
	
				#rzpts.append((X[ix]*np.cos(theta),Y[iy]*np.cos(theta),K_ES*(1.0-np.cos(theta))))
				
	xyzpts=np.asarray(xyzpts)
	EWDxyz=ES(xyzpts)
	#print rzpts.shape
	
	#EWDxyz=EWDxyz.reshape(xyzpts.shape[0]/D.shape[2],D.shape[2])
	rpts=np.sqrt(xyzpts[:,0]**2.0+xyzpts[:,1]**2.0)
	
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
	# for i in xrange(Hcount.shape[0]):
		# print i,np.amax(Hcount[i,:])
	# a=raw_input("input")
	# for i in xrange(Hcount.shape[1]):
		# print i,np.amax(Hcount[:,i])
	#exit()
	
	# print Hcount
	# print np.amax(Hcount)
	#exit()
	
	
	
	
	Hval,XEV,YEV=np.histogram2d(rpts,xyzpts[:,2],weights=EWDxyz,normed=True,bins=(X,Z))
	# a=raw_input("hvalx")
	# for i in xrange(Hval.shape[0]):
		# print i,np.amax(Hval[i,:])
	# a=raw_input("hvalz")
	# for i in xrange(Hval.shape[1]):
		# print i,np.amax(Hval[:,i])
	
	Hrz=Hval/Hcount
	
	#for ir in xrange(Hrz.shape[0]):
	for ir in xrange(0,Hrz.shape[0]/2):
		#Hrz[:,-iz+Hrz.shape[1]/2]=Hrz[:,iz+Hrz.shape[1]/2]
		Hrz[-ir-1+Hrz.shape[0]/2,:]=Hrz[ir+1+Hrz.shape[0]/2,:]
		#Hrz[ir+Hrz.shape[1]/2]=1.0
	# print Hcount
	# a=raw_input("input")
	# print Hval
	# a=raw_input("input")
	# print Hrz
	# exit()
	#ax = plt.add_subplot(111, title='pcolormesh: actual edges',
	#	aspect='equal')
	
	
	exev=(XEV[1]-XEV[0])/0.5
	eyev=(YEV[1]-YEV[0])/0.5
	XMG, YMG = np.meshgrid(XEV+eyev, YEV+eyev)
	# XMG, YMG = np.meshgrid(XEV, YEV)
	
	
	#tst=np.where(YMG[:-1,:-1]>0.0,np.log(Hrz),0.0)
	#plt.pcolormesh(XMG[:-1], YMG[:-1], np.log(Hrz))
		
	fig, ax = plt.subplots()
	cax=ax.pcolormesh(XMG[:-1,:], YMG[:-1,:], Hrz.T,vmin=0.0,vmax=np.amax(Hrz))
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
	cbar = fig.colorbar(cax)
	plt.savefig(path+"rzplot"+format,dpi=DPI)
	plt.clf()
	
	
	
	# fig, ax = plt.subplots()
	# cax=ax.pcolormesh(XMG[:-1,:], YMG[:-1,:], np.log(Hrz.T))
	#cbar = fig.colorbar(cax, ticks=[-1, 0, 1], orientation='horizontal')
	# cbar = fig.colorbar(cax)
	# plt.savefig(path+"rzplot_log"+format,dpi=DPI)
	# plt.clf()
	
	
	
	
	
	
	
	# print Hcount.shape
	# print "HC"
	# exit()
	
	
	# print EWDxyz.shape
	exit()
	
	
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
	

