import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import RegularGridInterpolator
import math

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
	
def Plot_Ewald_Sphere_Correction(D,wavelength_angstroms,**kwargs):  #pass full 3d data,SF,wavelength in angstroms


	if not os.path.exists(path):
		os.makedirs(path)
	
	X =D[:,0,0,0]
	Y =D[0,:,0,1]
	Z =D[0,0,:,2]
	SF=D[:,:,:,3]
	
	K_ES=2.0*math.pi/wavelength_angstroms  #calculate k for incident xrays in inverse angstroms
	
	ES = RegularGridInterpolator((X, Y, Z), SF)		
	
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
	

