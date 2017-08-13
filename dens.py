import numpy as np
import numpy.linalg
import math
try:
	import mayavi.mlab as mlab
except:
	print "mayavi not imported"
import matplotlib.pyplot as plt
import tqdm
from plot2d import pl

from tqdm import trange
#number of gridpoints

#Nspatialgrid=np.asarray([128,128,128])
#Nspatialgrid=128
# NSGx=128
# NSGy=128
# NSGZ=150

USE_BETTER_RESOLUTION=True  #ensure that resolution is never worse than requested
PRINT_DETAILS=True    #print spatial resolution details

RANDOM_NOISE=0		#this will override density with random noise

Nspatialgrid=np.asarray([0.0,0.0,0.0])

PRECISION=1.0E-12

dtyp=np.float32   #float64 is default.  If there are memory issues, this can be reduced

#returns a dictionary.  keys are atom symbols or types, and value is a list containing electron number and gaussian width
# for example print a["Sn"][0:1] gives you number and width for tin

# def set_Nspatialgrid(L,precision):
	# Nspatialgrid=L/precision
	

def load_radii(filename):
	out={}
	with open(filename) as f:
		lines=f.readlines()
		
	for line in lines:
		l=line.strip().split()
		out[l[1]]=(float(l[0]),float(l[2])/100.0)  
		#a list with 2 values, the number of electrons and the distribution sigma
		#divide by 100 to convert pm to angstroms
		
	return out
	
def get_borders(ad,dr,keys):  #compute width in each direction in spatial steps
	bdict={}
	for a in keys:
		bdict[a]=ad[a][1]*np.sqrt(np.log(ad[a][0]/PRECISION))/dr
	return bdict

#rescale coordinates so that cell dimensions are constant over the simulation
#returns rescaled coordinates, average length, and spatial step
def rescale(coords,dims):
	avgdims=np.average(dims,axis=0)
	a=avgdims/dims
	rc=coords

	# print rc.shape
	# exit()
	for it in xrange(rc.shape[0]):
		for i in xrange(3):
			rc[it,:,i]*=a[it,i]
			
	
	
	# dr=np.divide(avgdims,Nspatialgrid)
	# pl("avgdims",avgdims)
	# pl("Nsp",Nspatialgrid)
	# pl("dr",dr)
	# exit()
	return rc,avgdims
	

	
def remap_grid(d0,des,ori):

#remap borders onto periodic cell

	d1=np.copy(d0[ori[1]:ori[2],ori[1]:ori[2],ori[1]:ori[2]])
	for io1 in xrange(0,4,2):
		id1=2-io1
		d1[:,:,des[id1]:des[id1+1]]+=d0[ori[1]:ori[2],ori[1]:ori[2],ori[io1]:ori[io1+1]] #sides normal to z axis  (6 total)
		d1[:,des[id1]:des[id1+1],:]+=d0[ori[1]:ori[2],ori[io1]:ori[io1+1],ori[1]:ori[2]] #sides normal to y axis
		d1[des[id1]:des[id1+1],:,:]+=d0[ori[io1]:ori[io1+1],ori[1]:ori[2],ori[1]:ori[2]] #sides normal to x axis
		for io2 in xrange(0,4,2):
			id2=2-io2
			d1[:,des[id1]:des[id1+1],des[id2]:des[id2+1]]+=d0[ori[1]:ori[2],ori[io1]:ori[io1+1],ori[io2]:ori[io2+1]] #edges parallel to x axis (12 total)
			d1[des[id1]:des[id1+1],:,des[id2]:des[id2+1]]+=d0[ori[io1]:ori[io1+1],ori[1]:ori[2],ori[io2]:ori[io2+1]] #edges parallel to y axis
			d1[des[id1]:des[id1+1],des[id2]:des[id2+1],:]+=d0[ori[io1]:ori[io1+1],ori[io2]:ori[io2+1],ori[1]:ori[2]] #edges parallel to z axis
			for io3 in xrange(0,4,2):
				id3=2-io3
				d1[des[id1]:des[id1+1],des[id2]:des[id2+1],des[id2]:des[id2+1]]+=d0[ori[io1]:ori[io1+1],ori[io2]:ori[io2+1],ori[io3]:ori[io3+1]] #corners (8 total)
	return d1
	

def remap_grid_tcl(d0,des,ori):
#remap borders onto periodic cell

	# print "="*50
	# print "="*50
	# print "ori"
	# print ori
	# print "0,1"
	# print ori[0][1]

	d1=np.copy(d0[ori[0][1]:ori[0][2],ori[1][1]:ori[1][2],ori[2][1]:ori[2][2]])
	
	
	for io1 in xrange(0,4,2):
		id1=2-io1
		#Sides (6 total)
		d1[:,:,des[2][id1]:des[2][id1+1]]+=d0[ori[0][1]:ori[0][2],ori[1][1]:ori[1][2],ori[2][io1]:ori[2][io1+1]] #sides normal to z axis  
		d1[:,des[1][id1]:des[1][id1+1],:]+=d0[ori[0][1]:ori[0][2],ori[1][io1]:ori[1][io1+1],ori[2][1]:ori[2][2]] #sides normal to y axis
		d1[des[0][id1]:des[0][id1+1],:,:]+=d0[ori[0][io1]:ori[0][io1+1],ori[1][1]:ori[1][2],ori[2][1]:ori[2][2]] #sides normal to x axis
		for io2 in xrange(0,4,2):
			id2=2-io2
			#Edges (12 total)
			d1[:,des[1][id1]:des[1][id1+1],des[2][id2]:des[2][id2+1]]+=d0[ori[0][1]  :ori[0][2]    ,ori[1][io1]:ori[1][io1+1],ori[2][io2]:ori[2][io2+1]] #edges parallel to x axis 
			d1[des[0][id1]:des[0][id1+1],:,des[2][id2]:des[2][id2+1]]+=d0[ori[0][io1]:ori[0][io1+1],ori[1][1]  :ori[1][2]    ,ori[2][io2]:ori[2][io2+1]] #edges parallel to y axis
			d1[des[0][id1]:des[0][id1+1],des[1][id2]:des[1][id2+1],:]+=d0[ori[0][io1]:ori[0][io1+1],ori[1][io2]:ori[1][io2+1],ori[2][1]  :ori[2][2]] #edges parallel to z axis
			for io3 in xrange(0,4,2):
				id3=2-io3
				#corners (8 total)
				d1[des[0][id1]:des[0][id1+1],des[1][id2]:des[1][id2+1],des[2][id2]:des[2][id2+1]]+=d0[ori[0][io1]:ori[0][io1+1],ori[1][io2]:ori[1][io2+1],ori[2][io3]:ori[2][io3+1]] 
	return d1
		
	
	
	
def get_dplot(dmag):  #convert array to full 3d plottable array

	#assume dimensions of dmag are even
	
	print dmag.shape
	# exit()
	
	#dplot has odd dimensions so that zero position can be at center-removed
	dplot=np.empty([dmag.shape[0],dmag.shape[1],(dmag.shape[2])*2-1])
	ix1=dmag.shape[0]/2
	ix2=dmag.shape[0]
	iy1=dmag.shape[1]/2
	iy2=dmag.shape[1]
	iz1=dmag.shape[2]-1
	iz2=dmag.shape[2]*2-2
	
	
	
	# print "="*20
	# print "="*20
	# print ix1,ix2
	# print iy1,iy2
	# print iz1,iz2
	
	# exit()
	
	
	dplot[ix1:ix2,iy1:iy2,iz1:iz2]=dmag[  0:ix1,  0:iy1, 0:iz1]
	dplot[ 0:ix1, iy1:iy2,iz1:iz2]=dmag[ix1:ix2,  0:iy1, 0:iz1]
	dplot[ix1:ix2,  0:iy1,iz1:iz2]=dmag[  0:ix1,iy1:iy2, 0:iz1]
	dplot[ 0:ix1,   0:iy1,iz1:iz2]=dmag[ix1:ix2,iy1:iy2, 0:iz1]
		
	dplot[1:ix2-1,1:iy2-1,:iz1]=dplot[:1:-1,:1:-1,iz2:iz1:-1]
	dplot[dplot.shape[0]/2,dplot.shape[1]/2,dplot.shape[2]/2]=1.0/3.0*(dplot[dplot.shape[0]/2+1,dplot.shape[1]/2,dplot.shape[2]/2]+dplot[dplot.shape[0]/2,dplot.shape[1]/2+1,dplot.shape[2]/2]+dplot[dplot.shape[0]/2,dplot.shape[1]/2,dplot.shape[2]/2+1])
	return dplot[1:-1,1:-1,1:-1]

def compute_sf(r,L,typ,out_filename,rad,ucell,Sres):  
#compute 3d structure factor.
#arguments:
#r= 			coordinates
#L 				cell dimensions
#typ  			list of element types or names
#out_filename 	output filename to store numpy array
#rad: 			dictionary of atomic radii and numbers by type or element
#ucell:         unit cell as a 3x3 numpy array
#Sres:          Spatial resolution of the density grid in angstroms.  This setting determines the maximum scattering vector, which is 2*pi/Sres

	r,L=rescale(r,L)  #keep this inside compute_fft to allow coordinates to be rescaled to whichever numpy view is passed
	
	print "="*50
	print L
	print Sres
	# print L/Sres
	print type(Sres)
	
	
	Nspatialgrid=(L/Sres).astype(int)
	
	# exit()
	
	if USE_BETTER_RESOLUTION:  #ensure that resolution is never worse than requested  (should use ceil here instead)
		for i in xrange(3):
			if L[i]/Nspatialgrid[i]>Sres:
				Nspatialgrid[i]+=1	
	for i in xrange(3):  #ensure Nspatialgrid is even
		Nspatialgrid[i]+=Nspatialgrid[i]%2
	if PRINT_DETAILS:    #print spatial resolution details
		print "="*50
		print "="*50
		print "Unit cell has dimensions of ",L
		print "requested resolution: ",Sres, " Angstroms"
		print "Using a spatial grid of: ",Nspatialgrid
		print "actual resolution: ", L/Nspatialgrid, " Angstroms"
		print "This corresponds to maximum q vector of ",2*math.pi*Nspatialgrid/L, " inverse Angstroms"
		print "="*50
		print "="*50	
		# exit()
	
	dr=np.divide(L,Nspatialgrid)
	
	zv=[0.0,0.0,0.0] #zero vector
	
	BUFFSIZE=1000000
	print "remapping coordinates into periodic cell"
	#this must be done in case some atoms are outside the periodic cell (some MD programs save like this)
	if r.shape[1]<BUFFSIZE:	  #check if system is so large that it needs to be buffered.  
		for it in xrange(r.shape[0]):  								#looped to save memory
			r[it,...]=np.where(r[it,...]<L,r[it,...],r[it,...]-L)     								#get positions in periodic cell
			r[it,...]=np.where(r[it,...]>zv,r[it,...],r[it,...]+L)
	else:
		# for ic in trange(3):	#buffered to further save memory
			for it in trange(r.shape[0]):  								#looped to save memory
				for imin in xrange(0,r.shape[0],BUFFSIZE):
					imax=imin+BUFFSIZE
					if imax>r.shape[0]:
						imax=r.shape[0]
					
					r[it,imin:imax,:]=np.where(r[it,imin:imax,:]<L,r[it,imin:imax,:],r[it,imin:imax,:]-L)     								#get positions in periodic cell
					r[it,imin:imax,:]=np.where(r[it,imin:imax,:]>zv,r[it,imin:imax,:],r[it,imin:imax,:]+L)

	
	bdict=get_borders(rad,dr,set(typ))															#dictionary of borders by atom type
		
	Nborder=int(np.amax(bdict.values())+1)														#maximum border required to account for atoms near periodic boundaries
	
#	try:
	# print Nspatialgrid+2*Nborder,Nspatialgrid+2*Nborder,Nspatialgrid+2*Nborder
	
	print "allocating memory"
	print dtyp
	print (int(Nspatialgrid[0]),int(Nspatialgrid[1]),int(Nspatialgrid[2])/2+1)
	print (int(Nspatialgrid[0])+2*Nborder)*(int(Nspatialgrid[1])+2*Nborder)*(int(Nspatialgrid[2])+2*Nborder)*4/(2**20.0), "MB"
	d0=np.zeros((int(Nspatialgrid[0])+2*Nborder,int(Nspatialgrid[1])+2*Nborder,int(Nspatialgrid[2])+2*Nborder),dtype=dtyp)  		#density[x,y,z]

	print d0.nbytes/(2**20.0), "MB"
	
	
	sf=np.zeros((int(Nspatialgrid[0]),int(Nspatialgrid[1]),int(Nspatialgrid[2])/2+1),dtype=dtyp)  								#3d Structure factor

	#spatial grid.  grid[x,y,z,0] gives x coordinate of this grid position.  4th index value of 3 will ultimately store SF data
	grid=np.zeros((int(Nspatialgrid[0])+2*Nborder,int(Nspatialgrid[1])+2*Nborder,int(Nspatialgrid[2])+2*Nborder,4),dtype=dtyp)		#if this causes a MemoryError for your system, consider changing dtyp from np.float64 to np.float32 or smaller
	# except:
		# print "Memory Allocation Error.  Please reduce dtyp in file dens.py and try again"
		# print "current dtyp=",dtyp		
	
	for ix in xrange(grid.shape[0]):															
		grid[ix,:,:,0]+=(ix-Nborder)*dr[0]
	for iy in xrange(grid.shape[1]):
		grid[:,iy,:,1]+=(iy-Nborder)*dr[1]	
	for iz in xrange(grid.shape[2]):
		grid[:,:,iz,2]+=(iz-Nborder)*dr[2]

		# pl("dr",dr)
		# pl("ix",ix)
		# pl("iy",iy)
		# pl("iz",iz)
	# for i in xrange(4):
		# pl(i,np.amax(grid[:,:,:,i]))
	# exit()
	#declare variables outside of loop so they're (possibly) handled more efficiently
	buf=np.array([])
	sig=0
	A1=0
	ix0=0.0
	ix1=0.0
	iy0=0.0
	iy1=0.0
	iz0=0.0
	iz1=0.0
	Ax=1
	Ay=1
	Az=1
		
	#lists for remapping grid
	# des1=[0,Nborder,Nspatialgrid-Nborder,Nspatialgrid]
	# ori1=[0,Nborder,Nborder+Nspatialgrid,2*Nborder+Nspatialgrid]
	
	
	des_tcl=[]
	ori_tcl=[]
	
	for i in xrange(3):
		des_tcl.append([0,Nborder,Nspatialgrid[i]-Nborder,Nspatialgrid[i]])
		ori_tcl.append([0,Nborder,Nborder+Nspatialgrid[i],2*Nborder+Nspatialgrid[i]])
	
	
	print "Calculating Structure factor for ",r.shape[1]," atoms over ",r.shape[0], " timesteps. \n" ,"Progress: "
	
	Ntimesteps=r.shape[0]
	if RANDOM_NOISE>0:
		Ntimesteps*=1
		print "*"*80
		print "*"*80
		print "OVERRIDING DENSITY WITH RANDOM NOISE"
		print Ntimesteps, "timesteps"
		print "*"*80
		print "*"*80
		
	
	for it in tqdm.tqdm(xrange(Ntimesteps),unit='timesteps'):  #loop over timestep
	
			if RANDOM_NOISE>0:

				d1=np.random.rand(int(Nspatialgrid[0]),int(Nspatialgrid[1]),int(Nspatialgrid[2]))  		#density[x,y,z]
			else:
		
				for im in tqdm.tqdm(xrange(r.shape[1]),unit='atoms'): #loop over atoms

					ir=(r[it,im,:]/dr).astype(int)  #get integer coordinates of atoms
					
		
					Ax=np.int(bdict[typ[im]][0])  
					Ay=np.int(bdict[typ[im]][1])
					Az=np.int(bdict[typ[im]][2])
					
					ix0=ir[0]-Ax+Nborder
					ix1=ir[0]+Ax+Nborder
					iy0=ir[1]-Ay+Nborder
					iy1=ir[1]+Ay+Nborder
					iz0=ir[2]-Az+Nborder
					iz1=ir[2]+Az+Nborder
					
					buf=r[it,im,:]-grid[ix0:ix1,iy0:iy1,iz0:iz1,:3]
					
					#print "buf ",buf.shape
									
					
					#distSQ=buf[...,0]**2+buf[...,1]**2+buf[...,2]**2 
					#distSQ=np.einsum('ijkl,ml->ijkm',buf,ucell)
					
					buf2=np.einsum('ml,ijkm',ucell,buf)
					# buf2=buf
					
					distSQ=buf2[...,0]**2+buf2[...,1]**2+buf2[...,2]**2 
					
					
					sig=rad[typ[im]][1]	#width of gaussian
					Nel=rad[typ[im]][0]	#electron number
					
					d0[ix0:ix1,iy0:iy1,iz0:iz1]+=(Nel/np.power(sig,3.))*np.exp(-distSQ / (2 * np.power(sig, 2.)))  
					#print dr[0]*dr[1]*dr[2]*np.sum(d0)/np.power(2.0*math.pi,1.5)  #use this to verify normalization
				
				d1=remap_grid_tcl(d0,des_tcl,ori_tcl)			#remap density onto periodic cell
			
			
			# plt.contourf(d1[:,:,d1.shape[2]/2])
			# plt.show()
			# plt.clf()
			

			dfft=np.fft.rfftn(d1)				#compute 3d FFT
			
			dmag=dfft.conjugate()						
			dmag=np.real(np.multiply(dmag,dfft,out=dmag))	#compute magnitude of complex structure factor
			
			sf+=dmag	#add to time-averaged structure factor
			
			d1*=0.0   #reset arrays
			d0*=0.0
			

			# optionally plot each SF at each timestep (note that the corners correspond to the smallest k values)
			# plt.contourf(sf[:,:,sf.shape[2]/2])
			# plt.show()
			# plt.clf()
						
	
	
	sfplt=get_dplot(sf)
	
	# plt.plot(sfplt[sfplt.shape[0]/2,sfplt.shape[1]/2,:])
	# plt.show()
	# exit()
			
	
	kgrid=np.zeros((sf.shape[0],sf.shape[1],sf.shape[2],4))
	
	for ix in xrange(kgrid.shape[0]/2):
		kgrid[ix,:,:,0]=ix*2.0*math.pi/L[0]
		kgrid[kgrid.shape[0]-1-ix,:,:,0]=-(ix+0.5)*2.0*math.pi/L[0]
	for iy in xrange(kgrid.shape[1]/2):
		kgrid[:,iy,:,1]=iy*2.0*math.pi/L[1]
		kgrid[:,kgrid.shape[1]-1-iy,:,1]=-(iy+0.5)*2.0*math.pi/L[1]
	
	for iz in xrange(kgrid.shape[2]):
		kgrid[:,:,iz,2]=iz*2.0*math.pi/L[2]

	kgridplt=np.zeros((sf.shape[0]-2,sf.shape[1]-2,sf.shape[2]*2-3,4))
	for ix in xrange(kgridplt.shape[0]):
		kgridplt[ix,:,:,0]=(ix-kgridplt.shape[0]/2)*2.0*math.pi/L[0]
	for ix in xrange(kgridplt.shape[1]):	
		kgridplt[:,ix,:,1]=(ix-kgridplt.shape[1]/2)*2.0*math.pi/L[1]
	for ix in xrange(kgridplt.shape[2]):
		kgridplt[:,:,ix,2]=(ix-kgridplt.shape[2]/2)*2.0*math.pi/L[2]
	
	kgridplt[:,:,:,3]=sfplt
	
	# print "z: \n"
	# print kgridplt[0,0,0,2]
	# print kgridplt[0,0,kgridplt.shape[2]/2,2]
	# print kgridplt[0,0,kgridplt.shape[2]-1,2]
	
	# print 
	# plt.plot(kgridplt[kgridplt.shape[0]/2,kgridplt.shape[1]/2,:,3])
	# plt.show()
	
	
	
	np.savez_compressed(out_filename,sf=sf,sfplt=sfplt,L=L,N=Nspatialgrid,kgrid=kgrid,kgridplt=kgridplt)
	