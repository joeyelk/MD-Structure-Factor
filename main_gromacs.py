
import numpy as np

from scipy.interpolate import RegularGridInterpolator
import math
import os.path

import dens                
import plot2d as p2d  

import tqdm  #progress bar
from tqdm import trange

import platform
import matplotlib.pyplot as plt  #must be imported after anything that imports mayavi/mlab

import argparse

XRAY_WAVELENGTH=1.54
TRANSFORM_MONOCLINIC=True

parser = argparse.ArgumentParser(description='Calculate 3d Structure Factor')
#File names
parser.add_argument('-i', '--input', default='', type=str, help='Input topolgy and trajectory basename')
parser.add_argument('-top', '--topology', default='', type=str, help='Input topolgy filename')
parser.add_argument('-traj', '--trajectory', default='', type=str, help='Input trajectory filename')
#parser.add_argument('-o', '--output', default='', type=str, help='override output basename')  #there are issues with using this, so it is disabled for now

#other simulation parameters
parser.add_argument('-fi', '--first_frame', default=0, type=int, help='frame to start at')
parser.add_argument('-TRstride', default=1, type=int, help='stride to use for computing SF')

parser.add_argument("-tf","--traj_format",default="gromacs",type=str,help='format of MD trajectory file to load')


parser.add_argument('-fr', '--force_recompute', default=0, type=int, help='force recomputing SF (if >=1) or trajectory and SF(if >=2)')

#random trajectory parameters
parser.add_argument('-RC', '--random_counts', default=0, type=int, help='set this to specify number of random particles to use')
parser.add_argument('-RT', '--random_timesteps', default=1, type=int, help='set this to specify number of random timesteps to use')
parser.add_argument('-RL', '--random_label', default="R3", type=str, help='set this to specify the element type of the random particle')

#parser.add_argument('-HL', '--hex_label', default="R3", type=string, help='set this to specify the element type of the hexagonal lattice')


parser.add_argument('-LX', '--lattice_x', default=0, type=int, help='set this to specify the number of lattice points in the X direction')
parser.add_argument('-LY', '--lattice_y', default=1, type=int, help='set this to specify the number of lattice points in the Y direction')
parser.add_argument('-LZ', '--lattice_z', default=1, type=int, help='set this to specify the number of lattice points in the Z direction')

# parser.add_argument('-LN', '--lattice_x', default=0, type=int, help='set this to override the number of lattice points in each direction')


parser.add_argument('-LL', '--lattice_label', default="R3", type=str, help='set this to specify the element label of lattice particles')

parser.add_argument('-SR', '--spatial_resolution', default=1.0, type=float,help='set this to specify the spatial resolution for the density grid')
parser.add_argument('-RN', '--random_noise', default=0, type=int,help='set this to a positive value to use random noise for testing.  A conventional trajectory must still be loaded.  The number of timesteps to be used will be scaled by this integer')
parser.add_argument('-RS', '--random_seed', default=1, type=int,help='Set the random seed from the command line')
 
parser.add_argument('-NBR', '--number_bins_rad', default=0, type=int,help='Set this to a nonzero value to use that many radial bins.  These bins will be scaled such that they contain roughly the same number of points') 

# tdict={"orth":2.0,"hex":3.0}
# helpstring=''.join([i+" pi/"+str(j) for i,j in tdict.iteritems()])
# parser.add_argument('--cs','--cell_style' default="orth", type=str,help="choose cell type from "+helpstring) 
parser.add_argument('-ct','--cell_theta', default=90.0, type=float,help="choose cell theta (in degrees)") 

# theta=math.pi/3.0  #theta for monoclinic unit cell
# theta=math.pi/6.0  #theta for monoclinic unit cell
# theta=math.pi/2.0

# theta=math.pi/tdict[args.cell_style]
args=parser.parse_args()


theta=args.cell_theta*math.pi/180.0
ucell=np.array([[1,0,0],[np.cos(theta),np.sin(theta),0],[0,0,1]])



np.random.seed=args.random_seed #args.random_noise
p2d.NBINSRAD=args.number_bins_rad

if args.random_noise>0:
	
	dens.RANDOM_NOISE=args.random_noise

top_extension={"gromacs":".gro","namd":".psf"}
traj_extension={"gromacs":".trr","namd":".dcd"}

if len(args.input)>0:
	basename=args.input
	top_file=args.input+top_extension[args.traj_format.lower()]
	traj_file=args.input+traj_extension[args.traj_format.lower()]
else:
	top_file=args.topology
	traj_file=args.trajectory
	basename=args.topology.rsplit('.', 1)[0]


#if len(args.output)>0:
#	basename=args.output


print "running on",platform.system(),platform.release(),platform.version()

if platform.system()=="Windows":  #path separators are different between unix/windows
	fd="\\"
else:
	import load_traj as lt
	fd="/"


label="out_"+basename

tfname=label+"_traj"
sfname=label+"_sf"


#dens.Nspatialgrid=128 #number of points for spatial grid.  Since this is 3d, memory usage scales as N^3.  Although this can be reduced in the future




#dens.Nspatialgrid=np.asarray([128,128,200]) #number of points for spatial grid.  Since this is 3d, memory usage scales as N^3.  Although this can be reduced in the future
#dens.Nspatialgrid=64
#dens.Nspatialgrid=180

# if args.LN>0:
	# args.lattice_x=args.LN
	# args.lattice_y=args.LN
	# args.lattice_z=args.LN
	
Nlat=args.lattice_x*args.lattice_y*args.lattice_z

if Nlat>0:

	boxsize=100.0
	
	Nx=args.lattice_x
	Ny=args.lattice_y
	Nz=args.lattice_z
	
	Nsteps=1
	
	dims=np.ones((Nsteps,3))*boxsize
	coords=np.zeros((Nsteps,Nlat,3))
	
	cbufx=np.linspace(0,boxsize,Nx,endpoint=False)
	cbufy=np.linspace(0,boxsize,Ny,endpoint=False)
	cbufz=np.linspace(0,boxsize,Nz,endpoint=False)
	
	cbx,cby,cbz=np.meshgrid(cbufx,cbufy,cbufz)
	
	coords[0,:,0]=cbx.reshape((Nlat))
	coords[0,:,1]=cby.reshape((Nlat))
	coords[0,:,2]=cbz.reshape((Nlat))
	
	sfname='lattice_'+str(Nx)+"_"+str(Ny)+"_"+str(Nz)
	
	name=np.zeros(Nlat,dtype=object)
	mass=np.zeros(Nlat)
	typ=np.zeros(Nlat,dtype=object)

	typ[:]=args.lattice_label
	
	
	print "saving..."
	np.savez_compressed(sfname,dims=dims,coords=coords,name=name,typ=typ)
	rad=dens.load_radii("radii.txt")					#load radii definitions from file
	print "computing SF..."
	
	dens.compute_sf(coords,dims,typ,sfname,rad,ucell,args.spatial_resolution)		#compute time-averaged 3d structure factor and save to sfname.npz
	


elif args.random_counts>0:	#create a random trajectory
	spcheck=0  #check if simulating a single particle.  
	
	Rboxsize=100.0
	BUFFsize=10000000
	sfname='RND'
	
	print "generating random trajectory..."
	Rsteps=args.random_timesteps
	Ratoms=args.random_counts #100000#0
	
	dims=np.ones((Rsteps,3))*Rboxsize
	# dims[0,2]*=1.2
	coords=np.random.random((Rsteps,Ratoms,3))*dims[0,:]
	
	name=np.zeros(Ratoms,dtype=object)
	mass=np.zeros(Ratoms)
	typ=np.zeros(Ratoms,dtype=object)

	# for it in xrange(Ratoms):
		# typ[it]=args.random_label
	typ[:]=args.random_label
	
	print "saving..."
	np.savez_compressed("RAND",dims=dims,coords=coords,name=name,typ=typ)
	rad=dens.load_radii("radii.txt")					#load radii definitions from file
	print "computing SF..."
	
	dens.compute_sf(coords,dims,typ,sfname,rad,ucell,args.spatial_resolution)		#compute time-averaged 3d structure factor and save to sfname.npz

else:  #load trajectory or npz file

	if args.force_recompute>0 or not os.path.isfile(sfname+".npz"):			#check to see if SF needs to be calculated
		if args.force_recompute>1 or not os.path.isfile(tfname+".npz"): 		#check to see if trajectory needs to be processed
			print sfname
			print tfname
			if platform.system()=="Windows":  										#This part must be done in an environment that can import MDAnalysis
				print "Unable to process trajectory file on Windows"
				print "This part must be done in an environment that can import MDAnalysis"
				exit()
			else:
				print "processing trajectory file "+traj_file
				
				lt.process(top_file,traj_file,tfname,args.TRstride)   					#Process trajectory into numpy array.  
				# if args.traj_format=="gromacs":
					# lt.process_gro(top_file,traj_file,tfname,args.TRstride)   					#Process trajectory into numpy array.  
				# elif args.traj_format.lower()=="namd":
					# lt.process(top_file,traj_file,tfname,args.TRstride)
				print 'done'
				
		traj=np.load(tfname+".npz")							#load processed trajectory
		print "trajectory file",tfname+".npz", "loaded"
		
		print type(traj["coords"])
		T=traj['coords']		
		
		if TRANSFORM_MONOCLINIC and theta!=90.0:	
			print "transforming coordinates to monoclinic cell (theta={0:f} deg)".format(theta*180.0/np.pi)
			T[...,1]=T[...,1]/np.sin(theta)
			T[...,0]=T[...,0]-T[...,1]*np.cos(theta)
		
		# coords=traj["coords"]
		
		# print .shape
		# for i in traj:
			# print i #,val.shape
		# print traj.shape
		# exit()
		
		rad=dens.load_radii("radii.txt")					#load radii definitions from file

		
		dens.compute_sf(T[args.first_frame:,...],traj['dims'][args.first_frame:,...],traj['typ'],sfname,rad,ucell,args.spatial_resolution)		#compute time-averaged 3d structure factor and save to sfname.npz


print "reloading SF..."
dpl=np.load(sfname+".npz")					#load 3d SF

grid=dpl['kgridplt']  #grid contains Kx,Ky,Kz,S(Kx,ky,kz) information

#print grid[:,0,0,0]
#print grid.shape

#print grid[grid.shape[0]/2,grid.shape[1]/2,grid.shape[2]/2,3]
#print grid[grid.shape[0]/2+1,grid.shape[1]/2,grid.shape[2]/2,3]
# 



p2d.mainlabel=basename			#title for plots

dir=p2d.mainlabel+"_plots"+fd	#set up directories for saving plots
sfdir=dir+"structure_factor"+fd
sfsubdir=sfdir+"additional_plots"+fd
EWdir=dir+"Ewald_Corrected"+fd

print "making plots..."
if True:
	print "Ewald plots"
	p2d.path=EWdir
	#p2d.Plot_Ewald_Sphere_Correction(grid,1.54)		#compute Ewald-corrected SF cross sections in xy,xz,yz planes
	
	#print ucell.shape
	#print ucell
	# 
	p2d.Plot_Ewald_triclinic(grid,XRAY_WAVELENGTH,ucell)		#compute Ewald-corrected SF cross sections in xy,xz,yz planes
	print "EW done"
	#xy,yz,xz planes of SF
if True:
	print "xy,yz,xz plots"
	p2d.path=dir+sfdir
	# p2d.sfplot(grid[grid.shape[0]/2,:,:,:])		#plot yz plane
	# p2d.sfplot(grid[:,grid.shape[1]/2,:,:])		#plot xz plane
	# p2d.sfplot(grid[:,:,grid.shape[2]/2,:])		#plot xy plane
	p2d.radial_integrate(grid,300,dir+"radial.png")


if False:  #additional slices through SF
	print "additional plots"
	Nsteps=8
	p2d.path=sfsubdir+"xplots"+fd
	p2d.savelin=False
	for i in tqdm.tqdm(xrange(0,grid.shape[0],grid.shape[0]/Nsteps)):
		p2d.sfplot(grid[i,:,:,:])

	p2d.path=sfsubdir+"yplots"+fd		
	for i in tqdm.tqdm(xrange(0,grid.shape[1],grid.shape[1]/Nsteps)):
		p2d.sfplot(grid[:,i,:,:])

	p2d.path=sfsubdir+"zplots"+fd	
	for i in tqdm.tqdm(xrange(0,grid.shape[1],grid.shape[2]/Nsteps)):
		p2d.sfplot(grid[:,:,i,:])
