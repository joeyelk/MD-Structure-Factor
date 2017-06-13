
import numpy as np

from scipy.interpolate import RegularGridInterpolator
import math
import os.path


import dens                
import plot2d as p2d  

import tqdm  #progress bar

import platform
import matplotlib.pyplot as plt  #must be imported after anything that imports mayavi/mlab

import argparse


parser = argparse.ArgumentParser(description='Calculate 3d Structure Factor')

parser.add_argument('-i', '--input', default='', type=str, help='Input topolgy and trajectory basename')
parser.add_argument('-top', '--topology', default='', type=str, help='Input topolgy filename')
parser.add_argument('-traj', '--trajectory', default='', type=str, help='Input trajectory filename')
parser.add_argument('-fi', '--first_frame', default=0, type=int, help='frame to start at')
parser.add_argument('-fr', '--force_recompute', default=0, type=int, help='force recomputing SF (if >=1) or trajectory and SF(if >=2)')
parser.add_argument('-st', '--save_traj', default=0, type=int, help='set equal to 1 to save trajectory as npz file')
#parser.add_argument('-o', '--output', default='', type=str, help='override output basename')


args=parser.parse_args()

if len(args.input)>0:
	basename=args.input
	top_file=args.input+".gro"
	traj_file=args.input+".trr"
else:
	top_file=args.topology
	traj_file=args.trajectory
	basename=args.topology.rsplit('.', 1)[0]

#if len(args.output)>0:
#	basename=args.output
	


print "running on",platform.system(),platform.release(),platform.version()

if platform.system()=="Windows":  #path separators 
	fd="\\"
else:
	import load_traj as lt
	fd="/"

	
label="out_"+basename

tfname=label+"_traj"
sfname=label+"_sf"

dens.Nspatialgrid=128


if args.force_recompute>0 or not os.path.isfile(sfname+".npz"):					#check to see if SF needs to be calculated
	if args.force_recompute>1 or not os.path.isfile(tfname+".npz"):  				#check to see if trajectory needs to be processed
		if platform.system()=="Windows":
			print "Unable to process trajectory file on Windows"
			exit()
		if args.save_traj==1:
			lt.process_gro(top_file,traj_file,tfname)   					#Process trajectory into numpy array.  
			coords=traj['coords']
			dims=traj['dims']
			typ=traj['typ']
			mass=traj['mass']
			name=traj['name']			
		else:
			dims,coords,name,mass,typ=lt.process_gro(top_file,traj_file)   					#load trajectory
	else:
		traj=np.load(tfname+".npz")							#load processed trajectory
		coords=traj['coords']
		dims=traj['dims']
		typ=traj['typ']
	
	rad=dens.load_radii("radii.txt")					#load radii definitions from file
	dens.compute_sf(coords[args.first_frame:,...],dims[args.first_frame:,...],typ,sfname,rad)		#compute time-averaged 3d structure factor and save to sfname.npz


dpl=np.load(sfname+".npz")					#load 3d SF

grid=dpl['kgridplt']


p2d.mainlabel=basename			#title for plots

dir=p2d.mainlabel+"_plots"+fd	#set up directories for saving plots
sfdir=dir+"structure_factor"+fd
sfsubdir=sfdir+"additional_plots"+fd
EWdir=dir+"Ewald_Corrected"+fd


if True:
	print "Ewald plots"
	p2d.path=EWdir
	p2d.Plot_Ewald_Sphere_Correction(grid,1.54)		#compute Ewald-corrected SF cross sections in xy,xz,yz planes

#xy,yz,xz planes of SF
if True:
	print "xy,yz,xz plots"
	p2d.path=dir+sfdir
	p2d.sfplot(grid[grid.shape[0]/2,:,:,:])		#plot yz plane
	p2d.sfplot(grid[:,grid.shape[1]/2,:,:])		#plot xz plane
	p2d.sfplot(grid[:,:,grid.shape[2]/2,:])		#plot xy plane
	p2d.radial_integrate(grid,300,dir+"radial.png")


if True:  #additional slices through SF
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
