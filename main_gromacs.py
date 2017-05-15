
import numpy as np

from scipy.interpolate import RegularGridInterpolator
import math
import os.path

#import load_traj as lt
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

args=parser.parse_args()

if len(args.input)>0:
	basename=args.input
	top_file=args.input+".gro"
	traj_file=args.input+".trr"
else:
	top_file=args.topology
	traj_file=args.trajectory
	basename=args.topology.rsplit('.', 1)[0]
	
#print basename
#exit()

print "running on",platform.system(),platform.release(),platform.version()

if platform.system()=="Windows":  #path separators 
	fd="\\"
else:
	fd="/"



	
label="out_"+basename

tfname=label+"_traj"
sfname=label+"_sf"

dens.Nspatialgrid=128


if not os.path.isfile(sfname+".npz"):					#check to see if SF needs to be calculated
	if not os.path.isfile(tfname+".npz"):  				#check to see if trajectory needs to be processed
		print "processing trajectory file "+traj_file
		#lt.process_gro(top_file,traj_file,tfname)   					#Process trajectory into numpy array.  Commented to run on windows
		print 'done'
	traj=np.load(tfname+".npz")							#load processed trajectory
	rad=dens.load_radii("radii.txt")					#load radii definitions from file

	dens.compute_sf(traj['coords'],traj['dims'],traj['typ'],sfname,rad)		#compute time-averaged 3d structure factor and save to sfname.npz


dpl=np.load(sfname+".npz")					#load 3d SF
grid=dpl['kgridplt']

p2d.mainlabel=basename+"o"			#title for plots

dir=p2d.mainlabel+"_plots"+fd	#set up directories for saving plots
sfdir=dir+"structure_factor"+fd
sfsubdir=sfdir+"additional_plots"+fd
EWdir=dir+"Ewald_Corrected"+fd

p2d.path=EWdir
p2d.Plot_Ewald_Sphere_Correction(grid,1.54)		#compute Ewald-corrected SF cross sections in xy,xz,yz planes

#xy,yz,xz planes of SF
if True:
	p2d.path=dir+sfdir
	p2d.sfplot(grid[grid.shape[0]/2,:,:,:])		#plot yz plane
	p2d.sfplot(grid[:,grid.shape[1]/2,:,:])		#plot xz plane
	p2d.sfplot(grid[:,:,grid.shape[2]/2,:])		#plot xy plane


if True:  #additional slices through SF
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
