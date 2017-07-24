#! /usr/bin/env python

from __future__ import division
from __future__ import print_function
from builtins import range
# from past.utils import old_div
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import math
import os.path
import dens2
import plot2d as p2d
import tqdm  # progress bar
import platform
import matplotlib.pyplot as plt  # must be imported after anything that imports mayavi/mlab
import argparse
import warnings

parser = argparse.ArgumentParser(description='Calculate 3d Structure Factor')

parser.add_argument('-i', '--input', default='', type=str, help='Input topolgy and trajectory basename')
parser.add_argument('-top', '--topology', default='', type=str, help='Input topolgy filename')
parser.add_argument('-traj', '--trajectory', default='', type=str, help='Input trajectory filename')
parser.add_argument('--cscale', default=1, type=float, help='Scale color map on plots')
parser.add_argument('--lcscale', default=1, type=float, help='Scale color map on log plots')
parser.add_argument('-fi', '--first_frame', default=0, type=int, help='frame to start at')
parser.add_argument('-e', '--end_frame', default=-1, type=int, help='frame to end at')
parser.add_argument('-fr', '--force_recompute', default=0, type=int, help='force recomputing SF (if >=1) or trajectory and SF(if >=2)')
#parser.add_argument('-o', '--output', default='', type=str, help='override output basename')

args=parser.parse_args()

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in

warnings.simplefilter("ignore", RuntimeWarning)

if len(args.input) > 0:
    basename = args.input
    top_file = args.input+".gro"
    traj_file = args.input+".trr"
else:
    top_file = args.topology
    traj_file = args.trajectory
    basename = args.topology.rsplit('.', 1)[0]

# if len(args.output)>0:
# 	basename=args.output

print("running on", platform.system(), platform.release(), platform.version())

if platform.system() == "Windows":  # path separators
    fd = "\\"
else:
    import load_traj as lt
    fd = "/"

label = "out_"+basename

tfname = label+"_traj"
sfname = label+"_sf"

dens2.Nspatialgrid = 128

if args.force_recompute > 0 or not os.path.isfile(sfname+".npz"):						# check to see if SF needs to be calculated
    if args.force_recompute > 1 or not os.path.isfile(tfname+".npz"):  				# check to see if trajectory needs to be processed
        print("processing trajectory file "+traj_file)
        lt.process_gro_mdtraj(top_file, traj_file, tfname)  # Process trajectory into numpy array.  Commented to run on windows
        print('done')
    traj = np.load(tfname+".npz")							# load processed trajectory
    rad = dens2.load_radii("%s/radii.txt" % location)					# load radii definitions from file stored in git repo

    dens2.compute_sf(traj['coords'][args.first_frame:args.end_frame,...],traj['dims'],traj['typ'],sfname,rad)	 # compute time-averaged 3d structure factor and save to sfname.npz

dpl = np.load(sfname+".npz")  # load 3d SF
grid = dpl['kgridplt']

p2d.mainlabel = basename			# title for plots

dir = p2d.mainlabel+"_plots"+fd	 # set up directories for saving plots
sfdir = dir+"structure_factor"+fd
sfsubdir = sfdir+"additional_plots"+fd
EWdir = dir+"Ewald_Corrected"+fd

if True:
    print("Ewald plots")
    p2d.path = EWdir
    p2d.Plot_Ewald_Sphere_Correction(grid, 1.54, args.cscale, args.lcscale)		# compute Ewald-corrected SF cross sections in xy,xz,yz planes
    theta = math.pi / 3.0
    ucell = np.array([[1, 0, 0], [np.cos(theta), np.sin(theta), 0], [0, 0, 1]])
    p2d.Plot_Ewald_triclinic(grid, 1.54, ucell)
    p2d.radial_integrate(grid, 2000, dir + "radial.png")

exit()
# xy,yz,xz planes of SF
if True:
    # print "xy,yz,xz plots"
    # p2d.path = dir + sfdir
    # p2d.sfplot(grid[grid.shape[0]/2, :, :, :], args.lcscale)		# plot yz plane
    # p2d.sfplot(grid[:, grid.shape[1]/2, :, :], args.lcscale)		# plot xz plane
    # p2d.sfplot(grid[:, :, grid.shape[2]/2, :], args.lcscale)		# plot xy plane
    p2d.radial_integrate(grid, 1500, dir + "radial.png")

if True:  # additional slices through SF
    print("additional plots")
    Nsteps = 8
    p2d.path = sfsubdir+"xplots"+fd
    p2d.savelin = False
    for i in tqdm.tqdm(range(0, grid.shape[0], old_div(grid.shape[0],Nsteps))):
        p2d.sfplot(grid[i, :, :, :], args.lcscale)

    p2d.path = sfsubdir+"yplots"+fd
    for i in tqdm.tqdm(range(0, grid.shape[1], old_div(grid.shape[1],Nsteps))):
        p2d.sfplot(grid[:, i, :, :], args.lcscale)

    p2d.path = sfsubdir+"zplots"+fd
    for i in tqdm.tqdm(range(0, grid.shape[1], old_div(grid.shape[2],Nsteps))):
        p2d.sfplot(grid[:, :, i, :], args.lcscale)
