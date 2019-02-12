#!/usr/bin/env python

from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from past.utils import old_div
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import math
from tqdm import trange
import argparse
from PIL import Image
import os


def initialize():

    parser = argparse.ArgumentParser(description='Measure Membrane Thickness')

    parser.add_argument('-l', '--load', help='Compressed numpy arrays output from main_gromacs.py')
    parser.add_argument('-b', '--box', nargs='+', help='Crop box dimensions', default=[63, 73, 60, 40])
    parser.add_argument('-o', '--output', type=str, default='intensities.txt', help='Name of output file to write intensity to')
    parser.add_argument('-y', '--yes', action="store_true", help='Use if you know the box dimensions are correct')
    parser.add_argument('-f', '--full', action="store_true", help='See full XRD pattern before cropping')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    wavelength_angstroms = 1.54
    theta = math.pi/3.0
    ucell = np.array([[1, 0, 0], [np.cos(theta), np.sin(theta), 0], [0, 0, 1]])

    dpl = np.load(args.load)

    D = dpl['kgridplt']

    X = D[:, 0, 0, 0].copy()
    Y = D[0, :, 0, 1].copy()
    Z = D[0, 0, :, 2].copy()

    SF = D[:, :, :, 3]

    a1 = ucell[0]
    a2 = ucell[1]
    a3 = ucell[2]

    b1 = old_div((np.cross(a2, a3)), (np.dot(a1, np.cross(a2, a3))))
    b2 = old_div((np.cross(a3, a1)), (np.dot(a2, np.cross(a3, a1))))
    b3 = old_div((np.cross(a1, a2)), (np.dot(a3, np.cross(a1, a2))))

    Dnew = np.zeros_like(D)

    for ix in trange(D.shape[0]):
        Dnew[ix, :, :, 0:3] += X[ix]*b1
    for iy in trange(D.shape[1]):
        Dnew[:, iy, :, 0:3] += Y[iy]*b2
    for iz in trange(D.shape[2]):
        Dnew[:, :, iz, 0:3] += Z[iz]*b3

    D[..., :3] = Dnew[..., :3]

    K_ES = 2.0*math.pi/wavelength_angstroms  #calculate k for incident xrays in inverse angstroms

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html
    # scipy.interpolate.RegularGridInterpolator
    # Notes
    # Contrary to LinearNDInterpolator and NearestNDInterpolator, RegularGridInterpolator class avoids expensive
    # triangulation of the input data by taking advantage of the regular grid structure.
    # this is why this style of interpolation is so slow

    XGD = D[:, :, :, 0]  #X spatial grid view
    YGD = D[:, :, :, 1]
    ZGD = D[:, :, :, 2]
    VGD = D[:, :, :, 3]

    DC = D[:, :, :, 0:3]


    # check if fast interpolation can be used
    lbuf = True
    for i in range(3):
        for j in range(i + 1, 3):
            if ucell[i, j] != 0 or ucell[j, i] != 0:
                lbuf = False


    ES = RegularGridInterpolator((X, Y, Z), SF, bounds_error=False)

    xyzpts = []

    print("setting up points for radial integration")
    if True:
        for ix in trange(D.shape[0]):
            for iy in range(D.shape[1]):
                for iz in range(D.shape[2]):
                    # xyzpts.append((X[ix],Y[iy],Z[iz]))
                    xyzpts.append((D[ix, iy, iz, 0], D[ix, iy, iz, 1], D[ix, iy, iz, 2]))
    else:
        pass

    xyzpts = np.asarray(xyzpts)

    EWDxyz = ES(xyzpts)

    rpts = np.sqrt(xyzpts[:, 0]**2.0 + xyzpts[:, 1]**2.0)

    Hcount, XEC, YEC = np.histogram2d(rpts, xyzpts[:, 2], bins=(X, Z))

    Hval, XEV, YEV = np.histogram2d(rpts, xyzpts[:, 2], weights=EWDxyz, normed=True, bins=(X, Z))

    Hcount = np.where(Hcount == 0, 1, Hcount)

    Hrz = Hval / Hcount

    # if not switch1:
    #     Hrz = np.ma.masked_invalid(Hrz)

    for ir in range(old_div(Hrz.shape[0], 2)):
        Hrz[-ir + old_div(Hrz.shape[0], 2), :] = Hrz[ir + 1 + old_div(Hrz.shape[0], 2), :]

    XMG, YMG = np.meshgrid(XEV, YEV)

    if args.full:
        lx = 0
        by = 0
        rx = Hrz.shape[0]
        ty = Hrz.shape[1]
    else:
        lx = int(args.box[0])
        rx = int(args.box[1])
        ty = int(args.box[2])
        by = int(args.box[3])

    print("Please crop the picture to the region where you'd like to measure average and maximum intensity")

    fig = plt.figure()
    plt.pcolormesh(XMG[by:ty, lx:rx], YMG[by:ty, lx:rx], Hrz.T[by:ty, lx:rx], vmin=0.0, vmax=0.007*np.amax(Hrz))
    if not args.yes:
        plt.show(block=False)
        response = input("Is this good? ")
        while response != 'yes':
            new_box_info = input("Please enter new box dimensions: ")
            plt.close(fig)
            new = new_box_info.split()

            if 'left' in new:
                lx += int(new[new.index('left') + 1])
            if 'top' in new:
                ty += int(new[new.index('top') + 1])
            if 'right' in new:
                rx += int(new[new.index('right') + 1])
            if 'bottom' in new:
                by += int(new[new.index('bottom') + 1])
            if 'zoom' in new:
                factor = float(new[new.index('zoom') + 1])
                print("Zooming in by a factor of %s" % factor)
                lx += int((rx - lx) / (factor * 2))
                rx -= int((rx - lx) / (factor * 2))
                ty += int((by - ty) / (factor * 2))
                by -= int((by - ty) / (factor * 2))
            if 'max' in new:
                print("Maximum intensity in this region: %s" % np.amax(Hrz.T[by:ty, lx:rx]))
            if 'up' in new:
                ty += int(new[new.index('up') + 1])
                by += int(new[new.index('up') + 1])
            if 'full' in new:
                lx = 0
                by = 0
                rx = Hrz.shape[0]
                ty = Hrz.shape[1]

            fig = plt.figure()
            plt.pcolormesh(XMG[by:ty, lx:rx], YMG[by:ty, lx:rx], Hrz.T[by:ty, lx:rx], vmin=0.0, vmax=0.007*np.amax(Hrz))
            plt.show(block=False)
            response = input("How about now? ")

        print(lx, rx, ty, by)

    Iavg = np.mean(Hrz.T[by:ty, lx:rx])
    Imax = np.amax(Hrz.T[by:ty, lx:rx])

    print("Average Intensity in region = %s" % Iavg)
    print("Maximum Intensity in region = %s" % Imax)

    iter = int(str.split(str.split(args.load, '_')[3], '.')[0])  # very specific naming convention in main_gromacs.py

    if os.path.isfile(args.output):
        with open(args.output, 'a') as f:
            f.write("{:<5d}{:6.5f}\n".format(iter, Imax))
    else:
        with open(args.output, 'w') as f:
            f.write("{:<5d}{:6.5f}\n".format(iter, Imax))