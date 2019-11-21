from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from past.utils import old_div
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator
import math
from tqdm import trange
import tqdm
import time
import matplotlib
from matplotlib import ticker
from matplotlib.widgets import Slider, Button, RadioButtons
from scipy.optimize import curve_fit
import datetime


colorbar = True

#matplotlib.rc('axes', color_cycle=['r', 'g', 'b', '#004060'])

mainlabel = ""

units = "$\AA^{-1}$"

xunits = units
yunits = units
zunits = units

contours = 200
DPI = 300
format = ".png"
text = "Structure Factor"

PLOT_EWALDS = True  # enable ewald-corrected SF plots
savelog = True
savelin = True
NBINSRAD = 0
normplot = 1
FP_THRESHOLD = 1.0E-12
theta = np.pi / 2

title_fontsize = 9

path = ""


def make_flat_plot(D, xr, yr, zr):

    if len(xr) != 1 and len(yr) != 1 and len(zr) != 1:
        print("error in make_flat_plot!  one of these lengths must be 1")
        exit()

    for ix in xr:
        for iy in yr:
            for iz in zr:
                r = D[ix, iy, iz, :4]
                pts.append((r))


def pl(title, obj):
    delim = "="*20
    print(delim, title, delim)
    print(obj)


def pli(obj, title=""):
    pl(title, obj)
    buf = input("enter q to quit, anything else to continue")  # raw_input renamed to input() in python3
    if buf == 'q':
        exit()


def ple(title, obj):
    pl(title, obj)
    exit()


def csplot_wlog(X, Y, Z, contours, lab, xlab, ylab, **kwargs):

    csplot(X, Y, Z, contours, lab, xlab, ylab, **kwargs)
    csplot(X, Y, np.log(Z), contours, "log_"+lab, xlab, ylab, **kwargs)


def csplot(X, Y, Z, contours, lab, xlab, ylab,**kwargs):

    title = lab+" S("+xlab+","+ylab+")"
    fname = lab+"_"+xlab+"_"+ylab
    fig, ax = plt.subplots()
    plt.suptitle(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)

    if normplot == 1:
        cax = plt.contourf(X, Y, Z / np.amax(Z), contours, vmin=0.0, vmax=0.05, **kwargs)
    else:
        cax = plt.contourf(X, Y, Z, contours, vmax=0.01*np.amax(Z), **kwargs)

    # ax.set_aspect((np.amax(Y)-np.amin(Y))/(np.amax(X)-np.amin(X)))
    # ax.set_aspect('auto')
    cbar = fig.colorbar(cax)

    plt.savefig(path+fname+format, dpi=DPI)
    plt.clf()


def sfplot(data, lcscale, **kwargs):
    """ data: plot slice through structure factor"""

    if not os.path.exists(path):
        os.makedirs(path)

    cspos = 0.0
    la = []
    lb = 0
    an = ['x', 'y', 'z']  # axes names

    for i in range(data.shape[2] - 1):
        if np.unique(data[..., i]).size > 1:
            la.append(i)
        else:
            lb = i
            cspos = data[0, 0, i]

    title = mainlabel + "\n" + text + "\n" + an[lb] + "=" + str(round(cspos, 2)) + zunits
    ltitle = mainlabel + "\n" + "log " + text + "\n" + an[lb] + "=" + str(round(cspos, 2)) + zunits

    xlab = an[la[0]]
    ylab = an[la[1]]

    filename = path + an[lb] + "=" + str(round(cspos, 2))

    xlab += "(" + xunits + ")"
    ylab += "(" + yunits + ")"

    if savelog:
        plt.suptitle(ltitle, fontsize=title_fontsize)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        max_log = np.amax(np.log(data[..., 3]))
        plt.contourf(data[..., la[0]], data[..., la[1]], np.log(data[..., 3]), contours, vmax=lcscale*max_log, **kwargs)

        plt.savefig(filename+"_log"+format, dpi=DPI)
        plt.clf()

    if savelin:
        plt.suptitle(title, fontsize=title_fontsize)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.contourf(data[..., la[0]], data[..., la[1]], data[..., 3], contours, **kwargs)
        plt.savefig(filename+format, dpi=DPI)
        plt.clf()


def radial_integrate(D, Nbins, outputname):

    SF = D[:, :, :, 3]

    R = (D[:, :, :, 0]**2).astype(np.float16) + (D[:, :, :, 1]**2).astype(np.float16) + (D[:, :, :, 2]**2).astype(np.float16)
    H, E = np.histogram(R, bins=Nbins, weights=SF)
    Hc, E = np.histogram(R, bins=Nbins)
    Hc = np.where(Hc != 0, Hc, 1.0)
    H /= Hc
    H[:1] = 0.0
    H /= np.amax(H)
    plt.plot(E[:-1], H)
    plt.xlim(0, 5)
    plt.savefig(outputname, dpi=DPI)


def spherical_integrate(D):
    exit()


def Plot_Ewald_Sphere_Correction_old(D, wavelength_angstroms):
    """ pass full 3d data,SF,wavelength in angstroms """

    X = D[:, 0, 0, 0]
    Y = D[0, :, 0, 1]
    Z = D[0, 0, :, 2]
    SF = D[:, :, :, 3]

    K_ES = 2.0*math.pi/wavelength_angstroms  # calculate k for incident xrays in inverse angstroms

    ES = RegularGridInterpolator((X, Y, Z), SF)

    pts = []
    for ix in range(D.shape[0]):
        xsq = X[ix]**2.0
        for iy in range(D.shape[1]):
            R = np.sqrt(xsq+Y[iy]**2.0)
            theta = np.arctan(old_div(R,K_ES))
            xnew = X[ix]*np.cos(theta)
            ynew = Y[iy]*np.cos(theta)
            znew = K_ES*(1.0-np.cos(theta))
            pts.append((X[ix], Y[iy], xnew, ynew, znew))

    pts = np.asarray(pts)

    EWD = ES(pts[:, 2:])
    EWD = EWD.reshape(D.shape[0], D.shape[1])
    plt.contourf(D[:, :, 0, 0], D[:, :, 0, 1], EWD, 200, interpolation=interp)

    plt.savefig("EWxy.png",dpi=300)
    plt.clf()

    plt.contourf(D[:, :, 0, 0], D[:, :, 0, 1], np.log(EWD), 200, interpolation=interp)

    plt.savefig("EWxylog.png", dpi=300)
    plt.clf()


def Plot_Ewald_Sphere_Correction(D, wavelength_angstroms, ucell=[], cscale=1, lcscale=1, **kwargs):

    """ pass full 3d data,SF,wavelength in angstroms """
    # cscale : factor by which to scale the maximum value of the colorbar
    # lcscale : factor by which to scale the maximum value of the colorbar

    if not os.path.exists(path):
        os.makedirs(path)

    X = D[:, 0, 0, 0]
    Y = D[0, :, 0, 1]
    Z = D[0, 0, :, 2]
    SF = D[:, :, :, 3]

    K_ES = 2.0*math.pi/wavelength_angstroms  # calculate k for incident xrays in inverse angstroms

    ES = RegularGridInterpolator((X, Y, Z), SF, bounds_error=False)

    xypts = []
    for ix in range(D.shape[0]):
        xsq = X[ix]**2.0
        for iy in range(D.shape[1]):
            theta = np.arctan(old_div(np.sqrt(xsq + Y[iy]**2.0),K_ES))
            xypts.append((X[ix]*np.cos(theta), Y[iy]*np.cos(theta), K_ES*(1.0 - np.cos(theta))))

    xzpts = []
    for ix in range(D.shape[0]):
        xsq = X[ix]**2.0
        for iz in range(D.shape[2]):
            theta = np.arctan(old_div(np.sqrt(xsq + Z[iz]**2.0),K_ES))
            xzpts.append((X[ix]*np.cos(theta), K_ES*(1.0-np.cos(theta)), Z[iz]*np.cos(theta)))

    yzpts = []
    for iy in range(D.shape[1]):
        ysq = Y[iy]**2.0
        for iz in range(D.shape[2]):
            theta = np.arctan(old_div(np.sqrt(ysq+Z[iz]**2.0),K_ES))
            yzpts.append((K_ES*(1.0-np.cos(theta)), Y[iy]*np.cos(theta), Z[iz]*np.cos(theta)))

    xypts = np.asarray(xypts)
    xzpts = np.asarray(xzpts)
    yzpts = np.asarray(yzpts)

    EWDxy = ES(xypts)
    EWDxz = ES(xzpts)
    EWDyz = ES(yzpts)

    EWDxy = EWDxy.reshape(D.shape[0], D.shape[1])
    EWDxz = EWDxz.reshape(D.shape[0], D.shape[2])
    EWDyz = EWDyz.reshape(D.shape[1], D.shape[2])

    title = "Ewald Corrected Structure Factor \n $\lambda=$"+str(wavelength_angstroms)+" $\AA$   $k_{ew}=$"+str(round(K_ES,2))+" $\AA^{-1}$"
    ltitle = 'log ' + title

    xlab = 'x (' + units + ")"
    ylab = 'y (' + units + ")"
    zlab = 'z (' + units + ")"

    fname = "Ewald_"

    plt.figure(1)
    plt.suptitle(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    EWDmax_xy = np.amax(EWDxy)
    plt.contourf(D[:, :, 0, 0], D[:, :, 0, 1], EWDxy, contours, vmax=cscale*EWDmax_xy, **kwargs)
    plt.savefig(path + fname + "xy" + format, dpi=DPI)
    plt.clf()

    plt.figure(2)
    plt.suptitle(ltitle)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    EWDmax_xylog = np.amax(np.log(EWDxy))
    plt.contourf(D[:, :, 0, 0], D[:, :, 0, 1], np.log(EWDxy), contours, vmax=lcscale*EWDmax_xylog, **kwargs)
    plt.savefig(path + fname + "xylog" + format, dpi=DPI)
    plt.clf()

    plt.figure(3)
    plt.suptitle(title)
    plt.xlabel(xlab)
    plt.ylabel(zlab)
    EWDmax_xz = np.amax(EWDxz)
    plt.contourf(D[:, 0, :, 0], D[:, 0, :, 2], EWDxz, contours, vmax=cscale*EWDmax_xz, **kwargs)
    plt.savefig(path + fname + "xz" + format, dpi=DPI)
    plt.clf()

    plt.figure(4)
    plt.suptitle(ltitle)
    plt.xlabel(xlab)
    plt.ylabel(zlab)
    EWDmax_xzlog = np.amax(np.log(EWDxz))
    plt.contourf(D[:, 0, :, 0], D[:, 0, :, 2], np.log(EWDxz), contours, vmax=lcscale*EWDmax_xzlog, **kwargs)
    lims = [np.amax(D[:, 0, :, 0]), np.amax(D[:, 0, :, 2])]
    qmax = min(lims)
    plt.xlim([-qmax, qmax])
    plt.ylim([-qmax, qmax])
    plt.savefig(path + fname + "xzlog" + format, dpi=DPI)
    plt.clf()

    plt.figure(5)
    plt.suptitle(title)
    plt.xlabel(ylab)
    plt.ylabel(zlab)
    EWDmax_yz = np.amax(EWDyz)
    plt.contourf(D[0, :, :, 1], D[0, :, :, 2], EWDyz, contours, vmax=cscale*EWDmax_yz, **kwargs)
    plt.savefig(path + fname + "yz" + format, dpi=DPI)
    plt.clf()

    plt.figure(6)
    plt.suptitle(ltitle)
    plt.xlabel(ylab)
    plt.ylabel(zlab)
    EWDmax_yzlog = np.amax(np.log(EWDyz))
    plt.contourf(D[0, :, :, 1], D[0, :, :, 2], np.log(EWDyz), contours, vmax=lcscale*EWDmax_yzlog, **kwargs)
    plt.savefig(path + fname + "yzlog" + format, dpi=DPI)
    plt.clf()


def lorentz(points, a, b):
    """
    :param p: lorentzian parameters : [full width half max (FWHM), position of maximum]
    :param p: position
    :return:
    """

    w = np.pi / a

    x = (b - points) / (w/2)

    return 1 / (1 + x**2)


def inverse_ft(D, ucell):

    X = D[:, 0, 0, 0]
    Y = D[0, :, 0, 1]
    Z = D[0, 0, :, 2]
    Z += Z[np.argmin(abs(Z))]

    SF = D[..., 3]

    fbin_x = X[1] - X[0]  # size of x bins in fourier space
    fbin_y = Y[1] - Y[0]  # size of y bins in fourier space
    fbin_z = Z[1] - Z[0]  # size of z bins in fourier space

    real_x = 2 * np.pi / fbin_x  # largest x dimension in real space
    real_y = 2 * np.pi / fbin_y  # largest y dimension in real space
    real_z = 2 * np.pi / fbin_z  # largest z dimension in real space

    rbin_x = real_x / X.shape[0]
    rbin_y = real_y / Y.shape[0]
    rbin_z = real_z / Z.shape[0]

    X_real = np.linspace(-real_x / 2, real_x / 2, X.shape[0])
    Y_real = np.linspace(-real_y / 2, real_y / 2, Y.shape[0])
    Z_real = np.linspace(-real_z / 2, real_z / 2, Z.shape[0])

    # reorder lists so they conform to numpy (https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.fft.ifftn.html)
    start = list(X).index(0)
    X_reordered = np.concatenate((X[start:], X[:start]))
    ndx_x = [list(X).index(i) for i in X_reordered]

    start = list(Y).index(0)
    Y_reordered = np.concatenate((Y[start:], Y[:start]))
    ndx_y = [list(Y).index(i) for i in Y_reordered]

    start = list(Z).index(0)
    Z_reordered = np.concatenate((Z[start:], Z[:start]))
    ndx_z = [list(Z).index(i) for i in Z_reordered]

    SF_reordered = SF[ndx_x, :, :]
    SF_reordered = SF_reordered[:, ndx_y, :]
    SF_reordered = SF_reordered[:, :, ndx_z]

    # inverse fourier transform
    inverse_fft = np.fft.ifftn(SF_reordered)

    # reorder again
    inverse_fft = inverse_fft[ndx_x, :, :]
    inverse_fft = inverse_fft[:, ndx_y, :]
    inverse_fft = inverse_fft[:, :, ndx_z]

    # fourier transform of inversion as a test
    # ft = np.abs(np.fft.fftn(inverse_fft))**2
    # ft = ft[ndx_x, :]
    # ft = ft[:, ndx_y]
    # plt.imshow(ft)
    # plt.show()

    inverse_fft = inverse_fft.real / np.amax(inverse_fft.real)

    final, rfin, zfin = angle_average(X_real, Y_real, Z_real, inverse_fft, ucell=ucell)

    rbound1 = 0
    rbound2 = 0
    while rfin[rbound1] < -15:
        rbound1 += 1
    while rfin[rbound2] < 15:
        rbound2 += 1

    zbound1 = 0
    zbound2 = 0
    while zfin[0][zbound1] < -15:
        zbound1 += 1
    while zfin[0][zbound2] < 15:
        zbound2 += 1

    levels = np.linspace(np.amin(final), 0.001 * np.amax(final), 200)
    plt.contourf(rfin[rbound1:rbound2], zfin[0][zbound1:zbound2], final[rbound1:rbound2, zbound1:zbound2].T,
                 levels=levels, cmap='seismic', extend='max')
    plt.colorbar()
    plt.xlabel('r ($\AA$)')
    plt.ylabel('z ($\AA$)')
    plt.show()
    exit()


def angle_average(X, Y, Z, SF, ucell=None):

    ES = RegularGridInterpolator((X, Y, Z), SF, bounds_error=False)

    THETA_BINS_PER_INV_ANG = 20.
    MIN_THETA_BINS = 10  # minimum allowed bins
    RBINS = 100

    if ucell is not None:

        a1 = ucell[0]
        a2 = ucell[1]
        a3 = ucell[2]

        b1 = (np.cross(a2, a3)) / (np.dot(a1, np.cross(a2, a3)))
        b2 = (np.cross(a3, a1)) / (np.dot(a2, np.cross(a3, a1)))
        b3 = (np.cross(a1, a2)) / (np.dot(a3, np.cross(a1, a2)))

        b_inv = np.linalg.inv(np.vstack((b1, b2, b3)))

    ZBINS = Z.shape[0]  # 400

    XR = (X[-1] - X[0])
    YR = (Y[-1] - Y[0])

    Rmax = min(XR, YR) / 2.0
    Rmax *= 0.95

    rarr, rspace = np.linspace(0.0, Rmax, RBINS, retstep=True)
    zar = np.linspace(Z[0], Z[-1], ZBINS)

    oa = np.zeros((rarr.shape[0], zar.shape[0]))
    circ = 2.*np.pi*rarr  # circumference

    for ir in range(rarr.shape[0]):

        NTHETABINS = max(int(THETA_BINS_PER_INV_ANG*circ[ir]), MIN_THETA_BINS)  #calculate number of bins at this r
        thetas = np.linspace(0.0, np.pi*2.0, NTHETABINS, endpoint=False)  # generate theta array

        t, r, z = np.meshgrid(thetas, rarr[ir], zar)  # generate grid of cylindrical points

        xar = r*np.cos(t)  # set up x,y coords
        yar = r*np.sin(t)

        pts = np.vstack((xar.ravel(), yar.ravel(), z.ravel())).T  # reshape for interpolation

        if ucell is not None:
            # pts = mc_inv(pts, ucell)
            pts = np.matmul(pts, b_inv)

        oa[ir, :] = np.average(ES(pts).reshape(r.shape), axis=1)  # store average values in final array

    mn = np.nanmin(oa)
    oa = np.where(np.isnan(oa), mn, oa)

    rad_avg = np.average(oa)  # ???
    oa /= rad_avg  # normalize

    # set up data for contourf plot by making it symmetrical
    final = np.append(oa[::-1, :], oa[1:], axis=0)  # SF
    rfin = np.append(-rarr[::-1], rarr[1:])  # R
    zfin = np.append(z[:, 0, :], z[1:, 0, :], axis=0)  # Z

    return final, rfin, zfin


def Rspots(R, Z, waxs, theta=37, theta_sigma=(7, 5), bounds=(1.256, 1.57), cmap='jet'):

    """ Measure intensity of R-spots in specified region """

    spots = np.copy(waxs.T)
    inner = bounds[0]
    outer = bounds[1]
    I = []

    for i in range(R.shape[0]):
        for j in range(Z.shape[0]):
            if inner < np.linalg.norm([R[i], Z[j]]) < outer:
                angle = (180 / np.pi) * np.arctan(Z[j] / R[i])
                if (theta - theta_sigma[0]) < angle < (theta + theta_sigma[1]) or \
                        (theta - theta_sigma[0]) < (angle - 2*angle) < (theta + theta_sigma[1]):
                    spots[i, j] = 100
                    I.append(waxs[j, i])

    average_intensity = np.mean(I)

    plt.figure()
    levels = np.linspace(0, 3.1, 200)

    plt.contourf(R, Z, spots.T, cmap=cmap, levels=levels, extend='max')
    plt.xlim(-2.5, 2.5)
    plt.ylim(-2.5, 2.5)
    plt.figure()
    plt.hist(I, bins=25)
    plt.title('Average intensity of R-spots: %.2f' % average_intensity)
    # plt.show()

    return average_intensity


def gaussian(points, mean, sigma, amplitude, yshift):
    return yshift + (amplitude / np.sqrt(2 * np.pi * sigma ** 2)) * np.exp(
        -(points - mean) ** 2 / (2 * sigma ** 2))


def lorentz(points, a, b, c):
    """
    :param p: lorentzian parameters : [full width half max (FWHM), position of maximum, maximum heigth]
    :param p: position
    :return:
    """

    w = a / 2

    x = (b - points) / w

    return (c / (np.pi * w)) / (1 + x ** 2)


def triple_lorentz(x, a0, a1, a2, b0, b1, b2, c0, c1, c2):
    return lorentz(x, a0, b0, c0) + lorentz(x, a1, b1, c1) + lorentz(x, a2, b2, c2)


def PLOT_RAD_NEW(D, wavelength_angstroms, ucell, format=False, factor=3.1, **kwargs):

    """
    :param D: raw structure factor
    :param wavelength_angstroms: wavelength of X-ray (angstroms)
    :param ucell: 3 x 3 unitcell vectors
    :param factor: maximum colorbar value if using formatting from Coscia et al. manuscript
    :param format: plot simulated XRD patterns as they appear in Coscai et al. manuscript
    :return:
    """

    if not os.path.exists(path):
        os.makedirs(path)

    # inverse_ft(D, ucell)

    X = D[:, 0, 0, 0]
    Y = D[0, :, 0, 1]
    Z = D[0, 0, :, 2]
    SF = D[..., 3]

    ############## Plot z-slice down the middle of the raw structure factor ###################
    # plt.plot(Z, SF[len(X)//2, len(Y)//2, :])
    # plt.xlabel('q$_z$ ($\AA^{-1}$)')
    # plt.ylabel('Intensity')
    # plt.savefig('z_section.png')
    # plt.show()
    # exit()

    ES = RegularGridInterpolator((X, Y, Z), SF, bounds_error=False)

    THETA_BINS_PER_INV_ANG = 20.
    MIN_THETA_BINS = 1  # minimum allowed bins
    RBINS = 400
    NLEVELS = 200  # number of levels for contour plots

    a1 = ucell[0]
    a2 = ucell[1]
    a3 = ucell[2]

    b1 = (np.cross(a2, a3)) / (np.dot(a1, np.cross(a2, a3)))
    b2 = (np.cross(a3, a1)) / (np.dot(a2, np.cross(a3, a1)))
    b3 = (np.cross(a1, a2)) / (np.dot(a3, np.cross(a1, a2)))

    b_inv = np.linalg.inv(np.vstack((b1, b2, b3)))

    ZBINS = Z.shape[0]  # 400
    XR = (X[-1] - X[0])*ucell[0][0]
    YR = (Y[-1] - Y[0])*ucell[1][1]

    Rmax = min(XR, YR) / 2.0
    Rmax *= 0.95

    rarr, rspace = np.linspace(0.0, Rmax, RBINS, retstep=True)
    zar = np.linspace(Z[0], Z[-1], ZBINS)

    oa = np.zeros((rarr.shape[0], zar.shape[0]))
    circ = 2.*np.pi*rarr  # circumference

    for ir in trange(rarr.shape[0]):

        NTHETABINS = max(int(THETA_BINS_PER_INV_ANG*circ[ir]), MIN_THETA_BINS)  #calculate number of bins at this r
        thetas = np.linspace(0.0, np.pi*2.0, NTHETABINS, endpoint=False)  # generate theta array

        t, r, z = np.meshgrid(thetas, rarr[ir], zar)  # generate grid of cylindrical points

        xar = r*np.cos(t)  # set up x,y coords
        yar = r*np.sin(t)

        pts = np.vstack((xar.ravel(), yar.ravel(), z.ravel())).T  # reshape for interpolation

        MCpts = np.matmul(pts, b_inv)  # slower: MCpts = mc_inv(pts,ucell)

        oa[ir, :] = np.average(ES(MCpts).reshape(r.shape), axis=1)  # store average values in final array

    mn = np.nanmin(oa)
    oa = np.where(np.isnan(oa), mn, oa)

    if not format:
        rad_avg = np.average(oa)
        oa /= rad_avg  # normalize

    # set up data for contourf plot by making it symmetrical
    final = np.append(oa[::-1, :], oa[1:], axis=0)  # SF
    rfin = np.append(-rarr[::-1], rarr[1:])  # R
    zfin = np.append(z[:, 0, :], z[1:, 0, :], axis=0)  # Z

    unitlab = '($\AA^{-1}$)'  # Angstroms

    logfinal = np.log(final)

    MIN = np.amin(final)  # MIN = np.amin(np.ma.masked_invalid(final))
    MAX = np.amax(final)  # MAX = np.amax(np.ma.masked_invalid(final))

    lvls = np.linspace(MIN, MAX, NLEVELS)

    if format:
        alkane_intensity = normalize_alkanes(rfin, zfin[0], final, 1.4, 1.57, 120)  # 1.4, 1.57
        final /= alkane_intensity  # normalize according to R-alkanes

    # lvls = np.linspace(0, factor, NLEVELS)  # contour levels
    rlimits = [np.argmin(np.abs(rfin + 2.5)), np.argmin(np.abs(rfin - 2.5))]
    zlimits = [np.argmin(np.abs(zfin[0] + 2.5)), np.argmin(np.abs(zfin[0] - 2.5))]

    MIN = np.amin(final[rlimits[0]:rlimits[1], zlimits[0]:zlimits[1]])
    MIN = 0.4
    # MAX = 7.67

    #lvls = np.linspace(np.log10(MIN), np.log10(MAX), NLEVELS)

    if format:

        cmap = 'jet'
        print(factor)
        lvls = np.linspace(0, factor, NLEVELS)
        lvls_log = np.linspace(np.log10(final[-1, -1]), np.log10(np.amax(final)), NLEVELS)

        # plot 1D SAXS
        plt.figure()
        plt.plot(rfin, final[:, zfin[0].shape[0]//2], linewidth=2)
        plt.xlabel('$q_r\ (\AA^{-1})$', fontsize=14)
        plt.ylabel('Intensity', fontsize=14)
        plt.tight_layout()

        plt.figure()

        heatmap = plt.contourf(rfin, zfin[0], final.T, levels=lvls, cmap=cmap, extend='max')
        cbar = plt.colorbar(heatmap)
        plt.xlabel('$q_r\ (\AA^{-1}$)', fontsize=18)
        plt.ylabel('$q_z\ (\AA^{-1}$)', fontsize=18)
        plt.gcf().get_axes()[0].set_ylim(-2.5, 2.5)
        plt.gcf().get_axes()[0].set_xlim(-2.5, 2.5)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.gcf().get_axes()[0].set_aspect('equal')
        plt.tight_layout()
        plt.savefig('rzplot.png')
        print('rzplot.png saved')

        ################# Q_R and Q_Z CROSS_SECTIONS OF R_PI WITH GAUSSIAN AND LORENTZIAN FITS ##################
        ############################### FIT TO QR CROSS-SECTION OF R-PI #########################

        plt.figure()

        rpi_ndx = np.argmin(np.abs(zfin[0] - zfin[0][np.argmax(final[rfin.size // 2, :])]))

        plt.plot(rfin, final[:, rpi_ndx], linewidth=2, color='blue')  # its xkcd:blue in paper

        p = np.array([0, 0.3, 4, 1])
        solp, cov_x = curve_fit(gaussian, rfin, final[:, rpi_ndx], p,
                                bounds=((-np.inf, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf)))

        plt.plot(rfin, gaussian(rfin, solp[0], solp[1], solp[2], solp[3]), '--', color='blue', label='Gaussian Fit',
                 linewidth=2)

        print("Gaussian FWHM = %.3f +/- %.3f A^-1" % (2*np.sqrt(2*np.log(2))*solp[1],
                                               2 * np.sqrt(2 * np.log(2)) * cov_x[1, 1] ** 0.5))

        p = np.array([0.1, 0, 4])
        solp_lorentz, cov_x = curve_fit(lorentz, rfin, final[:, rpi_ndx], p,
                                bounds=[[0, -np.inf, 0], [np.inf, np.inf, np.inf]])

        plt.plot(rfin, lorentz(rfin, solp_lorentz[0], solp_lorentz[1], solp_lorentz[2]), '--', label='Lorentzian Fit',
                 linewidth=2, color='orange')  # its xkcd:orange in the paper

        print("Lorentzian FWHM = %.3f +/- %.3f A^-1" % (solp_lorentz[0], cov_x[0, 0] ** 0.5))

        plt.legend(fontsize=16)
        plt.xlabel('$q_r\ (\AA^{-1})$', fontsize=18)
        plt.ylabel('Intensity', fontsize=18)
        plt.gcf().get_axes()[0].tick_params(labelsize=18)
        plt.tight_layout()
        #plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/structure_paper/figures/sim_rsection_fit.pdf')

        # ######################## FIT TO QZ CROSS-SECTION OF R-PI #########################
        # plt.figure()
        #
        # rndx = rfin.size // 2
        # zstart = zfin[0].size // 2
        # plt.plot(zfin[0][zstart:], final[rndx, zstart:], linewidth=2, color='blue')
        #
        # p = np.array([1.4, 0.1, 7, 0])
        # solp, cov_x = curve_fit(gaussian, zfin[0][zstart:], final[rndx, zstart:], p,
        #                         bounds=([-np.inf, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
        #
        # fine_grid = np.linspace(zfin[0][zstart], zfin[0][-1], 1000)
        # plt.plot(fine_grid, gaussian(fine_grid, solp[0], solp[1], solp[2], solp[3]), '--', color='blue', label='Gaussian Fit',
        #          linewidth=2)
        #
        # print("Gaussian FWHM = %.3f +/- %.3f A^-1" % (2*np.sqrt(2*np.log(2))*solp[1],
        #                                        2 * np.sqrt(2 * np.log(2)) * cov_x[1, 1] ** 0.5))
        #
        # p = np.array([0.1, 0, 4])
        # solp_lorentz, cov_x = curve_fit(lorentz, zfin[0][zstart:], final[rndx, zstart:], p,
        #                         bounds=[[0, -np.inf, 0], [np.inf, np.inf, np.inf]])
        #
        # plt.plot(fine_grid, lorentz(fine_grid, solp_lorentz[0], solp_lorentz[1], solp_lorentz[2]), '--',
        #          label='Lorentzian Fit', linewidth=2, color='orange')
        #
        # print("Lorentzian FWHM = %.3f +/- %.3f A^-1" % (solp_lorentz[0], cov_x[0, 0] ** 0.5))
        #
        # plt.legend(fontsize=17)
        # plt.xlabel('$q_z\ (\AA^{-1})$', fontsize=18)
        # plt.ylabel('Intensity', fontsize=18)
        # plt.gcf().get_axes()[0].tick_params(labelsize=18)
        # plt.tight_layout()
        # #plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/structure_paper/figures/sim_zsection_fit.pdf')
        #
        # print('Average R-pi intensity: %.2f' % np.amax(final[rfin.size // 2, :]))
        #print('Average R-spots intensity : %.2f' % Rspots(rfin, zfin[0], final.T, theta=30, theta_sigma=(1, 1),
                                                          #bounds=(1.39, 1.49), cmap=cmap))

    else:

        plt.figure()
        plt.contourf(rfin, zfin[0], final.T, levels=lvls, cmap='jet')
        plt.colorbar()

        plt.title('S(r,z)')
        plt.xlabel('r ' + unitlab)
        plt.ylabel('z ' + unitlab)

        plt.savefig('new_rzplot.png')

        plt.figure()

        cs = plt.contourf(rfin, zfin[0], final.T, levels=lvls, cmap='jet', extend='both')
        cs.cmap.set_under('k')
        cs.set_clim(MIN, 0.1 * MAX)
        plt.title('S(r,z)')
        plt.xlabel('r ' + unitlab)
        plt.ylabel('z ' + unitlab)
        plt.colorbar()
        plt.savefig('cs.png')

        plt.figure()

        plt.contourf(rfin, zfin[0], final.T, levels=lvls, cmap='jet')
        plt.colorbar()

        plt.title('S(r,z)')
        plt.xlabel('r ' + unitlab)
        plt.ylabel('z ' + unitlab)
        plt.savefig('new_rzplot2.png')

        plt.figure()
        lglvls = np.linspace(np.amin(logfinal), np.amax(logfinal), NLEVELS)

        plt.contourf(rfin, zfin[0], logfinal.T, levels=lglvls, cmap='jet')
        plt.colorbar()

        plt.title('ln(S(r,z))')
        plt.xlabel('r ' + unitlab)
        plt.ylabel('z ' + unitlab)
        plt.savefig('new_log_rzplot.png')

        plt.figure()

        x2 = np.linspace(-Rmax, Rmax, RBINS * 2 - 1)
        z2 = np.linspace(Z[0], Z[-1], RBINS)

        xg2, yg2, zg2 = np.meshgrid(x2, np.asarray(0), z2)
        pts = np.vstack((xg2.ravel(), yg2.ravel(), zg2.ravel())).T
        out2 = ES(pts).reshape(xg2.shape[1], xg2.shape[2])

        o2n = out2[:, :] / rad_avg

        plt.contourf(xg2[0, :, :], zg2[0, :, :], o2n, levels=lvls, cmap='jet')

        plt.xlabel('x ' + unitlab)
        plt.ylabel('z ' + unitlab)
        plt.title('S(x,z)|$_{y=0}$')

        plt.colorbar()
        plt.savefig('new_xzplot.png')

        plt.figure()

        x2 = np.linspace(-Rmax, Rmax, RBINS * 2 - 1)
        y2 = np.linspace(-Rmax, Rmax, RBINS * 2 - 1)

        xg2, yg2, zg2 = np.meshgrid(x2, y2, np.asarray(0))
        pts = np.vstack((xg2.ravel(), yg2.ravel(), zg2.ravel())).T
        out2 = ES(pts).reshape(xg2.shape[0], xg2.shape[1])

        o2n = out2[:, :] / np.average(out2)

        lvlsxy = np.linspace(np.amin(o2n), np.amax(o2n), NLEVELS)  # contour levels

        plt.contourf(xg2[:, :, 0], yg2[:, :, 0], o2n, levels=lvlsxy, cmap='jet')

        plt.xlabel('x ' + unitlab)
        plt.ylabel('y ' + unitlab)
        plt.title('S(x,y)')  # |$_{y=0}$')

        plt.colorbar()
        plt.savefig('new_xyplot.png')

        if False:

            plt.figure()

            dif = o2n - final
            lvls2 = np.linspace(-0.4, 0.4, 100)

            plt.contourf(xg2[0, :, :], zg2[0, :, :], dif, levels=lvls2, cmap='seismic')
            plt.xlabel('x,r ' + unitlab)
            plt.ylabel('z ' + unitlab)
            plt.title('S(r,z)-S(x,z)|$_{y=0}$')

            plt.colorbar()
            plt.savefig('difference.png')

    plt.show()


def normalize_alkanes(R, Z, Raw_Intensity, inner, outer, angle):
    """
    Plot angular integration of 2D WAXS data bounded by a circle defined by radii 'inner' and 'outer'
    :param R: points in r direction
    :param Z: points in z direction
    :param Raw_Intensity: values at all (R, Z) points on grid
    :param inner: inside radius of region bounding alkane reflections
    :param outer: outside radius of region bounding alkane reflections
    :return: Intensity values normalized by average intensity inside alkane region
    """

    nbins = 90
    bins = np.linspace(-90, 90, nbins)

    bw = 180 / (nbins - 1)

    angles = []
    intensity = []
    for i in range(R.shape[0]):
        for j in range(Z.shape[0]):
            if inner < np.linalg.norm([R[i], Z[j]]) < outer:
                angles.append((180/np.pi)*np.arctan(Z[j]/R[i]))
                intensity.append(Raw_Intensity[i, j])

    inds = np.digitize(angles, bins)

    I = np.zeros([nbins])
    counts = np.zeros([nbins])
    for i in range(len(inds)):
        I[inds[i] - 1] += intensity[i]
        counts[inds[i] - 1] += 1

    #Get average intensity in ring excluding 60 degree slice around top and bottom #######

    bin_range = 180 / nbins  # degrees which a single bin covers

    start = int((angle/2) / bin_range)  # start at the bin which covers -60 degrees and above
    end = nbins - start  # because of symmetry

    total_intensity = np.sum(I[start:end])
    avg_intensity = total_intensity / np.sum(counts[start:end])

    print('Average Intensity in alkane chain region : %s' % avg_intensity)

    return avg_intensity


def tm2(D, ucell):

    a1 = ucell[0]
    a2 = ucell[1]
    a3 = ucell[2]

    V = np.dot(a1, np.cross(a2, a3))

    b1 = np.cross(a2, a3) / V
    b2 = np.cross(a3, a1) / V  # *2.0*math.pi
    b3 = np.cross(a1, a2) / V  #*2.0*math.pi

    Dnew = np.zeros_like(D)

    X = D[..., 0]
    Y = D[..., 1]
    Z = D[..., 2]

    for ix in range(D.shape[0]):
        Dnew[ix, 0:3] += X[ix]*b1  #(X[ix]-X[X.shape[0]/2])*b1

    for iy in range(D.shape[0]):
        Dnew[iy, 0:3] += Y[iy]*b2  #(Y[iy]-Y[Y.shape[0]/2])*b2

    for iz in range(D.shape[0]):
        Dnew[iz, 0:3] += Z[iz]*b3  #(Z[iz]-Z[Z.shape[0]/2])*b3

    return Dnew


def to_monoclinic(D, ucell):		#monoclinic for now

    a1 = ucell[0]
    a2 = ucell[1]
    a3 = ucell[2]

    b1 = (np.cross(a2, a3)) / (np.dot(a1, np.cross(a2, a3)))
    b2 = (np.cross(a3, a1)) / (np.dot(a2, np.cross(a3, a1)))#*2.0*math.pi
    b3 = (np.cross(a1, a2)) / (np.dot(a3, np.cross(a1, a2)))#*2.0*math.pi

    Dnew = np.zeros_like(D)

    X = D[..., 0]
    Y = D[..., 1]
    Z = D[..., 2]

    for ix in range(D.shape[0]):
        Dnew[ix, 0:3] += X[ix]*b1  #(X[ix]-X[X.shape[0]/2])*b1

    for iy in range(D.shape[0]):
        Dnew[iy, 0:3] += Y[iy]*b2  #(Y[iy]-Y[Y.shape[0]/2])*b2

    for iz in range(D.shape[0]):
        Dnew[iz, 0:3] += Z[iz]*b3  #(Z[iz]-Z[Z.shape[0]/2])*b3

    return Dnew


def mc_inv(D, ucell):

    a1 = ucell[0]
    a2 = ucell[1]
    a3 = ucell[2]

    b1 = (np.cross(a2, a3))/(np.dot(a1, np.cross(a2, a3)))
    b2 = (np.cross(a3, a1))/(np.dot(a2, np.cross(a3, a1)))
    b3 = (np.cross(a1, a2))/(np.dot(a3, np.cross(a1, a2)))

    b_inv = np.linalg.inv(np.vstack((b1, b2, b3)))
    Dnew = np.zeros_like(D)

    X = D[..., 0]
    Y = D[..., 1]
    Z = D[..., 2]

    for ix in range(D.shape[0]):
        Dnew[ix, 0:3] += X[ix]*b_inv[0]

    for iy in range(D.shape[0]):
        Dnew[iy, 0:3] += Y[iy]*b_inv[1]

    for iz in range(D.shape[0]):
        Dnew[iz, 0:3] += Z[iz]*b_inv[2]

    return Dnew


def Plot_Ewald_triclinic(D, wavelength_angstroms, ucell, factor=3.1, format=True, **kwargs):  # pass full 3d data,SF,wavelength in angstroms

    PLOT_RAD_NEW(D, wavelength_angstroms, ucell, factor=factor, format=format, **kwargs)
    exit()

    if not os.path.exists(path):
        os.makedirs(path)

    X = D[:, 0, 0, 0].copy()
    Y = D[0, :, 0, 1].copy()
    Z = D[0, 0, :, 2].copy()

    NBINSZ = 1 * D[0, 0, :, 2].size
    ZBNS = np.linspace(Z[0], Z[-1], NBINSZ)

    if NBINSRAD > 0:
        XBNSRD = np.linspace(-NBINSRAD, NBINSRAD, num=NBINSRAD*2)
        XBNSRD = np.sqrt(np.abs(XBNSRD))*np.sign(XBNSRD)
        XBNSRD *= (X[-1]/XBNSRD[-1])
    else:
        XBNSRD = X
        print("setting XBNSRD=", X)

    dx1 = X[1 + int(X.shape[0]/2)] - X[int(X.shape[0]/2)]

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

    K_ES = 2.0*math.pi/wavelength_angstroms  # calculate k for incident xrays in inverse angstroms

    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RegularGridInterpolator.html#scipy.interpolate.RegularGridInterpolator
    # Notes
    # Contrary to LinearNDInterpolator and NearestNDInterpolator, RegularGridInterpolator class avoids expensive triangulation of the input data by taking advantage of the regular grid structure.
    # this is why this style of interpolation is so slow

    XGD = D[:, :, :, 0]  # X spatial grid view
    YGD = D[:, :, :, 1]
    ZGD = D[:, :, :, 2]
    VGD = D[:, :, :, 3]

    DC = D[:, :, :, 0:3]

    DR = DC.reshape(DC.size/3, 3)

    # check if fast interpolation can be used
    lbuf = True
    for i in range(3):
        for j in range(i + 1, 3):
            if ucell[i, j] != 0 or ucell[j, i] != 0:
                lbuf = False

    print("Interpolating grid...")

    if ucell[0, 0] == ucell[1, 1] and ucell[0, 0] == ucell[2, 2] and lbuf:

        print("using fast interpolation for orthorhombic cell")
        ES = RegularGridInterpolator((X, Y, Z), SF, bounds_error=False)

    else:

        print("Interpolating non-orthorhombic cell")
        dtime = 480.0 * XGD.size / (98 * 98 * 99)  # empirical time estimate

        print("interpolation time estimate: ", round(dtime / 60, 1), " minutes, finishing around ", (
                    datetime.datetime.now() + datetime.timedelta(seconds=dtime)).strftime('%I:%M %p'))

        start = time.time()
        coords = list(zip(XGD.ravel(), YGD.ravel(), ZGD.ravel()))

        if False:
            ES = LinearNDInterpolator(coords, VGD.ravel())
        else:
            ES = NearestNDInterpolator(coords, VGD.ravel())
        end = time.time()

        print("interpolation finished, taking %4.2f seconds" % (end-start))

    xyzpts = np.asarray([])
    print("setting up points for radial integration")
    Scale=1

    if False:
        for ix in trange(D.shape[0]):
            for iy in range(D.shape[1]):
                for iz in range(D.shape[2]):
                    xyzpts.append((D[ix, iy, iz, 0], D[ix, iy, iz, 1], D[ix, iy, iz, 2]))
    else:
        XPTS = np.linspace(D[0, 0, 0, 0], D[-1, 0, 0, 0], Scale * D.shape[0], dtype=np.float16)
        YPTS = np.linspace(D[0, 0, 0, 1], D[0, -1, 0, 1], Scale * D.shape[1], dtype=np.float16)
        ZPTS = np.linspace(D[0, 0, 0, 2], D[0, 0, -1, 2], Scale * D.shape[2], dtype=np.float16)
        print("mesh")
        xyzpts = np.meshgrid(XPTS, YPTS, ZPTS)
        print("stack")
        xyzpts = np.stack(xyzpts, -1).reshape(-1, 3)
        print("done")

    xyzpts = np.reshape(D[:, :, :, :3], (D.shape[0]*D.shape[1]*D.shape[2], 3))  # 5000x faster than above loop

    NSP = 20
    NSP = np.minimum(NSP, xyzpts.shape[0])  # split into at most 20 chunks before processing to limit memory usage

    xyzpieces = np.array_split(xyzpts, NSP)
    EWDxyz = np.asarray([])
    print("interpolating")
    for i in tqdm.tqdm(xyzpieces):
        buf = ES(i)
        EWDxyz = np.append(EWDxyz, buf, axis=0)

    print("EWD done")

    rpts = np.sqrt(xyzpts[:, 0]**2.0 + xyzpts[:, 1]**2.0)

    Hcount, XEC, YEC = np.histogram2d(rpts, xyzpts[:, 2], bins=(XBNSRD, ZBNS))

    Hval, XEV, YEV = np.histogram2d(rpts, xyzpts[:, 2], weights=EWDxyz, normed=False, bins=(XBNSRD, ZBNS))

    switch1 = True

    if switch1:
        Hcount = np.where(Hcount == 0, 1, Hcount)

    Hrz = Hval / Hcount

    if not switch1:
        Hrz = np.ma.masked_invalid(Hrz)

    S1 = np.sum(Hrz)
    S3 = np.sum(Hrz[Hrz.shape[0]/2, :])

    Condition1 = False # Need to figure this out-when should this be true?

    if Condition1:
        for ir in range(1, Hrz.shape[0] / 2 - 1):
            Hrz[-ir + Hrz.shape[0] / 2, :] = Hrz[ir + Hrz.shape[0] / 2,
                                             :]  # this needs to be tested for both even and odd numbers of bins
    else:
        for ir in range(1, Hrz.shape[0] / 2 - 1):
            Hrz[-ir + 2 + Hrz.shape[0] / 2, :] = Hrz[ir + Hrz.shape[0] / 2,
                                                 :]  # this needs to be tested for both even and odd numbers of bins


    S2 = np.sum(Hrz)

    XMG, YMG = np.meshgrid(XEV, YEV)

    plt.pcolormesh(XMG[:-1, :], YMG[:-1, :], np.log10(Hrz.T), vmin=np.amin(np.log10(Hrz)), vmax=np.amax(np.log10(Hrz)))
    plt.savefig(path+"_log_rzplot"+format, dpi=DPI)
    plt.clf()
    print("_log_rzplot saved")

    mn = np.amin(Hrz[np.nonzero(Hrz)])
    Hbuf = np.where(Hrz > 0.0, Hrz, mn)
    Log_HRZ = np.log10(Hbuf)

    plt.pcolormesh(XMG[:-1, :] - dx1 / 2.0, YMG[:-1, :], Log_HRZ.T, vmin=np.amin(Log_HRZ), vmax=np.amax(Log_HRZ),
                   cmap='nipy_spectral')
    plt.colorbar()
    plt.savefig(path + "_log_rzplot" + format, dpi=DPI)
    plt.clf()

    Nx = D.shape[0]
    Ny = D.shape[1]
    Nz = D.shape[2]

    #==============flat and Ewald-corrected plots=================

    xypts = []
    xyflat = []
    for ix in range(D.shape[0]):
        for iy in range(D.shape[1]):
            xp = D[ix, iy, int(Nz/2), 0]
            yp = D[ix, iy, int(Nz/2), 1]

            theta = np.arctan(np.sqrt(xp**2.0 + yp**2.0)/K_ES)
            xypts.append((xp*np.cos(theta), yp*np.cos(theta), K_ES*(1.0 - np.cos(theta))))
            xyflat.append((xp, yp, 0.0))

    xzpts = []
    xzflat = []

    for ix in range(D.shape[0]):
        for iz in range(D.shape[2]):
            xp = D[ix, int(Ny/2), iz, 0]
            zp = D[ix, int(Ny/2), iz, 2]
            theta = np.arctan(np.sqrt(xp**2.0 + yp**2.0)/K_ES)
            xzpts.append((xp*np.cos(theta), K_ES*(1.0-np.cos(theta)), zp*np.cos(theta)))
            xzflat.append((xp, 0.0, zp))

    yzpts = []
    yzflat = []
    for iy in range(D.shape[1]):
        for iz in range(D.shape[2]):
            yp = D[int(Nz/2), iy, iz, 1]
            zp = D[int(Nz/2), iy, iz, 2]
            theta = np.arctan(np.sqrt(yp**2.0 + zp**2.0)/K_ES)
            yzpts.append((K_ES*(1.0-np.cos(theta)), yp*np.cos(theta), zp*np.cos(theta)))
            yzflat.append((0.0, yp, zp))

    xypts = np.asarray(xypts)
    xzpts = np.asarray(xzpts)
    yzpts = np.asarray(yzpts)

    xyflat = np.asarray(xyflat)
    xzflat = np.asarray(xzflat)
    yzflat = np.asarray(yzflat)

    EWDxy = ES(xypts)
    EWDxz = ES(xzpts)
    EWDyz = ES(yzpts)

    EWDxyflat = ES(xyflat)
    EWDxzflat = ES(xzflat)
    EWDyzflat = ES(yzflat)

    EWDxy = EWDxy.reshape(D.shape[0], D.shape[1])
    EWDxz = EWDxz.reshape(D.shape[0], D.shape[2])
    EWDyz = EWDyz.reshape(D.shape[1], D.shape[2])

    EWDxyflat = EWDxyflat.reshape(D.shape[0], D.shape[1])
    EWDxzflat = EWDxzflat.reshape(D.shape[0], D.shape[2])
    EWDyzflat = EWDyzflat.reshape(D.shape[1], D.shape[2])

    title = "Ewald Corrected Structure Factor \n $\lambda=$"+str(wavelength_angstroms)+" $\AA$   $k_{ew}=$"+str(round(K_ES,2))+" $\AA^{-1}$"
    ltitle = 'log ' + title

    xlab = 'x ('+units + ")"
    ylab = 'y ('+units + ")"
    zlab = 'z ('+units + ")"

    fname = "Ewald_"

    iz = 0
    plt.suptitle(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.contourf(D[:, :, iz, 0], D[:, :, iz, 1], EWDxy, contours, **kwargs)
    plt.savefig(path+fname+"xy"+str(iz)+format,dpi=DPI)
    plt.clf()

    lax = ['x', 'y', 'z']

    ewlab = "Ewald"
    flab = "Flat"

    iax1 = 0
    iax2 = 1

    EWDxy = np.ma.masked_invalid(EWDxy)
    EWDxyflat = np.ma.masked_invalid(EWDxyflat)

    EWDxz = np.ma.masked_invalid(EWDxz)
    EWDxzflat = np.ma.masked_invalid(EWDxzflat)

    EWDyz = np.ma.masked_invalid(EWDyz)
    EWDyzflat = np.ma.masked_invalid(EWDyzflat)

    if PLOT_EWALDS:
        csplot_wlog(D[:, :, int(Nz / 2) + 1, iax1], D[:, :, int(Nz / 2) + 1, iax2], EWDxy, contours, ewlab, lax[iax1],
                    lax[iax2], **kwargs)

    csplot_wlog(D[:,:,int(Nz/2)+1,iax1],D[:,:,int(Nz/2)+1,iax2],EWDxyflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)

    iax1 = 0
    iax2 = 2
    if PLOT_EWALDS:
        csplot_wlog(D[:, int(Ny / 2), :, iax1], D[:, int(Ny / 2), :, iax2], EWDxz, contours, ewlab, lax[iax1],
                    lax[iax2], **kwargs)

    csplot_wlog(D[:,int(Ny/2),:,iax1],D[:,int(Ny/2),:,iax2],EWDxzflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)

    iax1 = 1
    iax2 = 2
    if PLOT_EWALDS:
        csplot_wlog(D[int(Nx / 2), :, :, iax1], D[int(Nx / 2), :, :, iax2], EWDyz, contours, ewlab, lax[iax1],
                    lax[iax2], **kwargs)

    csplot_wlog(D[int(Nx/2),:,:,iax1],D[int(Nx/2),:,:,iax2],EWDyzflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)

