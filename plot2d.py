from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from past.utils import old_div
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import NearestNDInterpolator
import math
from tqdm import trange
import time
import datetime
import matplotlib
print('Fuck you')
matplotlib.rc('axes', color_cycle=['r', 'g', 'b', '#004060'])

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

title_fontsize = 9

path = ""


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
    # cbar = fig.colorbar(cax)

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
    plt.ylim(0, 0.5)
    plt.xlim(0.1, 0.5)
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
    EWD = EWD.reshape(D.shape[0],D.shape[1])
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
    # plt.xlim([-2.5, 2.5])
    # plt.ylim([-2.5, 2.5])
    plt.savefig(path + fname + "yzlog" + format, dpi=DPI)
    plt.clf()


def PLOT_RAD_NEW(D, wavelength_angstroms, ucell, **kwargs):

    if not os.path.exists(path):
        os.makedirs(path)

    X = D[:, 0, 0, 0]
    Y = D[0, :, 0, 1]
    Z = D[0, 0, :, 2]
    SF = D[..., 3]

    ES = RegularGridInterpolator((X, Y, Z), SF, bounds_error=False)

    THETA_BINS_PER_INV_ANG = 20.
    MIN_THETA_BINS = 10  # minimum allowed bins
    RBINS = 400
    NLEVELS = 200  # number of levels for contour plots

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

        MCpts = to_monoclinic(pts, ucell)  # transform to monoclinic cell
        # MCpts = pts

        oa[ir, :] = np.average(ES(MCpts).reshape(r.shape), axis=1)  # store average values in final array
        # if ir == 373:
        #     print(ES(MCpts).reshape(r.shape).shape)
            # print(np.mean(ES(MCpts).reshape(r.shape), axis=1))

    # rad_avg = np.average(oa)
    rad_avg = np.mean(np.ma.masked_invalid(oa))

    oa /= rad_avg  # normalize

    # set up data for contourf plot

    final = np.append(oa[::-1, :], oa[1:], axis=0)  # SF
    rfin = np.append(-rarr[::-1], rarr[1:])  # R
    zfin = np.append(z[:, 0, :], z[1:, 0, :], axis=0)  # Z

    # plt.plot(rfin, final.T[int(zfin.shape[1] / 2), :])
    # plt.show()
    # exit()
    unitlab = '($\AA^{-1}$)'  # Angstroms

    MIN = np.amin(np.ma.masked_invalid(final))
    MAX = np.amax(np.ma.masked_invalid(final))

    print(MIN, MAX, NLEVELS)
    lvls = np.linspace(MIN, MAX, NLEVELS)  # contour levels

    cs = plt.contourf(rfin, zfin[0], final.T, levels=lvls, cmap='jet')
    cs.cmap.set_under('k')
    cs.set_clim(MIN, 0.1*MAX)
    # fig, ax = plt.subplots()
    # heatmap = ax.pcolormesh(rfin, zfin[0], final.T, cmap='jet', vmax=0.1)  # jet matches experiment
    #
    # cbar = plt.colorbar(heatmap)
    # from matplotlib import ticker
    # tick_locator = ticker.MaxNLocator(nbins=5)
    # cbar.locator = tick_locator
    # cbar.update_ticks()
    # cbar.ax.set_yticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1'])
    # plt.gcf().get_axes()[0].set_ylim(-2.5, 2.5)
    # plt.gcf().get_axes()[0].set_xlim(-2.5, 2.5)
    # plt.gcf().get_axes()[0].set_xlabel('$q_r$' + unitlab, fontsize=14)
    # plt.gcf().get_axes()[0].set_ylabel('$q_z$' + unitlab, fontsize=14)
    # plt.gcf().get_axes()[0].set_aspect('equal')
    # plt.gcf().get_axes()[0].tick_params(labelsize=14)
    # plt.tight_layout()
    # plt.savefig('new_rzplot.png')
    # fig.clf()

    plt.colorbar()

    plt.title('S(r,z)')
    plt.xlabel('r ' + unitlab)
    plt.ylabel('z ' + unitlab)
    plt.gcf().get_axes()[0].set_ylim(-2.5, 2.5)
    plt.gcf().get_axes()[0].set_xlim(-2.5, 2.5)
    plt.savefig('new_rzplot.png')
    plt.clf()
    print('new_rzplot.png saved')

    x2 = np.linspace(-Rmax, Rmax, RBINS*2 - 1)
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
    plt.clf()

    if False:
        dif = o2n - final
        lvls2 = np.linspace(-0.4, 0.4, 100)

        plt.contourf(xg2[0, :, :], zg2[0, :, :], dif, levels=lvls2, cmap='seismic')
        plt.xlabel('x,r ' + unitlab)
        plt.ylabel('z ' + unitlab)
        plt.title('S(r,z)-S(x,z)|$_{y=0}$')

        plt.colorbar()
        plt.savefig('difference.png')


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


def Plot_Ewald_triclinic(D, wavelength_angstroms, ucell, load, **kwargs):  #pass full 3d data,SF,wavelength in angstroms

    PLOT_RAD_NEW(D, wavelength_angstroms, ucell, **kwargs)

    rzscale = kwargs["rzscale"]

    if not os.path.exists(path):
        os.makedirs(path)

    X = D[:, 0, 0, 0].copy()
    Y = D[0, :, 0, 1].copy()
    Z = D[0, 0, :, 2].copy()

    if NBINSRAD > 0:
        XBNSRD = np.linspace(-NBINSRAD, NBINSRAD, num = NBINSRAD*2)
        XBNSRD = np.sqrt(np.abs(XBNSRD))*np.sign(XBNSRD)
        XBNSRD *= (X[-1]/XBNSRD[-1])
    else:
        XBNSRD=X

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

    print("Interpolating grid...")
    if not load:

        if ucell[0, 0] == ucell[1, 1] and ucell[0, 0] == ucell[2, 2] and lbuf:

            ES = RegularGridInterpolator((X, Y, Z), SF, bounds_error=False)

        else:

            # ES = np.load("ES.npy")
            dtime = 480.0*XGD.size / (98 * 98 * 99)  # empirical time estimate

            start = time.time()
            coords = list(zip(XGD.ravel(), YGD.ravel(), ZGD.ravel()))

            if False:
                ES = LinearNDInterpolator(coords, VGD.ravel())
            else:
                ES = NearestNDInterpolator(coords, VGD.ravel())
            end = time.time()

            print("interpolation finished, taking %4.2f seconds" % (end-start))

    xyzpts = []
    print("setting up points for radial integration")

    # if True:
    #     for ix in trange(D.shape[0]):
    #         for iy in range(D.shape[1]):
    #             for iz in range(D.shape[2]):
    #                 # xyzpts.append((X[ix],Y[iy],Z[iz]))
    #                 xyzpts.append((D[ix, iy, iz, 0], D[ix, iy, iz, 1], D[ix, iy, iz, 2]))
    # else:
    #     pass
    #
    # xyzpts = np.asarray(xyzpts)

    xyzpts = np.reshape(D[:, :, :, :3], (D.shape[0]*D.shape[1]*D.shape[2], 3))  # 5000x faster than above loop

    if not load:
        EWDxyz = ES(xyzpts)
        np.save("EWDxyz", EWDxyz)
    else:
        EWDxyz = np.load("EWDxyz.npy")

    rpts = np.sqrt(xyzpts[:, 0]**2.0 + xyzpts[:, 1]**2.0)

    # plt.hist(rpts, bins=80)
    # plt.show()
    # exit()

    Hcount, XEC, YEC = np.histogram2d(rpts, xyzpts[:, 2], bins=(XBNSRD, Z))

    Hval, XEV, YEV = np.histogram2d(rpts, xyzpts[:, 2], weights=EWDxyz, normed=True, bins=(XBNSRD, Z))

    switch1=True

    if switch1:
        Hcount = np.where(Hcount == 0, 1, Hcount)

    Hrz = Hval / Hcount

    if not switch1:
        Hrz = np.ma.masked_invalid(Hrz)

    for ir in range(old_div(Hrz.shape[0], 2)):
        Hrz[-ir + old_div(Hrz.shape[0], 2), :] = Hrz[ir + 1 + old_div(Hrz.shape[0], 2), :]

    np.save("Hrz", Hrz)

    # with open('Intensities.txt', 'w') as f:
    #
    #     for i in range(np.shape(Hrz)[0]):
    #         line = ''
    #         for j in range(np.shape(Hrz)[1]):
    #             line += '{:3.2f}'.format(Hrz[i, j])
    #
    #         f.write(line + '\n')

    XMG, YMG = np.meshgrid(XEV, YEV)

    ###################################################################################################################
    # # Plot angular integration of 2D WAXS data bounded by a circle defined by radii 'lower' and 'upper'

    lower = 1.4  #
    upper = 1.57

    nbins = 45
    bins = np.linspace(-90, 90, nbins)

    bw = 180 / (nbins - 1)

    angles = []
    intensity = []
    for i in range(XMG.shape[0]):
        for j in range(XMG.shape[1]):
            if lower < np.linalg.norm([XMG[i, j], YMG[i, j]]) < upper:
                angles.append((180/np.pi)*np.arctan(YMG[i, j]/XMG[i, j]))
                intensity.append(Hrz[j, i])

    inds = np.digitize(angles, bins)
    I = np.zeros([nbins])
    counts = np.zeros([nbins])
    for i in range(len(inds)):
        I[inds[i]] += intensity[i]
        counts[inds[i]] += 1

    #Get average intensity in ring excluding 60 degree slice around top and bottom #######

    bin_range = 180 / nbins  # degrees which a single bin covers

    start = int(30 / bin_range)  # start at the bin which covers -60 degrees and above
    end = nbins - start  # because of symmetry

    total_intensity = np.sum(I[start:end])
    avg_intensity = total_intensity / np.sum(counts[start:end])
    avg_intensity = 0.0042

    print('Average Intensity in alkane chain region : %s' % avg_intensity)
    #
    # I /= (counts*np.amax(intensity))
    #
    # plt.bar(bins, I, bw, color='#1f77b4')
    #
    # plt.xlabel('Angle with respect to $q_z=0$', fontsize=14)
    # plt.ylabel('Normalized integrated intensity', fontsize=14)
    # plt.gcf().get_axes()[0].tick_params(labelsize=14)
    #
    # plt.tight_layout()
    # plt.savefig('angular_integration.png')
    # plt.show()
    ###################################################################################################################

    # m = np.amax(Hrz[:, :-49])  # 49 for full sized thing a majig, 23 for half size
    # # exit()
    # # m = 0.00832122075357  # offset
    # # m = 0.00935587721282  # layered
    # m = 0.0036316001  # alkanes (average of PD and sandwiched average intensity in alkane region)
    # # m *= factor
    factor = 2.5  # this needs to match the experimental factor (see WAXS.py)
    m = avg_intensity*factor

    fig, ax = plt.subplots()
    #plt.pcolormesh(XMG[:-1, :], YMG[:-1, :], Hrz.T, vmin=0.0, vmax=rzscale*np.amax(Hrz), cmap='viridis')
    heatmap = ax.pcolormesh(XMG[:-1, :], YMG[:-1, :], Hrz.T/m, vmin=0.0, vmax=1.0, cmap='jet')  # jet matches experiment
    # heatmap = ax.pcolormesh(XMG[:-1, :], YMG[:-1, :], Hrz.T, cmap='jet')  # jet matches experiment

    # heatmap = ax.pcolormesh(XMG[:-50, :-49], YMG[:-50, :-49], Hrz.T[:-49, :-49]/m, vmin=0.0, vmax=1.0, cmap='jet')  # jet matches experiment

    ## Use this block (and comment out previous line to plot raw waxs data ##
    # waxs = np.load('/home/bcoscia/PycharmProjects/MD-Structure-Factor/waxs.npy')
    # qmax = 2.5
    # heatmap = ax.imshow(waxs, cmap='jet', extent=[-qmax, qmax, -qmax, qmax])
    ## end ##

    cbar = plt.colorbar(heatmap)
    from matplotlib import ticker
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.ax.set_yticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1'])

    # x = XMG[0, :]
    # y = YMG[:, 0]
    #
    # bins = 40
    # angles = np.zeros([bins])
    # norm = np.zeros([bins])
    #
    # area_to_integrate = np.zeros([x.shape[0], y.shape[0]])
    # r_inner = 1.1
    # r_outer = 1.7
    # for i in range(x.shape[0]):
    #     for k in range(y.shape[0]):
    #         if r_inner <= np.linalg.norm([x[i], y[k]]) <= r_outer:
    #             area_to_integrate[i, k] = Hrz[i, k]
    #             angle = np.arctan(y[k]/x[i]) * (180/np.pi)
    #
    #             bin = int((bins/2) + ((angle/90)*(bins/2)))
    #             if bin == bins:
    #                 bin -= 1
    #
    #             angles[bin] += Hrz[i, k]
    #             norm[bin] += 1
    #
    # avg = angles / norm  # normalize intensities so it is on a per count basis
    # avg /= np.sum(np.ma.masked_invalid(avg))
    #
    # angles = np.linspace(-90, 90, bins + 1)  # We will only see angles in the range of -90 to 90 since we use np.arctan
    # bin_angles = np.array([(angles[i] + angles[i + 1])/2 for i in range(bins)])  # bars will be placed in the middle of the bins
    # width = angles[1] - angles[0]  # width of bins
    #
    # keep = np.array([i for i in range(len(avg)) if not np.isnan(avg[i])])
    #
    # plt.figure()
    # plt.bar(bin_angles[keep], avg[keep], width=width, color='#1f77b4')
    # plt.xlabel('Angle with x axis', fontsize=14)
    # plt.ylabel('Normalized integrated intensity', fontsize=14)
    # plt.gcf().get_axes()[0].tick_params(labelsize=14)
    # plt.tight_layout()
    # plt.savefig('angle_v_I.png')
    # # plt.imshow(area_to_integrate.T)
    # plt.show()
    # exit()

    # fig.ylim(np.amin(YMG), np.amax(YMG))
    # fig.xlim(np.amin(XMG), np.amax(XMG))
    #

    plt.gcf().get_axes()[0].set_ylim(-2.5, 2.5)
    plt.gcf().get_axes()[0].set_xlim(-2.5, 2.5)
    plt.gcf().get_axes()[0].set_xlabel('$q_r$', fontsize=14)
    plt.gcf().get_axes()[0].set_ylabel('$q_z$', fontsize=14)
    plt.gcf().get_axes()[0].set_aspect('equal')
    plt.gcf().get_axes()[0].tick_params(labelsize=14)
    plt.tight_layout()
    fig.savefig(path+"rzplot"+format, dpi=DPI)
    fig.clf()

    measure_intensity = False

    if measure_intensity:

        from PIL import Image
        im = Image.open(path+"rzplot"+format)
        pix = im.load()
        left = 0
        upper = 0
        right = im.size[0]
        lower = im.size[1]
        box = (left, upper, right, lower)
        box = (428, 171, 1542, 1283)
        region = im.crop(box)
        region.show()
        print("Please crop the image so only the diffraction pattern is showing")
        response = input("Is this good? ")

        while response != 'yes':
            new_box_info = input("Please enter new box dimensions: ")
            new = new_box_info.split()
            if 'left' in new:
                left += int(new[new.index('left') + 1])
            if 'upper' in new:
                upper += int(new[new.index('upper') + 1])
            if 'right' in new:
                right += int(new[new.index('right') + 1])
            if 'lower' in new:
                lower += int(new[new.index('lower') + 1])
            if 'zoom' in new:
                factor = float(new[new.index('zoom') + 1])
                print("Zooming in by a factor of %s" % factor)
                left += (right - left) / (factor * 2)
                right -= (right - left) / (factor * 2)
                upper += (lower - upper) / (factor * 2)
                lower -= (lower - upper) / (factor * 2)

            box = (left, upper, right, lower)
            region = im.crop(box)
            region.show()
            response = input("Is this good? ")

        full_plot = box
        yrange = full_plot[2] - full_plot[0]  # right - left
        xrange = full_plot[3] - full_plot[1]  # lower - upper

        print(full_plot)

        print("Now crop the picture to the region where you'd like to measure average intensity")

        response = input("Is this good? ")
        while response != 'yes':
            new_box_info = input("Please enter new box dimensions: ")
            new = new_box_info.split()
            if 'left' in new:
                left += int(new[new.index('left') + 1])
            if 'upper' in new:
                upper += int(new[new.index('upper') + 1])
            if 'right' in new:
                right += int(new[new.index('right') + 1])
            if 'lower' in new:
                lower += int(new[new.index('lower') + 1])
            if 'zoom' in new:
                factor = float(new[new.index('zoom') + 1])
                print("Zooming in by a factor of %s" % factor)
                left += (right - left) / (factor * 2)
                right -= (right - left) / (factor * 2)
                upper += (lower - upper) / (factor * 2)
                lower -= (lower - upper) / (factor * 2)

            box = (left, upper, right, lower)
            region = im.crop(box)
            region.show()
            response = input("Is this good? ")

        region = box
        print(region)
        print(Hrz.shape)
        print(full_plot)
        print(xrange, yrange)
        lx = int(Hrz.shape[1]*(region[0] - full_plot[0]) / xrange)
        rx = int(Hrz.shape[1]*(region[2] - full_plot[0]) / xrange)
        ty = Hrz.shape[0] - int(Hrz.shape[0]*(region[1] - full_plot[1]) / yrange)
        by = Hrz.shape[0] - int(Hrz.shape[0]*(region[3] - full_plot[1]) / yrange)
        print(lx,rx,ty,by)

        #plt.pcolormesh(XMG[lx:rx, by:ty], YMG[lx:rx, by:ty], Hrz.T[lx:rx, by:ty], vmin=0.0, vmax=0.007*np.amax(Hrz))
        plt.pcolormesh(XMG[by:ty, lx:rx], YMG[by:ty, lx:rx], Hrz.T[by:ty, lx:rx], vmin=0.0, vmax=0.007*np.amax(Hrz))
        plt.show()
        print(np.mean(Hrz.T[by:ty, lx:rx]))

    plt.pcolormesh(XMG[:-1, :], YMG[:-1, :], np.log10(Hrz.T), vmin=np.amin(np.log10(Hrz)), vmax=np.amax(np.log10(Hrz)))
    plt.savefig(path+"_log_rzplot"+format, dpi=DPI)
    plt.clf()
    print("_log_rzplot saved")

    Nx = D.shape[0]
    Ny = D.shape[1]
    Nz = D.shape[2]

#==============flat and Ewald-corrected plots=================

    xypts=[]
    xyflat=[]
    for ix in range(D.shape[0]):
        for iy in range(D.shape[1]):
            xp = D[ix, iy, int(Nz/2), 0]
            yp = D[ix, iy, int(Nz/2), 1]

            theta = np.arctan(np.sqrt(xp**2.0 + yp**2.0)/K_ES)
            xypts.append((xp*np.cos(theta), yp*np.cos(theta), K_ES*(1.0 - np.cos(theta))))
            xyflat.append((xp, yp, 0.0))

    xzpts =[]
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

    lax=['x','y','z']

    ewlab = "Ewald"
    flab = "Flat"

    iax1=0
    iax2=1

    EWDxy = np.ma.masked_invalid(EWDxy)
    EWDxyflat = np.ma.masked_invalid(EWDxyflat)

    EWDxz = np.ma.masked_invalid(EWDxz)
    EWDxzflat = np.ma.masked_invalid(EWDxzflat)

    EWDyz = np.ma.masked_invalid(EWDyz)
    EWDyzflat = np.ma.masked_invalid(EWDyzflat)

    if PLOT_EWALDS:
        csplot_wlog(D[:, :, int(Nz/2)+1, iax1],D[:,:,int(Nz/2)+1,iax2],EWDxy, contours,ewlab,lax[iax1],lax[iax2],**kwargs)

    csplot_wlog(D[:,:,int(Nz/2)+1,iax1],D[:,:,int(Nz/2)+1,iax2],EWDxyflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)

    iax1=0
    iax2=2
    if PLOT_EWALDS:
        csplot_wlog(D[:,int(Ny/2),:,iax1],D[:,int(Ny/2),:,iax2],EWDxz,    contours,ewlab,lax[iax1],lax[iax2],**kwargs)

    csplot_wlog(D[:,int(Ny/2),:,iax1],D[:,int(Ny/2),:,iax2],EWDxzflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)

    iax1=1
    iax2=2
    if PLOT_EWALDS:
        csplot_wlog(D[int(Nx/2),:,:,iax1],D[int(Nx/2),:,:,iax2],EWDyz,    contours,ewlab,lax[iax1],lax[iax2],**kwargs)

    csplot_wlog(D[int(Nx/2),:,:,iax1],D[int(Nx/2),:,:,iax2],EWDyzflat,contours,flab ,lax[iax1],lax[iax2],**kwargs)
