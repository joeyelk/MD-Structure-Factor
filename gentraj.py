#!/usr/bin/env python

import numpy as np
import mdtraj as md
import argparse
import scipy.stats
import tqdm

"""
Create an arbitrarily shapped triclinic box uniformly filled with particles. The particles are sodium ions by default.
Change the particle type with the -particle flag. Define the box angles and vector lengths using the -box and -angles
flags.

One can also create features in the box, however this requires a monoclinic unitcell (alpha = beta = 90 degrees).
You can create layers, pores and disks which are stacked in layers. (These might work okay for arbitrary triclinic boxes
but it has not been tested). Create the features with the appropriate --pores, --layers or --disks flags.
Edit feature parameters (pore radius, distance between layers etc.) using the relevant flags, which are grouped together,
in the initialize() function.

By default, all features are created by not allowing particles in the defined regions. Using the --gaussian flag will
allow the features to be filled with points with a gaussian probability which is a function of pore/disk radius or
layer width. The lowest probability is in the middle of those features. NOTE: Gaussian disks are not implemented yet.
Using the -invert flag will allow particles in the defined regions only, thus inverting the configurations.

The output is a .trr GROMACS trajectory file and .gro coordinate file of the last frame. Specify the number of frames
with the -f flag and the name of the output files with the -o flag. Specify the number of particles in each frame with
the -n flag.
"""


def initialize():

    parser = argparse.ArgumentParser(description='Create random configurations with certain features while the rest of '
                                                 'space is filled uniformly with particles')

    parser.add_argument('-o', '--output', default='out', type=str, help='Name of output .gro and .trr files')
    parser.add_argument('-box', '--box_lengths', nargs='+', default=[8, 8, 8], type=float, help='box vector lengths '
                                                                                                '[x y z] ')
    parser.add_argument('-angles', '--angles', nargs='+', default=[90, 90, 60], type=float, help='angles between box'
                        'vectors (yz, xz, xy) -- the order of these is inspired by the gromacs tool "editconf"')
    parser.add_argument('-n', '--npoints', default=10000, type=int, help='number of points to place in each frame')
    parser.add_argument('-f', '--nframes', default=100, type=int, help='number of frames')
    parser.add_argument('-particle', default='NA', type=str, help='name of atom/particle to place in space')
    parser.add_argument('--gaussian', action="store_true", help='Allow points inside features with gaussian probability')
    # if you are making pores
    parser.add_argument('--pores', action="store_true", help='build a system with pores')
    parser.add_argument('-r', '--radius', default=0.5, type=float, help='radius of pores')
    parser.add_argument('-rows', default=2, type=int, help='Number of rows of pores. This will also be number of columns')
    # if you are making layers
    parser.add_argument('--layers', action="store_true", help='build a system with layers')
    parser.add_argument('-dbwl', default=0.37, type=float, help='Distance between layers (nm)')
    parser.add_argument('-lw', '--layer_width', default=0.1, type=float, help='width of layers')
    # if you want to make disks. You'll want to specify -dbwl, -r and -lw as well unless you want the defaults
    parser.add_argument('--disks', action="store_true", help='Create a system with disks stacked in layers in the z '
                        'direction. The disks are arranged around pore centers with their centers a distance specified'
                        'by the -r flag. The number of disks in each layer is specified with the -dpl flag')
    parser.add_argument('-dpl', '--disks_per_layer', default=5, help='Number of disks in each layer')
    parser.add_argument('-rdisk', '--disk_radius', default=0.1, help='Radius of disks')
    parser.add_argument('--offset', action="store_true", help='Offset disks with respect to layers above/below')
    parser.add_argument('-invert', '--invert', action="store_true", help='Fill in features and leave non-features empty')

    args = parser.parse_args()

    return args


def translate(pt, translation):
    """
    :param pt: 3D coordinates of points to translate
    :param translation: How far to translate the point in each direction
    :return: translated point
    """
    translation = np.array(translation)  # in case the array is not a numpy array already
    t = np.append(translation, [1])  # append a 1 at the end of the translation vector so we can do matrix multiplication
    T = np.zeros([4, 4])  # initialize translation matrix
    for i in range(4):
        T[i, i] = 1  # diagonals are all 1

    T[:3, 3] = pt  # last column contains the point to be translated

    return np.dot(T, t.T)[:3]


def Rz(pt, theta):
    """
    :param pt: 3D coordinates of point to be rotated
    :param theta: angle to rotate with respect z axis
    :return: rotated point
    """

    R = np.zeros([3, 3])
    R[2, 2] = 1
    R[0, 0] = np.cos(theta)
    R[0, 1] = -np.sin(theta)
    R[1, 0] = np.sin(theta)
    R[1, 1] = np.cos(theta)

    return np.dot(R, pt)


def check_pores(pt, pores, r):
    """
    :param pt: point which we are checking
    :param pores: coordinates of pore locations
    :param r: radius of pores
    :return: True/False - whether pt is contained in one of the pore regions
    """

    contained = False  # start off assuming the point is not contained in the pore region
    npores = pores.shape[0]

    for i in range(npores):
        radius = np.linalg.norm(pores[i] - pt[:2])  # xy distance between point and pore axis
        if radius <= r:
            if args.gaussian:
               # allow points inside the pores based on a gaussian probability
               # The full cdf sums to 1. We are using half the pdf so the max the cdf will get is 0.5. So multiply by 2
                probability = 2*scipy.stats.norm(r, r/2).cdf(radius)
                if np.random.rand() >= probability:  # generate random number between 0 and 1 as a test
                    contained = True  # Will only be excluded if it meets the condition
            else:
                contained = True
            break

    return contained


def check_layers(pt, layers, width):
    """
    :param pt: point which we are checking
    :param layers: coordinates of layer locations
    :param width: layer width
    :return: True/False - whether pt is contained in one of the pore regions
    """

    contained = False  # start off assuming the point is not contained in the layer region
    nlayers = layers.shape[0]
    r = 0.5*width

    for i in range(nlayers):
        radius = abs(pt[2] - layers[i])  # distance between point and center (based on z) of layer
        if radius <= r:
            if args.gaussian:
                # allow points inside the pores based on a gaussian probability
                # The full cdf sums to 1. We are using half the pdf so the max the cdf will get is 0.5. So multiply by 2
                probability = 2*scipy.stats.norm(r, r/2).cdf(radius)
                if np.random.rand() >= probability:  # generate random number between 0 and 1 as a test
                    contained = True  # Will only be excluded if it meets the condition
            else:
                contained = True
            break

    return contained


def check_disks(pt, disk_locations, layer_width, disk_radius):
    """
    :param pt: point which we are checking
    :param disk_locations: coordinates of disk center locations
    :param layer_width: width of disk in z direction. Same as width of layers
    :param disk_radius: radius of disk
    :return: True/False - whether pt is contained in one of the disk regions
    """

    contained = False  # start off assuming the point is not contained in the disk region
    ndisks = disk_locations.shape[0]

    for i in range(ndisks):
        radius = np.linalg.norm(disk_locations[i, :2] - pt[:2])  # xy distance between point and disk center
        width = abs(disk_locations[i, 2] - pt[2])  # z distance between point and disk center
        if radius <= disk_radius and width <= 0.5*layer_width:
            # if args.gaussian:
            #     # allow points inside the pores based on a gaussian probability
            #     # The full cdf sums to 1. We are using half the pdf so the max the cdf will get is 0.5. So multiply by 2
            #     probability = 2*scipy.stats.norm(r, r/2).cdf(radius)
            #     if np.random.rand() >= probability:  # generate random number between 0 and 1 as a test
            #         contained = True  # Will only be excluded if it meets the condition
            # else:
            contained = True
            break

    return contained


def write_gro_pos(pos, out, name='NA', box=[0, 0, 0], ids=[]):
    """
    write a .gro file from positions
    :param pos: xyz coordinates in
    :param out: name of output .gro file
    :param name: name to give atoms being put in the .gro
    :return: A .gro file
    """

    with open(out, 'w') as f:

        f.write('This is a .gro file\n')
        f.write('%s\n' % pos.shape[0])

        for i in range(pos.shape[0]):

            if not ids:
                f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format((i + 1) % 100000, '%s' % name, '%s' % name,
                                                                        (i + 1) % 100000, pos[i, 0], pos[i, 1], pos[i, 2]))
            else:
                f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format((i + 1) % 100000, '%s' % ids[i], '%s' % ids[i],
                                                                        (i + 1) % 100000, pos[i, 0], pos[i, 1], pos[i, 2]))
        for i in range(len(box)):
            f.write('{:10f}'.format(box[i]))

        f.write('\n')


if __name__ == "__main__":

    # load arguments
    args = initialize()

    # box dimensions, following crystallographic conventions
    a, b, c = args.box_lengths
    alpha, beta, gamma = [angle*(np.pi / 180) for angle in args.angles]  # convert input to radians

    # define number of points in each frame and number of total frames
    npts = args.npoints
    frames = args.nframes

    if args.layers or args.disks:

        nlayers = c / args.dbwl  # number of layers
        c = nlayers * args.dbwl  # recalculate z for consistency
        layer_locations = np.linspace(0, c, nlayers)  # centers of disks or layers

    # Define unit cell. See here: https://en.wikipedia.org/wiki/Fractional_coordinates#In_Crystallography for equations
    # Volume of unit cell
    V = a*b*c*np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2*np.cos(alpha)*np.cos(beta)*np.sin(gamma))
    # vectors defining unit cell
    O = np.array([0, 0, 0])  # origin (needed later when randomly placing points)
    A = np.array([a, 0, 0])  # vector in x direction
    B = np.array([b*np.cos(gamma), b*np.sin(gamma), 0])  # vector in y direction
    C = np.array([c*np.cos(beta), c*((np.cos(alpha) - np.cos(gamma)*np.cos(beta))/np.sin(gamma)), V / (a*b*np.sin(gamma))])  # vector in z direction

    if args.pores or args.disks:
        pore_radius = args.radius
        rows = args.rows  # number of rows of pores
        npores = rows**2  # assumes that x and y are the same dimension
        pore_locations = np.zeros([npores, 2])  # x, y positions of pore centers
        p2p = a / rows  # distance between pores

        for i in range(rows):
            for j in range(rows):
                pore_locations[i*rows + j, :] = [p2p*(0.5 + j) + p2p*(0.5 + i)*np.cos(gamma), p2p/2*np.sin(gamma) + i*p2p*np.sin(gamma)]

    if args.disks:
        # disks are placed with reference to pore centers. A single layer of disks consists of disks_per_layer disks
        # rotated about the z axis of each pore with angular spacing of 360 / disks_per_layer
        disks_per_layer = args.disks_per_layer
        nlayers = int(c / args.dbwl)  # number of layers based on z dimension and distance between layers
        pore_radius = args.radius  # This will say how far to place disks away from pore centers
        disk_radius = args.disk_radius
        disk_locations = np.zeros([npores*nlayers*disks_per_layer, 3])
        for i in range(npores):
            for j in range(nlayers):
                for k in range(disks_per_layer):
                    theta = 2*np.pi*k / disks_per_layer  # angle by which to rotate about pore axis
                    if args.offset:
                        theta += (j % 2) * (np.pi / disks_per_layer)
                    pt = np.array([pore_radius, 0, 0])  # create a point an x distance away from the origin
                    pt = Rz(pt, theta)  # rotate the point about the origin
                    pore = np.append(pore_locations[i, :], layer_locations[j])  # tack on the z component
                    disk_locations[i*nlayers*disks_per_layer + j*disks_per_layer + k, :] = translate(pt, pore)  # translate the point to the pore

    points = np.zeros([frames, npts, 3])  # will contain all positions for all frames

    if args.invert:
        # generate random points inside box
        for t in tqdm.tqdm(range(frames)):
            for i in tqdm.tqdm(range(npts)):
                u, v, w = np.random.rand(3)  # generate 3 random numbers between 0 and 1
                pt = O + u * A + v * B + w * C  # places point inside 3D box defined by box vector A, B and C
                # if --pores, --layers or --disks is specified, check to make the random point is not in the region. If it
                # is, keep generating random points until the point is outside the region
                if args.pores:
                    while not check_pores(pt, pore_locations, pore_radius):  # force all points outside the pore region
                        u, v, w = np.random.rand(3)
                        pt = O + u * A + v * B + w * C
                if args.layers:  # force all points outside layer region
                    while not check_layers(pt, layer_locations, args.layer_width):
                        u, v, w = np.random.rand(3)
                        pt = O + u * A + v * B + w * C
                if args.disks:
                    while not check_disks(pt, disk_locations, args.layer_width, disk_radius):
                        u, v, w = np.random.rand(3)
                        pt = O + u * A + v * B + w * C

                points[t, i, :] = pt

    else:
        # generate random points inside box
        for t in tqdm.tqdm(range(frames)):
            for i in tqdm.tqdm(range(npts)):
                u, v, w = np.random.rand(3)  # generate 3 random numbers between 0 and 1
                pt = O + u * A + v * B + w * C  # places point inside 3D box defined by box vector A, B and C
                # if --pores, --layers or --disks is specified, check to make the random point is not in the region. If it
                # is, keep generating random points until the point is outside the region
                if args.pores:
                    while check_pores(pt, pore_locations, pore_radius):  # force all points outside the pore region
                        u, v, w = np.random.rand(3)
                        pt = O + u * A + v * B + w * C
                if args.layers:  # force all points outside layer region
                    while check_layers(pt, layer_locations, args.layer_width):
                        u, v, w = np.random.rand(3)
                        pt = O + u * A + v * B + w * C
                if args.disks:
                    while check_disks(pt, disk_locations, args.layer_width, disk_radius):
                        u, v, w = np.random.rand(3)
                        pt = O + u * A + v * B + w * C

                points[t, i, :] = pt

    # Now write the trajectory in GROMACS format

    box = [A[0], B[1], C[2], A[1], A[2], B[0], B[2], C[0], C[1]]  # how its written in a .gro file

    unitcell_vectors = np.zeros([frames, 3, 3])
    for i in range(frames):
        # vectors don't change but need them as a trajectory
        unitcell_vectors[i, 0, :] = A
        unitcell_vectors[i, 1, :] = B
        unitcell_vectors[i, 2, :] = C

    write_gro_pos(points[-1, :, :], '%s.gro' % args.output, name=args.particle, box=box)  # write .gro of last frame
    traj = md.formats.TRRTrajectoryFile('%s.trr' % args.output, mode='w', force_overwrite=True)  # creat mdtraj TRR trajectory object
    time = np.linspace(0, 1000, frames)  # arbitrary times. Times are required by mdtraj
    traj.write(points, time=time, box=unitcell_vectors)  # write the trajectory in .trr format