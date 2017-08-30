import MDAnalysis
import numpy as np
import tqdm
#process trajectory and topology file and store in output file as np array
def process(topology_filename,trajectory_filename,output_filename,stride):

	print "processing ",trajectory_filename
	u=MDAnalysis.Universe(topology_filename,trajectory_filename)

	dims=np.zeros((len(u.trajectory)/stride,3))
	coords=np.zeros((len(u.trajectory)/stride,u.trajectory.n_atoms,3))
	name=np.zeros(u.trajectory.n_atoms,dtype=object)
	mass=np.zeros(u.trajectory.n_atoms)
	typ=np.zeros(u.trajectory.n_atoms,dtype=object)
	charge=np.zeros(u.trajectory.n_atoms)

	im=0
	for ia in tqdm.tqdm(u.atoms):
		name[im]  =u.atoms[im].name
		mass[im]  =u.atoms[im].mass
		charge[im]=u.atoms[im].charge
		typ[im] =u.atoms[im].type
		im+=1
	it=0

	for icount,ts in tqdm.tqdm(enumerate(u.trajectory)):
		if icount%stride==0:
			dims[it,:]=ts.dimensions[0:3]
			coords[it,:,:]=ts.positions
			it+=1
	print "saving ",output_filename
	#print typ
	np.savez_compressed(output_filename,dims=dims,coords=coords,name=name,mass=mass,typ=typ,charge=charge)
	print 'done saving'


def process_gro(topology_filename,trajectory_filename,output_filename,stride):

	print "processing ",trajectory_filename
	u=MDAnalysis.Universe(topology_filename,trajectory_filename)

	dims=np.zeros((len(u.trajectory)/stride,3))
	coords=np.zeros((len(u.trajectory)/stride,u.trajectory.n_atoms,3))
	name=np.zeros(u.trajectory.n_atoms,dtype=object)
	mass=np.zeros(u.trajectory.n_atoms)
	typ=np.zeros(u.trajectory.n_atoms,dtype=object)


	im=0
	for ia in tqdm.tqdm(u.atoms):
		name[im]  =u.atoms[im].name
		mass[im]  =u.atoms[im].mass
		typ[im] =u.atoms[im].type
		im+=1
	it=0

	for icount,ts in tqdm.tqdm(enumerate(u.trajectory)):
		if icount%stride==0:
			dims[it,:]=ts.dimensions[0:3]
			coords[it,:,:]=ts.positions
			it+=1
	print "saving ",output_filename
	#print typ
	np.savez_compressed(output_filename,dims=dims,coords=coords,name=name,mass=mass,typ=typ)
	print 'done saving'





#process("W11_large.psf","DH5_423_0.dcd","o2")
