import MDAnalysis
import numpy as np
import tqdm
#process trajectory and topology file and store in output file as np array
def process(topology_filename,trajectory_filename,output_filename):

	print "processing ",trajectory_filename
	u=MDAnalysis.Universe(topology_filename,trajectory_filename)

	dims=np.zeros((len(u.trajectory),3))
	coords=np.zeros((len(u.trajectory),u.trajectory.n_atoms,3))
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

	for ts in tqdm.tqdm(u.trajectory):
		dims[it,:]=ts.dimensions[0:3]
		coords[it,:,:]=ts.positions
		it+=1
	print "saving ",output_filename
	#print typ
	np.savez_compressed(output_filename,dims=dims,coords=coords,name=name,mass=mass,typ=typ,charge=charge)
	print 'done saving'

#output_filename is optional.  If it's not included, no npz trajectory file will be saved
def process_gro(topology_filename,trajectory_filename,output_filename=""):

	print "processing ",trajectory_filename
	u=MDAnalysis.Universe(topology_filename,trajectory_filename)

	dims=np.zeros((len(u.trajectory),3))
	coords=np.zeros((len(u.trajectory),u.trajectory.n_atoms,3))
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

	for ts in tqdm.tqdm(u.trajectory):
		dims[it,:]=ts.dimensions[0:3]
		coords[it,:,:]=ts.positions
		it+=1
	
	
	if len(output_filename)>0:
		print "saving ",output_filename	
		np.savez_compressed(output_filename,dims=dims,coords=coords,name=name,mass=mass,typ=typ)
		print 'done saving'
		
	return dims,coords,name,mass,typ





#process("W11_large.psf","DH5_423_0.dcd","o2")
