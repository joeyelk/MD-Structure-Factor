# MD-Structure-Factor
Calculates 3d structure factor and simulated x-ray diffraction pattern from molecular dynamics trajectories



example program usage:

main_gromacs.py -i foo    
(loads foo.gro and foo.trr)

main_gromacs.py -top foo.gro -traj foo.trr

The program will then create two intermediate files, out_foo_traj.npz and out_foo_sf.npz.  On subsequent runs of the program with the same inputs, these files will be automatically reloaded.  

Output:

The program will create a folder foo_plots, which will contain subfolders with various cross sections through the 3d structure factor.  Additionally, Ewald-corrected plots will be produced.  These more accurately show what would be seen in an x-ray diffraction experiment.  


