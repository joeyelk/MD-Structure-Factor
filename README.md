# MD-Structure-Factor
Calculates 3d structure factor and simulated x-ray diffraction pattern from molecular dynamics trajectories

## Dependencies

 * [scipy](https://www.scipy.org/)
 * [mayavi](https://docs.enthought.com/mayavi/mayavi/)
 * [tqdm](https://pypi.org/project/tqdm/)
 * [MDAnalysis](https://www.mdanalysis.org/)
 * [duecredit](https://github.com/duecredit/duecredit)
 * future
 * setuptools
 * wheel
 * vtk
 
```
python pip install setuptools wheel vtk scipy tqdm MDAnalsysis duecredit mayavi future
```

## Example

example program usage:
```
main_gromacs.py -i foo    
```
(loads foo.gro and foo.trr)

```
main_gromacs.py -top foo.gro -traj foo.trr
```

The program will then create two intermediate files, out_foo_traj.npz and out_foo_sf.npz.  On subsequent runs of the program with the same inputs, these files will be automatically reloaded.  

Output:

The program will create a folder foo_plots, which will contain subfolders with various cross sections through the 3d structure factor.  Additionally, Ewald-corrected plots will be produced.  These more accurately show what would be seen in an x-ray diffraction experiment.  


