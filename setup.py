
from distutils.core import setup

setup(
    name='MD-Structure-Factor',
    version='1.0',
    packages=['Md-Structure-Factor',],
    license=open('LICENSE').read(),
    long_description=open('README.md').read(),
    install_requires=['wheel','setuptools','duecredit','tqdm','MDAnalysis','mdtraj','future','vtk','scipy','mayavi'],
        )
