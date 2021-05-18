from ase import Atoms, units

from asap3.io.trajectory import Trajectory
from asap3.md.langevin import Langevin
from asap3.md.velocitydistribution import MaxwellBoltzmannDistribution
from asap3 import LennardJones
from ase.io import read, write
import io, os

def simulate_asap(atoms, temp, time_step, length, friction, fn = None, writeout=None):
    """
    Sets up and runs atoms cell using the Maxwell Boltzmann Distribution and
    Langevin Dynamics for the given inputs using the Asap-C++ optimized method. 
    Atoms object must have a calculator attached.
    
    Parameters
    Accepts:
    atoms          (Atoms Object) ASE atoms object. Must have calculator and 
                               periodic boundary conditions attached
    temp           (float) Temperature of simulation, in Kelvin
    time_step      (int) Time step, in femptoseconds
    length         (int) Number of iterations to run
    friction       (float) Friction coefficient. Usually between 0.01 - 0.001
    fn             (str) When not None, sets the file path to write the output to
    writeout       (int) When not None, Writes out positions every x iterations
    
    Returns:
    atoms          (Atoms Object) Returned atoms object will be the last frame
                               of the simulation
    
    """
    MaxwellBoltzmannDistribution(atoms, temperature_K=temp)
    dyn = Langevin(atoms, time_step*units.fs, temperature_K=temp,
                   friction=friction )
    if writeout is not None:
        if fn is not None:
            traj = Trajectory(fn, 'w', atoms)
            dyn.attach( traj.write, writeout)
            dyn.run(length)
    return atoms
    
    
def convert_traj(fn_read, fn_write):
    """
    Converts ASE trajectory file to .xyz format for analysis in other workflows,
    including VMD. Can be used for analysis in other python packages, including
    MDAnalysis.
    
    Parameters:
    Accepts:
    fn_read        (str) File name to read
    fn_write       (str) File name to write
    
    Returns:
    Nothing. Will write out a .xyz trajectory in filepath
    """

    traj = Trajectory(fn_read)
    string = 'structure'

    #get current working directory and make a scratch directory
    path = os.getcwd()
    path = path + '/scratch'
    if not os.path.exists(path): os.makedirs(path)

    #write each structure from the .traj file in .xyz format
    for i in range(len(traj)):
        atoms = traj[i]
        string = 'structure%03d' % (i,) +'.xyz'
        outStruct = os.path.join(path, string)
        write(outStruct, atoms)
    #combines all optimization structures in one trajectory file
        inFile = open(os.path.join(path, 'structure%03d' % 
                  (i,)  +'.xyz'), 'r')
        fileStr = inFile.read()
        outFile = open(fn_write, 'a')
        outFile.write(fileStr)
        
    