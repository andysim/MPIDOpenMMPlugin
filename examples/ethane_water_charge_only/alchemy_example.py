from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem
from sys import stdout, argv
import mpidplugin
import numpy as np
from mdtraj.reporters import NetCDFReporter

pdb = PDBFile('solvated_ethane_from_openmm_setup.pdb')
forcefield = ForceField('ethane_water.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, polarization='direct',
                                 nonbondedCutoff=8*angstrom, constraints=HBonds, defaultTholeWidth=8)
temperature = 298*kelvin
system.addForce(MonteCarloBarostat(1*atmosphere, temperature, 25))
integrator = LangevinIntegrator(temperature, 1/picosecond, 2*femtoseconds)

try:
    myplatform = Platform.getPlatformByName('CUDA')
    # Figure out which GPU to run on, i.e. did the user tell us?
    deviceid = argv[1] if len(argv) > 1 else '0'
    #myproperties = {'DeviceIndex': deviceid, 'Precision': 'double'}
    myproperties = {'DeviceIndex': deviceid, 'Precision': 'mixed'}
except:
    print("CUDA NOT FOUND!!!!!!!!!!")
    myplatform = None
    deviceid = "N/A"

if myplatform:
    simulation = Simulation(pdb.topology, system, integrator, myplatform, myproperties)
else:
    simulation = Simulation(pdb.topology, system, integrator)

context = simulation.context
if pdb.topology.getPeriodicBoxVectors():
    context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

print("Running on ", context.getPlatform().getName(), " Device ID:", deviceid)

# Initialize
context.setPositions(pdb.positions)
#simulation.loadState('equilibrated.xml')

# Get a list of the forces.  We use a dictionary here because there is only one force present
# that identifies itself as just a "Force".  If multiple plugins are used, this needs to change.
forces = {system.getForce(i).__class__.__name__: system.getForce(i) for i in range(system.getNumForces())}
if mpidplugin.MPIDForce.isinstance(forces['Force']):
    mpidforce = mpidplugin.MPIDForce.cast(forces['Force'])
else:
    raise "Unknown custom force detected!"

print(f'initial energy = ',context.getState(getEnergy=True).getPotentialEnergy())
# Turn on the electrostatic terms in increments, after storing the initial values
lambdas = np.arange(0.0, 1.05, 0.1)
selection = range(8)  # We know the first 8 atoms are ethane
elec_params = []
for i in selection:
    #charge, dipole, qpole, opole, axis, atomZ, atomX, atomY, thole, alphas
    elec_params.append(mpidforce.getMultipoleParameters(i))

for lam in lambdas:
    for i in selection:
        params = elec_params[i].copy()
        # charges
        params[0] *= lam
        # polarizabilities
        params[9] = [a * lam for a in params[9]]
        mpidforce.setMultipoleParameters(i, *params)
    mpidforce.updateParametersInContext(context)
    print(f'elec lambda = {lam:.2f}, energy = ',context.getState(getEnergy=True).getPotentialEnergy())

# Turn on the VDW terms in increments, after storing the initial values
vdw_params = []
for i in selection:
    #charge, sigma, epsilon
    vdw_params.append(forces['NonbondedForce'].getParticleParameters(i))

for lam in lambdas:
    for i in selection:
        params = [lam*p for p in vdw_params[i]]
        forces['NonbondedForce'].setParticleParameters(i, *params)
    forces['NonbondedForce'].updateParametersInContext(context)
    print(f'vdw lambda = {lam:.2f}, energy = ',context.getState(getEnergy=True).getPotentialEnergy())
print(f'final energy = ',context.getState(getEnergy=True).getPotentialEnergy())
