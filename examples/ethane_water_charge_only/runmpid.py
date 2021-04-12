from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem
from sys import stdout, argv
import mpidplugin
import numpy as np

pdb = PDBFile('waterbox_31ang.pdb')
forcefield = ForceField('mpidwater.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=LJPME, nonbondedCutoff=8*angstrom, constraints=HBonds, defaultTholeWidth=8)
temperature = 298*kelvin
system.addForce(MonteCarloBarostat(1*atmosphere, temperature, 25))
integrator = LangevinIntegrator(temperature, 10000/picosecond, 0.002*femtoseconds)
#integrator = VerletIntegrator(1*femtoseconds)

# Make sure all the forces we expect are present
for force in range(system.getNumForces()):
    print(system.getForce(force))

try:
    myplatform = Platform.getPlatformByName('CUDA')
    # Figure out which GPU to run on, i.e. did the user tell us?
    deviceid = argv[1] if len(argv) > 1 else '0'
    myproperties = {'DeviceIndex': deviceid, 'Precision': 'double'}
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
simulation.loadState('equilibrated.xml')

# Dump trajectory info every 10ps
#simulation.reporters.append(DCDReporter('output.dcd', 5000))
# Dump simulation info every 1ps
nsteps = 1000000
simulation.reporters.append(StateDataReporter(stdout, 500, totalSteps=nsteps, speed=True, volume=True, density=True,
                            step=True, potentialEnergy=True, totalEnergy=True, temperature=True))
simulation.reporters.append(PDBReporter('equilibrated.pdb', nsteps))
# Run 100ps of simulation
simulation.step(nsteps)
simulation.saveState('equilibrated.xml')
