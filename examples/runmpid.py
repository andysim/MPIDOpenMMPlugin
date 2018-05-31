from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem
import simtk.openmm.app.forcefield as forcefield
from sys import stdout, argv

pdb = PDBFile('mpidwater.pdb')
forcefield = ForceField('mpidwater.xml')
system = forcefield.createSystem(pdb.topology, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
context = simulation.context
# Set up platform properties
platformName = context.getPlatform().getName()
# Make sure there weren't any issues loading plugins
if context.getPlatform().getPluginLoadFailures():
    print("FAILED PLUGINS:", context.getPlatform().getPluginLoadFailures())
if platformName == 'CUDA':
    platform = Platform.getPlatformByName('CUDA')
    # Figure out which GPU to run on, i.e. did the user tell us?
    deviceid = argv[1] if len(argv) > 1 else '0'
    properties = {'DeviceIndex': deviceid, 'Precision': 'mixed'}
    platform.setPropertyValue(context, 'Precision', 'mixed')
    platform.setPropertyValue(context, 'DeviceIndex', deviceid)
    #simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    print("Running on ", context.getPlatform().getName(), " Device ID:", deviceid)
else:
    print("Running on the %s platform" % context.getPlatform().getName())
# Initialize
context.setPositions(pdb.positions)
if pdb.topology.getPeriodicBoxVectors():
    context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())


# Dump trajectory info every 1ps
simulation.reporters.append(DCDReporter('output.dcd', 500))
# Dump simulation info every 10fs
simulation.reporters.append(StateDataReporter(stdout, 5,
                            step=True, potentialEnergy=True, temperature=True))
# Run 100ps of simulation
simulation.step(50000)
