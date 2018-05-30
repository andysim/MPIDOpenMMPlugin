from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import simtk.openmm.app.element as elem
import simtk.openmm.app.forcefield as forcefield
from sys import stdout

pdb = PDBFile('mpidwater.pdb')
forcefield = ForceField('mpidwater.xml')
system = forcefield.createSystem(pdb.topology, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
if pdb.topology.getPeriodicBoxVectors():
    simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Dump trajectory info every 1ps
simulation.reporters.append(DCDReporter('output.dcd', 500))
# Dump simulation info every 10fs
simulation.reporters.append(StateDataReporter(stdout, 5,
                            step=True, potentialEnergy=True, temperature=True))
# Run 100ps of simulation
simulation.step(50000)
