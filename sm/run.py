#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import importlib
import os
import sys

import numpy as np
import simtk.openmm as mm
import simtk.unit as u
from simtk.openmm.app import PDBFile, ForceField, Simulation, PDBReporter, StateDataReporter

from utils import sizeof_fmt, plot_data

parser = argparse.ArgumentParser()
parser.add_argument('config', help="config file (python module)")
args = parser.parse_args()

settings_file = args.config
if settings_file.endswith('.py'):
    settings_file = settings_file[:-3]
if '.' in settings_file:
    sys.exit('settings file must not contain dot in its file name (except .py extension)')
cfg = importlib.import_module(settings_file)

print("Initial setup...")
total_time = cfg.N_STEPS * cfg.TIME_STEP
if cfg.RANDOM_SEED == 0:
    random_seed = np.random.randint(2147483647)
else:
    random_seed = cfg.RANDOM_SEED
reporting_to_screen_freq = max(1, int(round(cfg.N_STEPS / cfg.N_OF_STATE_READS_REPORTED_TO_SCREEN)))
reporting_to_file_freq = max(1, int(round(cfg.N_STEPS / cfg.N_OF_STATE_READS_REPORTED_TO_FILE)))
trajectory_freq = max(1, int(round(cfg.N_STEPS / cfg.TRAJECTORY_FRAMES)))

print("   OpenMM version:                  {}".format(mm.__version__))
print("   Number of steps:                 {} steps".format(cfg.N_STEPS))
print("   Time step:                       {}".format(cfg.TIME_STEP))
print("   Temperature:                     {}".format(cfg.TEMPERATURE))
print("   Total simulation time:           {}".format(total_time.in_units_of(u.nanoseconds)))
print("   Number of state reads:           {} reads".format(cfg.N_OF_STATE_READS_REPORTED_TO_SCREEN))
print("   State reporting to screen every: {} step".format(reporting_to_screen_freq))
print("   State reporting to file every:   {} step".format(reporting_to_file_freq))
print("   Number of trajectory frames:     {} frames".format(cfg.TRAJECTORY_FRAMES))
print("   Trajectory frame every:          {} step".format(trajectory_freq))
print('   Random seed:', random_seed)
print()

print("Loading initial structure:\n  {}".format(cfg.INITIAL_STRUCTURE_FILENAME))
pdb = PDBFile(cfg.INITIAL_STRUCTURE_FILENAME)
print("Loading forcefield file:\n  {}".format(cfg.FORCEFIELD_FILE))
forcefield = ForceField(cfg.FORCEFIELD_FILE)
print("Building system...")
system = forcefield.createSystem(pdb.topology)

if cfg.EXTERNAL_FIELD_FILE:
    print('Loading external field...')
    size = os.stat(cfg.EXTERNAL_FIELD_FILE).st_size
    print("   Reading {} file ({})...".format(cfg.EXTERNAL_FIELD_FILE, sizeof_fmt(size)))
    img = np.load(cfg.EXTERNAL_FIELD_FILE)
    print("   Array of shape {} loaded.".format(img.shape))
    print("   Number of values: {}".format(img.size))
    print("   Min: {}".format(np.min(img)))
    print("   Max: {}".format(np.max(img)))
    if cfg.F_SWAP_AXES:
        print('   [INFO] X and Z axis will be swapped! ', end='')
        img = img.swapaxes(0, 2)
        print("Shape after swapping: ", img.shape)
    if cfg.F_NORMALIZE:
        print('   [INFO] Field will be normalized to [0, 1]')
        a = np.min(img)
        b = np.max(img)
        img = (img + abs(a)) / abs(b-a)
    print("  Creating a force based on density...")
    shape = img.shape
    img = img.flatten().astype(np.float64)
    field_function = mm.Continuous3DFunction(shape[2], shape[1], shape[0], img,
                                             cfg.F_XMIN, cfg.F_XMAX, cfg.F_YMIN, cfg.F_YMAX, cfg.F_ZMIN, cfg.F_ZMAX)
    field_force = mm.CustomCompoundBondForce(1, 'ksi*fi(x1,y1,z1)')
    field_force.addTabulatedFunction('fi', field_function)
    field_force.addGlobalParameter('ksi', cfg.F_SCALING_FACTOR)
    print("  Adding force to the system...")
    for i in range(system.getNumParticles()):
        field_force.addBond([i], [])
    system.addForce(field_force)
else:
    print("External force (density) file was not provided.")

if cfg.RESTRAINTS:
    print("  Adding restraints...")
    flat_bottom_force = mm.CustomBondForce('step(r-r0) * (k/2) * (r-r0)^2')
    flat_bottom_force.addPerBondParameter('r0')
    flat_bottom_force.addPerBondParameter('k')
    system.addForce(flat_bottom_force)
    with open(cfg.CONTACT_MAP_FILE) as input_file:
        for line in input_file:
            columns = line.split()
            atom_index_i = int(columns[0][1:]) - 1
            atom_index_j = int(columns[1][1:]) - 1
            r0 = cfg.SPRING_EQ_DIST
            k =  cfg.SPRING_CONST
            flat_bottom_force.addBond(atom_index_i, atom_index_j, [r0, k])
            print('    {} - {}'.format(atom_index_i, atom_index_j))

print("Integrator initialization...")
if cfg.INTEGRATOR_TYPE == "langevin":
    integrator = mm.LangevinIntegrator(cfg.TEMPERATURE, cfg.FRICTION_COEFF, cfg.TIME_STEP)
    integrator.setRandomNumberSeed(random_seed)
elif cfg.INTEGRATOR_TYPE == "verlet":
    integrator = mm.VerletIntegrator(cfg.TIME_STEP)
else:
    sys.exit("Integrator initialization error!")

print("Setting up simulation...")
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(cfg.TEMPERATURE)
if cfg.MINIMIZE:
    print('  Energy minimizing...')
    simulation.minimizeEnergy()
    if cfg.MINIMIZED_FILE:
        print('  Saving minimized structure in {}'.format(cfg.MINIMIZED_FILE))
        state = simulation.context.getState(getPositions=True)
        PDBFile.writeFile(pdb.topology, state.getPositions(), open(cfg.MINIMIZED_FILE, 'w'))

if cfg.RUN_SIMULATION:
    print('Setting up reporters...')
    simulation.reporters.append(PDBReporter(cfg.TRAJECTORY_FILENAME, trajectory_freq))
    simulation.reporters.append(StateDataReporter(sys.stdout, reporting_to_screen_freq,
                                                  step=True, potentialEnergy=True, kineticEnergy=False,
                                                  totalEnergy=False, temperature=True))
    simulation.reporters.append(StateDataReporter(cfg.STATE_FILE_NAME, reporting_to_file_freq,
                                                  step=True, potentialEnergy=True, kineticEnergy=False,
                                                  totalEnergy=False, temperature=True))
    print('Running simulation...')
    simulation.step(cfg.N_STEPS)
    if cfg.PLOT_DATA:
        plot_data(cfg.STATE_FILE_NAME, cfg.PLOT_FILE_NAME)
print()
print("Everything is done")

