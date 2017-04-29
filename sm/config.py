# coding=utf-8
import simtk.unit as unit

# Basic input
INITIAL_STRUCTURE_FILENAME = 'initial_structure.pdb'
FORCEFIELD_FILE            = 'chr_ff.xml'                   # OpenMM forcefield XML file

# Restraints based on contact map
RESTRAINTS        = True
CONTACT_MAP_FILE = 'mySprings.rst'
SPRING_EQ_DIST   = 0.2
SPRING_CONST     = 300000.0

# External Field parameters
EXTERNAL_FIELD_FILE = None   # .npy filename or None
F_UNIT = unit.angstrom
F_XMIN = -3.5 * F_UNIT
F_XMAX =  3.5 * F_UNIT
F_YMIN = -3.0 * F_UNIT
F_YMAX =  3.0 * F_UNIT
F_ZMIN = -3.0 * F_UNIT
F_ZMAX =  3.0 * F_UNIT
F_SWAP_AXES = True                # Should be X swapped with Z? Usually True for fields from TIFF files
F_NORMALIZE = True                # Should the field be normalized to [0;1]?
F_SCALING_FACTOR = 500

# Energy minimization
MINIMIZE = True
MINIMIZED_FILE = "minimized.pdb"            # Set to None if you don't want to save minimized file

# integrator
INTEGRATOR_TYPE = 'langevin'                # Alternative: langevin, verlet
FRICTION_COEFF = 100 / unit.picosecond      # Used only with langevin integrator

# Simulation parameters
RUN_SIMULATION = True
N_STEPS   = 100000
TIME_STEP = 200 * unit.femtoseconds
TEMPERATURE = 309.75 * unit.kelvin    # 309.75K = 36.6°C
RANDOM_SEED = 0                       # 0 to make it random

# Trajectory settings
TRAJECTORY_FRAMES   = 10000
TRAJECTORY_FILENAME = 'trajectory.pdb'

# State reporting
N_OF_STATE_READS_REPORTED_TO_SCREEN = 20
N_OF_STATE_READS_REPORTED_TO_FILE   = 1000
STATE_FILE_NAME                     = 'state.csv'
PLOT_DATA                           = True
PLOT_FILE_NAME                      = 'energy.pdf'
