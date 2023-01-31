# CellModeller script for two bacterial colonies dividing in a microchannel
# To achieve 1200 min of real time, 240,000 simulation steps must be made
# with 0.005 min timestep

import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
import numpy
import math

max_cells = 10**3

# Cell parameters
# With a fixed simualtion time step dt = 0.005 min,
# real growth rate corresponds to 0.035 um/min
# which corresponds to a doubling time of 20 min
dt = 0.005 # [min]
base_growthRate = 0.035 # [um/min^{-1}]
growthRate_scale = 20

# Imaging frequency in real time
imaging_freq = 3 # [min]
pickle_freq = math.floor((imaging_freq/dt)/growthRate_scale)

# Microchannel dimensions
chan_wid = 1.5 # [um]
chan_len = 30.0 # [um]

# Uniform variation of cell parameters during division
growthRate_var = 0.1
length_var = 0.1

cell_cols = {0:[0,1.0,0], 1:[1.0,0,0]} # RGB cell colours ... green and red
cell_lens = {0:2.0, 1:2.0} # target cell lengths ... 2 um
cell_growr = {0:growthRate_scale * base_growthRate, 
    1:growthRate_scale * base_growthRate}

def get_x_coord(cell):
    # Finds x-coordinate of a cell
    return cell.pos[0]

def apply_variation(x, variation=0.1):
    """ Samples a number from a uniform realtive variation 
        defined by the passed number.
    """
    return random.uniform((1-variation) * x, (1+variation) * x)

def setup(sim):
    # Set biophysics and regulation models
    biophys = CLBacterium(sim,
        gamma=10, 
        jitter_z=False, 
        max_planes=2)

    # add the planes to set physical boundaries of cell growth
    # 4 micrometer width
    biophys.addPlane((0, chan_wid/2, 0), (0,-1,0), 1.0)
    biophys.addPlane((0, -chan_wid/2, 0), (0,1,0), 1.0)

    regul = ModuleRegulator(sim, sim.moduleName)

    # Only biophys and regulation
    sim.init(biophys, regul, None, None)

    # Specify the initial cells and their location in the simulation
    # For each simulation put inital cells in the same position
    # 5 um away from the center of the channel
    sim.addCell(cellType=0, pos=(5,0,0), dir=(-1,0,0), len = cell_lens[0]) 
    sim.addCell(cellType=1, pos=(-5,0,0), dir=(1,0,0), len = cell_lens[1])

    sim.pickleSteps = pickle_freq

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetLen = 2 * apply_variation(cell_lens[cell.cellType], 
        variation=length_var)
    cell.length = apply_variation(cell_lens[cell.cellType], 
        variation=length_var)

    # Specify growth rate of cells
    cell.growthRate = apply_variation(cell_growr[cell.cellType], 
        variation=growthRate_var)

    # Specify colour of cells
    cell.color = cell_cols[cell.cellType]

    # Initialize killFlag
    cell.killFlag = False

def update(cells):
    # Iterate through each cell and flag cells that reach target size for division
    # or remove cells from the simualtion that are out of the channel
    for (id, cell) in cells.items():
        if cell.length > cell.targetLen:
            cell.divideFlag = True
        elif abs(get_x_coord(cell)) > (chan_len / 2):
            cell.killFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetLen = 2 * apply_variation(cell_lens[parent.cellType], 
        variation=length_var)
    d2.targetLen = 2 * apply_variation(cell_lens[parent.cellType], 
        variation=length_var)

    # Specify daugther growth rate
    d1.growthRate = apply_variation(cell_growr[parent.cellType], 
        variation=growthRate_var)
    d2.growthRate = apply_variation(cell_growr[parent.cellType], 
        variation=growthRate_var)
