# CellModeller script for two bacterail collonies dividing in a channel
import random
from CellModeller.Regulation.ModuleRegulator import ModuleRegulator
from CellModeller.Biophysics.BacterialModels.CLBacterium import CLBacterium
import numpy
import math

max_cells = 10**3

cell_cols = {0:[0,1.0,0], 1:[1.0,0,0]} #RGB cell colours
cell_lens = {0:2.5, 1:2.5} #target cell lengths
cell_growr = {0:1.0, 1:1.0} #growth rates

def setup(sim):
    # Set biophysics, signalling, and regulation models
    biophys = CLBacterium(sim, 
        jitter_z=False, 
        max_planes=2)

    # add the planes to set physical boundaries of cell growth
    biophys.addPlane((0,-2,0), (0,1,0), 1)
    biophys.addPlane((0,2,0), (0,-1,0), 1)

    regul = ModuleRegulator(sim)
    # Only biophys and regulation
    sim.init(biophys, regul, None, None)


    # Specify the initial cell and its location in the simulation
    sim.addCell(cellType=0, pos=(1,0,0), dir=(1,0,0)) 
    sim.addCell(cellType=1, pos=(-1,0,0), dir=(1,0,0))

    if sim.is_gui:
        # Add some objects to draw the models
        from CellModeller.GUI import Renderers
        therenderer = Renderers.GLBacteriumRenderer(sim)
        sim.addRenderer(therenderer)

    sim.pickleSteps = 10

def init(cell):
    # Specify mean and distribution of initial cell size
    cell.targetVol = cell_lens[cell.cellType] + random.uniform(0.0,0.5)
    # Specify growth rate of cells
    cell.growthRate = cell_growr[cell.cellType]
    
    cell.color = cell_cols[cell.cellType]

def update(cells):
    #Iterate through each cell and flag cells that reach target size for division
    for (id, cell) in cells.items():
        cell.color = [cell.cellType*0.6+0.1, 1.0-cell.cellType*0.6, 0.3]
        if cell.volume > cell.targetVol:
            cell.divideFlag = True

def divide(parent, d1, d2):
    # Specify target cell size that triggers cell division
    d1.targetVol = 2.5 + random.uniform(0.0,0.5)
    d2.targetVol = 2.5 + random.uniform(0.0,0.5)

