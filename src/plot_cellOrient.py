# %% [markdown]
# # Analysis of bacterial cell oritetation in microchannels
# Analyse cell oritentation distribution in microchannels simulated with `CellModeller` [1].
#
# ## Purpose
# Compare simualtion and experimental cell oritentation data to calibrate simulation parameters.
#
# ## Methodology
# 1. Read `.pickle` file into dictionary with `pickle` module.
# 2. Plot cell oritentation histogram for each simulation step.
# 3. Coallete all histograms into one interactive figure.
#
# ## WIP - improvements
# Use this section only if the notebook is not final.
#
# ## Notable TODOs:
#
# - Load multiple `.pickle` files. Store in a dictionary `{simulation_step: cell_state}`
# - Plot oritentation histogram for each simulation step.
# - Use `plotly` to visualize multiple steps in one figure.
# - Copy model file to output directory.
#
# ## Results
# `norm` - cell orientation vector. **Note** Already normalized.
#
# ## Suggested next steps
# - Load, process and analyse experimetal data
# - Compare experimetal and simualtion data plots

# %% [markdown]
# # Setup
# ## Library import

# %%
import plotly.express as px
import os
from os.path import join as join_paths, basename as get_basename
import glob
import sys

import pickle
import numpy as np
from scipy.stats import gaussian_kde

from matplotlib import pyplot as plt
# %matplotlib inline

if len(sys.argv) != 2:
    print("Directory path with pickle files must be passed")
    print("Example: pyhton src/plot_cellOrient.py data/microch_30x4-22-11-29-09-57")
    exit
elif not os.path.isdir(sys.argv[1]):
    print("Check if you enetered a right data dir path")
    exit    


# %%
os.chdir("/home/i/igors-dubanevics/projects/bacteria-microchannel")
os.getcwd()

# %% [markdown]
# ## Parameter definition

# %%
data_dir_path = sys.argv[1]

out_dir_path = join_paths("scratch", get_basename(data_dir_path))

# Time-step
dt = 0.005  # [min]

# %% [markdown]
# ## Pre-data import

# %%
# Create output dir
os.makedirs(out_dir_path, exist_ok=True)

# %% [markdown]
# ## Data import

# %%
# Store data in a dictionary of form {step_number: cell_state}
cell_states = {}

# Get data from saved pickle files
filenames = sorted(glob.glob(join_paths(data_dir_path, 'step-*.pickle')))
for filename in filenames:
    # Extract simualtion step number from a filename
    step_num = int(get_basename(filename).replace(
        ".pickle", "").replace("step-", ""))
    data = pickle.load(open(filename, 'rb'))
    cs = data['cellStates']
    cell_states[step_num] = cs

# %% [markdown]
# ## Data processing

# %%
lengths = {}
postions = {}
orients = {}

orient_per_cell = {}

for step_num, cell_state in cell_states.items():
    # Make some convenient data arrays from cell properties
    length = np.array([cell.length for (id, cell)
                      in cell_state.items()])  # [um]
    # Postion in respect to the center of simulation
    # Assume it's center of the channel (0,0)
    pos = np.array([cell.pos for (id, cell) in cell_state.items()])
    # Normalized orientation vector
    orient = np.array([cell.dir for (id, cell) in cell_state.items()])

    # Orientation per cell
    for (id, cell) in cell_state.items():
        if id in orient_per_cell.keys():
            orient_per_cell[id] = np.append(
                orient_per_cell[id], [cell.dir], axis=0)
        else:
            orient_per_cell[id] = np.array([cell.dir])

    lengths[step_num] = length
    postions[step_num] = pos
    orients[step_num] = orient

# %%
orients[0][0]

# %%
# Calculate angle in respect to a channel
angles_per_step = {}
# Unit-vector parallel to the channel (pointing left-to-right)
y_norm = np.array((0, 1, 0), dtype=float)
for step_num, orient in orients.items():
    # Dot product with a unit-vector parallel to channel
    inner_prods = np.inner(y_norm, orient)
    # Dot product to radians
    rads = np.arccos(inner_prods)
    # Convert radians to degrees
    angles_per_step[step_num] = np.rad2deg(rads)

# %%
angles = np.array([angles for (step_num, angles)
                  in angles_per_step.items()], dtype=object)
angles = np.concatenate(angles).reshape(-1)
angles = angles - 90

# %%
# Calculate angle in respect to a channel
angles_per_cell = {}

# Unit-vector parallel to the channel (pointing left-to-right)
y_norm = np.array((0, 1, 0), dtype=float)

for id, orient in orient_per_cell.items():
    # Dot product with a unit-vector parallel to channel
    inner_prods = np.inner(y_norm, orient)
    # Dot product to radians
    rads = np.arccos(inner_prods)
    # Convert radians to degrees
    angles_per_cell[id] = np.rad2deg(rads)

# %%
# Average angle per cell for the whole simulation time
avg_angles = np.array([np.average(angles)
                      for (id, angles) in angles_per_cell.items()])
avg_angles = avg_angles - 90

# %%
# Plot histograms of some cell properties

# Orientation
fig, ax = plt.subplots(figsize=(5, 3.5))
plt.title("filename: {}".format(get_basename(data_dir_path)), loc='left')
bin_num = 90
ax.hist(avg_angles, bin_num, fc='gray', histtype='stepfilled', alpha=0.3,
        density=True, label='COMP ({:d} bins)'.format(bin_num))
# ax.plot(x_grid, angles_kde[key](x_grid), linewidth=2, label='KDE')

plt.xlabel(r"cell orientation ($\arccos{\hat{r}\cdot\hat{x}}$) [deg]")
plt.ylabel("Density (N={:d})".format(len(avg_angles)))
plt.tight_layout()

# Convert simualtion time into real time
# plt.text(.1, .99, "REAL time: {:.1f} min".format(
    # 500), ha='left', va='top', transform=ax.transAxes)

plt.legend(loc='upper right', frameon=False)

# Save figure
for ext in ['png']:
    fig_outname = "orient.{ext}".format(ext=ext)
    fig.savefig(join_paths(out_dir_path, fig_outname), dpi=300)



# %% [markdown]
# ## References
# We report here relevant references:
# 1. https://github.com/cellmodeller/CellModeller
# 2.
