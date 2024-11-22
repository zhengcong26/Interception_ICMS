# -*- coding: utf-8 -*-
"""
@author: zheng

"""

#%% Figure 3a

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

# Load the MAT file data (Choose the appropriate file by uncommenting)
mat_data = loadmat(r'G_PSTH_ICMS.mat')
# mat_data = loadmat(r'L_PSTH_ICMS.mat')

selected_key = 'ans'
data_1 = mat_data[selected_key]

# Specify columns to delete (from 12th to 21st columns
columns_to_delete = list(range(12, 22))  

# Remove the specified columns from the data array
data = np.delete(data_1, columns_to_delete, axis=1)

# Get the shape of the data (rows and columns)
rows, cols = data.shape

# Create a figure with a custom size (height = 8, width = 5)
plt.figure(figsize=(5, 8)) 


plt.imshow(data, cmap='RdBu_r', interpolation='nearest', extent=[0, cols, 0, rows], aspect='auto')

# 隐藏坐标轴
plt.axis('off')

# 添加 colorbar 并设置 scale
cbar = plt.colorbar()
plt.clim(-10, 10) 
cbar.set_ticks([-10, 0, 10])  # 设置刻度位置

cbar.outline.set_visible(False)  # 隐藏 colorbar 的轮廓线

plt.savefig(r'ICMS_heatmap.png', dpi=800)

plt.show() 

#%% Figure 4a

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

# Load the data
mat_data = loadmat(r'G_unit_heatmap.mat')
# mat_data = loadmat(r'L_unit_heatmap.mat')

selected_key = 'ZZ'
data = mat_data[selected_key]

rows, cols = data.shape

# fig, ax = plt.subplots()
fig, ax = plt.subplots(figsize=(6,12))  # Adjust the figure size here

fig.dpi = 600

# Plot the heatmap
im = ax.imshow(data, cmap='RdBu_r', interpolation='nearest', extent=[0, cols, 0, 280], aspect=1)

plt.axis('off')

# Mark vertical lines
ax.axvline(x=11, color='blue', linestyle='-', linewidth=0.6) # TO
ax.axvline(x=61, color='green', linestyle='-', linewidth=0.6) # GO
ax.axvline(x=71, color='red', linestyle='-', linewidth=0.6) # MO
ax.axvline(x=86, color='black', linestyle='-', linewidth=0.6) # TA

# Set color limits
im.set_clim(-1.5, 1.5)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, orientation='vertical', fraction=0.05, pad=0.05)
cbar.set_ticks([-1, 0, 1])  # Set tick marks

# Adjust colorbar position
cbar.ax.set_position([0.83, 0.13, 0.6, 0.1])
cbar.outline.set_visible(False)  # Hide outline of the colorbar

# Save the figure
plt.savefig(r'units_heatmap.png', dpi=600)

# Show the plot
plt.show() 