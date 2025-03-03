# -*- coding: utf-8 -*-
"""
@author: zheng

"""

#%% Figure 3a

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

# Load the MAT file data (Choose the appropriate file by uncommenting)
# mat_data = loadmat(r'G_PSTH_ICMS.mat')
mat_data = loadmat(r'L_PSTH_ICMS.mat')

selected_key = 'ans'
data_1 = mat_data[selected_key]

# Specify columns to delete (from 12th to 21st columns
columns_to_delete = list(range(12, 22))  

# Remove the specified columns from the data array
data = np.delete(data_1, columns_to_delete, axis=1)

# Get the shape of the data (rows and columns)
rows, cols = data.shape

plt.figure(figsize=(3.5, 9.5)) 
plt.pcolormesh(data, cmap='RdBu_r', vmin=-10, vmax=10)  # shading='auto' to handle the spacing between cells

plt.gca().invert_yaxis()  # 反转y轴

# 隐藏坐标轴
plt.axis('off')

# cbar = plt.colorbar()
# plt.clim(-10, 10) 
# cbar.set_ticks([-10, 0, 10])  # 设置刻度位置
# cbar.outline.set_visible(False)  # 隐藏 colorbar 的轮廓线

plt.savefig(r'L_ICMS_heatmap.svg', dpi=600, format="svg")

plt.show() 

#%% Figure 4a

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

# Load the data
# mat_data = loadmat(r'G_unit_heatmap.mat')
mat_data = loadmat(r'L_unit_heatmap.mat')

selected_key = 'ZZ'
data = mat_data[selected_key]

rows, cols = data.shape

# 创建画布
fig, ax = plt.subplots(figsize=(3, 9.8))

# 画热图
pcmesh = ax.pcolormesh(data, cmap='RdBu_r', vmin=-1.5, vmax=1.5)

# 反转y轴，使得数据按正确方向显示
ax.invert_yaxis()

# 隐藏坐标轴
ax.axis('off')

# 添加垂直线
ax.axvline(x=11, color='blue', linestyle='-', linewidth=1)  # TO
ax.axvline(x=61, color='green', linestyle='-', linewidth=1)  # GO
ax.axvline(x=71, color='red', linestyle='-', linewidth=1)  # MO
ax.axvline(x=86, color='black', linestyle='-', linewidth=1)  # TA

# 设置 colorbar
# cbar = fig.colorbar(pcmesh, ax=ax, orientation='vertical', fraction=0.05, pad=0.05)
# cbar.set_ticks([-10, 0, 10])  # 设置刻度
# cbar.outline.set_visible(False)  # 隐藏 colorbar 轮廓

# 保存矢量图
plt.savefig(r'L_COINT_heatmap.svg', format='svg')

# 显示图像
plt.show()
