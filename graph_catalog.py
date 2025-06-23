import sys
from galaxy_catalog import GalaxyCatalog
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.preprocessing import PowerTransformer

import scienceplots


sys.path.append("../")
import os
os.chdir("/global/homes/d/dkololgi/TNG/Illustris/")
# from TNG.Illustris.Network_stats import network
from Network_stats import network
# from TNG.Illustris import Utilities

import torch
from torch_geometric.data import Data
from torch_geometric.utils import from_networkx

DESI_NETWORK = network(
    masscut=9.,
    from_DESI=True
)

print('DESI network object created')
G=DESI_NETWORK.subhalo_delauany_network(xyzplot=False)
print('DESI delaunay graph created')
DESI_NETWORK.network_stats_delaunay()
print('DESI delaunay network stats calculated')
DESI_geom = from_networkx(G, group_edge_attrs='all')
print('DESI delaunay graph converted to torch_geometric Data object')
scaler = PowerTransformer(method = 'box-cox')
DESI_features = pd.DataFrame(scaler.fit_transform(DESI_NETWORK.data), index=DESI_NETWORK.data.index, columns=DESI_NETWORK.data.columns)
DESI_geom.x = torch.tensor(DESI_features.values, dtype=torch.float32)
print('DESI features scaled and converted to torch tensor')
# Declare model for inference
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch_geometric.nn import GATv2Conv
class SimpleGAT(nn.Module):
    '''simple GAT model using the GATLayer
    
    Args:
        input_dim (int): Dimensionality of the input feature vectors
        output_dim (int): Dimensionality of the output softmax distribution
        num_heads (int): Number of attention heads
    
    '''
    def __init__(self, input_dim, output_dim, num_heads = 1): # keep a default of 1 head
        super(SimpleGAT, self).__init__()

        hidden_dim = 20  # per-head output size
        total_hidden = hidden_dim * num_heads  # total output size if concat=True


        # self.gat_layer1 = GATv2Conv(input_dim, 15, edge_dim=1, heads=num_heads, concat=True)
        # self.gat_layer2 = GATv2Conv(15, 15, edge_dim=1, heads=num_heads, concat=True)
        # self.gat_layer3 = GATv2Conv(15, output_dim, edge_dim=1, heads=num_heads, concat=False)
        # # self.gat_layer4 = GATv2Conv(10, output_dim)

        self.gat_layer1 = GATv2Conv(input_dim, hidden_dim, edge_dim=1, heads=num_heads, concat=True)
        self.gat_layer2 = GATv2Conv(total_hidden, hidden_dim, edge_dim=1, heads=num_heads, concat=True)
        self.gat_layer3 = GATv2Conv(total_hidden, hidden_dim, edge_dim=1, heads=num_heads, concat=True)
        self.gat_layer4 = GATv2Conv(hidden_dim * num_heads, output_dim, edge_dim=1, heads=1, concat=False)

        self.gat_dropout = nn.Dropout(p=0.1)
        self.norm1 = nn.LayerNorm(total_hidden)
        self.norm2 = nn.LayerNorm(total_hidden)
        self.norm3 = nn.LayerNorm(total_hidden)


    def forward(self, x, edge_index, edge_weight=None, return_attention=False):
        """Forward pass of the GNN model"""

        if return_attention:
            x1, (ei1, attn1) = self.gat_layer1(x, edge_index, edge_attr=edge_weight, return_attention_weights = True)
            x1 = self.gat_dropout(F.relu(self.norm1(x1)))

            x2, (ei2, attn2) = self.gat_layer2(x1, edge_index, edge_attr=edge_weight, return_attention_weights = True)
            x2 = self.gat_dropout(F.relu(self.norm2(x2) + x1))

            x3, (ei3, attn3) = self.gat_layer3(x2, edge_index, edge_attr=edge_weight, return_attention_weights = True)
            x3 = self.gat_dropout(F.relu(self.norm3(x3) + x2))

            y_hat, (ei4, attn4) = self.gat_layer4(x3, edge_index, edge_attr=edge_weight, return_attention_weights = True)


            # y_hat, (ei3, attn3) = self.gat_layer3(x2, edge_index, edge_attr=edge_weight, return_attention_weights = True)

            return y_hat, [(ei1, attn1), (ei2, attn2), (ei3, attn3), (ei4, attn4)]
        
        else:

            x1 = self.gat_dropout(F.relu(self.norm1(self.gat_layer1(x, edge_index, edge_attr=edge_weight))))
            x2 = self.gat_dropout(F.relu(self.norm2(self.gat_layer2(x1, edge_index, edge_attr=edge_weight)) + x1))
            x3 = self.gat_dropout(F.relu(self.norm3(self.gat_layer3(x2, edge_index, edge_attr=edge_weight)) + x2))
            y_hat = self.gat_layer4(x3, edge_index, edge_attr=edge_weight)
        # # x = self.gat_dropout(x)
        # x = F.relu(self.gat_layer1(x, edge_index, edge_attr=edge_weight))
        # # x = self.gat_dropout(x)
        # x = F.relu(self.gat_layer2(x, edge_index, edge_attr=edge_weight))

        # x = F.relu(self.gat_layer3(x, edge_index, edge_attr=edge_weight))

        # y_hat = self.gat_layer4(x, edge_index, edge_attr=edge_weight)
        # # y_hat = self.gat_layer3(x, edge_index, edge_attr=edge_weight)
        # # y_hat = self.gat_layer4(x, edge_index, edge_weight=edge_weight)
        return y_hat

# Load model for inference
inference_model = SimpleGAT(
    10, 4, num_heads=4
)
inference_model.load_state_dict(
    torch.load("/global/homes/d/dkololgi/TNG/Illustris/trained_gat_simulation.pth", map_location='cpu')
)

inference_model.eval()
print('Inference model loaded and set to eval mode')
with torch.no_grad():
    DESI_out = inference_model(
        DESI_geom.x, 
        DESI_geom.edge_index, 
        edge_weight=DESI_geom.edge_attr
    )
    DESI_pred = DESI_out.argmax(dim=1).numpy()
    DESI_probs = DESI_out.numpy()
print('Inference completed, predictions and probabilities obtained')
# Define environment labels and custom palette
environ_dicts = {
    0: 'Void',
    1: 'Wall',
    2: 'Filament',
    3: 'Cluster'}

custom_palette = {
    0: 'blue',     # void
    1: 'green',     # wall
    2: 'orange', # filament
    3: 'red'   # clusters
}

cmap4 = plt.get_cmap('magma', 4)
custom_palette2 = cmap4(np.arange(4))

# Historgram of counts with labels giving percentage of each environment
plt.figure(figsize=(10, 6))
plt.hist(DESI_pred, bins=np.arange(5)-0.5, rwidth=0.8, color='skyblue', edgecolor='black')
plt.xticks(np.arange(4), [environ_dicts[i] for i in range(4)])
plt.xlabel('Cosmic Web Environment')
plt.ylabel('Count')
plt.title('Distribution of Cosmic Web Environments in DESI Network')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()

# plt.hist(DESI_pred)
# plt.xlabel('CW Environment')
# plt.ylabel('Count')

import plotly.io as pio
pio.renderers.default = "vscode"
import plotly.graph_objects as go
import plotly.express as px

# Map labels and colors
labels = [environ_dicts[int(label)] for label in DESI_pred]
colors = [custom_palette[int(label)] for label in DESI_pred]

# 2D projection plot of DESI galaxies with cosmic web predictions
ax.set_facecolor('none')    # makes axes transparent

zlims = (-10, 10)  # Set z slab limits in Mpc
fig = plt.figure(figsize=(10, 8), facecolor='white')
ax = fig.add_subplot(facecolor='white')
# set z slab between -10 and 10 Mpc
proj_x = DESI_geom.pos[:, 0][(DESI_geom.pos[:,2]<zlims[1])&(DESI_geom.pos[:,2]>zlims[0])]
proj_y = DESI_geom.pos[:, 1][(DESI_geom.pos[:,2]<zlims[1])&(DESI_geom.pos[:,2]>zlims[0])]
ax.scatter(
    proj_x,
    proj_y,
    c=[custom_palette[int(label)] for label in DESI_pred[(DESI_geom.pos[:,2]<zlims[1])&(DESI_geom.pos[:,2]>zlims[0])]],
    s=2,  # Size of points
    # alpha=0,  # Transparency
    edgecolor='none'  # No edge color
)
ax.grid(False)
ax.legend(handles=[
    plt.Line2D([0], [0], marker='o', color='w', label=environ_dicts[i],
               markerfacecolor=custom_palette[i], markersize=5) for i in range(4)
], title='Cosmic Web Environments', loc='best')
ax.set_xlim(0, 300)  # Set x limits in Mpc
ax.set_ylim(0, 300)  # Set y limits in Mpc
ax.set_xlabel('X (Mpc)')
ax.set_ylabel('Y (Mpc)')
ax.set_title('DESI Galaxy Network with Cosmic Web Predictions (2D Projection)')
# Set aspect ratio to equal for better visualization
ax.set_aspect('equal', adjustable='box')
# Show the plot
plt.show()


# Build the interactive plot
fig = go.Figure(data=[go.Scatter3d(
    x=DESI_geom.pos[:, 0],
    y=DESI_geom.pos[:, 1],
    z=DESI_geom.pos[:, 2],
    mode='markers',
    marker=dict(
        size=2,
        color=colors,  # Custom hex colors
        opacity=0.8
    ),
    text=labels,  # Hover shows label name
    hovertemplate='<b>Environment</b>: %{text}<extra></extra>',
)])

fig.update_layout(
    title='DESI Galaxy Network with Cosmic Web Predictions',
    scene=dict(
        xaxis_title='X (Mpc)',
        yaxis_title='Y (Mpc)',
        zaxis_title='Z (Mpc)'
    ),
    legend_title='CW Environments',
    width=900,
    height=800
)

fig.show()

# Histograms of stellar mass and colour for each environment
LOGMSTAR = DESI_NETWORK.DESI_GAL_CAT.zcat['LOGMSTAR']

fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharex=False, sharey=False, layout='constrained')
axs = axs.flatten()
for i in range(4):
    mask = DESI_pred == i
    axs[i].hist(LOGMSTAR[mask], bins=50, alpha=0.7, label=f'{environ_dicts[i]}', density=True, color=custom_palette[i])
    axs[i].set_xlim([9.,12.5])
    axs[i].set_ylim([0., 1.])
    axs[i].set_xlabel('Log Stellar Mass')
    axs[i].set_ylabel('Frequency')
    axs[i].legend()
# This code is for creating a galaxy network from the DESI BGS catalog and predicting the cosmic web environment using a GAT model.

FLUXG = DESI_NETWORK.DESI_GAL_CAT.zcat['FLUX_G']
FLUXR = DESI_NETWORK.DESI_GAL_CAT.zcat['FLUX_R']
FLUXZ = DESI_NETWORK.DESI_GAL_CAT.zcat['FLUX_Z']

G_MAG_mask = FLUXG > 0
R_MAG_mask = FLUXR > 0
Z_MAG_mask = FLUXZ > 0

G_MAG = 22.5 - 2.5*np.log10(FLUXG[np.where(G_MAG_mask)[0]])
R_MAG = 22.5 - 2.5*np.log10(FLUXR[np.where(R_MAG_mask)[0]])
Z_MAG = 22.5 - 2.5*np.log10(FLUXZ[np.where(Z_MAG_mask)[0]])


fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharex=False, sharey=False, layout='constrained')
axs = axs.flatten()
for i in range(4):
    mask = (DESI_pred == i)[G_MAG_mask]
    axs[i].hist(G_MAG[mask], bins=50, alpha=0.7, label=f'{environ_dicts[i]}', density=True, color=custom_palette[i])
    axs[i].set_xlabel('g Mag.')
    axs[i].set_ylabel('Frequency')
    axs[i].legend()


fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharex=False, sharey=False, layout='constrained')
axs = axs.flatten()
for i in range(4):
    mask = (DESI_pred == i)[R_MAG_mask]
    axs[i].hist(R_MAG[mask], bins=50, alpha=0.7, label=f'{environ_dicts[i]}', density=True, color=custom_palette[i])
    axs[i].set_xlabel('r Mag.')
    axs[i].set_ylabel('Frequency')
    axs[i].legend()


fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharex=False, sharey=False, layout='constrained')
axs = axs.flatten()
for i in range(4):
    mask = (DESI_pred == i)[Z_MAG_mask]
    axs[i].hist(Z_MAG[mask], bins=50, alpha=0.7, label=f'{environ_dicts[i]}', density=True, color=custom_palette[i])
    axs[i].set_xlabel('z Mag.')
    axs[i].set_ylabel('Frequency')
    axs[i].legend()

# g-r color, removing all galaxies with negative or zero flux in either band
BAD_GAL = np.unique(list(np.where((FLUXG<=0))[0])+list(np.where((FLUXR<=0))[0])+list(np.where((FLUXZ<=0))[0]))
GOOD_FLUXG = np.delete(FLUXG, BAD_GAL)
GOOD_FLUXR = np.delete(FLUXR, BAD_GAL)
GOOD_FLUXZ = np.delete(FLUXZ, BAD_GAL)

GOOD_G_MAG = 22.5 - 2.5*np.log10(GOOD_FLUXG)
GOOD_R_MAG = 22.5 - 2.5*np.log10(GOOD_FLUXR)
GOOD_Z_MAG = 22.5 - 2.5*np.log10(GOOD_FLUXZ)

g_r = GOOD_G_MAG-GOOD_R_MAG
#Plotting g-r color for each environment
fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharex=False, sharey=False, layout='constrained')
axs = axs.flatten()
for i in range(4):
    mask = np.delete((DESI_pred == i), BAD_GAL)
    axs[i].hist(g_r[mask], bins=50, alpha=0.7, label=f'{environ_dicts[i]}', density=True, color=custom_palette[i])
    axs[i].set_xlabel('g-r Color')
    axs[i].set_ylabel('Frequency')
    axs[i].legend()

#Plotting g-r color for each environment
fig, axs = plt.subplots(2,2, figsize=(12, 10), sharex=False, sharey=False, layout='constrained')
axs = axs.flatten()
for i in range(4):
    mask = np.delete((DESI_pred == i), BAD_GAL)
    axs[i].scatter(GOOD_R_MAG[mask],g_r[mask],marker='.',  alpha=0.5, label=f'{environ_dicts[i]}', color=custom_palette[i])
    axs[i].set_yscale('log')
    axs[i].set_ylim([1e-2,1e1])
    axs[i].set_xlabel('r Mag.')
    axs[i].set_ylabel('g-r Colour')
    axs[i].legend()

# animation
from matplotlib import animation

fig, ax = plt.subplots(figsize=(10, 8))
sc = ax.scatter([], [], s=1, alpha=0.8)
ax.set_xlim(-310, 310)
ax.set_ylim(-310, 310)
ax.set_xlabel('X (Mpc)')
ax.set_ylabel('Y (Mpc)')
ax.set_aspect('equal')

def update(frame):
    z0, z1 = -300 + frame * 10, -300 + frame * 10 + 10
    mask = (DESI_geom.pos[:, 2] > z0) & (DESI_geom.pos[:, 2] < z1)
    x, y = DESI_geom.pos[mask, 0], DESI_geom.pos[mask, 1]
    colors = [custom_palette[int(l)] for l in DESI_pred[mask]]
    sc.set_offsets(np.c_[x, y])
    sc.set_color(colors)
    ax.set_title(f'z in [{z0}, {z1}] Mpc')
    return sc,
from IPython.display import HTML

ani = animation.FuncAnimation(fig, update, frames=60, interval=200, blit=True)
HTML(ani.to_jshtml())

# ani.save("cosmic_web_z_slab.gif", writer="pillow", fps=5)

# DESI_GAL_CAT = GalaxyCatalogue(
#     PATH="/global/homes/d/dkololgi/GraphWeb_DESI/loa-combined-lowz.fits" # Path to reduced fastspecfit BGS catalog
# )

# DESI_GAL_CAT.cartesian_coord() # Creating cartesian coordinates (x, y, z) for DESI BGS galaxies

# torch_data = Data(
#     pos=torch.tensor(
#         [DESI_GAL_CAT.X, DESI_GAL_CAT.Y, DESI_GAL_CAT.Z], 
#         dtype=torch.float32
#     ).T
# )