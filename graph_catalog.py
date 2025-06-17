import sys
from galaxy_catalog import GalaxyCatalog
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.preprocessing import PowerTransformer

import scienceplots
plt.style.use(['science', 'no-latex']) 


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

with torch.no_grad():
    DESI_out = inference_model(
        DESI_geom.x, 
        DESI_geom.edge_index, 
        edge_weight=DESI_geom.edge_attr
    )
    DESI_pred = DESI_out.argmax(dim=1).numpy()
    DESI_probs = DESI_out.numpy()

plt.hist(DESI_pred)
plt.xlabel('CW Environment')
plt.ylabel('Count')

environ_dicts = {
    0: 'Void',
    1: 'Wall',
    2: 'Filament',
    3: 'Cluster'}

custom_palette = {
    0: '#4c78a8',     # deep teal
    1: '#a05eb5',     # violet
    2: '#76b7b2', # sky teal
    3: '#e17c9a'   # plum pink
}

custom_palette = {
    0: 'red',     # void
    1: 'green',     # wall
    2: 'blue', # filament
    3: 'yellow'   # plum pink
}

import plotly.io as pio
pio.renderers.default = "vscode"
import plotly.graph_objects as go
import plotly.express as px

# Map labels and colors
labels = [environ_dicts[int(label)] for label in DESI_pred]
colors = [custom_palette[int(label)] for label in DESI_pred]

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

fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharex=False, sharey=False, layout='constrained')
axs = axs.flatten()
for i in range(4):
    mask = DESI_pred == i
    axs[i].hist(FLUXG[mask], bins=50, alpha=0.7, label=f'{environ_dicts[i]}', density=True, color=custom_palette[i])
    axs[i].set_xlabel('Total g Band Flux')
    axs[i].set_ylabel('Frequency')
    axs[i].legend()

FLUXR = DESI_NETWORK.DESI_GAL_CAT.zcat['FLUX_R']
fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharex=False, sharey=False, layout='constrained')
axs = axs.flatten()
for i in range(4):
    mask = DESI_pred == i
    axs[i].hist(FLUXR[mask], bins=50, alpha=0.7, label=f'{environ_dicts[i]}', density=True, color=custom_palette[i])
    axs[i].set_xlabel('Total r Band Flux')
    axs[i].set_ylabel('Frequency')
    axs[i].legend()

FLUXZ = DESI_NETWORK.DESI_GAL_CAT.zcat['FLUX_Z']
fig, axs = plt.subplots(2, 2, figsize=(12, 10), sharex=False, sharey=False, layout='constrained')
axs = axs.flatten()
for i in range(4):
    mask = DESI_pred == i
    axs[i].hist(FLUXZ[mask], bins=50, alpha=0.7, label=f'{environ_dicts[i]}', density=True, color=custom_palette[i])
    axs[i].set_xlabel('Total z Band Flux')
    axs[i].set_ylabel('Frequency')
    axs[i].legend()


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