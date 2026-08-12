# ==============================================================
# 1. Import libraries
# ==============================================================

import numpy as np
import os
import networkx as nx
from numpy import genfromtxt
from gtda.graphs import GraphGeodesicDistance
from gtda.homology import FlagserPersistence
from gtda.diagrams import PersistenceLandscape
from gtda.plotting import plot_diagram
from IPython.display import display

# ==============================================================
# 2. Define paths 
# ==============================================================

input_dir = "../data/GRN_adjacency_matrices"       # folder containing GRN adjacency matrices
output_dir = "../output/PH"           # folder to store PH 

# ==============================================================
# 3. Define function to compute PH
# ==============================================================

def compute_persistence(adjacency_matrix):
    """Computes persistent homology from adjacency matrix."""
    G = nx.convert.to_networkx_graph(adjacency_matrix)
    adj = nx.adjacency_matrix(G)

    X_ggd = GraphGeodesicDistance(directed=False, unweighted=False).fit_transform([adj])
    PH = FlagserPersistence(directed=False).fit_transform(X_ggd)

    return PH[0]  

# ==============================================================
# 4. Batch run over all matrices
# ==============================================================

files = sorted([f for f in os.listdir(input_dir) if f.endswith(".csv")])
print(f"Found {len(files)} input matrices.")

for i, fname in enumerate(files):
    print(f"\nProcessing file {i+1}/{len(files)}: {fname}")

    # Read matrix
    mat = genfromtxt(os.path.join(input_dir, fname), delimiter=',')
    
    # Compute PH
    PH = compute_persistence(mat)

    # Save PH output
    out_path = os.path.join(output_dir, fname.replace(".csv", "_PH.csv"))
    np.savetxt(out_path, PH, delimiter=',', fmt="%.4f")

    print(f"Saved PH results at: {out_path}")












