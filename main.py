import functions
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolfiles
import os
import sys

"""
Whole script idea: https://stackoverflow.com/questions/35709562/how-to-calculate-clustering-entropy-a-working-example-or-software-code%E2%80%8B

Script to calculate the entropy of an array of clusters. 

by XFF
"""

# assign directories
directory = "data"  # PUT YOUR TEST FILES HERE
ref_directory = "ref"  # PUT YOUR XRAY STRUCT HERE. ALIGN PREVIOUSLY WITH TEST FILES.
threshold = 3.2  # CHANGE THE THRESHOLD VALUE
total_entropy = []
N_total = 0
xray_positions = functions.get_xray_positions(ref_directory)


channel_pdbs = []
# iterate over files in
# that directory
print("Loading water channels in pdb format")
for filename in os.listdir(directory):
    f = os.path.join(directory, filename)
    # checking if it is a file
    if os.path.isfile(f):

        channel_pdbs.append(
            rdmolfiles.MolFromPDBFile(
                f,
                sanitize=True,
                removeHs=False,
                proximityBonding=False,
            )
        )
print("Loaded ", len(channel_pdbs), " PDBs")
print("Begin cluster entropy calculation")
channel_entropy = []
for mol in channel_pdbs:
    # Define arrays for X and O points within channel
    o_ = []
    x_ = []

    for i, at in enumerate(mol.GetAtoms()):  # for each atom in the given molecule

        position = mol.GetConformer().GetAtomPosition(i)  # calculate its position and
        for xyz in xray_positions:  # given the xray positions
            val = functions.compute_geom_distance(
                position, xyz
            )  # check the distance between all of them
            if val <= threshold:  # append 1 to O if the condition is met
                o_.append(1)
                break
            elif val > threshold:  # otherwise continue
                continue

        else:  # and if it fails to append to O, append to X
            x_.append(1)

    N_w = len(x_) + len(o_)  # define total number of points in channel
    total_entropy.append(
        (N_w, functions.calculate_entropy_channel(N_w, len(x_), len(o_)))
    )
print("Done.\nCalculating entropy of system")
# calculate total number of channels. Convert tuple into dict for easier iteration
dict_total_entropy = dict(total_entropy)
for key in dict_total_entropy.keys():
    N_total += int(key)

result = functions.total_entropy(N_total, dict_total_entropy)
print("Done.\nTotal entropy of system: ", result)

"""Pseudocode:
channel_entropies = []

Get reference PDB (previously aligned with the channels)

for each xyz in reference_struct:
    store position in xray_positions_array
    
For each PDB in data folder:
    channel = Channel([],0) #create channel object to calculate entropy
    convert pdb in rdkit mol
    for at in mol.GetAtoms():
        for i in xray_positions_array:
            if distance(at, i) <= thresh:
                channel.append("O")
            else:
                channel.append("X")
    channel_entropies.append(calculate_entropy_of_channel(channel)      
    """
