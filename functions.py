import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolfiles
import os
import math


def get_xray_positions(ref_directory):
    xray_positions = []

    # assign directories

    for filename in os.listdir(ref_directory):
        f = os.path.join(ref_directory, filename)
        # checking if it is a file
        if os.path.isfile(f):
            # get mol of xray structure
            print("Getting PDB: ", f)
            xray_mol = rdmolfiles.MolFromPDBFile(
                f,
                sanitize=True,
                removeHs=False,
                proximityBonding=False,
            )

    # get 3Dpoints of atoms in xray structure, store in xray_positions
    for i, atom in enumerate(xray_mol.GetAtoms()):
        positions = xray_mol.GetConformer().GetAtomPosition(i)
        xray_positions.append(positions)
    print("XRAY positions appended")
    return xray_positions


def compute_geom_distance(centroid_1, point_2):

    return np.sqrt(
        ((point_2.x) - (centroid_1.x)) ** 2
        + ((point_2.y) - (centroid_1.y)) ** 2
        + ((point_2.z) - (centroid_1.z)) ** 2
    )


def total_entropy(N_total, dict_total_entropy):
    result = 0

    for key, value in dict_total_entropy.items():
        result += (int(key) / (N_total)) * float(value)

    return result


def calculate_entropy_channel(N_w, X, O):
    x_i = X / N_w
    o_i = O / N_w

    try:
        H_w = -((x_i) * math.log2(x_i) + (o_i) * math.log2(o_i))
        return H_w
    except ValueError:
        if x_i == 0:
            H_w = -((o_i) * math.log2(o_i))
            return H_w

        elif o_i == 0:
            H_w = -((x_i) * math.log2(x_i))
            return H_w
        else:
            return 0
