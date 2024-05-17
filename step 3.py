from Bio.PDB import PDBParser
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import pandas as pd

def calculate_energy_vector(row):
    vector_components = np.array([row['vector_x'], row['vector_y'], row['vector_z']])
    bond_energy = row['bond_energy']
    return vector_components * bond_energy

def calculate_flux_differentials(residue_vectors, coords):
    flux_differentials = {}
    for i in range(1, len(coords) - 1):
        vec1 = coords[i] - coords[i - 1]
        vec2 = coords[i + 1] - coords[i]
        area_vec = np.cross(vec1, vec2)
        norm_area_vec = area_vec / np.linalg.norm(area_vec)

        # we take only interactions. not directions(부호)
        total_flux = sum(np.dot(vec, norm_area_vec) for vec in residue_vectors.get(i, []))
        flux_differentials[i] = np.abs(total_flux)

    return flux_differentials

def process_intra_interactions(file_path, atom_id_cols):
    intra_interaction_df = pd.read_csv(file_path)
    print("Column names in intra_interaction_df:", intra_interaction_df.columns)
    
    intra_interaction_vectors = {}
    for index, row in intra_interaction_df.iterrows():
        vector = np.array([row['vector_x'], row['vector_y'], row['vector_z']])
        bond_energy = row['bond_energy']
        for atom_id in atom_id_cols:  # Use the provided column names
            residue_id = row[atom_id]
            if residue_id not in intra_interaction_vectors:
                intra_interaction_vectors[residue_id] = []
            intra_interaction_vectors[residue_id].append((vector, bond_energy))
    return intra_interaction_vectors

# Process intra-protein interactions for 1QLZ and its mutant with correct column names
intra_vectors_1qlz = process_intra_interactions('/Users/myunghyunjeong/Desktop/1qlz_intra_noncovalent_interactions.csv', ['1qlz_atom_id1', '1qlz_atom_id2'])
intra_vectors_1qlz_mutant = process_intra_interactions('/Users/myunghyunjeong/Desktop/1qlz_mutant_intra_noncovalent_interactions.csv', ['1qlz_atom_id', '1qlz_atom_id.1'])

# Combine intra_vectors for 1QLZ and its mutant
combined_intra_vectors = {**intra_vectors_1qlz, **intra_vectors_1qlz_mutant}

def calculate_combined_energy_vector(row, intra_interaction_vectors):
    inter_vector = np.array([row['vector_x'], row['vector_y'], row['vector_z']])
    inter_energy = row['bond_energy']

    residue_id = row['1qlz_atom_id']
    intra_vectors = intra_interaction_vectors.get(residue_id, [])
    
    combined_vector = inter_vector * inter_energy
    for intra_vector, intra_energy in intra_vectors:
        combined_vector += intra_vector * intra_energy

    return combined_vector

def parse_pdb_for_ribbon(filepath):
    parser = PDBParser()
    structure = parser.get_structure('protein', filepath)
    backbone_atoms = []
    residue_indices = []

    for model in structure:
        if model.get_id() == 0:
            for chain in model:
                if chain.get_id() == 'A':
                    for residue in chain:
                        if residue.get_id()[0] == ' ':
                            for atom in residue:
                                if atom.get_name() in ['N', 'CA', 'C', 'O']:
                                    backbone_atoms.append(atom.get_coord())
                                    residue_indices.append(residue.get_id()[1])
            break

    return backbone_atoms, residue_indices

# Function to plot ribbon diagram with transparency and color gradient
def plot_ribbon_and_linear_sequence_0(coords, diff_changes, title, residue_indices):
    fig = plt.figure(figsize=(10, 8))  # Adjusted for a single plot
    ax = fig.add_subplot(111, projection='3d')

    # Plot the 3D ribbon with transparency
    values = list(diff_changes.values())
    color_norm = Normalize(vmin=min(values), vmax=max(values))
    color_mapper = ScalarMappable(norm=color_norm, cmap='rainbow')

    for i in range(len(coords) - 1):
        start, end = coords[i], coords[i + 1]
        diff_change = diff_changes.get(i, 0)
        color = color_mapper.to_rgba(diff_change)
        ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], color=color, alpha=0.5, linewidth=2)

    # Add transparent text annotations
    for i, (coord, index) in enumerate(zip(coords, residue_indices)):
        ax.text(coord[0], coord[1], coord[2], str(index), color='black', alpha=0.5, fontsize=8)  # Adjust fontsize as needed

    # Adding a color bar
    color_mapper.set_array(values)
    cbar = fig.colorbar(color_mapper, ax=ax, orientation='vertical')
    cbar.set_label('Energy Gradient')
    
    ax.set_title(title)
    plt.show()

def plot_ribbon_relative_color(coords, flux_diff_changes, title):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Normalize color based on the flux changes for this protein
    flux_values = list(flux_diff_changes.values())
    color_norm = Normalize(vmin=min(flux_values), vmax=max(flux_values))
    color_mapper = ScalarMappable(norm=color_norm, cmap='rainbow')

    for i in range(len(coords) - 1):
        start, end = coords[i], coords[i + 1]
        flux_change = flux_diff_changes.get(i, 0)
        color = color_mapper.to_rgba(flux_change)
        ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], color=color, alpha=0.5, linewidth=2)

    # Create the color bar
    color_mapper.set_array(flux_values)
    cbar = fig.colorbar(color_mapper, ax=ax)
    cbar.set_label('Relative Flux Differential')

    ax.set_title(title)
    plt.show()


def plot_ribbon_and_linear_sequence_1(coords, diff_changes, title, residue_indices):
    fig = plt.figure(figsize=(10, 8))  # Adjusted for a single plot
    ax = fig.add_subplot(111, projection='3d')

    # Plot the 3D ribbon with transparency
    values = list(diff_changes.values())
    color_norm = Normalize(vmin=min(values), vmax=max(values))
    color_mapper = ScalarMappable(norm=color_norm, cmap='rainbow')

    for i in range(len(coords) - 1):
        start, end = coords[i], coords[i + 1]
        diff_change = diff_changes.get(i, 0)
        color = color_mapper.to_rgba(diff_change)
        ax.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], color=color, alpha=0.5, linewidth=2)

    # Adding a color bar
    color_mapper.set_array(values)
    cbar = fig.colorbar(color_mapper, ax=ax, orientation='vertical')
    cbar.set_label('Energy Gradient')
    
    # Add transparent text annotations
    ax.set_title(title)
    plt.show()

  
# Step 1: Parse PDB files for ribbon diagrams
ribbon_coords_1qlz, residue_indices_1qlz = parse_pdb_for_ribbon('/Users/myunghyunjeong/Desktop/1qlz.pdb')
ribbon_coords_etc, residue_indices_etc = parse_pdb_for_ribbon('/Users/myunghyunjeong/Desktop/1qlz_mutant.pdb')


# Step 2: Read CSV files and calculate energy vectors for normal and mutant
residue_energy_vectors_normal = {}
residue_energy_vectors_etc = {}

for i in range(17):  # Adjust the range as needed
    df_normal = pd.read_csv(f'/Users/myunghyunjeong/Desktop/results/normal_iteration_{i}_output_vectors.csv')
    df_etc = pd.read_csv(f'/Users/myunghyunjeong/Desktop/results/mutant_iteration_{i}_output_vectors.csv')
    # Print column names for debugging
    print("Column names in df_normal:", df_normal.columns)

    # Use calculate_combined_energy_vector function
    df_normal['combined_energy_vector'] = df_normal.apply(
        lambda row: calculate_combined_energy_vector(row, combined_intra_vectors), axis=1)
    df_etc['combined_energy_vector'] = df_etc.apply(
        lambda row: calculate_combined_energy_vector(row, combined_intra_vectors), axis=1)

    # Use combined_energy_vector for flux calculations
    for _, row in df_normal.iterrows():
        residue_id = row['1qlz_atom_id']
        combined_vector = row['combined_energy_vector']
        residue_energy_vectors_normal.setdefault(residue_id, []).append(combined_vector)

    for _, row in df_etc.iterrows():
        residue_id = row['1qlz_atom_id']
        combined_vector = row['combined_energy_vector']
        residue_energy_vectors_etc.setdefault(residue_id, []).append(combined_vector)

# Debugging: Print out a few entries from the dictionaries
print("Sample from normal:", list(residue_energy_vectors_normal.items())[:3])
print("Sample from etc:", list(residue_energy_vectors_etc.items())[:3])


# Step 3: Calculate differential changes
flux_diff_changes_normal = calculate_flux_differentials(residue_energy_vectors_normal, ribbon_coords_1qlz)
flux_diff_changes_etc = calculate_flux_differentials(residue_energy_vectors_etc, ribbon_coords_etc)


# Step 4: Plot ribbon diagrams
plot_ribbon_and_linear_sequence_1(ribbon_coords_1qlz, flux_diff_changes_normal, '1QLZ Prion Protein: Flux Differential Changes Visualization', residue_indices_1qlz)
plot_ribbon_and_linear_sequence_1(ribbon_coords_etc, flux_diff_changes_etc, '1QLZ_G127V mutant Prion Protein: Flux Differential Changes Visualization', residue_indices_etc)

plot_ribbon_and_linear_sequence_0(ribbon_coords_1qlz, flux_diff_changes_normal, '1QLZ Prion Protein: Flux Differential Changes Visualization', residue_indices_1qlz)
plot_ribbon_and_linear_sequence_0(ribbon_coords_etc, flux_diff_changes_etc, '1QLZ_G127V mutant Prion Protein: Flux Differential Changes Visualization', residue_indices_etc)

plot_ribbon_relative_color(ribbon_coords_1qlz, flux_diff_changes_normal, '1QLZ Prion Protein: Flux Differential Changes Visualization')
plot_ribbon_relative_color(ribbon_coords_etc, flux_diff_changes_etc, '1QLZ_G127V mutant Prion Protein: Flux Differential Changes Visualization')
