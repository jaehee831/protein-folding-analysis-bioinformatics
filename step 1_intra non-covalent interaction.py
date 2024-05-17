zimport numpy as np
from pymol import cmd, stored
import csv

# Define 5 types of non-covalent bond cutoff distance
hbond_cutoff_distance = 3.0         # 1 Ångströms
salt_bridge_cutoff_distance = 5.0   # 2 Ångströms
#ss_cutoff_distance =               # Define disulfide bond cutoff-- covalent bond
pi_pi_cutoff_distance = 4.0         # 3 Ångströms
pi_cation_cutoff_distance = 4.0     # 4 Define pi-cation interaction cutoff
van_der_waals_cutoff_distance = 5   # 5 Define Van der Waals interaction cutoff

# intra-residue non-covalent interactions
def calculate_bond_energy(bond_type, bond_length):
    if bond_length == 0:
        return 0  # or you could return None and handle it later

    # Constants for calculations
    k_coulomb = 8.9875e9  # Coulomb's constant in Nm^2/C^2
    charge_e = 1.602e-19  # Elementary charge in C
    k_pi = 1e-21          # Arbitrary constant for Pi interactions
    k_vdw = 1e-24         # Arbitrary constant for Van der Waals interactions

    energy = 0
    if bond_type == 'Salt Bridge' or bond_type == 'HBond':
        energy = k_coulomb * (charge_e ** 2) / bond_length
    elif bond_type == 'Pi-Stacking' or bond_type == 'Pi-Cation':
        energy = k_pi / (bond_length ** 2)
    elif bond_type == 'Van der Waals':
        energy = k_vdw / (bond_length ** 6)

    return energy


def get_alpha_carbon_coordinates(model, atom_id):
    stored.residue_number = []
    # Find the residue number for the given atom ID
    cmd.iterate(f"{model} and id {atom_id}", "stored.residue_number.append(resi)")

    if not stored.residue_number:
        print(f"No residue found for atom ID {atom_id} in model {model}.")
        return None

    # Now that we have the residue number, fetch the alpha carbon coordinates for that residue
    stored.alpha_coords = []
    cmd.iterate_state(1, f"{model} and chain A and resi {stored.residue_number[0]} and name CA", "stored.alpha_coords.append((x,y,z))")

    if stored.alpha_coords:
        print(f"Alpha carbon coordinates for {model}, Chain A, Residue {stored.residue_number[0]}: {stored.alpha_coords[0]}")
        return stored.alpha_coords[0]
    else:
        print(f"Alpha carbon not found for {model}, Chain A, Residue {stored.residue_number[0]}")
        return None
    
def analyze_interactions(selection1, selection2, cutoff_distance, interaction_type):
    interaction_details = []
    selection1_list = cmd.index(selection1)
    selection2_list = cmd.index(selection2)

    # Check if the selections are valid
    if not selection1_list or not selection2_list:
        print(f"Invalid selections: {selection1} or {selection2}")
        return interaction_details

    # Iterate through each pair of atoms from the two selections
    for item1 in selection1_list:
        for item2 in selection2_list:
            try:
                # Calculate the distance between the two atoms
                distance = cmd.get_distance(item1, item2)

                # Append interaction details if within the cutoff distance
                if distance <= cutoff_distance:
                    interaction_details.append([str(item1[1]), str(item2[1]), str(distance), interaction_type])
            except pymol.CmdException as e:
                # Print error message if distance calculation fails
                #print(f"Error with atoms {item1[1]} in {selection1} and {item2[1]} in {selection2}: {e}")
                continue

    return interaction_details

# hydrogen bond within 1qlz
cmd.select("1qlz_donors", "1qlz and donor")
cmd.select("1qlz_acceptors", "1qlz and acceptor")
hbond_details_1qlz = analyze_interactions("1qlz_donors", "1qlz_acceptors", hbond_cutoff_distance, "HBond")

# salt bridges within 1qlz
cmd.select("1qlz_salt_bridges", "1qlz and (resn LYS+ARG+HIS+ASP+GLU)")
salt_bridge_details_1qlz = analyze_interactions("1qlz_salt_bridges", "1qlz_salt_bridges", salt_bridge_cutoff_distance, "Salt Bridge")

# pi-pi stacking interactions within 1qlz
cmd.select("1qlz_pi_pi", "1qlz and resn PHE+TYR+TRP+HIS")
pi_stacking_details_1qlz = analyze_interactions("1qlz_pi_pi", "1qlz_pi_pi", pi_pi_cutoff_distance, "Pi-Stacking")

# pi-cation interactions within 1qlz
cmd.select("1qlz_pi_cation", "1qlz and (resn PHE+TYR+TRP+HIS+LYS+ARG)")
pi_cation_details_1qlz = analyze_interactions("1qlz_pi_cation", "1qlz_pi_cation", pi_cation_cutoff_distance, "Pi-Cation")

# vdw interactions within 1qlz
cmd.select("1qlz_vdw", "1qlz within " + str(van_der_waals_cutoff_distance) + " of 1qlz")
van_der_waals_details_1qlz = analyze_interactions("1qlz_vdw", "1qlz_vdw", van_der_waals_cutoff_distance, "Van der Waals")

# interpoint response
print("Intra-protein interaction analysis for '1qlz' completed.")

# Intra-protein non-covalent interactions for 1qlz
hbond_details_1qlz_intra = analyze_interactions("1qlz_donors", "1qlz_acceptors", hbond_cutoff_distance, "HBond")
salt_bridge_details_1qlz_intra = analyze_interactions("1qlz_salt_bridges", "1qlz_salt_bridges", salt_bridge_cutoff_distance, "Salt Bridge")
pi_stacking_details_1qlz_intra = analyze_interactions("1qlz_pi_pi", "1qlz_pi_pi", pi_pi_cutoff_distance, "Pi-Stacking")
pi_cation_details_1qlz_intra = analyze_interactions("1qlz_pi_cation", "1qlz_pi_cation", pi_cation_cutoff_distance, "Pi-Cation")
van_der_waals_details_1qlz_intra = analyze_interactions("1qlz_vdw", "1qlz_vdw", van_der_waals_cutoff_distance, "Van der Waals")

print("Intra-protein interaction analysis for '1qlz' completed.")

all_interactions_1qlz_intra = (hbond_details_1qlz_intra + salt_bridge_details_1qlz_intra +
                               pi_stacking_details_1qlz_intra + pi_cation_details_1qlz_intra +
                               van_der_waals_details_1qlz_intra)

def calculate_unit_vector(coord1, coord2):
    vector = np.array(coord2) - np.array(coord1)
    norm = np.linalg.norm(vector)
    if norm == 0:
        return None
    return vector / norm

# Define the CSV header
header = ["1qlz_atom_id", "1qlz_atom_id", "vector_x", "vector_y", "vector_z", "bond_type", "bond_length", "bond_energy"]

# Define the file path for the new CSV file
csv_file_path = "/users/myunghyunjeong/desktop/results/intra_noncovalent_interactions.csv"

# Write the interaction details to a CSV file
with open(csv_file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(header)

    for interaction in all_interactions_1qlz_intra:
        # Extract the interaction details
        atom_id_1, atom_id_2, distance, interaction_type = interaction
        
        # Calculate the bond energy
        bond_energy = calculate_bond_energy(interaction_type, float(distance))
        
        # Retrieve coordinates for both atoms involved in the interaction
        coord_1 = get_alpha_carbon_coordinates('1qlz', atom_id_1)
        coord_2 = get_alpha_carbon_coordinates('1qlz', atom_id_2)

        if coord_1 and coord_2:
            # Calculate the unit vector
            unit_vector = calculate_unit_vector(coord_1, coord_2)
            if unit_vector is not None:
                # Write the row to the CSV file only if unit vector is not None
                writer.writerow([
                    atom_id_1,
                    atom_id_2,
                    unit_vector[0],
                    unit_vector[1],
                    unit_vector[2],
                    interaction_type,
                    distance,
                    bond_energy
                ])
