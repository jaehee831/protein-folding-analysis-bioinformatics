import csv
import numpy as np
from pymol import cmd, stored

# Define 5 types of non-covalent bond cutoff distance
hbond_cutoff_distance = 3.0         # 1 Ångströms
salt_bridge_cutoff_distance = 5.0   # 2 Ångströms
#ss_cutoff_distance =               # Define disulfide bond cutoff-- covalent bond
pi_pi_cutoff_distance = 4.0         # 3 Ångströms
pi_cation_cutoff_distance = 4.0     # 4 Define pi-cation interaction cutoff
van_der_waals_cutoff_distance = 5   # 5 Define Van der Waals interaction cutoff


for i in range(100):  # Adjust the range as necessary
    # Translate and rotate the objects
    cmd.translate([-.1, .1, .1], "1qlz")
    cmd.translate([-.1, .1, .1], "1qlz_mutant")
    
    cmd.rotate('x', .1, "1qlz")
    cmd.rotate('x', .1, "1qlz_mutant")
    
    cmd.rotate('y', 1, "1qlz")
    cmd.rotate('y', .1, "1qlz_mutant")
    
    cmd.rotate('z', .1, "1qlz")
    cmd.rotate('z', .1, "1qlz_mutant")

    print(f"retrieving iteration {i}")
    
    # hydrogen bond (donors from chain A of 6lni and all acceptors from 1qlz)
    cmd.select("6lni_chainA_donors", "chain A and 6lni and donor")
    cmd.select("1qlz_acceptors", "1qlz and acceptor")
    cmd.select("1qlz_mutant_acceptors", "1qlz_mutant and acceptor")


    # salt bridges (positively and negatively charged residues)
    cmd.select("6lni_salt_bridges", "chain A and 6lni and (resn LYS+ARG+HIS+ASP+GLU)")
    cmd.select("1qlz_salt_bridges", "1qlz and (resn LYS+ARG+HIS+ASP+GLU)")
    cmd.select("1qlz_mutant_salt_bridges", "1qlz_mutant and (resn LYS+ARG+HIS+ASP+GLU)")


    #  pi-pi stacking interactions (aromatic residues)
    cmd.select("6lni_pi_pi", "chain A and 6lni and resn PHE+TYR+TRP+HIS")
    cmd.select("1qlz_pi_pi", "1qlz and resn PHE+TYR+TRP+HIS")
    cmd.select("1qlz_mutant_pi_pi", "1qlz_mutant and resn PHE+TYR+TRP+HIS")


    # pi-cation interactions (aromatic residues and positively charged residues)
    cmd.select("6lni_pi_cation", "chain A and 6lni and (resn PHE+TYR+TRP+HIS+LYS+ARG)")
    cmd.select("1qlz_pi_cation", "1qlz and (resn PHE+TYR+TRP+HIS+LYS+ARG)")
    cmd.select("1qlz_mutant_pi_cation", "1qlz_mutant and (resn PHE+TYR+TRP+HIS+LYS+ARG)")


    # vdw interactions
    cmd.select("6lni_vdw", "chain A and 6lni within " + str(van_der_waals_cutoff_distance) + " of 1qlz")
    cmd.select("1qlz_vdw", "1qlz within " + str(van_der_waals_cutoff_distance) + " of 6lni")
    cmd.select("1qlz_mutant_vdw", "1qlz_mutant within " + str(van_der_waals_cutoff_distance) + " of 6lni")


    # interpoint response
    print("Analysis completed. Interaction details saved.")


    # Function to analyze interactions and populate details list
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



    hbond_details_1qlz = analyze_interactions("6lni_chainA_donors", "1qlz_acceptors", hbond_cutoff_distance, "HBond")
    pi_stacking_details_1qlz = analyze_interactions("6lni_pi_pi", "1qlz_pi_pi", pi_pi_cutoff_distance, "Pi-Stacking")
    salt_bridge_details_1qlz = analyze_interactions("6lni_salt_bridges", "1qlz_salt_bridges", salt_bridge_cutoff_distance, "Salt Bridge")
    pi_cation_details_1qlz = analyze_interactions("6lni_pi_cation", "1qlz_pi_cation", pi_cation_cutoff_distance, "Pi-Cation")
    van_der_waals_details_1qlz = analyze_interactions("6lni_vdw", "1qlz_vdw", van_der_waals_cutoff_distance, "Van der Waals")


    hbond_details_1qlz_mutant = analyze_interactions("6lni_chainA_donors", "1qlz_mutant_acceptors", hbond_cutoff_distance, "HBond")
    pi_stacking_details_1qlz_mutant = analyze_interactions("6lni_pi_pi", "1qlz_mutant_pi_pi", pi_pi_cutoff_distance, "Pi-Stacking")
    salt_bridge_details_1qlz_mutant = analyze_interactions("6lni_salt_bridges", "1qlz_mutant_salt_bridges", salt_bridge_cutoff_distance, "Salt Bridge")
    pi_cation_details_1qlz_mutant = analyze_interactions("6lni_pi_cation", "1qlz_mutant_pi_cation", pi_cation_cutoff_distance, "Pi-Cation")
    van_der_waals_details_1qlz_mutant = analyze_interactions("6lni_vdw", "1qlz_mutant_vdw", van_der_waals_cutoff_distance, "Van der Waals")


    # Combine interaction details for '1qlz', '1qlz_mutant'
    all_interactions_1qlz = hbond_details_1qlz + pi_stacking_details_1qlz + salt_bridge_details_1qlz + pi_cation_details_1qlz + van_der_waals_details_1qlz
    all_interactions_1qlz_mutant = hbond_details_1qlz_mutant + pi_stacking_details_1qlz_mutant + salt_bridge_details_1qlz_mutant + pi_cation_details_1qlz_mutant + van_der_waals_details_1qlz_mutant

    # Retrieve alpha carbon coordinates for a given atom ID at pymol -- not coord in PDB
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


    def calculate_unit_vector(coord1, coord2):
        vector = np.array(coord2) - np.array(coord1)
        norm = np.linalg.norm(vector)
        if norm == 0:
            return None
        return vector / norm


    def calculate_bond_energy(bond_type, bond_length):
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
    
    # Process interactions and calculate unit vectors
    normal_unit_vectors = []

    for interaction in all_interactions_1qlz:
        atom_id_6lni = interaction[0]
        atom_id_1qlz = interaction[1]

        coord_6lni = get_alpha_carbon_coordinates('6lni', atom_id_6lni)
        coord_1qlz = get_alpha_carbon_coordinates('1qlz', atom_id_1qlz)
        
        if coord_6lni and coord_1qlz:
            normal_unit_vector = calculate_unit_vector(coord_6lni, coord_1qlz)
            if normal_unit_vector is not None:
                normal_unit_vectors.append(normal_unit_vector)
            else:
                print(f"Could not calculate unit vector normal 1qlz & 6lni for interaction {interaction}")

    mutant_unit_vectors = []
    
    for interaction in all_interactions_1qlz_mutant:
        atom_id_6lni = interaction[0]
        atom_id_1qlz_mutant = interaction[1]

        coord_6lni = get_alpha_carbon_coordinates('6lni', atom_id_6lni)
        coord_1qlz_mutant = get_alpha_carbon_coordinates('1qlz_mutant', atom_id_1qlz_mutant)
        
        if coord_6lni and coord_1qlz_mutant:
            mutant_unit_vector = calculate_unit_vector(coord_6lni, coord_1qlz_mutant)
            if mutant_unit_vector is not None:
                mutant_unit_vectors.append(mutant_unit_vector)
            else:
                print(f"Could not calculate unit vector for mutant 1qlz & 6lni interaction {interaction}")
                
    # Save the unit vectors to a CSV file
    normal_output_file_path = f'/Users/myunghyunjeong/Desktop/results/normal_iteration_{i}_output_vectors.csv'
    with open(normal_output_file_path, "w", newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["6lni_atom_id", "1qlz_atom_id", "vector_x", "vector_y", "vector_z", "bond_type", "bond_length", "bond_energy"])
        
        for interaction, vector in zip(all_interactions_1qlz, normal_unit_vectors):
            bond_energy = calculate_bond_energy(interaction[3], float(interaction[2]))
            row = [
                interaction[0],  # 6lni_atom_id
                interaction[1],  # 1qlz_atom_id
                vector[0],       # vector_x
                vector[1],       # vector_y
                vector[2],       # vector_z
                interaction[3],  # bond_type
                interaction[2],  # bond_length
                bond_energy      # bond_energy
            ]
            csv_writer.writerow(row)

    mutant_output_file_path = f'/Users/myunghyunjeong/Desktop/results/mutant_iteration_{i}_output_vectors.csv'
    with open(mutant_output_file_path, "w", newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["6lni_atom_id", "1qlz_atom_id", "vector_x", "vector_y", "vector_z", "bond_type", "bond_length", "bond_energy"])
        
        for interaction, vector in zip(all_interactions_1qlz_mutant, mutant_unit_vectors):
            bond_energy = calculate_bond_energy(interaction[3], float(interaction[2]))
            row = [
                interaction[0],  # 6lni_atom_id
                interaction[1],  # 1qlz_atom_id
                vector[0],       # vector_x
                vector[1],       # vector_y
                vector[2],       # vector_z
                interaction[3],  # bond_type
                interaction[2],  # bond_length
                bond_energy      # bond_energy
            ]
            csv_writer.writerow(row)
    
    print(f"Saved all unit vectors and bond energies to the CSV file.")
