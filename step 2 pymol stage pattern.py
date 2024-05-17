load 6lni.pdb
load 1qlz.pdb
load 1qlz_mutant.pdb


# Align 1QLZ to 6LNI based on the residues 185-195
align 1qlz and i. 185-195, 6lni and i. 185-195

# Align 1QLZ_G127V mutant to 6LNI based on the residues 185-195
align 1qlz_mutant and i. 185-195, 6lni and i. 185-195

center 6lni

translate [-15,5,5], object = 1qlz
translate [-15,5,5], object = 1qlz_mutant

load inter_automate.py
