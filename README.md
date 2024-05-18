## 🧬Dynamic Temporal Energy Vectors for Protein Folding Analysis

### Overview
This project focuses on the **analysis of dynamic temporal energy vectors to identify residues in protein structures that are susceptible to structural changes**. By examining both intra- and inter-residue non-covalent interactions, this research aims to pinpoint residues that experience intense non-covalent interactions, making them critical points for protein folding and misfolding.

### Objective
The primary objective of this project is to calculate energy vectors and flux differentials from non-covalent interactions in protein structures. This analysis is crucial for understanding how certain residues contribute to structural stability or instability, particularly **in the context of protein misfolding diseases such as prion diseases**.

### Significance
Calculating flux differentials from energy vectors is essential as it highlights regions of the protein that are under significant stress due to non-covalent interactions. These stressed regions are often the sites where structural changes occur, leading to misfolding or aggregation. Identifying these critical residues can provide insights into the mechanisms of protein folding diseases and guide potential therapeutic interventions.

### Methodology
The project employs a combination of computational tools and visualization techniques to analyze protein structures. The workflow involves the following steps:

1. **Loading and Aligning Protein Structures**: Protein structures are loaded into PyMOL and aligned for consistent analysis.
2. **Calculating Intra- and Inter-Residue Interactions**: Scripts calculate non-covalent interactions within and between residues.
3. **Computing Energy Vectors**: Energy vectors are computed based on interaction types and distances.
4. **Calculating Flux Differentials**: Flux differentials are calculated to identify residues under intense non-covalent interactions.
5. **Visualizing Results**: Results are visualized using 3D ribbon diagrams to highlight the stressed regions in the protein.

### Data and Software Requirements
- **Proteins Analyzed**: The code is designed for analyzing the prion protein structures 6LNI (PDB ID) and 1QLZ (PDB ID).
- **Software Required**: 
  - **PyMOL**: A molecular visualization system used for preparing and aligning protein structures.
  - **Python Libraries**: numpy, pandas, Biopython, matplotlib for data analysis and visualization.
- **PDB Files**: Protein structures can be obtained from the Protein Data Bank (PDB) or from AlphaFold predictions if available.

### AlphaFold Reference
If you use an AlphaFold DB prediction in your work, please cite the following papers:
- Jumper, J et al. Highly accurate protein structure prediction with AlphaFold. Nature (2021).
- Varadi, M et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic Acids Research (2021).

### Scripts and Their Functions

1. **`step 1_intra non-covalent interaction.py`**:
    - Calculates and analyzes intra-residue non-covalent interactions within a single protein structure.
    - Identifies key interactions contributing to the protein's internal stability and folding mechanisms.

2. **`step 1_inter non-covalent interaction.py.py`**:
    - Calculates and analyzes inter-residue non-covalent interactions between different residues in the protein.
    - Focuses on interactions such as hydrogen bonds, salt bridges, pi-pi stacking, pi-cation interactions, and Van der Waals forces.

3. **`step 2_pymol stage pattern.py`**:
    - Prepares and aligns protein structures for analysis by loading PDB files, centering, and translating structures.
    - Ensures that structures are correctly positioned for accurate interaction calculations.
    - Valid protein trajectories are subject to protein conformational and structural change. i.e. 6LNI sides are not susceptible for elongation of 1QLZ. Hence, pymol stage was postulated as a static xy plane. 

4. **`step 3.py`**:
    - Calculates energy vectors and flux differentials for residues based on interaction data.
    - Identifies residues susceptible to structural changes due to intense non-covalent interactions.
    - Utilizes visualization techniques to create 3D ribbon diagrams, highlighting critical residues and their energy gradients.

### Usage Instructions

1. **Set Up Environment**:
    - Ensure that PyMOL and necessary Python libraries (e.g., numpy, pandas, Biopython, matplotlib) are installed.

2. **Load and Align Structures**:
    - Use `step 2_pymol stage pattern.py` to load, align, and prepare protein structures (6LNI and 1QLZ) for analysis.

3. **Analyze Interactions**:
    - Run `step 1_intra non-covalent interaction.py` to analyze intra-residue interactions.
    - Execute `step 1_inter non-covalent interaction.py.py` to analyze inter-residue interactions.

4. **Calculate Energy Vectors and Flux**:
    - Run `step 3.py` to calculate energy vectors, flux differentials, and identify key residues contributing to protein stability and folding mechanisms.

5. **Visualize Results**:
    - Use the visualization functions in `step 3.py` to generate 3D ribbon diagrams that highlight regions of the protein under significant stress.

### Conclusion
This project provides a comprehensive methodology for analyzing protein folding dynamics through the calculation of dynamic temporal energy vectors and flux differentials. By identifying residues susceptible to structural changes, this analysis contributes to a deeper understanding of protein misfolding diseases and offers potential pathways for therapeutic development.

### Contribution
Myunghyun Jeong - mhjonathan at gist.ac.kr
Jaehee Lee - jhlee.ug@gm.gist.ac.kr
