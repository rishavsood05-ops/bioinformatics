# HER2-Trastuzumab Epitope Mapping and Thermodynamic Affinities Workspace

This repository houses the computational workflows, localized force-field automation steps, and data-processing scripts used to map structural hotspots at the HER2-Trastuzumab interface. This study explicitly characterizes the thermodynamic and biophysical destabilization penalties of the single-point missense mutation E558F.

## Repository Contents
* **scripts/affinity_calculator.py**: Custom Python calculation script executed in Thonny to model physiological dissociation constant (Kd) shifts.
* **scripts/pymol_visualization.pml**: Command script tracking localized non-covalent networks and structural clashing vectors.

---

## Technical Workflow Execution Guide

### Phase 1: Structural Input Retrieval and Isolation
1. Launch PyMOL and access the command terminal input lines.
2. Retrieve the high-resolution crystalline coordinate files from the RCSB Protein Data Bank repository by executing: `fetch 1n8z`
3. Clean the biological workspace to maintain the core interacting variable complex by removing water molecules and extra chains: `remove solvent`
4. Verify your chain assignments using the Sequence Viewer window under **Display -> Sequence** to map the Trastuzumab light chain to Chain A, the heavy chain to Chain B, and the HER2 extracellular domain to Chain C.

### Phase 2: In Silico Mutagenesis and Minimization
1. Submit your pristine wild-type complex to the MutaBind2 architecture server.
2. Input the chain variable groupings, assigning the antibody loops (Chains A and B) to Partner Box 1, and the HER2 receptor (Chain C) to Partner Box 2.
3. Formulate the single-point substitution mutation string exactly as: `EC558F`
4. Execute the structural minimization run to calculate the total binding free energy shift, yielding a deleterious destabilization penalty of +1.94 kcal/mol.

### Phase 3: Interactive Macromolecular Visualization in PyMOL
Because the international server environment restricts coordinate file downloads, reconstruct the optimized side-chain rotamers manually using a hybrid GUI and command workflow:

1. In the top PyMOL menu bar, click **Wizard -> Mutagenesis** and ensure the mode is configured to **Protein**.
2. Click directly on any atom of Glutamate 558 on Chain C to activate the selection disk.
3. Change the active amino acid type to **PHE** inside the lower-right control panel.
4. Use the frame arrows (`<-` and `->`) to cycle through rotamers, picking the conformation that faces Trastuzumab's paratope loops without breaking standard atomic valency. Click **Apply** and click **Done**.
5. Change your visual scene styles to clean, publication-ready parameters by executing the command lines inside `scripts/pymol_visualization.pml`.

### Phase 4: Drawing Steric Clashing Vectors via Measurement Wizard
1. In the top file menu bar of PyMOL, navigate to **Wizard -> Measurement**.
2. Ensure the lower-right dropdown menu is set to **Distances**.
3. Move your mouse to the viewport and click directly on an outermost tip carbon atom of your new Phenylalanine-558 ring.
4. Click directly on the closest facing side-chain atom situated on Trastuzumab's cyan loop (Chain A).
5. PyMOL will instantly snap a dashed line between the atoms and render your exact sub-Van der Waals spatial distance label, measuring a severe steric clash of 2.8 Ångströms.

### Phase 5: Thermodynamic Calculations via Thonny Python
1. Launch the Thonny Integrated Development Environment.
2. Open a blank script tab and copy the complete processing logic located in `scripts/affinity_calculator.py`.
3. Click the green **Run** button to process the Gibbs free energy relationship at a fixed human physiological body temperature of 310.15 K.
4. The calculation script will automatically output your finalized, publication-quality bar visualization and prove a ~23.3-fold drop in antibody binding affinity.
