# Charge Susceptibility Calculations for Hybridized Materials

Contains the Python and C scripts used to calculate the wavefunction-dependent charge susceptibility (Lindhard function) of materials. The values associated with parameters and functions in these scripts are tailored for calculations associated with hydrogen and lithium chains, which led to the publication of the following paper showing the suppression of expected Peierls instabilities in highly hybridized materials.  
**[Link to paper](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.064102)**

If you use or are inspired by the contents of this repository, please cite the associated paper like so:\
**Nassim Derriche, Ilya Elfimov, and George Sawatzky. Suppression of Peierls-like nesting-based instabilities in solids. Physical Review B, 106(6):064102, August 2022**

However, these scripts can be easily adapted to perform similar calculations for different materials; one must simply follow the following series of steps while making sure to modify parameter values (such as imported DFT data file locations) appropriately.

# Steps for a Complete Calculation



## 1. Tight Binding Model and DFT Calculations

After selecting a material to study, construct a tight binding model capturing the relevant low energy physics. The relevant states can be figured out by identifying the bands near the Fermi energy from a Density Functional Theory (DFT) band structure calculation. The numerical values of the parameters of the derived Hamiltonian can be chosen through literature values, experiment fits or by fitting its energy eigenvalues to the DFT bands. In the provided code, the tight binding model used is for a one-dimensional chain of atomic sites which host *s* and *p* orbitals, and it includes all nearest neighbor hopping parameters *t<sub>ss</sub>*, *t<sub>pp</sub>* and *t<sub>sp</sub>*. It also includes an onsite energy difference *Δ* between the two orbitals.




## 2. Orbital Overlap Calculations

Using the script [orbitals_integrals.py](orbital_integrals.py), calculate the overlap integrals up to the desired n<sup>th</sup> neighbor level (the code provided only goes up to the second nearest neighbor range). The overlap integrals must be calculated for many k-points (refer to the paper), so this is an expensive process. Consequently, all of the integrals are computed through the C script [integrals.c](integrals.c), where the actual functions representing the position and momentum-dependent functional forms of the matrix elements that enter into the susceptibility are defined and called. The Python script then exports a file containing all of that data for later use. Make sure to modify the calculation parameters such as the number of k-points, the dimensionality of the integrals and the lattice constants to those of the material under study.




## 3. Susceptibility Calculations

The Python script [susceptibility.py](susceptibility.py) is separated into thoroughly commented code cells that each serve a different purpose that allow one to calculate and visualize the susceptibility for a chosen system. To start, the first cell imports and plots a DFT-calculated orbital-weighted band structure (in this case from the DFT software FPLO). Then, a simple script importing DFT band structure data and fitting the selected bands to a tight binding cosine function in order to extract the Hamiltonian parameters is presented. Then, the script contains code that calculates the eigenvalues and eigenstates of a tight binding model of interest directly for many k-points which are then exported for later use. The full susceptibility which utilizes these parameters, along with the matrix elements previously numerically calculated with C, is then computed (in this case again for a 1D chain composed of *s* and *p* orbitals). The other cells allow for plotting the susceptibility and its different orbital-dependent components, analyzing the temperature dependence of the susceptibility through a Fermi-Dirac distribution function that is also defined in the script, plotting the overlap integrals in real space and plotting the chosen bands' orbital eigenstates.





