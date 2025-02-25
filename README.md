# Charge Susceptibility Calculations for Hybridized Materials

Contains the Python and C scripts used to calculate the wavefunction-dependent charge susceptibility (Lindhard function) of materials. The values associated with parameters and functions in these scripts are tailored for calculations associated with hydrogen and lithium chains, which led to the publication of the following paper showing the suppression of expected Peierls instabilities in highly hybridized materials:   
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.064102

However, these scripts can be easily adapted to perform similar calculations for different materials; one must simply follow the following series of steps while making sure to modify parameter values (such as imported DFT data file locations) appropriately.

# Steps for a Complete Calculation



## 1. Tight Binding Model and DFT Calculations

After selecting a material to study, construct a tight binding model capturing the relevant low energy physics. The relevant states can be figured out by identifying the bands near the Fermi energy from a Density Functional Theory (DFT) band structure calculation. THe numerical values of the parameters of the derived Hamiltonian can be chosen through literature values, experiment fits or by fitting its energy eigenvalues to the DFT bands. In the provided code, the tight binding model used is for a one dimensional chain of atomic sites which host *s* and *p* orbitals, and it includes all nearest neighbor hopping parameters *t<sub>ss</sub>*, *t<sub>pp</sub>* and *t<sub>sp</sub>*. It also includes an onsite energy difference between the two orbitals *Δ*




## 2. Orbital Overlap Calculations

Using the script [orbitals_integrals.py](orbital_integrals.py), calculate the overlap integrals up to the desired n<sup>th</sup> neighbor level (the code provided only goes up to the second nearest neighbor range). The overlap intergrals must be calculated for many k-points (refer to the paper), so this is an expensive process. Consequently, all of the integrals are computed through the C script [integrals.c](integrals.c), where the actual functions representing the position and momentum-dependent functional forms of the matrix elements that enter into the susceptibility are defined and called. The python script then exports a file containign all of that data for later use. Make sure to modify the calculation parameters such as the number of k-points, the dimensionality of the integrals and the lattice constants to those of material under study.




## 3. 





