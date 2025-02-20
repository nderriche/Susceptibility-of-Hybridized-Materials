# Charge Susceptibility Calculations for Hybridized Materials

Contains the Python and C scripts used to calculate the wavefunction-dependent charge susceptibility (Lindhard function) of materials. The values associated with parameters and functions in these scripts are tailored for calculations associated with hydrogen and lithium chains, which led to the publication of the following paper showing the suppression of expected Peierls instabilities in highly hybridized materials: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.106.064102

However, these scripts can be easily adapted to perform similar calculations for different materials; one must simply follow the following series of steps while making sure to modify parameter values (such as imported DFT data file locations) appropriately.

# Steps for a Complete Calculation

## 1. Orbital Overlap Calculations


