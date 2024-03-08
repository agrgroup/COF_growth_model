# COF Growth Model

### Generalized Theory of 2D Polymer Growth to Model the Synthesis of Covalent Organic Frameworks 
> Shubhani Paliwal, Wei Li, Pingwei Liu, and Ananth Govind Rajan

## Contents
This repository outlines a global kinetic model (using the files OptimiseParameters.m, Rate.m, Residual.m, F3calc.m, and F45calc.m) and a local kinetic Monte Carlo (KMC) simulation (using the files KMC_RPT.m and COF_growth.m) to describe the growth of a two-dimensional (2D) covalent organic framework (COF). The code includes the following functions/files:

## OptimiseParameters.m
This function optimizes the parameters involved in the kinetic model using nonlinear least squares to reduce the error between the predicted and the experimental COF yield at various times.
  
## Rate.m
This function calculates the rates of the different types of reactions involved in the process of COF growth.

## Residual.m
This function calculates the COF yield at various times and determines the error between the experimental yield and the yield predicted from the kinetic model.
 
## F3calc.m
This function calculates the number of free ends of the monomers involved in the dimerisation reaction and the number of neighbors this dimer will attach to in the COF polymer. It is applicable to reaction type 3 in the kinetic model.
	
## F45calc.m
This function calculates the number of free ends of the monomer that is attaching to the COF and the number of neighbors this monomer will attach to in the COF polymer. It is applicable to reaction types 4 and 5 in the kinetic model.
		
## KMC_RPT.m
This function performs a single kinetic Monte Carlo (KMC) simulation for modeling COF growth in a 2D hexagonal lattice following the reversible polycondensation termination (RPT) mechanism involving 2 trifunctional monomers and 2 monofunctional inhibitor species. 

## COF_growth.m
This program repeatedly calls the KMC routine and stochastically adds/removes monomers/inhibitors on lattice sites.
