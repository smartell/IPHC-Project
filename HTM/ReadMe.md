**Halibut Toy Model (HTM)**
---
Author: Steven Martell

This is a spatially implicit 8 area model for the Pacific halibut stock where 
delay-difference accounting is used to represent the biomass dynamics in each
of the regulatory areas.  The delay-difference model is fit to CPUE data in each
of the regulatory areas, and movement is based on a simple gravity model.  The 
stock-recruitment relationship is based on the total coast-wide stock, and 
recruitment to each regulatory area is based on the eigenvector of the dispersal
matrix.

	R=a*S/(1-b*S)
