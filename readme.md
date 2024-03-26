# Fit function for spr data with variable supplied motor concentrations

This project attempts to fit data created by an SPR experiment using the lmfit library and trying to find the origind of the brick state.

The base assumption is the existence of two washing stages, the first one, where motors and ligands are supplied and the reponse increases and the second stage where the motor concentration is removed but any ligands are still supplied. 

The response is given by how many motors bind to the DNA that is on the chip.

It fits the data using two approaches:

* assuming that the brick state creation rate depends linearly on the motor concentration added $k_{23} = a_{23} C_\text{motor}$. Further, as that corresponds to e.g. a capture of another motor, this might mean that the response increases further for this brick state (double the mass on the chip)

* The brick state creation rate does not depend on the added motor concentration $k_{23} = a_{23}$. This might mean e.g. conformational change, i.e. the response stays the same (mass on chip does not change)


The 'Fit_program_collective_fit.ipynb' corresponds to the program and is heavily commented to explain how to use the data of an xlsx file.


# Needed dependencies:
scipy
matplotlib
openpyxl
lmfit