
# File Organization

In each molecule directory there are separate tarballs for the dimer and monomer 
calculations. Each tarball contains a series of subdirectories labeled 
either `{molecule}_{dimer/trimer}_{distance}` or 
`{molecule}_{dimer/trimer}_{x}_{y}_{z}` where distance is the distance 
separating the monomers and x,y,z are the translation distances between the 
monomers. 

Note: All trimers are in an equilateral triangle position so there is only one 
distance.

Inside each subdirectory is all of the nwchem input and output files to 
calculate the BSSE and total interaction energy for that geometry. The files are
labeled as follows:
`{input/output}_E_{basis}_{monomers}.txt`

where basis is the full basis used including ghost atoms and monomers are just 
the monomers the calculation is being run on.
