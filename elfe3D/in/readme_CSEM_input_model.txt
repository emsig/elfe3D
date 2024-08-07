###############################################
'CSEM_input_model' is an example model.
You can find a figure showing the model in 
elfe3D/docs/PhD_thesis_Rulff_Kappa.pdf
on page 55.
You can the example to understand how to run
simulations and refinement with elfe3D.
The original elfe3D_input file contains all 
specifications to run this example.
Only the mesh files have to be generated from
CSEM_input_model.poly via calling tetgen
in the elfe3D/in directory.
e.g.: tetgen -pq1.3kAaen CSEM_input_model.poly
###############################################

Model specifications:

modelling domain: 30 x 30 x 30 km
z positive upwards

CSEM_input_model.poly file is the input file 
for the mesh generator tetgen.
It contains the starting nodes, facets and regions, 
that are the basis for mesh generation and model
parameter assignment.

receivers: 200 m to 2000 m in x-direction, every 25 m
I put my receivers in small tetrahedra to ensure accurate
responses at the receiver sites. 
Thats why I allocate triangle facets for the receiver sites
at the air-earth interface
and put the z-coordinates of the receivers slightly below 
the air-earth interface.

source: -50 m to 50 m in x-direction on the surface
The source is along 5 m long edge segments.
source moment = 1

anomaly coordinates: 
x y z
900 -100 -600                            
1100 -100 -600                             
1100 -100 -200                           
900 -100 -200                                                           
900 1100 -600                              
1100 1100 -600                              
1100 1100 -200                              
900 1100 -200  

region attributes and resistivities:
3 regions            
1  air      10E-8 Ohmm              
2  Earth    1000 Ohmm             
3  anomaly  10 Ohmm   


To run refinement:

generate initial mesh with 
e.g. tetgen -pq1.6kAaen CSEM_input_model.poly 

choose one frequency in elfe3D_input.txt
e.g. 10 Hz

adjust refinement parameters in elfe3D_input.txt
e.g.
maxRefSteps             10
maxUnknowns             10000000  
betaRef                 0.85
accuracyTol             0.00003
vtk                     1
errorEst_method         4
refStrategy             1


