# elfe3D
Modelling with the total electric field approach using finite elements in 3-D

About:


Getting started:
• gfortran should be part of your system in Gnu compiler collection
(gcc). Intel Fortran options are included, but not fully tested.
• OpenBLAS and make packages are required.
• tetgen: The open source mesh generator tetgen must be installed. It
can be downloaded from
http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1#Download
• MUMPS is an open source direct solver available at https://mumps-
solver.org. Read the MUMPS documentation for compiling instruc-
tions. Link the MUMPS routines in your makefile.
• Modify your makefile as appropriate.
• Compile with ’make all’.
• Adjust your input files.
• Run elfe3d.
