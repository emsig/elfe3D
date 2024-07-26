# elfe3D
Modelling with the total **el**ectric field approach using **f**inite **e**lements in **3D**

_About:_

**elfe3D** is a 3D forward modelling code that can be used to simulate electric and magnetic field responses from frequency-domain controlled-source electromagnetic setups. It uses tetrahedral meshes and first-order finite-element approximations.

_Contributions:_

A first version of **elfe3D** was developed by Paula Rulff with contributions from Laura Maria Buntin and Thomas Kalscheuer at Uppsala University from 2018-2023.

A modified version of **elfe3D** is implemented in the inversion software emilia (cite).

Further developments of **elfe3D** by Paula Rulff, now at Delft University of Technology, are ongoing. Suggestions for improvements are welcome!

_Contact_: p.rulff@tudelft.nl

_Getting started:_

- gfortran should be part of your system via the gnu compiler collection (gcc). Intel Fortran compilation is possible, but not fully tested.
- OpenBLAS and make packages are required.
- tetgen: The open source mesh generator tetgen must be installed. It can be downloaded from http://wias-berlin.de/software/index.jsp?id=TetGen&lang=1#Download or via `sudo apt install tetgen` on Ubuntu.
- MUMPS is an open source direct solver available at https://mumps-solver.org. Read the MUMPS documentation for compiling instructions. Link the MUMPS routines in your makefile. In the code documentation you find detailed instructions to install and link MUMPS in an Ubuntu system.
- Modify your makefile as appropriate.
- Compile with `make all`.
- Adjust your input file `in/elfe3D_input_txt`.
- Run tetgen to generate your mesh input files.
- Run **elfe3d**.

_Credits:_

If you publish results generated with **elfe3D**, please give credit to the **elfe3D** developers by citing:

Paula Rulff, Laura M Buntin, Thomas Kalscheuer, Efficient goal-oriented  mesh refinement in 3-D finite-element modelling adapted for controlled source electromagnetic surveys, Geophysical Journal International, Volume 227, Issue
3, December 2021, Pages 1624–1645, https://doi.org/10.1093/gji/ggab264

Do not forget to aslo acknowledge MUMPS and tetgen developers!

and refer to the **elfe3D** the version you used with the ZENODO DOI: Link to Zenodo DOI
