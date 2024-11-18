# elfe3D
Modelling with the total **el**ectric field approach using **f**inite **e**lements in **3D**

_About:_

`elfe3D` is a 3D forward modelling code that can simulate electric and magnetic field responses from frequency-domain controlled-source electromagnetic geophysical setups. It uses tetrahedral meshes and first-order finite-element approximations. In addition, adaptive mesh refinement approaches are implemented.

_Statement of need:_

_Contributions:_

An earlier version of the code that `elfe3D` is based on was developed by Paula Rulff with contributions from Laura Maria Buntin and Thomas Kalscheuer at Uppsala University from 2018-2023 financed by the Smart Exploration project (European Union’s Horizon 2020 funding, grant agreement No. 775971).

The present version of `elfe3D` was released in 2024 under the Apache License, Version 2.0. Further developments of `elfe3D` by Paula Rulff, now at Delft University of Technology, are ongoing. Suggestions for improvements are welcome!

_Contact_: p.rulff@tudelft.nl

_Getting started:_

You find the `elfe3D` manual including instalation instructions in `elfe3D/elfe3D/README.md`.
`elfe3D` can be compiled with the provided Makefile.
Note that, the open source mesh generator `tetgen` and the direct solver `MUMPS` must be installed additionally.

_Credits:_

If you publish results generated with `elfe3D`, please give credit to the `elfe3D` developers by citing:

Paula Rulff, Laura M Buntin, Thomas Kalscheuer, Efficient goal-oriented  mesh refinement in 3-D finite-element modelling adapted for controlled source electromagnetic surveys, Geophysical Journal International, Volume 227, Issue
3, December 2021, Pages 1624–1645, https://doi.org/10.1093/gji/ggab264

and refer to the `elfe3D` version you used via the ZENODO DOI: https://doi.org/10.5281/zenodo.13309721

Do not forget to acknowledge `MUMPS` and `tetgen` developers!

