# elfe3D
Modelling with the total **el**ectric field approach using **f**inite **e**lements in **3D**

_About:_

`elfe3D` is a 3D forward modelling code that can simulate electric and magnetic field responses from frequency-domain controlled-source electromagnetic geophysical setups. It uses tetrahedral meshes and first-order finite-element approximations. In addition, adaptive mesh refinement approaches are implemented.

_Statement of need:_

`elfe3D`  solves forward problems arising from the curl-curl equation in terms of the total electric field using a direct forward solver. The code is designed for Earth Scientists who want to simulate electric and magnetic field responses originating from a transmitter and the interaction of its transmitted signal with the 3D Earth. This so-called controlled-source electromagnetic method is used to search for resources and environmental applications, such as geothermal energy, minerals or groundwater. The air and the Earth’s subsurface consist of cells hosting variable model parameters: isotropic electric resistivities and magnetic permeabilities. Compared to standard electromagnetic geophysical simulation software, `elfe3D` excels in flexibility regarding subsurface geometries and survey settings, i.e. receivers can be arbitrarily placed in the modelling domain and the electrical properties can be flexibly distributed in the subsurface upon model design. Implemented adaptive mesh refinement approaches can automatically design problem-specific meshes and optimise computational load and solution accuracy.

_Contributions:_

An earlier version of the code that `elfe3D` is based on was developed by Paula Rulff with contributions from Laura Maria Buntin and Thomas Kalscheuer at Uppsala University from 2018-2023 financed by the Smart Exploration project (European Union’s Horizon 2020 funding, grant agreement No. 775971).

The present version of `elfe3D` was released in 2024 under the Apache License, Version 2.0. Further developments of `elfe3D` by Paula Rulff, now at Delft University of Technology, are ongoing. Suggestions for improvements are welcome!

_Contact_: 

If you would like to contribute to `elfe3D`, report issues with `elfe3D` or seek support, please send an email to p.rulff@tudelft.nl.

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

