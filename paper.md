---
title: 'elfe3D: Modelling with the total **el**ectric field approach using **f**inite **e**lements in **3D**'
tags:
  - Simulations
  - Geophysics
  - Electromagnetics
  - Controlled-sources
  - Frequency Domain
  - Fortran
authors:
  - name: Paula Rulff
    orcid: 0000-0001-6714-8008
    affiliation:  "1, 2"
affiliations:
 - name: Department of Geoscience and Engineering, TU Delft, Netherlands
   index: 1
 - name: Department of Earth Sciences, Uppsala University, Sweden
   index: 2

date: 18 November 2024
bibliography: paper.bib
---

# Summary

The controlled-source electromagnetic method is a geophysical technique that detects variations in electrical and magnetic material properties in the subsurface. The method is valuable for both environmental assessments and resource exploration. To evaluate electromagnetic data obtained over complex subsurface structures, precise three-dimensional numerical modeling software is required. This is particularly important for tasks related to survey design and for inversion routines that estimate subsurface models based on field data inputs. The forward modelling program `elfe3D` simulates electric and magnetic field responses for three-dimensional subsurface models, i.e. the data of frequency-domain controlled-source electromagnetic configurations. To optimise the balance between the size of the problems and the accuracy of the solutions, adaptive mesh refinement strategies are employed. The `elfe3D` program allows the user to define an arbitrary survey setup: The sensor locations, where the electric and magnetic field responses are calculated, can be placed at or below the subsurface, or in the air. This flexibility facilitates realistic modelling of surface-based, borehole and airborne controlled-source electromagnetic surveys. 


# Statement of need

Contemporary electromagnetic geophysical field investigations employ advanced measurement configurations designed to image complex subsurface structures. This context necessitates the development of flexible and manageable three-dimensional (3D) simulation software tailored for the evaluation of controlled-source electromagnetic data [@Rochlitz2019]. To address this task, `elfe3D` was developed to numerically calculate forward responses, i.e. electric and magnetic field components, for 3D subsurface models studied with frequency-domain controlled-source electromagnetic surveys [Figure \ref{fig:example}]. 
Compared to standard electromagnetic geophysical simulation software, `elfe3D` excels in flexibility regarding subsurface properties and geometries as well as survey settings, i.e. receivers can be arbitrarily placed in the modelling domain. Implemented adaptive mesh refinement approaches can automatically design problem-specific meshes and balance computational load and solution accuracy.

![Key steps of the forward modelling procedure including the choice of a subsurface model and source-receiver setup (Step I; note that the subsurface anomaly and survey setup are enlarged in the image for better visibility), the meshing of the modelling domain (Step II; note that a slice through the inner modelling domain is displayed) and the calculation of field responses with `elfe3D` (Step III; note that only the electric field component in x-direction (Ex) is displayed). Figure adapted from PhD thesis: @Rulff2023.\label{fig:example}](modelling-procedure-elfe3D.png){ width=90% }

# State of the field

Common strategies in electromagnetic modelling often avoid 3D simulations due to their substantial computational demands and the predominant focus tends to be on electrical resistivity models, as this material property is dominant at typical frequencies used in controlled-source electromagnetic surveys. The numerical techniques required to facilitate 3D simuations and their integration into controlled-source electromagnetic modelling code is currently a prominent area of academic research [@Boerner2010] and some related open-source codes became available in recent years [@Heagy2017; @Castillo-Reyes2018; @Rochlitz2019; @Werthmueller2019]. Objectives for designing these and other new electromagnetic forward modelling codes are to (i) enable modelling of complex subsurface settings [@Rochlitz2019], (ii) enhance the accuracy of forward responses while minimising computational costs [@Grayver2015], (iii) optimise and understand the effects of model discretisation [@Oldenburg2020; Werthmuller2021; @Castillo-Reyes2023] and (iv) to develop efficient solvers [@Weiss2023].


# Special code capabilities

In `elfe3D`, three special functionalities are combined. These features address the above-mentioned objectives i-iii and, in summary, distinguish `elfe3D` from other available codes:

1. **Variable Magnetic Permeability**: In traditional 3D electromagnetic forward modeling algorithms, magnetic permeability and dielectric permittivity are typically not considered variable model parameters. This oversight is often due to the relatively minor contrasts in these properties compared to conductivity contrasts prevalent in most geological environments at frequencies in the Hz- to kHz-range. However, it is well documented that mineralised zones and metallic infrastructure [@Heagy2023] may exhibit pronounced contrasts in magnetic permeability relative to their host rocks. Dielectric permittivities influence only high-frequency data in very resistive environments [@kals08]. Consequently, in `elfe3D`, magnetic permeability is treated as a variable model parameter, in conjunction with electrical resistivity. @Rulff2021 present a numerical test example with magnetic anomalies.

2. **Automatic Mesh Refinement**: The refinement approach implemented in `elfe3D` is inspired by @Ren2013 and @Grayver2015. It aims to achieve an optimal balance between the accuracy of the solutions at specified areas of interest and the overall size of the computational problem. This is accomplished through a goal-oriented mesh refinement strategy, specifically tailored for models characterised by variable electrical conductivity and magnetic permeability. The goal-oriented mesh refinement strategy including an adjoint problem formulation and the utilised error estimation approaches are detailed in @Rulff2021.

3. **Perfect Electric Conductors**: `elfe3D` incorporates a methodology for approximating highly conductive infrastructures within the model, such as steel-cased wells that are frequently encountered in exploration contexts. These structures are modeled as perfect electric conductors [@Um2020], utilising an approach similar to boundary conditions. Within the computational domain, this formulation enforces the electric field to zero at mesh edges corresponding to the perfect electric conductors. @Castillo-Reyes2023 provide details on the perfect electric conductor implementation and application.

# Target audience

The program `elfe3D` is designed for geophysicists who want to:

- Plan controlled-source electromagnetic experiments and require a tool for conducting numerical feasibility studies.
- Generate problem-specific tetrahedral meshes for controlled-source electromagnetic settings, which can be utilised in other simulation environments.
- Integrate the code into their inversion frameworks.
- Further develop the code.

# Workflow overview

Isotropic electric resistivities and magnetic permeabilities are variable, but element-wise constant, model parameters arranged in a tetrahedral mesh generated with `TetGen` [@Si2015]. The mesh files must be provided as input to `elfe3D` along with a start file specifying survey parameters such as frequencies, receiver and source positions, etc. Line or loop sources are represented along the edges of the mesh elements.

The physical behavior of electromagnetic fields is described by Maxwell's equations, which serve as the foundation for the curl-curl equation in terms of the total electric field $\bf{E}$ that `elfe3D` is solving:

$$\nabla \times \frac{1}{\mu} \nabla \times \mathbf{E} - i\omega \frac{1}{\rho} \mathbf{E} - \omega^2 \epsilon \mathbf{E}  = \color{black} i\omega \mathbf{J}_{p} \color{black} \quad \text{in} \quad \Omega,$$

where $\Omega$ is the computational domain, $\omega$ is the angular frequency, $\rho$ the electrical resistivity, $\varepsilon$ the dielectric permittivity and $\mu$ the magnetic permeability. A time dependence $e^{-i\omega t}$ is assumed and $\textbf{J}_{p}$ describes the source term.  Dirichlet boundary conditions are imposed on the outer domain boundaries.

This governing equation is discretised with linear vector finite-element shape functions [@Jin2014]. Developed using modern Fortran, `elfe3D` leverages vectorisation and shared-memory parallelism and employs a direct method, the `MUMPS` solver [@Amestoy2001], for solving the system of equations. At the end of the workflow [Figure \ref{fig:example}], the synthetic data, that is the electric and magnetic field components in Cartesian space, are calculated from the solution $\bf{E}$ and written to output files.


# Development history and ongoing research
   
The underlying code of `elfe3D` was developed between 2018-2023 and validated in @Rulff2021. This initial code is implemented in the inversion software `emilia` [@kals08; @kals10;@kals15] to enable 3D controlled-source electromagnetic inversion [@Rulff2023; @Rulff2024].
Adaptively refined meshes and parts of the synthetic data reported in @Castillo-Reyes2023 were generated with `elfe3D`.
The program is currently used and extended to design surface-to-borehole controlled-source electromagnetic surveys for geothermal applications [@Rulff2024emiw].

# Acknowledgements

I acknowledge contributions from Laura Maria Buntin and Thomas Kalscheuer to an earlier version of the code that `elfe3D` is based on.
The orginal code development was financed by the Smart Exploration project (European Union’s Horizon 2020 funding, grant agreement No. 775971).
My acknoledgements also go to the `MUMPS` and `TetGen` developers as well as to Dieter Werthmüller, who helped with making `elfe3D` open-source and kindly offered to host `elfe3D` in the `emsig` project.

# Availability
Version 1.0.0 of `elfe3D` is freely available under the Apache License, Version 2.0. The source code, along with documentation, an example and reference solutions, is hosted at https://github.com/emsig/elfe3D/tree/main. Further developments of `elfe3D` are ongoing. Collaboration and community feedback are welcome.



# References
