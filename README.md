# ϵlastσPlasm.jl 👻
[![Build Status](https://github.com/ewyser/ElastoPlasm.jl/workflows/CI/badge.svg)](https://github.com/ewyser/ElastoPlasm.jl/actions)
[![Documentation](https://github.com/ewyser/ElastoPlasm.jl/actions/workflows/docs.yaml/badge.svg)](https://github.com/ewyser/ElastoPlasm.jl/actions/workflows/docs.yaml)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ewyser.github.io/ElastoPlasm.jl/)
[![](https://img.shields.io/badge/NVIDIA-CUDA-green.svg?logo=nvidia)](https://developer.nvidia.com/cuda-toolkit)
[![](https://img.shields.io/badge/AMD-ROCm-red.svg?logo=amd)](https://www.amd.com/en/products/software/rocm.html)
<!--
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaci.github.io/PkgTemplates.jl/stable)
[![](https://img.shields.io/badge/docs-stable-blue.svg?logo=quicklook)](https://github.com/LandslideSIM/MaterialPointSolver.jl/wiki)
[![](https://img.shields.io/badge/version-v0.3.0-926116)]()

[![](https://img.shields.io/badge/Intel-oneAPI-blue.svg?logo=intel)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
[![](https://img.shields.io/badge/Apple-Metal-purple.svg?logo=apple)](https://developer.apple.com/metal/)
--->

## Overview
This package is an evolution of the somewhat cumbersome-to-use [`ep2-3De v1.0`](https://github.com/ewyser/ep2-3De), and is entirely written in **Julia**. It is designed for **rapid prototyping** while maintaining **reasonable production capabilities**. It addresses the following key aspects:

- **Updated Lagrangian explicit formulation** for elastoplastic simulations.
- Supports both **finite** and **infinitesimal deformation** frameworks:
  - **Finite deformation**: employs logarithmic strains and Kirchhoff stresses.
  - **Infinitesimal deformation**: based on a **Jaumann rate** formulation.
- Compatible with multiple **shape function bases**:
    - Standard linear shape function $N_n(\boldsymbol{x}_p)$
    - GIMP shape function $S_n(\boldsymbol{x}_p)$
    - Boundary-modified cubic B-spline shape function $\phi_n(\boldsymbol{x}_p)$
- Provides mappings between nodes (denoted $n$ or $v$) and material points (denoted $p$), using:
    - FLIP with augmented mUSL procedure
    - TPIC with standard USL procedure

## How to ```plasming``` ?  

0. (opt.) Get Julia [here](https://julialang.org/downloads/) and follow instructions for installation

1. Clone [```ElastoPlasm.jl```](https://github.com/ewyser/elastoPlasm.jl/tree/main)  and ```cd``` to your local repo 

2. Launch Julia (on macOS, drag & drop ```start_macOS.sh``` in the terminal) and enter pkg mode ``` ] ```, then ```activate .``` the project ```ElastoPlasm``` and ```instantiate``` its environment and related packages.

4. Once ```ElastoPlasm``` has been correctly instantiated, you can ```using ElastoPlasm```

