# ϵlastσPlasm.jl 👻

[![Build Status](https://github.com/ewyser/ElastoPlasm.jl/workflows/CI/badge.svg)](https://github.com/ewyser/ElastoPlasm.jl/actions)
[![Documentation](https://github.com/ewyser/ElastoPlasm.jl/actions/workflows/docs.yaml/badge.svg)](https://github.com/ewyser/ElastoPlasm.jl/actions/workflows/docs.yaml)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ewyser.github.io/ElastoPlasm.jl/stable)

<!--
[![](https://img.shields.io/badge/NVIDIA-CUDA-green.svg?logo=nvidia)](https://developer.nvidia.com/cuda-toolkit)
[![](https://img.shields.io/badge/AMD-ROCm-red.svg?logo=amd)](https://www.amd.com/en/products/software/rocm.html)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaci.github.io/PkgTemplates.jl/stable)
[![](https://img.shields.io/badge/docs-stable-blue.svg?logo=quicklook)](https://github.com/LandslideSIM/MaterialPointSolver.jl/wiki)
[![](https://img.shields.io/badge/version-v0.3.0-926116)]()

[![](https://img.shields.io/badge/Intel-oneAPI-blue.svg?logo=intel)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
[![](https://img.shields.io/badge/Apple-Metal-purple.svg?logo=apple)](https://developer.apple.com/metal/)
--->

---

## ✨ Overview

**ϵlastσPlasm.jl** is the spiritual successor to the somewhat cumbersome [`ep2-3De v1.0`](https://github.com/ewyser/ep2-3De), now reborn in pure **Julia** 🦑. Designed for **rapid prototyping** and **serious production**, it brings:

- **Updated Lagrangian explicit formulation** for elastoplastic simulations.
- **Finite** and **infinitesimal deformation** frameworks:
  - **Finite deformation**: logarithmic strains & Kirchhoff stresses.
  - **Infinitesimal deformation**: Jaumann rate formulation.
- Multiple **shape function bases**:
    - Linear: $N_n(\boldsymbol{x}_p)$
    - GIMP: $S_n(\boldsymbol{x}_p)$
    - Cubic B-spline: $\phi_n(\boldsymbol{x}_p)$
- Mapping scheme for Node ($n$ or $v$) ↔ Material point ($p$):
    - PIC-FLIP
    - TPIC-TFLIP
    - APIC-AFLIP
- Update scheme
    - standard USL
    - augmented mUSL

---

## 🛠️ Installation

1. **Install Julia:** [Download here](https://julialang.org/downloads/) and follow the instructions.
2. **Clone the repo:**
   ```sh
   git clone https://github.com/ewyser/ElastoPlasm.jl.git
   cd ElastoPlasm.jl
   ```
3. **Start Julia** (on macOS, drag & drop `start_macOS.sh` in the terminal).
   ```sh
   julia --project=. 
   ```
4. **Quick Start:**
    ```julia
    using Pkg
    Pkg.instantiate()
    julia> using ElastoPlasm
    ┌ Welcome to ϵlastσPlasm 👻 v0.4.2
    │ New comer ? Try slumping out
    │   L,nel  = [64.1584,64.1584/4.0],[40,10];
    │   jld2   = slump_problem(L,nel;cli()...);
    └   out    = elastoplasm(jld2; workflow = [elastodynamic!,elastoplastic!]);

    julia> 
    # plasming begins here!
    ```

---

## ✍️ Contributing

Contributions, issues, and feature requests are welcome  
Feel free to check the [issues page](https://github.com/ewyser/ElastoPlasm.jl/issues), the [project page](https://github.com/users/ewyser/projects/2) or open a [pull request](https://github.com/ewyser/ElastoPlasm.jl/pulls).

---

## 🧠 Philosophy

ϵlastσPlasm.jl aims to be a fun sandbox to use. Happy plasming! 🎉

---

