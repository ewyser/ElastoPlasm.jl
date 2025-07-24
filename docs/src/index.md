```@meta
CurrentModule = ElastoPlasm
```

# Documentation 

```@index
```

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

The solver can generate initial fields $f(\boldsymbol{x})$—such as the cohesion $c(\boldsymbol{x}_p)$ or internal friction angle $\varphi(\boldsymbol{x}_p)$—using random Gaussian fields, with $\boldsymbol{x}_p$ representing a material point’s coordinate.

<p align="center">
  <img src="./assets/img/epII.png" alt="Plastic strain" style="width:80%;"/>>
  <br/>
  <em>Figure: Slumping dynamics (without volumetric locking corrections) showing the accumulated plastic strain $\epsilon_p^{\mathrm{acc}}$ after an elastic load of 8 s and an additional elasto-plastic load of ≈ 7 s.</em>
</p>

<p align="center">
  <img src="./assets/img/c0.png" alt="Initial cohesion field" style="width:80%;"/>
  <br/>
  <em>Figure: Initial cohesion field \( c_0(\boldsymbol{x}_p) \) with average \( \mu = 20 \,\text{kPa} \) and variance \( \sigma = \pm 5 \,\text{kPa} \).</em>
</p>

## **Content**
1. [Usage](#id-section2)
<div id='id-section2'/>
