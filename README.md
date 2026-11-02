# FDM-cell-based-AMR-Cartesian
A third-order finite difference method on cell-based adaptive Cartesian grids for compressible flows

This repository provides a streamlined two-dimensional Fortran 90 implementation
(8 source files, ~2800 lines) of the algorithm presented in the associated paper.
The default configuration reproduces the **2D isentropic vortex evolution accuracy
test**, including uniform and center-refined
AMR meshes.

## Features

- Cell-based adaptive Cartesian mesh on a quadtree (fully threaded tree with
  cue-pointing neighbor pointers)
- Strict 2:1 mesh balancing with a two-cell transition band; mesh quality
  correction loop after each refinement/coarsening step
- Hanging-cell treatment by directional stencil selection
  and a multi-resolution WENO interpolation, applied in both fine-to-coarse and
  coarse-to-fine directions
- Third-order WCNS spatial discretization with local characteristic decomposition
- Explicit third-order TVD Runge--Kutta time integration
- OpenMP shared-memory parallelization at the cell level (dynamic scheduling)

## Code structure

| File | Contents |
|---|---|
| `src/MAIN.F90` | Program driver; the parameter routine `set_parameter` |
| `src/module.F90` | Precision constants and global data types (quadtree cell type) |
| `src/mesh.F90` | Mesh generation, block/solution-based refinement and coarsening, mesh quality correction |
| `src/neighbor.F90` | Recursive establishment and update of neighbor pointers |
| `src/XuanGua.F90` | Hanging-cell interpolation (directional stencils, multi-resolution WENO) |
| `src/nnd_center_RK.F90` | Third-order WCNS flux evaluation |
| `src/solution.F90` | Initialization, main time-integration loop, AMR driver, OpenMP loops |
| `src/output_ADT.F90` | Residual history and error-norm output |

## Compilation

Any standard Fortran 90 compiler with OpenMP support works, e.g.:

```bash
# Intel Fortran
ifort -qopenmp -O2 src/*.F90 -o amr2d

# GNU Fortran
gfortran -fopenmp -O2 src/*.F90 -o amr2d
```

Double precision is used throughout (set in `module.F90`).


## Running

```bash
./amr2d
```

Field results are written in Tecplot format, and the density error norms for the
accuracy test are appended to `accuracy_error.dat`.

## Parameter input

There is no external input file: all case settings are specified in the
parameter routine `set_parameter` in `src/MAIN.F90`. Edit the values there and
recompile. The meaning of each parameter is as follows.

### Domain and base grid

| Parameter | Meaning | Default |
|---|---|---|
| `Xmin, Xmax, Ymin, Ymax` | Computational domain | `[0,10] x [0,10]` |
| `h` | Base cell size (level 0) | `0.5` |
| `gama` | Ratio of specific heats | `1.4` |

### Mesh adaptation

| Parameter | Meaning | Default |
|---|---|---|
| `iflag_AMR` | Solution-adaptive refinement switch (`0`=off, `1`=on) | `0` |
| `lvlmax` | Global maximum refinement level | `0` |
| `lvlAMR` | Maximum level allowed for solution-adaptive refinement | `0` |
| `lvlblock` | Initial block-refinement level applied to `Block_size` at mesh creation | `0` |
| `Block_size(1:4)` | Refined block `[x1,x2] x [y1,y2]` | `[2.5,7.5] x [2.5,7.5]` |
| `Refine_goal_rot/div` | Refinement thresholds `r1, r2` (relative to the RMS values) | `1.4` |
| `Coarse_goal_rot/div` | Coarsening thresholds `r3, r4` | `0.4` |
| `Nstep_Adapt` | First step at which adaptation is performed | `0` |
| `Pstep_Adapt` | Adaptation interval (time steps) | `30` |

To reproduce the center-refined meshes of the paper's accuracy test, set
`lvlblock = 1` (or `2`) and `lvlmax = 1` (or `2`) with the default `Block_size`.

### Boundary conditions

| Parameter | Meaning | Default |
|---|---|---|
| `kind_xpbc, kind_xmbc, kind_ypbc, kind_ymbc` | Boundary types at x+/x-/y+/y- (`6` = periodic) | `6` |

### Time integration

| Parameter | Meaning | Default |
|---|---|---|
| `t_end` | Final physical time | `10.0` |
| `iflag_localdt` | Local time stepping switch (`0` = global step) | `0` |
| `tstep` | Fixed time step (used when `iflag_localdt = 0`) | `1.0e-4` |
| `CFL` | CFL number (used for adaptive stepping) | `0.45` |

### Parallelization and output

| Parameter | Meaning | Default |
|---|---|---|
| `Omp_num_threads` | Number of OpenMP threads | `4` |
| `Pstep_Res` | Output interval for the residual history | `100` |

## References

[1] X. Li, M. Yang, L. Bi, R. Xu, C. Luo, S. Yuan, X. Yuan, Z. Tang, An efficient Cartesian mesh generation strategy for complex geometries, *Comput. Methods Appl. Mech. Eng.* 418 (2024) 116564.

[2] S. Meng, D. Zhou, X. Yuan, X. Li, H. Chen, L. Bi, Enhanced strategy for adaptive Cartesian grid generation with arbitrarily complex 3D geometry, *Adv. Eng. Softw.* 174 (2022) 103304.

[3] S. Meng, D. Zhou, L. Bi, X. Li, C. Luo, H. Du, Anisotropic Cartesian grid generation strategy for arbitrarily complex geometry based on a fully threaded tree, *Comput.-Aided Des.* 155 (2023) 103447.

[4] C. Luo, L. Bi, S. Meng, H. Gao, X. Yuan, Z. Tang, A novel algorithm for fluid dynamics simulation with cell-adaptive Cartesian grids and immersed boundary method, *Phys. Fluids* 37 (2025) 046122.

## Citation

@article{mi2026third,
  title   = {A third-order finite difference method on cell-based adaptive
             Cartesian grids for compressible flows},
  author  = {Mi, Jian and Tang, Qinglin and Ma, Qiang and Yuan, Xianxu and Tang, Zhigong and Bi, Lin},
  journal = {Computer Physics Communications},
  volume  = {329},
  pages   = {110413},
  year    = {2026}
}
