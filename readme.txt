__     _________     ______ _                    
\ \   / /  ___\ \   / / ___| |__   ___ _ __ ___  
 \ \ / /| |_   \ \ / / |   | '_ \ / _ \ '_ ` _ \ 
  \ V / |  _|   \ V /| |___| | | |  __/ | | | | |
   \_/  |_|      \_/  \____|_| |_|\___|_| |_| |_|

Vertex-Centered Finite Volume Method Solver for Chemical Engineering Applications
1D Vertex-Centered Finite Volume Method Module
Copyright (c) 2026 VFVChem Development Team

Overview
--------
VFVChem is a vertex-centered finite volume method (FVM) solver for heat, mass, and/or momentum transfer problems in chemical engineering.
It supports steady-state and transient simulations with non-constant properties and multiple material domains.
It is a pure Rust library. Simulations are configured as Rust code, and no separate installation or build system is required.

Note: The project is under active development. The public API and numerical models may change between releases.

Features
--------
- Steady or transient simulations of heat, mass, or momentum transport and couplings thereof
- Fixed or adaptive time stepping for transient simulations
- Non-constant inputs, such as temperature and concentration-dependent properties and reactions
- Multicomponent mass transfer with cross-diffusion
- Multiple domains, boundaries, and interfaces along a line
- Automatic placement of vertices with spacing graded from the boundaries and interfaces; no separate mesh generation is required
- Control volumes centered on the vertices with faces at the midpoints between neighboring vertices
- Power law scheme for advective transport
- Velocity from the continuity equation with an inlet velocity or mass flux; no momentum equation or pressure in 1D
- Direct LU linear solver and an iterative Jacobi-preconditioned GMRES solver
- Outputs simulation results as CSV or VTU files for post-processing

Requirements
------------
- A Rust toolchain with Rust 2024 edition support
- ParaView for visualizing output (optional)

Quick Start
-----------
Clone or download the project, change to the vfvchem_vfvm1d directory, and run `cargo run --release --example quickstart`.
This example solves a steady heat transfer problem on a line and writes the temperature field to `examples/output_quickstart/temp_0000.csv`.
Other examples can be found in the `examples` folder. An overview of the examples is provided in `examples/readme.txt`.

Adaptive Time Stepping
----------------------
Every transient solver can be run with adaptive time steps through `solve_adaptive` instead of `solve`:

    let mut opt = AdaptiveOptions::new(dt_init, t_end, "examples/output_name")?;
    opt.rtol = 1e-3;                   // relative error tolerance (default 1e-3)
    opt.out_times = vec![1.0, 2.5];    // optional; default is 11 evenly spaced times from 0 to t_end
    opt.scl_exclude = vec![vel];       // leave velocity out of the error estimate on flow problems
    phys.solve_adaptive(&mut vars, &mut lin, &opt)?;

Other fields of AdaptiveOptions set dt_min, dt_max, atol, the step controller, and the Picard iteration (max_iter, tol, damp).
Output files are written exactly at t = 0, at each time in out_times, and at t_end; steps are shortened to land on them.
Two files are written to the folder given to AdaptiveOptions::new:
- adaptive_times.csv : output file index and the corresponding time (e.g., temp_0003.csv is the row with file_index 3)
- adaptive_stats.csv : time, dt, and error estimate of every accepted step, with the number of accepted steps,
                       rejected steps, and Picard iterations repeated on each row
Velocity has no time derivative in 1D (it follows from continuity), so it should be excluded from the error estimate
with scl_exclude; it is still solved.
Time-dependent functions are evaluated at the end of each step. If a source or boundary condition changes suddenly
(e.g., a heat pulse from t0 to t1), add t0 and t1 to out_times; otherwise a change that falls inside one step is missed.

Convergence Scaling
-------------------
The Picard tolerance (tol) and the adaptive error tolerance (rtol) are measured per unknown field against a reference scale,
and the largest field is used. The reference scale is max|T|, max|c|, or max|u| of the field.
The reference scales can be given lower bounds, e.g., a reference velocity at very low flow rates:

    vars.set_scldom_scale(vel, 0.02)?;    // reference scale of a scalar

Setting `vars.scale_off = true` restores the unscaled measures (|x_new - x_old| / |x_old| for the Picard iteration).

Project Structure
-----------------
A guide to the source code is provided in `src/readme.txt`.

Contributing
------------
Due to limited maintenance capacity, we are not currently accepting external contributions.
However, users are welcome to fork and modify the repository in accordance with the MIT License.

License
-------
VFVChem is available under the MIT License. See license.txt for details.
