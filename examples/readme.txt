__     _________     ______ _                    
\ \   / /  ___\ \   / / ___| |__   ___ _ __ ___  
 \ \ / /| |_   \ \ / / |   | '_ \ / _ \ '_ ` _ \ 
  \ V / |  _|   \ V /| |___| | | |  __/ | | | | |
   \_/  |_|      \_/  \____|_| |_|\___|_| |_| |_|

Vertex-Centered Finite Volume Method Solver for Chemical Engineering Applications
1D Vertex-Centered Finite Volume Method Module - Examples

How to Use Examples
-------------------
First-time users are encouraged to run quickstart.rs

Change the directory to vfvchem_vfvm1d then run quickstart.rs:
`cargo run --release --example quickstart`

The other examples are registered in `Cargo.toml` and can be run in place:
`cargo run --release --example name_of_file`

Output is written to `examples/output_name_of_file/` as CSV files (one file per output index).

List of Examples
----------------
Below is a list of examples and the concepts they introduce:

Heat Transfer
- quickstart.rs     : introduction to VFVChem code
- heat_transient.rs : transient heat transfer simulation
- heat_tworeg.rs    : multiple domains joined by an interface; time-dependent boundary conditions
- heat_analytic.rs  : steady-state conduction in a uniformly heated slab
- heat_rod.rs       : cooling fin with convective losses as production and source terms

Mass Transfer
- mass_slab.rs     : steady-state multicomponent diffusion with cross-diffusion
- mass_analytic.rs : steady-state diffusion with a first-order reaction using a production coefficient

Momentum Transfer
- flow_analytic.rs : steady flow with a prescribed density profile; u = G / rho(x)
- flow_tworeg.rs   : flow through two domains joined by an interface
- flow_heated.rs   : transient flow of an ideal gas that expands as it is heated; mass flux inlet

Heat and Momentum Transfer
- heatflow_tube.rs : air heated by a hot wall with an ideal-gas density

Mass and Momentum Transfer
- massflow_reactor.rs : axial-dispersion reactor with a second-order reaction

Heat, Mass, and Momentum Transfer
- heatmassflow_reactor.rs : transient start-up of an adiabatic exothermic reactor with an Arrhenius rate

Features
- feat_prod.rs : production coefficient versus mass source for a fast first-order reaction
- feat_itf.rs  : heat and mass transfer through regions with different properties joined by an interface
- feat_adt.rs  : adaptive time stepping
