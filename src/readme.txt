__     _________     ______ _                    
\ \   / /  ___\ \   / / ___| |__   ___ _ __ ___  
 \ \ / /| |_   \ \ / / |   | '_ \ / _ \ '_ ` _ \ 
  \ V / |  _|   \ V /| |___| | | |  __/ | | | | |
   \_/  |_|      \_/  \____|_| |_|\___|_| |_| |_|

Vertex-Centered Finite Volume Method Solver for Chemical Engineering Applications
1D Vertex-Centered Finite Volume Method Module - Source Code

Overview
--------
VFVChem is a vertex-centered finite volume method (FVM) solver for heat, mass, and/or momentum transfer problems in chemical engineering.
It supports steady-state and transient simulations with non-constant properties and multiple material domains.

Nodes are placed at the ends of domains (boundaries and interfaces), and points are placed inside domains
with a spacing that grows from h_min at the nodes to h_max at a distance ref_len.
A control volume is built around every vertex, with faces at the midpoints between neighboring vertices.
Currently, only backward Euler time stepping is implemented, with constant or adaptive time steps.
Advective transport is discretized via the power law scheme.
There is no momentum equation in 1D; the velocity follows from the continuity equation and an inlet condition.

The source code is organized into the following folders:
- base      : core data structures (e.g., geometry, vertices, scalars, writers)
- linear    : linear matrix solvers (e.g., LU, GMRES)
- operator  : applies discretized conservation equations into global matrix equation
- steady    : steady-state solvers
- transient : transient solvers

The operator, steady, and transient folders are split into two subfolders:
- core   : heat, mass, and momentum transfer in single-phase materials
- extend : additional physics built on the core (none yet in 1D)
The base classes (operator_base, steady_base, transient_base, transient_adaptive)
and the preludes stay at the top of each folder. The preludes re-export both
core and extend types, so `use vfvchem_vfvm1d::*` covers both.

Base
----
- error    : defines errors returned by the library
- geom_bnd : stores 0D boundary geometry
- geom_dom : stores 1D domain geometry
- geom_itf : stores 0D interface geometry; interfaces join two domains
- geom_seg : stores the line segments and control volume data of a domain
- nds      : stores nodes on boundaries and interfaces
- pts      : stores points inside domains; placed with graded spacing
- scaling  : per-field reference scales for the Picard residual and the adaptive error norm
- scl_bnd  : stores a scalar field defined on a boundary
- scl_dom  : stores a scalar field defined on a domain
- scl_flw  : stores the face mass flux shared by the transport operators; filled by the solvers
- scl_itf  : stores a scalar field defined on an interface
- vars     : owns all geometry and scalar data
- write_csv: writes scalar values to CSV files
- write_vtu: writes scalar values to VTU files

Vertices are either points or nodes. Within a domain, vertex ids are ordered with
points first, then nodes, i.e., vid = pid for points and vid = num_pts + nid for nodes.
Points are numbered in order of x, and the nodes of a domain are its left (nid = 0) and right (nid = 1) ends.

Geometries are created through Variables helpers:
- new_bounds   : one domain "dom" from x_min to x_max with boundaries "bnd_l" and "bnd_r"
- new_nodes    : one domain "dom" with given vertex positions and boundaries "bnd_l" and "bnd_r"
- new_segments : domains "dom_0", "dom_1", ... between given node positions, with boundaries "bnd_l" and "bnd_r"
                 at the ends and interfaces "itf_0", "itf_1", ... between consecutive domains

Field types are registered through Variables helpers:
- *_cnst    : constant value
- *_func    : function of time, position x, and unknown scalar values;
              interface functions receive time and unknown values from both domains
- *_unkn    : unknown field (solved for)

Flow data is registered through add_sclflw and is used by the transport operators so that
all of them use the same face mass flux. SteadyBase and TransientBase call update_func and
update_flow before each matrix assembly so that functions and face mass fluxes stay consistent
with the unknowns. update_flow is a required method of each solver: flow solvers call
flow_helper::update_flow, and solvers without flow leave it empty. Unlike 2D, update_flow
takes the time step, since the face mass flux of a transient flow includes the change of density.

Operator
--------
Operator names use the form op + field type + location + term. For example,
opscl_dom_diff is a scalar domain diffusion operator. Files ending in unity
are variants whose weighting coefficient is fixed at one.
Domain operators apply to the points of the domain, boundary operators apply to the
control volumes of boundary nodes, and interface operators apply to the
combined control volumes of interface nodes.
Every boundary of a domain needs a boundary condition (e.g., add_hflx_bnd with zero for an insulated end).

Helpers hold code shared by several solvers. They are not re-exported by the
prelude and are called by path (e.g., flow_helper::update_flow).

- operator_base       : defines the common operator trait and global matrix assembly helpers
- prelude             : re-exports operator types for convenient internal imports

Core operators (operator/core)
- flow_helper         : builds the flow chains from the inlets and computes the face mass flux of a free fluid
- opscl_bnd_adv       : adds weighted advection to the control volumes of boundary nodes
- opscl_bnd_advunity  : adds advection with unit weight to the control volumes of boundary nodes
- opscl_bnd_diff      : adds diffusion to the control volumes of boundary nodes
- opscl_bnd_dir       : enforces a scalar Dirichlet boundary condition
- opscl_bnd_mfr       : enforces an inlet mass flux
- opscl_bnd_neu       : applies a prescribed scalar outward flux
- opscl_bnd_prod      : adds a production term (s * c) to the control volumes of boundary nodes
- opscl_bnd_src       : adds a source to the control volumes of boundary nodes
- opscl_bnd_time      : adds a weighted time derivative to the control volumes of boundary nodes
- opscl_bnd_timeprev  : adds a time derivative weighted by the previous coefficient to the control volumes of boundary nodes
- opscl_bnd_timeunity : adds a time derivative with unit weight to the control volumes of boundary nodes
- opscl_dom_adv       : adds weighted advection to a scalar transport equation
- opscl_dom_advunity  : adds advection with unit weight to a scalar transport equation
- opscl_dom_diff      : adds diffusion to a scalar transport equation
- opscl_dom_div       : adds the continuity equation of each segment to the row of its downstream vertex
- opscl_dom_prod      : adds a production term (s * c) to a scalar transport equation
- opscl_dom_src       : adds a source to a scalar transport equation
- opscl_dom_time      : adds a weighted time derivative to a scalar transport equation
- opscl_dom_timeprev  : adds a time derivative weighted by the previous coefficient to a scalar transport equation
- opscl_dom_timeunity : adds a time derivative with unit weight to a scalar transport equation
- opscl_itf_adv       : adds weighted advection to the combined control volumes of interface nodes
- opscl_itf_advunity  : adds advection with unit weight to the combined control volumes of interface nodes
- opscl_itf_cont      : enforces scalar continuity between two domains
- opscl_itf_diff      : adds diffusion to the combined control volumes of interface nodes
- opscl_itf_div       : enforces mass flux continuity between two domains with flow
- opscl_itf_prod      : adds a production term (s * c) to the combined control volumes of interface nodes
- opscl_itf_src       : adds a source to the combined control volumes of interface nodes
- opscl_itf_time      : adds a weighted time derivative to the combined control volumes of interface nodes
- opscl_itf_timeprev  : adds a time derivative weighted by the previous coefficient to the combined control volumes of interface nodes
- opscl_itf_timeunity : adds a time derivative with unit weight to the combined control volumes of interface nodes

Extended operators (operator/extend)
There are no extended operators in 1D yet.

Linear
------
- linear_base  : defines the common linear solver trait
- linear_gmres : solves sparse systems with the Jacobi-preconditioned restarted GMRES method
- linear_lu    : solves sparse systems with an LU factorization

Steady
------
- steady_base         : implements the common steady nonlinear solution loop

Core solvers (steady/core)
- steady_flow         : solves steady momentum transfer (continuity)
- steady_heat         : solves steady heat transfer
- steady_heatflow     : solves coupled steady heat and momentum transfer
- steady_heatmassflow : solves coupled steady heat, mass, and momentum transfer
- steady_mass         : solves steady mass transfer
- steady_massflow     : solves coupled steady mass and momentum transfer

Extended solvers (steady/extend)
There are no extended solvers in 1D yet.

Transient
---------
- transient_base         : implements the common time-stepping and nonlinear solution loops
- transient_adaptive     : implements adaptive time stepping for every transient solver

Core solvers (transient/core)
- transient_flow         : solves transient momentum transfer (continuity)
- transient_heat         : solves transient heat transfer
- transient_heatflow     : solves coupled transient heat and momentum transfer
- transient_heatmassflow : solves coupled transient heat, mass, and momentum transfer
- transient_mass         : solves transient mass transfer
- transient_massflow     : solves coupled transient mass and momentum transfer

Extended solvers (transient/extend)
There are no extended solvers in 1D yet.

Flow solvers
------------
The flow solvers solve the continuity equation d(rho)/dt + d(rho * u)/dx = 0 for the velocity u;
the density is given (constant, or a function of the other unknowns, e.g., the ideal gas law).
Each domain is registered through add_flow_dom (u, rho).

Flow goes along chains of domains joined end to end, starting at an inlet boundary:
- add_vel_bnd  : prescribed inlet velocity (+x component)
- add_mfr_bnd  : prescribed inlet mass flux (per unit area); positive into the domain
- add_fcnt_itf : flow continuity (rho * u continuous; u may jump if rho does)
The outlet (the other end of a chain) needs no flow condition. Every flow domain must be reached
from exactly one inlet. The continuity equation of each segment goes to the row of its downstream
vertex, so the inlet row holds the inlet condition and the interface row of the downstream domain
holds the mass flux continuity.

The face mass flux used by the transport operators is marched from the inlet so that every control
volume satisfies continuity exactly: F_out = F_in - V * (rho - rho_prev) / dt.
The transport operators use the non-conservative (advective) form with this mass flux, so the
boundary half of a control volume is upwinded from the node; heat and mass leave an outlet by
advection only (add_hout_bnd and add_mout_bnd).
