use vfvchem_vfvm1d::*;
use std::fs::create_dir_all;

/// Steady-state continuity equation for flow through two regions with different densities joined by an interface.
/// Only continuity is solved: d(rho u)/dx = 0 (no momentum equation).
/// Run with: `cargo run --release --example flow_tworeg`
///
/// Geometry:
/// Line split into two domains joined by an interface.
/// The flow goes from right to left; the velocity jumps at the interface so that rho u is continuous.
///
///   bnd_l +---- dom_0 ----+ itf_0 +---- dom_1 ----+ bnd_r
///
/// Domains:
/// dom_0 - left region (0.5 m)
/// dom_1 - right region (0.5 m)
///
/// Boundaries:
/// bnd_l - outlet
/// bnd_r - inlet
///
/// Interfaces:
/// itf_0 - interface between dom_0 and dom_1
///
/// Properties:
/// dom_0 - density (rho = 800.0 kg m-3)
/// dom_1 - density (rho = 1000.0 kg m-3)
///
/// Boundary conditions:
/// bnd_l - outlet (continuity only)
/// bnd_r - mass flow rate per unit area (G = 2.0 kg m-2 s-1; positive into the domain)
///
/// Interface conditions:
/// itf_0 - mass flux continuity
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_flow_tworeg").unwrap();

    // nodes, points, and variables
    // new_segments - create nodes at the given positions and internal points between them
    // arguments: node_positions, h_min, h_max, ref_len
    // h_min - point spacing at the nodes; h_max - point spacing far from the nodes
    // ref_len - distance from the nodes over which the spacing grows from h_min to h_max
    let mut vars = Variables::new_segments(&[0.0, 0.5, 1.0], 0.02, 0.02, 0.02)?;

    // geometry
    // if using new_segments, the domains are "dom_0", "dom_1", ... from left to right
    // the boundaries are "bnd_l" and "bnd_r"; the interfaces are "itf_0", "itf_1", ... from left to right
    // add_itf - arguments: domain_a, domain_b, group_name
    let dom_0 = vars.add_dom("dom_0")?;
    let dom_1 = vars.add_dom("dom_1")?;
    let _bnd_l = vars.add_bnd(dom_0, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom_1, "bnd_r")?;
    let itf_0 = vars.add_itf(dom_0, dom_1, "itf_0")?;

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., vel_0_0000.csv)
    let vel_0 = vars.add_scldom_unkn(dom_0, 0.0, "examples/output_flow_tworeg/vel_0.csv")?;  // velocity (+x component)
    let vel_1 = vars.add_scldom_unkn(dom_1, 0.0, "examples/output_flow_tworeg/vel_1.csv")?;  // velocity (+x component)

    // constant domain values
    // arguments: domain, value, output_file
    let den_0 = vars.add_scldom_cnst(dom_0, 800.0, "")?;  // density
    let den_1 = vars.add_scldom_cnst(dom_1, 1000.0, "")?;  // density

    // constant boundary values
    // arguments: boundary, value, output_file
    let mfr_r = vars.add_sclbnd_cnst(bnd_r, 2.0, "")?;  // mass flow rate per unit area (positive into the domain)

    // steady-state flow solver
    // add_flow_dom - register domain with flow
    // add_mfr_bnd - register inlet boundary with mass flow rate per unit area
    // add_fcnt_itf - register interface with mass flux continuity
    // the flow goes from the inlet to the other end of the chain of domains; the outlet needs no condition
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom_0, vel_0, den_0)?;  // arguments: domain, u, rho
    phys.add_flow_dom(dom_1, vel_1, den_1)?;  // arguments: domain, u, rho
    phys.add_mfr_bnd(bnd_r, mfr_r)?;  // arguments: boundary, G
    phys.add_fcnt_itf(itf_0)?;  // arguments: interface

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-8, 1.0)?;
    Ok(())
}
