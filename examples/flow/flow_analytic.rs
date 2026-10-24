use vfvchem_vfvm1d::*;
use std::fs::create_dir_all;

/// Steady-state continuity equation for a gas whose density falls along a pipe.
/// Only continuity is solved: d(rho u)/dx = 0 (no momentum equation).
/// Run with: `cargo run --release --example flow_analytic`
///
/// Geometry:
/// Line domain defined by its end points.
/// The analytic solution is u(x) = rho_in u_in / rho(x).
///
///   bnd_l +--------- dom ---------+ bnd_r
///
/// Domains:
/// dom - pipe (1.0 m)
///
/// Boundaries:
/// bnd_l - inlet
/// bnd_r - outlet
///
/// Properties:
/// [non-constant properties in SI units]
/// dom - density (kg m-3)
/// rho = 1.2 / (1 + x)
///
/// Boundary conditions:
/// bnd_l - velocity (u = 2.0 m s-1)
/// bnd_r - outlet (continuity only)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_flow_analytic").unwrap();

    // nodes, points, and variables
    // new_bounds - create boundary nodes and internal points between two end points
    // arguments: x_min, x_max, h_min, h_max, ref_len
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    let mut vars = Variables::new_bounds(0.0, 1.0, 0.01, 0.05, 0.2)?;

    // geometry
    // if using new_bounds, the domain is "dom"
    // the left and right boundaries are "bnd_l" and "bnd_r"
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let _bnd_r = vars.add_bnd(dom, "bnd_r")?;

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., vel_0000.csv)
    let vel = vars.add_scldom_unkn(dom, 2.0, "examples/output_flow_analytic/vel.csv")?;  // velocity (+x component)

    // non-constant domain values
    // arguments: domain, value_func, unknown_ids, output_file
    // value_func - returns value as a function of time, position x, and *unknown* values
    // - unknown values are given in the same order as in unknown_ids
    let den_func: SclDomFunc = Box::new(|_t, x, _unkn| { 1.2 / (1.0 + x) });
    let den = vars.add_scldom_func(dom, den_func, vec![], "examples/output_flow_analytic/den.csv")?;  // density

    // constant boundary values
    // arguments: boundary, value, output_file
    let vel_l = vars.add_sclbnd_cnst(bnd_l, 2.0, "")?;  // velocity (+x component)

    // steady-state flow solver
    // add_flow_dom - register domain with flow
    // add_vel_bnd - register inlet boundary with velocity
    // the flow goes from the inlet to the other end of the domain; the outlet needs no condition
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom, vel, den)?;  // arguments: domain, u, rho
    phys.add_vel_bnd(bnd_l, vel_l)?;  // arguments: boundary, u

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-8, 1.0)?;
    Ok(())
}
