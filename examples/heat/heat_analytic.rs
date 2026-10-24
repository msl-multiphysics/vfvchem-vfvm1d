use vfvchem_vfvm1d::*;
use std::fs::create_dir_all;

/// Steady-state heat conduction in a uniformly heated water slab.
/// Run with: `cargo run --release --example heat_analytic`
///
/// Geometry:
/// Line domain defined by its end points.
/// The analytic solution is T(x) = 350 - 500 x + Q x (L - x) / (2 k).
///
///   bnd_l +--------- dom ---------+ bnd_r
///
/// Domains:
/// dom - main domain (0.1 m)
///
/// Boundaries:
/// bnd_l - left end
/// bnd_r - right end
///
/// Properties:
/// dom - thermal conductivity (k = 0.6 W m-1 K-1)
/// dom - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom - heat source (Q = 1.0e4 W m-3)
///
/// Boundary conditions:
/// bnd_l - temperature (T = 350 K)
/// bnd_r - temperature (T = 300 K)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heat_analytic").unwrap();

    // nodes, points, and variables
    // new_bounds - create boundary nodes and internal points between two end points
    // arguments: x_min, x_max, h_min, h_max, ref_len
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    let mut vars = Variables::new_bounds(0.0, 0.1, 0.004, 0.004, 0.01)?;

    // geometry
    // if using new_bounds, the domain is "dom"
    // the left and right boundaries are "bnd_l" and "bnd_r"
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., temp_0000.csv)
    let temp = vars.add_scldom_unkn(dom, 300.0, "examples/output_heat_analytic/temp.csv")?;

    // constant domain values
    // arguments: domain, value, output_file
    // heat production is multiplied by the temperature; heat source adds the source directly
    // use heat source by default; use heat production (or both) if heat source does not converge
    let cond = vars.add_scldom_cnst(dom, 0.6, "")?;  // thermal conductivity
    let hprd = vars.add_scldom_cnst(dom, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hsrc = vars.add_scldom_cnst(dom, 1.0e4, "")?;  // heat source (Q; positive if source; negative if sink)

    // constant boundary values
    // arguments: boundary, value, output_file
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 350.0, "")?;  // temperature
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0, "")?;  // temperature

    // steady-state heat transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_temp_bnd - register boundary with temperature
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom, temp, cond, hprd, hsrc)?;  // arguments: domain, T, k, s, Q
    phys.add_temp_bnd(bnd_l, temp_l)?;  // arguments: boundary, T
    phys.add_temp_bnd(bnd_r, temp_r)?;  // arguments: boundary, T

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-6, 1.0)?;
    Ok(())
}
