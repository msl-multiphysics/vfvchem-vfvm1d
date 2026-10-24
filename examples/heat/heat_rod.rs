use vfvchem_vfvm1d::*;
use std::fs::create_dir_all;

/// Steady-state heat conduction in a cooling fin (steel rod) with convective losses along its length.
/// Run with: `cargo run --release --example heat_rod`
///
/// Geometry:
/// Line domain defined by its end points.
/// The analytic solution is T = T_inf + (T_b - T_inf) cosh(m (L - x)) / cosh(m L), m = sqrt(h P / (k A)).
///
///   bnd_l +--------- dom ---------+ bnd_r
///
/// Domains:
/// dom - steel rod (L = 0.2 m, diameter d = 0.01 m)
///
/// Boundaries:
/// bnd_l - base of the fin
/// bnd_r - tip of the fin
///
/// Properties:
/// dom - thermal conductivity (k = 16 W m-1 K-1)
/// dom - convective loss h P / A (T_inf - T) with h = 25 W m-2 K-1, P / A = 4 / d, T_inf = 300 K
///       written as heat production (s = -h P / A) and heat source (Q = h P / A T_inf)
///
/// Boundary conditions:
/// bnd_l - temperature (T = 400 K)
/// bnd_r - no flux (q = 0.0 W m-2; insulated tip)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heat_rod").unwrap();

    // nodes, points, and variables
    // new_bounds - arguments: x_min, x_max, h_min, h_max, ref_len
    let mut vars = Variables::new_bounds(0.0, 0.2, 0.004, 0.004, 0.004)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;

    // unknown domain values
    let temp = vars.add_scldom_unkn(dom, 300.0, "examples/output_heat_rod/temp.csv")?;

    // constant domain values
    let hpa = 25.0 * 4.0 / 0.01;  // h P / A
    let cond = vars.add_scldom_cnst(dom, 16.0, "")?;  // thermal conductivity
    let hprd = vars.add_scldom_cnst(dom, -hpa, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hsrc = vars.add_scldom_cnst(dom, hpa * 300.0, "")?;  // heat source (Q; positive if source; negative if sink)

    // constant boundary values
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 400.0, "")?;  // temperature
    let hflx_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;  // heat flux (positive if outward; negative if inward)

    // steady-state heat transfer solver
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom, temp, cond, hprd, hsrc)?;  // arguments: domain, T, k, s, Q
    phys.add_temp_bnd(bnd_l, temp_l)?;  // arguments: boundary, T
    phys.add_hflx_bnd(bnd_r, hflx_r)?;  // arguments: boundary, q

    // physics solver
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-6, 1.0)?;
    Ok(())
}
