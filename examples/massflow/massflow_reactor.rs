use vfvchem_vfvm1d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Steady-state isothermal plug-flow reactor with axial dispersion: 2 A -> B (second order in A).
/// Run with: `cargo run --release --example massflow_reactor`
///
/// Geometry:
///
/// Line domain defined by its end points.
///
///   bnd_l +--------- dom ---------+ bnd_r
///
/// Domains:
/// dom - reactor (L = 2 m)
///
/// Boundaries:
/// bnd_l - inlet
/// bnd_r - outlet
///
/// Properties:
/// dom - liquid density (rho = 1000 kg m-3)
/// dom - dispersion rho D_ax (D_ax = 1e-3 m2 s-1)
/// dom - reaction rate r = k rho w_A^2 (kg m-3 s-1), k = 0.5 s-1
///       A: mass production coefficient (-k rho w_A); B: mass source (+r)
///
/// Boundary conditions:
/// bnd_l - inlet mass fractions (w_A = 0.2, w_B = 0) and velocity (u = 0.05 m s-1)
/// bnd_r - outlet (zero diffusive mass flux)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_massflow_reactor").unwrap();

    // nodes, points, and variables
    // new_bounds - arguments: x_min, x_max, h_min, h_max, ref_len
    let mut vars = Variables::new_bounds(0.0, 2.0, 0.02, 0.02, 0.02)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;

    // unknown domain values
    let vel = vars.add_scldom_unkn(dom, 0.05, "")?;
    let frac_a = vars.add_scldom_unkn(dom, 0.2, "examples/output_massflow_reactor/frac_a.csv")?;
    let frac_b = vars.add_scldom_unkn(dom, 0.0, "examples/output_massflow_reactor/frac_b.csv")?;

    // domain values
    let den = vars.add_scldom_cnst(dom, 1000.0, "")?;
    let diff = vars.add_scldom_cnst(dom, 1000.0 * 1e-3, "")?;  // rho D (per unit mass with flow)
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;
    let mprd_a = vars.add_scldom_func(dom, Box::new(|_t, _x, u| -0.5 * 1000.0 * u[0]), vec![frac_a], "")?;
    let msrc_b = vars.add_scldom_func(dom, Box::new(|_t, _x, u| 0.5 * 1000.0 * u[0] * u[0]), vec![frac_a], "")?;

    // boundary values
    let vel_l = vars.add_sclbnd_cnst(bnd_l, 0.05, "")?;
    let frac_a_l = vars.add_sclbnd_cnst(bnd_l, 0.2, "")?;
    let frac_b_l = vars.add_sclbnd_cnst(bnd_l, 0.0, "")?;

    // steady-state mass transfer and flow solver
    let mut phys = SteadyMassFlow::new(2)?;
    phys.add_mass_dom(0, dom, frac_a, HashMap::from([(0, diff)]), mprd_a, zero)?;
    phys.add_mass_dom(1, dom, frac_b, HashMap::from([(1, diff)]), zero, msrc_b)?;
    phys.add_conc_bnd(0, bnd_l, frac_a_l)?;
    phys.add_conc_bnd(1, bnd_l, frac_b_l)?;
    phys.add_mout_bnd(0, bnd_r)?;
    phys.add_mout_bnd(1, bnd_r)?;
    phys.add_flow_dom(dom, vel, den)?;
    phys.add_vel_bnd(bnd_l, vel_l)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-8, 1.0)?;
    Ok(())
}
