use vfvchem_vfvm1d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Start-up of an adiabatic gas-phase tubular reactor with an exothermic first-order reaction A -> B.
/// Run with: `cargo run --release --example heatmassflow_reactor`
///
/// Geometry:
///
/// Line domain defined by its end points.
///
///   bnd_l +--------- dom ---------+ bnd_r
///
/// Domains:
/// dom - reactor (L = 1 m), initially filled with inert gas at 600 K
///
/// Boundaries:
/// bnd_l - inlet
/// bnd_r - outlet
///
/// Properties:
/// dom - ideal gas density rho = p M / (R T) (p = 1 atm, M = 0.029 kg mol-1); cp = 1100 J kg-1 K-1
/// dom - k = 0.05 W m-1 K-1 (with axial dispersion, effective 1.0 W m-1 K-1); rho D = 1e-3 kg m-1 s-1
/// dom - rate r = A0 exp(-Ea / (R T)) rho w_A (kg m-3 s-1), A0 = 1e8 s-1, Ea = 80 kJ mol-1
/// dom - heat of reaction dH = -5e5 J kg-1 (exothermic)
///
/// Boundary conditions:
/// bnd_l - inlet T = 600 K, w_A = 0.05, mass flux G = 0.5 kg m-2 s-1
/// bnd_r - outlet (zero diffusive flux)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heatmassflow_reactor").unwrap();

    // nodes, points, and variables
    // new_bounds - arguments: x_min, x_max, h_min, h_max, ref_len
    let mut vars = Variables::new_bounds(0.0, 1.0, 0.01, 0.01, 0.01)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;

    // unknown domain values
    let temp = vars.add_scldom_unkn(dom, 600.0, "examples/output_heatmassflow_reactor/temp.csv")?;
    let vel = vars.add_scldom_unkn(dom, 0.85, "examples/output_heatmassflow_reactor/vel.csv")?;
    let frac_a = vars.add_scldom_unkn(dom, 0.0, "examples/output_heatmassflow_reactor/frac_a.csv")?;

    // function domain values
    // rate constant times density; used by the mass production (implicit in w_A) and the heat source
    fn rho(t: f64) -> f64 { 101325.0 * 0.029 / (8.314 * t) }
    fn krho(t: f64) -> f64 { 1e8 * (-8.0e4 / (8.314 * t)).exp() * rho(t) }
    let den = vars.add_scldom_func(dom, Box::new(|_t, _x, u| rho(u[0])), vec![temp], "")?;
    let vlcp = vars.add_scldom_func(dom, Box::new(|_t, _x, u| rho(u[0]) * 1100.0), vec![temp], "")?;
    let mprd_a = vars.add_scldom_func(dom, Box::new(|_t, _x, u| -krho(u[0])), vec![temp], "")?;
    let hsrc = vars.add_scldom_func(dom, Box::new(|_t, _x, u| 5e5 * krho(u[0]) * u[1]), vec![temp, frac_a], "")?;

    // constant domain values
    let cond = vars.add_scldom_cnst(dom, 1.0, "")?;
    let diff = vars.add_scldom_cnst(dom, 1e-3, "")?;
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;

    // boundary values
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 600.0, "")?;
    let frac_a_l = vars.add_sclbnd_cnst(bnd_l, 0.05, "")?;
    let mfr_l = vars.add_sclbnd_cnst(bnd_l, 0.5, "")?;

    // transient heat transfer, mass transfer, and flow solver
    let mut phys = TransientHeatMassFlow::new(1)?;
    phys.add_heat_dom(dom, temp, vlcp, cond, zero, hsrc)?;
    phys.add_temp_bnd(bnd_l, temp_l)?;
    phys.add_hout_bnd(bnd_r)?;
    phys.add_mass_dom(0, dom, frac_a, HashMap::from([(0, diff)]), mprd_a, zero)?;
    phys.add_conc_bnd(0, bnd_l, frac_a_l)?;
    phys.add_mout_bnd(0, bnd_r)?;
    phys.add_flow_dom(dom, vel, den)?;
    phys.add_mfr_bnd(bnd_l, mfr_l)?;

    // adaptive physics solver
    // AdaptiveOptions::new - arguments: dt_init, t_end, folder for adaptive_times.csv and adaptive_stats.csv
    let mut opt = AdaptiveOptions::new(0.01, 10.0, "examples/output_heatmassflow_reactor")?;
    opt.rtol = 1e-3;
    opt.dt_max = 0.5;
    opt.tol = 1e-8;
    opt.out_times = vec![0.5, 1.0, 2.0, 5.0, 10.0];
    let mut lin = LinearLu::new(1)?;
    phys.solve_adaptive(&mut vars, &mut lin, &opt)?;
    Ok(())
}
