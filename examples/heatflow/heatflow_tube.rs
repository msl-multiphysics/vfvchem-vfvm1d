use vfvchem_vfvm1d::*;
use std::fs::create_dir_all;

/// Steady-state heating of air flowing through a tube with a hot wall.
/// Run with: `cargo run --release --example heatflow_tube`
///
/// Geometry:
///
/// Line domain defined by its end points.
///
///   bnd_l +--------- dom ---------+ bnd_r
///
/// Domains:
/// dom - tube (L = 1 m, diameter d = 0.02 m)
///
/// Boundaries:
/// bnd_l - inlet
/// bnd_r - outlet
///
/// Properties:
/// dom - air, ideal gas density rho = p M / (R T); cp = 1005 J kg-1 K-1; k = 0.03 W m-1 K-1
/// dom - wall heat transfer h P / A (T_w - T) with h = 20 W m-2 K-1, P / A = 4 / d, T_w = 500 K
///
/// Boundary conditions:
/// bnd_l - inlet temperature (T = 300 K) and mass flux (G = 2 kg m-2 s-1)
/// bnd_r - outlet (zero diffusive heat flux)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heatflow_tube").unwrap();

    // nodes, points, and variables
    // new_bounds - arguments: x_min, x_max, h_min, h_max, ref_len
    let mut vars = Variables::new_bounds(0.0, 1.0, 0.01, 0.01, 0.01)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;

    // unknown domain values
    let temp = vars.add_scldom_unkn(dom, 300.0, "examples/output_heatflow_tube/temp.csv")?;
    let vel = vars.add_scldom_unkn(dom, 1.0, "examples/output_heatflow_tube/vel.csv")?;

    // function domain values
    let den = vars.add_scldom_func(dom, Box::new(|_t, _x, u| 101325.0 * 0.029 / (8.314 * u[0])), vec![temp], "")?;
    let vlcp = vars.add_scldom_func(dom, Box::new(|_t, _x, u| 101325.0 * 0.029 / (8.314 * u[0]) * 1005.0), vec![temp], "")?;

    // constant domain values
    let hpa = 20.0 * 4.0 / 0.02;  // h P / A
    let cond = vars.add_scldom_cnst(dom, 0.03, "")?;
    let hprd = vars.add_scldom_cnst(dom, -hpa, "")?;
    let hsrc = vars.add_scldom_cnst(dom, hpa * 500.0, "")?;

    // boundary values
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 300.0, "")?;
    let mfr_l = vars.add_sclbnd_cnst(bnd_l, 2.0, "")?;

    // steady-state heat transfer and flow solver
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom, temp, vlcp, cond, hprd, hsrc)?;  // arguments: domain, T, rho cp, k, s, Q
    phys.add_temp_bnd(bnd_l, temp_l)?;
    phys.add_hout_bnd(bnd_r)?;
    phys.add_flow_dom(dom, vel, den)?;
    phys.add_mfr_bnd(bnd_l, mfr_l)?;

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-8, 1.0)?;
    Ok(())
}
