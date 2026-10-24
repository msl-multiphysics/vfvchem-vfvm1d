use vfvchem_vfvm1d::*;
use std::fs::create_dir_all;

/// Adaptive time stepping: a steel plate heated by a short, strong volumetric pulse and then cooling.
/// The step shrinks during the pulse and grows while the plate relaxes.
/// Run with: `cargo run --release --example feat_adt`
///
/// Geometry:
///
/// Line domain defined by its end points.
///
///   bnd_l +--------- dom ---------+ bnd_r
///
/// Domains:
/// dom - steel plate (L = 0.02 m)
///
/// Boundaries:
/// bnd_l - left end
/// bnd_r - right end
///
/// Properties:
/// dom - k = 45 W m-1 K-1; rho cp = 3.6e6 J m-3 K-1
/// dom - heat source Q = 5e8 W m-3 for 5 s < t < 6 s, centered profile
///
/// Boundary conditions:
/// bnd_l - temperature (T = 300 K)
/// bnd_r - temperature (T = 300 K)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_feat_adt").unwrap();

    // nodes, points, and variables
    // new_bounds - arguments: x_min, x_max, h_min, h_max, ref_len
    let mut vars = Variables::new_bounds(0.0, 0.02, 5e-4, 5e-4, 5e-4)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;

    // domain values
    let temp = vars.add_scldom_unkn(dom, 300.0, "examples/output_feat_adt/temp.csv")?;
    let vlcp = vars.add_scldom_cnst(dom, 3.6e6, "")?;
    let cond = vars.add_scldom_cnst(dom, 45.0, "")?;
    let hprd = vars.add_scldom_cnst(dom, 0.0, "")?;
    let hsrc = vars.add_scldom_func(dom, Box::new(|t, x, _u| {
        if t > 5.0 && t <= 6.0 { 5e8 * (std::f64::consts::PI * x / 0.02).sin() } else { 0.0 }
    }), vec![], "")?;

    // boundary values
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 300.0, "")?;
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0, "")?;

    // transient heat transfer solver
    let mut phys = TransientHeat::new()?;
    phys.add_heat_dom(dom, temp, vlcp, cond, hprd, hsrc)?;
    phys.add_temp_bnd(bnd_l, temp_l)?;
    phys.add_temp_bnd(bnd_r, temp_r)?;

    // adaptive physics solver
    // the pulse edges (5 s and 6 s) are listed as output times so that steps land on them
    let mut opt = AdaptiveOptions::new(0.1, 20.0, "examples/output_feat_adt")?;
    opt.rtol = 1e-3;  // relative local error per step
    opt.dt_max = 5.0;  // largest time step
    opt.tol = 1e-8;  // tolerance of the iteration in each time step
    opt.out_times = vec![5.0, 6.0, 10.0, 20.0];
    let mut lin = LinearLu::new(1)?;
    phys.solve_adaptive(&mut vars, &mut lin, &opt)?;
    Ok(())
}
