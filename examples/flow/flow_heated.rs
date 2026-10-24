use vfvchem_vfvm1d::*;
use std::fs::create_dir_all;

/// Transient continuity in a pipe whose gas is heated by a prescribed, time-dependent temperature field.
/// Only continuity is solved: d(rho)/dt + d(rho u)/dx = 0 (no momentum equation).
/// Run with: `cargo run --release --example flow_heated`
///
/// Geometry:
/// Line domain defined by its end points.
///
///   bnd_l +--------- dom ---------+ bnd_r
///
/// Domains:
/// dom - pipe (L = 1 m)
///
/// Boundaries:
/// bnd_l - inlet
/// bnd_r - outlet
///
/// Properties:
/// [non-constant properties in SI units]
/// dom - ideal gas density rho = p M / (R T)
/// T(x, t) = 300 + 300 x (1 - exp(-t / 2))
///
/// Boundary conditions:
/// bnd_l - mass flow rate per unit area (G = 0.5 kg m-2 s-1; positive into the domain)
/// bnd_r - outlet (continuity only)
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_flow_heated").unwrap();

    // nodes, points, and variables
    // new_bounds - arguments: x_min, x_max, h_min, h_max, ref_len
    let mut vars = Variables::new_bounds(0.0, 1.0, 0.02, 0.02, 0.02)?;

    // geometry
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let _bnd_r = vars.add_bnd(dom, "bnd_r")?;

    // unknown domain values
    let vel = vars.add_scldom_unkn(dom, 0.4, "examples/output_flow_heated/vel.csv")?;  // velocity (+x component)

    // non-constant domain values
    // arguments: domain, value_func, unknown_ids, output_file
    let den_func: SclDomFunc = Box::new(|t, x, _unkn| {
        let temp = 300.0 + 300.0 * x * (1.0 - (-t / 2.0).exp());
        101325.0 * 0.029 / (8.314 * temp)
    });
    let den = vars.add_scldom_func(dom, den_func, vec![], "examples/output_flow_heated/den.csv")?;  // density

    // constant boundary values
    let mfr_l = vars.add_sclbnd_cnst(bnd_l, 0.5, "")?;  // mass flow rate per unit area (positive into the domain)

    // transient flow solver
    // add_flow_dom - register domain with flow
    // add_mfr_bnd - register inlet boundary with mass flow rate per unit area
    let mut phys = TransientFlow::new()?;
    phys.add_flow_dom(dom, vel, den)?;  // arguments: domain, u, rho
    phys.add_mfr_bnd(bnd_l, mfr_l)?;  // arguments: boundary, G

    // physics solver
    // arguments: vars, linear_solver, dt, num_ts, num_ts_write, max_iter, tol, damping_factor
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 0.1, 100, 10, 10, 1e-8, 1.0)?;
    Ok(())
}
