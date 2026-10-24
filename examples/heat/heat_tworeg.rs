use vfvchem_vfvm1d::*;
use std::fs::create_dir_all;

/// Transient heat conduction through a two-layer wall (brick and insulation) heated on one side.
/// Run with: `cargo run --release --example heat_tworeg`
///
/// Geometry:
/// Line split into two domains joined by an interface.
///
///   bnd_l +---- dom_0 ----+ itf_0 +---- dom_1 ----+ bnd_r
///
/// Domains:
/// dom_0 - brick (0.10 m)
/// dom_1 - insulation (0.05 m)
///
/// Boundaries:
/// bnd_l - inner surface
/// bnd_r - outer surface
///
/// Interfaces:
/// itf_0 - brick-insulation interface
///
/// Properties:
/// dom_0 - volumetric heat capacity (rho*cp = 1.6e6 J m-3 K-1)
/// dom_0 - thermal conductivity (k = 0.7 W m-1 K-1)
/// dom_0 - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_0 - heat source (Q = 0.0 W m-3)
/// dom_1 - volumetric heat capacity (rho*cp = 3.0e4 J m-3 K-1)
/// dom_1 - thermal conductivity (k = 0.04 W m-1 K-1)
/// dom_1 - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_1 - heat source (Q = 0.0 W m-3)
///
/// Initial conditions:
/// dom_0 - temperature (T = 290 K)
/// dom_1 - temperature (T = 290 K)
///
/// Boundary conditions:
/// [non-constant properties in SI units]
/// bnd_l - inward heat flux ramping up to 200 W m-2 over the first hour
/// q = -200 * min(t / 3600, 1)
/// bnd_r - temperature (T = 290 K)
///
/// Interface conditions:
/// itf_0 - temperature and flux continuity
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_heat_tworeg").unwrap();

    // nodes, points, and variables
    // new_segments - create nodes at the given positions and internal points between them
    // arguments: node_positions, h_min, h_max, ref_len
    // h_min - point spacing at the nodes; h_max - point spacing far from the nodes
    // ref_len - distance from the nodes over which the spacing grows from h_min to h_max
    let mut vars = Variables::new_segments(&[0.0, 0.10, 0.15], 0.0025, 0.0025, 0.0025)?;

    // geometry
    // if using new_segments, the domains are "dom_0", "dom_1", ... from left to right
    // the boundaries are "bnd_l" and "bnd_r"; the interfaces are "itf_0", "itf_1", ... from left to right
    // add_itf - arguments: domain_a, domain_b, group_name
    let dom_0 = vars.add_dom("dom_0")?;
    let dom_1 = vars.add_dom("dom_1")?;
    let bnd_l = vars.add_bnd(dom_0, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom_1, "bnd_r")?;
    let itf_0 = vars.add_itf(dom_0, dom_1, "itf_0")?;

    // unknown domain values
    // arguments: domain, initial_value, output_file
    // initial_value - initial guess for steady-state problems; initial value for transient problems
    // output_file - can be .csv or .vtu; if empty string, no file is written
    // output index is automatically appended to end of output file name (e.g., temp_0_0000.csv)
    let temp_0 = vars.add_scldom_unkn(dom_0, 290.0, "examples/output_heat_tworeg/temp_0.csv")?;
    let temp_1 = vars.add_scldom_unkn(dom_1, 290.0, "examples/output_heat_tworeg/temp_1.csv")?;

    // constant domain values
    // arguments: domain, value, output_file
    let vlcp_0 = vars.add_scldom_cnst(dom_0, 1.6e6, "")?;  // volumetric heat capacity (rho * cp)
    let vlcp_1 = vars.add_scldom_cnst(dom_1, 3.0e4, "")?;  // volumetric heat capacity (rho * cp)
    let cond_0 = vars.add_scldom_cnst(dom_0, 0.7, "")?;  // thermal conductivity (brick)
    let cond_1 = vars.add_scldom_cnst(dom_1, 0.04, "")?;  // thermal conductivity (insulation)
    let hprd_0 = vars.add_scldom_cnst(dom_0, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hprd_1 = vars.add_scldom_cnst(dom_1, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hsrc_0 = vars.add_scldom_cnst(dom_0, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)
    let hsrc_1 = vars.add_scldom_cnst(dom_1, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)

    // non-constant boundary values
    // arguments: boundary, value_func, unknown_ids, output_file
    // value_func - returns value as a function of time, position x, and *unknown* values
    // - unknown values are given in the same order as in unknown_ids
    let hflx_l_func: SclBndFunc = Box::new(|t, _x, _unkn| { -200.0 * (t / 3600.0).min(1.0) });
    let hflx_l = vars.add_sclbnd_func(bnd_l, hflx_l_func, vec![], "")?;  // heat flux (positive if outward; negative if inward)

    // constant boundary values
    // arguments: boundary, value, output_file
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 290.0, "")?;  // temperature

    // transient heat transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_hflx_bnd - register boundary with heat flux
    // add_temp_bnd - register boundary with temperature
    // add_hcnt_itf - register interface with temperature and flux continuity
    let mut phys = TransientHeat::new()?;
    phys.add_heat_dom(dom_0, temp_0, vlcp_0, cond_0, hprd_0, hsrc_0)?;  // arguments: domain, T, rho*cp, k, s, Q
    phys.add_heat_dom(dom_1, temp_1, vlcp_1, cond_1, hprd_1, hsrc_1)?;  // arguments: domain, T, rho*cp, k, s, Q
    phys.add_hflx_bnd(bnd_l, hflx_l)?;  // arguments: boundary, q
    phys.add_temp_bnd(bnd_r, temp_r)?;  // arguments: boundary, T
    phys.add_hcnt_itf(itf_0)?;  // arguments: interface

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, dt, num_ts, num_ts_write, max_iter, tol, damping_factor
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 300.0, 96, 12, 10, 1e-6, 1.0)?;  // 8 h
    Ok(())
}
