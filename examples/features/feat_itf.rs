use vfvchem_vfvm1d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Heat and mass transfer through two regions with different properties joined by an interface.
/// Run with: `cargo run --release --example feat_itf`
///
/// The interface conditions make the unknown and its flux continuous across the interface.
/// Heat: water on the left and glass on the right, joined with add_hcnt_itf.
/// Mass: two diffusion coefficients, joined with add_mcnt_itf.
///
/// Geometry:
/// Line split into two domains joined by an interface.
///
///   bnd_l +---- dom_0 ----+ itf_0 +---- dom_1 ----+ bnd_r
///
/// Domains:
/// dom_0 - left region (0.5 m)
/// dom_1 - right region (0.5 m)
///
/// Boundaries:
/// bnd_l - left end of dom_0
/// bnd_r - right end of dom_1
///
/// Interfaces:
/// itf_0 - interface between dom_0 and dom_1
///
/// Heat Transfer
///
/// Properties:
/// dom_0 - thermal conductivity (k = 0.6 W m-1 K-1)
/// dom_0 - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_0 - heat source (Q = 0.0 W m-3)
/// dom_1 - thermal conductivity (k = 1.4 W m-1 K-1)
/// dom_1 - heat production coefficient (s = 0.0 W m-3 K-1)
/// dom_1 - heat source (Q = 0.0 W m-3)
///
/// Boundary conditions:
/// bnd_l - temperature (T = 350 K)
/// bnd_r - temperature (T = 300 K)
///
/// Interface conditions:
/// itf_0 - temperature and flux continuity
///
/// Mass Transfer
///
/// Properties:
/// dom_0 - diffusion coefficient (D_00 = 1.0e-9 m2 s-1)
/// dom_0 - production coefficient (s_0 = 0.0 s-1)
/// dom_0 - reaction rate (R_0 = 0.0 mol m-3 s-1)
/// dom_1 - diffusion coefficient (D_00 = 4.0e-9 m2 s-1)
/// dom_1 - production coefficient (s_0 = 0.0 s-1)
/// dom_1 - reaction rate (R_0 = 0.0 mol m-3 s-1)
///
/// Boundary conditions:
/// bnd_l - concentration (c_0 = 1.0 mol m-3)
/// bnd_r - concentration (c_0 = 0.0 mol m-3)
///
/// Interface conditions:
/// itf_0 - concentration and flux continuity
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_feat_itf").unwrap();

    // heat transfer

    // nodes, points, and variables
    // new_segments - create nodes at the given positions and internal points between them
    // arguments: node_positions, h_min, h_max, ref_len
    // h_min - point spacing at the nodes; h_max - point spacing far from the nodes
    // ref_len - distance from the nodes over which the spacing grows from h_min to h_max
    let mut vars = Variables::new_segments(&[0.0, 0.5, 1.0], 0.02, 0.04, 0.2)?;

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
    let temp_0 = vars.add_scldom_unkn(dom_0, 300.0, "examples/output_feat_itf/temp_0.csv")?;
    let temp_1 = vars.add_scldom_unkn(dom_1, 300.0, "examples/output_feat_itf/temp_1.csv")?;

    // constant domain values
    // arguments: domain, value, output_file
    // heat production is multiplied by the temperature; heat source adds the source directly
    // use heat source by default; use heat production (or both) if heat source does not converge
    let cond_0 = vars.add_scldom_cnst(dom_0, 0.6, "")?;  // thermal conductivity (water)
    let cond_1 = vars.add_scldom_cnst(dom_1, 1.4, "")?;  // thermal conductivity (glass)
    let hprd_0 = vars.add_scldom_cnst(dom_0, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hprd_1 = vars.add_scldom_cnst(dom_1, 0.0, "")?;  // heat production coefficient (s; positive if source; negative if sink)
    let hsrc_0 = vars.add_scldom_cnst(dom_0, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)
    let hsrc_1 = vars.add_scldom_cnst(dom_1, 0.0, "")?;  // heat source (Q; positive if source; negative if sink)

    // constant boundary values
    // arguments: boundary, value, output_file
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 350.0, "")?;  // temperature
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0, "")?;  // temperature

    // steady-state heat transfer solver
    // add_heat_dom - register domain with heat transfer
    // add_temp_bnd - register boundary with temperature
    // add_hcnt_itf - register interface with temperature and flux continuity
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom_0, temp_0, cond_0, hprd_0, hsrc_0)?;  // arguments: domain, T, k, s, Q
    phys.add_heat_dom(dom_1, temp_1, cond_1, hprd_1, hsrc_1)?;  // arguments: domain, T, k, s, Q
    phys.add_temp_bnd(bnd_l, temp_l)?;  // arguments: boundary, T
    phys.add_temp_bnd(bnd_r, temp_r)?;  // arguments: boundary, T
    phys.add_hcnt_itf(itf_0)?;  // arguments: interface

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-8, 1.0)?;

    // mass transfer
    // each problem uses its own variables

    // nodes, points, and variables
    let mut vars = Variables::new_segments(&[0.0, 0.5, 1.0], 0.02, 0.04, 0.2)?;

    // geometry
    let dom_0 = vars.add_dom("dom_0")?;
    let dom_1 = vars.add_dom("dom_1")?;
    let bnd_l = vars.add_bnd(dom_0, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom_1, "bnd_r")?;
    let itf_0 = vars.add_itf(dom_0, dom_1, "itf_0")?;

    // unknown domain values
    let conc_0 = vars.add_scldom_unkn(dom_0, 0.0, "examples/output_feat_itf/conc_0.csv")?;
    let conc_1 = vars.add_scldom_unkn(dom_1, 0.0, "examples/output_feat_itf/conc_1.csv")?;

    // constant domain values
    // mass production is multiplied by the concentration; mass source adds the source directly
    // use mass source by default; use mass production (or both) if mass source does not converge
    let diff_0 = vars.add_scldom_cnst(dom_0, 1.0e-9, "")?;  // D_00
    let diff_1 = vars.add_scldom_cnst(dom_1, 4.0e-9, "")?;  // D_00
    let mprd_0 = vars.add_scldom_cnst(dom_0, 0.0, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
    let mprd_1 = vars.add_scldom_cnst(dom_1, 0.0, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
    let msrc_0 = vars.add_scldom_cnst(dom_0, 0.0, "")?;  // mass source (R_0; positive if source; negative if sink)
    let msrc_1 = vars.add_scldom_cnst(dom_1, 0.0, "")?;  // mass source (R_0; positive if source; negative if sink)

    // constant boundary values
    let conc_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "")?;  // concentration
    let conc_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;  // concentration

    // steady-state mass transfer solver
    // SteadyMass::new - arguments: number of components
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mcnt_itf - register interface with concentration and flux continuity
    // diffusion coefficients must be collected in a hashmap with driving concentration index as key
    let mut phys = SteadyMass::new(1)?;
    phys.add_mass_dom(0, dom_0, conc_0, HashMap::from([(0, diff_0)]), mprd_0, msrc_0)?;  // arguments: component, domain, c, D, s, R
    phys.add_mass_dom(0, dom_1, conc_1, HashMap::from([(0, diff_1)]), mprd_1, msrc_1)?;  // arguments: component, domain, c, D, s, R
    phys.add_conc_bnd(0, bnd_l, conc_l)?;  // arguments: component, boundary, c
    phys.add_conc_bnd(0, bnd_r, conc_r)?;  // arguments: component, boundary, c
    phys.add_mcnt_itf(0, itf_0)?;  // arguments: component, interface

    // physics solver
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 10, 1e-8, 1.0)?;

    Ok(())
}
