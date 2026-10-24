use vfvchem_vfvm1d::*;
use std::collections::HashMap;
use std::fs::create_dir_all;

/// Steady-state multicomponent diffusion equation with cross-diffusion of component 0.
/// Run with: `cargo run --release --example mass_slab`
///
/// Geometry:
/// Line domain defined by its end points.
///
///   bnd_l +--------- dom ---------+ bnd_r
///
/// Domains:
/// dom - main domain (1.0 m)
///
/// Boundaries:
/// bnd_l - left end
/// bnd_r - right end
///
/// Properties:
/// dom - diffusion coefficients (m2 s-1)
/// [D_00 D_01] = [1.0e-4 0.5e-4]
/// [D_10 D_11]   [0.0    1.0e-4]
/// dom - production coefficients (s-1; positive if production)
/// [s_0] = [0.0]
/// [s_1]   [0.0]
/// dom - reaction rates (mol m-3 s-1; positive if production)
/// [R_0] = [-0.2]
/// [R_1]   [+0.1]
///
/// Boundary conditions:
/// bnd_l - outward flux 0; concentration 1 (mol m-2 s-1; mol m-3)
/// [N_0] = [0.5]
/// [c_1]   [1.0e-4]
/// bnd_r - concentration (mol m-3)
/// [c_0] = [3.0e-4]
/// [c_1]   [0.0]
///
fn main() -> Result<(), VFVChemError> {
    // output directory
    create_dir_all("examples/output_mass_slab").unwrap();

    // nodes, points, and variables
    // new_bounds - create boundary nodes and internal points between two end points
    // arguments: x_min, x_max, h_min, h_max, ref_len
    // h_min - point spacing at the boundary; h_max - point spacing far from the boundary
    // ref_len - distance from the boundary over which the spacing grows from h_min to h_max
    let mut vars = Variables::new_bounds(0.0, 1.0, 0.05, 0.05, 0.05)?;

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
    // output index is automatically appended to end of output file name (e.g., conc_0_0000.csv)
    let conc_0 = vars.add_scldom_unkn(dom, 0.0, "examples/output_mass_slab/conc_0.csv")?;
    let conc_1 = vars.add_scldom_unkn(dom, 0.0, "examples/output_mass_slab/conc_1.csv")?;

    // constant domain values
    // arguments: domain, value, output_file
    // no need to declare zero diffusion coefficients
    // component 0 diffuses down its own gradient and the gradient of component 1
    // mass production is multiplied by the concentration; mass source adds the source directly
    // use mass source by default; use mass production (or both) if mass source does not converge
    let diff_00 = vars.add_scldom_cnst(dom, 1.0e-4, "")?;  // D_00
    let diff_01 = vars.add_scldom_cnst(dom, 0.5e-4, "")?;  // D_01
    let diff_11 = vars.add_scldom_cnst(dom, 1.0e-4, "")?;  // D_11
    let mprd_0 = vars.add_scldom_cnst(dom, 0.0, "")?;  // mass production coefficient (s_0; positive if source; negative if sink)
    let mprd_1 = vars.add_scldom_cnst(dom, 0.0, "")?;  // mass production coefficient (s_1; positive if source; negative if sink)
    let msrc_0 = vars.add_scldom_cnst(dom, -0.2, "")?;  // mass source (R_0; positive if source; negative if sink)
    let msrc_1 = vars.add_scldom_cnst(dom,  0.1, "")?;  // mass source (R_1; positive if source; negative if sink)

    // constant boundary values
    // arguments: boundary, value, output_file
    let mflx_0_l = vars.add_sclbnd_cnst(bnd_l, 0.5, "")?;  // molar flux (positive if outward; negative if inward)
    let conc_0_r = vars.add_sclbnd_cnst(bnd_r, 3e-4, "")?;  // concentration
    let conc_1_l = vars.add_sclbnd_cnst(bnd_l, 1e-4, "")?;  // concentration
    let conc_1_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;  // concentration

    // steady-state mass transfer solver
    // SteadyMass::new - arguments: number of components
    // add_mass_dom - register domain with mass transfer
    // add_conc_bnd - register boundary with concentration
    // add_mflx_bnd - register boundary with molar flux
    // diffusion coefficients must be collected in a hashmap with driving concentration index as key
    let mut phys = SteadyMass::new(2)?;
    phys.add_mass_dom(0, dom, conc_0, HashMap::from([(0, diff_00), (1, diff_01)]), mprd_0, msrc_0)?;  // arguments: component, domain, c, D, s, R
    phys.add_mass_dom(1, dom, conc_1, HashMap::from([(1, diff_11)]), mprd_1, msrc_1)?;  // arguments: component, domain, c, D, s, R
    phys.add_mflx_bnd(0, bnd_l, mflx_0_l)?;  // arguments: component, boundary, N
    phys.add_conc_bnd(0, bnd_r, conc_0_r)?;  // arguments: component, boundary, c
    phys.add_conc_bnd(1, bnd_l, conc_1_l)?;  // arguments: component, boundary, c
    phys.add_conc_bnd(1, bnd_r, conc_1_r)?;  // arguments: component, boundary, c

    // physics solver
    // LinearLu::new - sparse LU linear solver; arguments: num_thread
    // arguments: vars, linear_solver, max_iter, tol, damping_factor
    // damping_factor - between 0.0 and 1.0; lower for stability and higher for speed (if linear or nearly linear)
    // for highly non-linear problems, using a lower damping factor (e.g., 0.8-0.9) may be faster
    let mut lin = LinearLu::new(1)?;
    phys.solve(&mut vars, &mut lin, 50, 1e-3, 0.8)?;

    Ok(())
}
