// Analytic and consistency tests of the 1D solvers
// Run with: cargo test --release
use vfvchem_vfvm1d::*;
use std::collections::HashMap;
use std::f64::consts::PI;

fn fmt(e: &[f64]) -> String { e.iter().map(|x| format!("{x:.3e}")).collect::<Vec<_>>().join(", ") }

fn lin() -> Box<dyn LinearBase> { LinearLu::new(1).unwrap() }

fn pos(vars: &Variables, dom_id: usize) -> Vec<f64> {
    // vid -> x; the points come first, then the nodes
    let dom = &vars.dom[dom_id];
    let mut pos: Vec<f64> = dom.pts_gpt.iter().map(|&gpt| vars.pts.pts_pos[gpt]).collect();
    pos.extend(dom.nds_gnd.iter().map(|&gnd| vars.nds.nds_pos[gnd]));
    pos
}

fn order(vars: &Variables, dom_id: usize) -> Vec<usize> {
    // vid in order of x
    let seg = &vars.seg[dom_id];
    let mut vid_ord = vec![seg.seg_vid[0][0]];
    vid_ord.extend(seg.seg_vid.iter().map(|&[_, vid_1]| vid_1));
    vid_ord
}

fn bnd_vid(vars: &Variables, bnd_id: usize) -> usize {
    let bnd = &vars.bnd[bnd_id];
    let dom = &vars.dom[bnd.dom_id];
    dom.num_pts() + dom.nds_nid[&bnd.nds_gnd[0]]
}

fn max_err(vars: &Variables, id: usize, exact: &dyn Fn(f64) -> f64) -> f64 {
    let pos = pos(vars, vars.scldom[id].dom_id);
    pos.iter().enumerate().map(|(vid, &x)| (vars.scldom[id].value[vid] - exact(x)).abs()).fold(0.0, f64::max)
}

fn tmp_dir(name: &str) -> String {
    let p = std::env::temp_dir().join(format!("vfvm1d_test_{name}"));
    std::fs::create_dir_all(&p).unwrap();
    p.display().to_string()
}

// heat

#[test]
fn heat_conduction_source() -> Result<(), VFVChemError> {
    // k T'' + Q = 0, T(0) = 350, T(L) = 300; FV is exact for the quadratic profile at the vertices
    let (l, k, q) = (0.1, 0.6, 1.0e4);
    let exact = |x: f64| 350.0 - 500.0 * x + q * x * (l - x) / (2.0 * k);
    let mut vars = Variables::new_bounds(0.0, l, 0.005, 0.005, 0.005)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let temp = vars.add_scldom_unkn(dom, 300.0, "")?;
    let cond = vars.add_scldom_cnst(dom, k, "")?;
    let hprd = vars.add_scldom_cnst(dom, 0.0, "")?;
    let hsrc = vars.add_scldom_cnst(dom, q, "")?;
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 350.0, "")?;
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0, "")?;
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom, temp, cond, hprd, hsrc)?;
    phys.add_temp_bnd(bnd_l, temp_l)?;
    phys.add_temp_bnd(bnd_r, temp_r)?;
    phys.solve(&mut vars, &mut lin(), 10, 1e-10, 1.0)?;
    let err = max_err(&vars, temp, &exact);
    println!("heat_conduction_source: max error = {err:.3e}");
    assert!(err < 1e-8);
    Ok(())
}

#[test]
fn heat_conduction_graded() -> Result<(), VFVChemError> {
    // as heat_conduction_source on a graded mesh; still exact for the quadratic profile
    let (l, k, q) = (0.1, 0.6, 1.0e4);
    let exact = |x: f64| 350.0 - 500.0 * x + q * x * (l - x) / (2.0 * k);
    let mut vars = Variables::new_bounds(0.0, l, 0.001, 0.01, 0.03)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let temp = vars.add_scldom_unkn(dom, 300.0, "")?;
    let cond = vars.add_scldom_cnst(dom, k, "")?;
    let hprd = vars.add_scldom_cnst(dom, 0.0, "")?;
    let hsrc = vars.add_scldom_cnst(dom, q, "")?;
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 350.0, "")?;
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 300.0, "")?;
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom, temp, cond, hprd, hsrc)?;
    phys.add_temp_bnd(bnd_l, temp_l)?;
    phys.add_temp_bnd(bnd_r, temp_r)?;
    phys.solve(&mut vars, &mut lin(), 10, 1e-10, 1.0)?;
    let x = pos(&vars, dom);
    let h: Vec<f64> = vars.seg[dom].ngh_dist.iter().flatten().copied().collect();
    let (h_min, h_max) = (h.iter().cloned().fold(f64::MAX, f64::min), h.iter().cloned().fold(0.0, f64::max));
    let err = max_err(&vars, temp, &exact);
    println!("heat_conduction_graded: {} vertices; h = {h_min:.3e} to {h_max:.3e}; max error = {err:.3e}", x.len());
    assert!(h_max > 5.0 * h_min);
    assert!(err < 1e-8);
    Ok(())
}

#[test]
fn heat_two_region_nonuniform_flux() -> Result<(), VFVChemError> {
    // two materials with an interface; inward heat flux q0 on the left, T(L) = T_r
    // piecewise linear profile with flux continuity: T = T_r + q0 * (L - x) / k2 (x > a); slope q0 / k1 (x < a)
    let (a, l, k1, k2, q0, t_r) = (0.3, 1.0, 2.0, 0.5, 100.0, 300.0);
    let exact = |x: f64| if x >= a { t_r + q0 * (l - x) / k2 } else { t_r + q0 * (l - a) / k2 + q0 * (a - x) / k1 };
    let mut vars = Variables::new_segments(&[0.0, a, l], 0.01, 0.08, 0.2)?;
    let dom_0 = vars.add_dom("dom_0")?;
    let dom_1 = vars.add_dom("dom_1")?;
    let bnd_l = vars.add_bnd(dom_0, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom_1, "bnd_r")?;
    let itf = vars.add_itf(dom_0, dom_1, "itf_0")?;
    let temp_0 = vars.add_scldom_unkn(dom_0, 300.0, "")?;
    let temp_1 = vars.add_scldom_unkn(dom_1, 300.0, "")?;
    let cond_0 = vars.add_scldom_cnst(dom_0, k1, "")?;
    let cond_1 = vars.add_scldom_cnst(dom_1, k2, "")?;
    let zero_0 = vars.add_scldom_cnst(dom_0, 0.0, "")?;
    let zero_1 = vars.add_scldom_cnst(dom_1, 0.0, "")?;
    let hflx_l = vars.add_sclbnd_cnst(bnd_l, -q0, "")?;  // inward
    let temp_r = vars.add_sclbnd_cnst(bnd_r, t_r, "")?;
    let mut phys = SteadyHeat::new()?;
    phys.add_heat_dom(dom_0, temp_0, cond_0, zero_0, zero_0)?;
    phys.add_heat_dom(dom_1, temp_1, cond_1, zero_1, zero_1)?;
    phys.add_hflx_bnd(bnd_l, hflx_l)?;
    phys.add_temp_bnd(bnd_r, temp_r)?;
    phys.add_hcnt_itf(itf)?;
    phys.solve(&mut vars, &mut lin(), 10, 1e-10, 1.0)?;
    let err = max_err(&vars, temp_0, &exact).max(max_err(&vars, temp_1, &exact));
    println!("heat_two_region_nonuniform_flux: max error = {err:.3e}");
    assert!(err < 1e-8);
    Ok(())
}

fn heat_transient_err(n: usize) -> Result<f64, VFVChemError> {
    // T_t = alpha T_xx, T(0) = T(1) = 0, T(x, 0) = sin(pi x); T = sin(pi x) exp(-alpha pi^2 t); dt ~ h^2
    let alpha = 0.1;
    let h = 1.0 / n as f64;
    let dt = 0.5 * h * h;
    let t_end = 0.5;
    let num_ts = (t_end / dt).round() as usize;
    let mut vars = Variables::new_bounds(0.0, 1.0, h, h, h)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let temp = vars.add_scldom_unkn(dom, 0.0, "")?;
    for (vid, x) in pos(&vars, dom).into_iter().enumerate() { vars.scldom[temp].value[vid] = (PI * x).sin(); }
    let vlcp = vars.add_scldom_cnst(dom, 1.0, "")?;
    let cond = vars.add_scldom_cnst(dom, alpha, "")?;
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 0.0, "")?;
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let mut phys = TransientHeat::new()?;
    phys.add_heat_dom(dom, temp, vlcp, cond, zero, zero)?;
    phys.add_temp_bnd(bnd_l, temp_l)?;
    phys.add_temp_bnd(bnd_r, temp_r)?;
    phys.solve(&mut vars, &mut lin(), dt, num_ts, num_ts, 10, 1e-10, 1.0)?;
    let t = dt * num_ts as f64;
    Ok(max_err(&vars, temp, &|x| (PI * x).sin() * (-alpha * PI * PI * t).exp()))
}

#[test]
fn heat_transient_refinement() -> Result<(), VFVChemError> {
    let errs: Vec<f64> = [10, 20, 40].iter().map(|&n| heat_transient_err(n)).collect::<Result<_, _>>()?;
    println!("heat_transient_refinement: errors (N = 10, 20, 40) = {}; ratios = {:.2}, {:.2}", fmt(&errs), errs[0] / errs[1], errs[1] / errs[2]);
    assert!(errs[0] / errs[1] > 3.0 && errs[1] / errs[2] > 3.0);  // second order with dt ~ h^2
    assert!(errs[2] < 1e-3);
    Ok(())
}

// mass

fn mass_reaction_err(n: usize) -> Result<f64, VFVChemError> {
    // D c'' - k c = 0, c(0) = 1, zero flux at L; c = cosh(m (L - x)) / cosh(m L), m = sqrt(k / D)
    let (l, d, k): (f64, f64, f64) = (1.0, 1.0e-2, 0.4);
    let m = (k / d).sqrt();
    let h = l / n as f64;
    let mut vars = Variables::new_bounds(0.0, l, h, h, h)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let conc = vars.add_scldom_unkn(dom, 0.5, "")?;
    let diff = vars.add_scldom_cnst(dom, d, "")?;
    let mprd = vars.add_scldom_cnst(dom, -k, "")?;
    let msrc = vars.add_scldom_cnst(dom, 0.0, "")?;
    let conc_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "")?;
    let mflx_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let mut phys = SteadyMass::new(1)?;
    phys.add_mass_dom(0, dom, conc, HashMap::from([(0, diff)]), mprd, msrc)?;
    phys.add_conc_bnd(0, bnd_l, conc_l)?;
    phys.add_mflx_bnd(0, bnd_r, mflx_r)?;
    phys.solve(&mut vars, &mut lin(), 10, 1e-10, 1.0)?;
    Ok(max_err(&vars, conc, &|x| (m * (l - x)).cosh() / (m * l).cosh()))
}

#[test]
fn mass_reaction_refinement() -> Result<(), VFVChemError> {
    let errs: Vec<f64> = [20, 40, 80].iter().map(|&n| mass_reaction_err(n)).collect::<Result<_, _>>()?;
    println!("mass_reaction_refinement: errors (N = 20, 40, 80) = {}; ratios = {:.2}, {:.2}", fmt(&errs), errs[0] / errs[1], errs[1] / errs[2]);
    assert!(errs[0] / errs[1] > 3.5 && errs[1] / errs[2] > 3.5);
    Ok(())
}

#[test]
fn mass_transient_source_conservation() -> Result<(), VFVChemError> {
    // closed slab with a nonuniform, concentration-dependent sink and a constant source
    // sum(c V) must change by exactly dt * sum((mprd c + msrc) V) per step
    let mut vars = Variables::new_bounds(0.0, 1.0, 0.01, 0.05, 0.2)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let conc = vars.add_scldom_unkn(dom, 1.0, "")?;
    let diff = vars.add_scldom_cnst(dom, 1e-2, "")?;
    let mprd = vars.add_scldom_func(dom, Box::new(|_t, x, _u| -2.0 * x), vec![], "")?;
    let msrc = vars.add_scldom_cnst(dom, 0.5, "")?;
    let mflx_l = vars.add_sclbnd_cnst(bnd_l, 0.0, "")?;
    let mflx_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let mut phys = TransientMass::new(1)?;
    phys.add_mass_dom(0, dom, conc, HashMap::from([(0, diff)]), mprd, msrc)?;
    phys.add_mflx_bnd(0, bnd_l, mflx_l)?;
    phys.add_mflx_bnd(0, bnd_r, mflx_r)?;
    let mut opt = AdaptiveOptions::new(0.05, 1.0, &tmp_dir("mass_cons"))?;
    opt.fixed = true;
    opt.tol = 1e-12;
    let mut worst: f64 = 0.0;
    phys.solve_adaptive_with(&mut vars, &mut lin(), &opt, &mut |v: &Variables, _t, dt| {
        let cvol = &v.seg[dom].cvol_len;
        let (c, cp, s) = (&v.scldom[conc].value, &v.scldom[conc].value_prev, &v.scldom[mprd].value);
        let bal: f64 = (0..cvol.len()).map(|i| ((c[i] - cp[i]) / dt - s[i] * c[i] - 0.5) * cvol[i]).sum();
        worst = worst.max(bal.abs());
    })?;
    println!("mass_transient_source_conservation: worst balance = {worst:.3e}");
    assert!(worst < 1e-10);
    Ok(())
}

fn advdiff_err(n: usize) -> Result<f64, VFVChemError> {
    // rho u w' = D w'', w(0) = 1, w(L) = 0 with variable density (rho u constant = G)
    // w = (exp(Pe x / L) - exp(Pe)) / (1 - exp(Pe)), Pe = G L / D
    let (l, g, d) = (1.0, 0.2, 0.02);
    let pe = g * l / d;
    let h = l / n as f64;
    let mut vars = Variables::new_bounds(0.0, l, h, h, h)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let vel = vars.add_scldom_unkn(dom, 0.0, "")?;
    let den = vars.add_scldom_func(dom, Box::new(|_t, x, _u| 1.0 + x), vec![], "")?;
    let conc = vars.add_scldom_unkn(dom, 0.5, "")?;
    let diff = vars.add_scldom_cnst(dom, d, "")?;
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;
    let mfr_l = vars.add_sclbnd_cnst(bnd_l, g, "")?;
    let conc_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "")?;
    let conc_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let mut phys = SteadyMassFlow::new(1)?;
    phys.add_mass_dom(0, dom, conc, HashMap::from([(0, diff)]), zero, zero)?;
    phys.add_conc_bnd(0, bnd_l, conc_l)?;
    phys.add_conc_bnd(0, bnd_r, conc_r)?;
    phys.add_flow_dom(dom, vel, den)?;
    phys.add_mfr_bnd(bnd_l, mfr_l)?;
    phys.solve(&mut vars, &mut lin(), 20, 1e-10, 1.0)?;
    Ok(max_err(&vars, conc, &|x| ((pe * x / l).exp() - pe.exp()) / (1.0 - pe.exp())))
}

#[test]
fn massflow_advection_diffusion_refinement() -> Result<(), VFVChemError> {
    let errs: Vec<f64> = [20, 40, 80, 160].iter().map(|&n| advdiff_err(n)).collect::<Result<_, _>>()?;
    println!("massflow_advection_diffusion_refinement: errors (N = 20, 40, 80, 160) = {}", fmt(&errs));
    for k in 0..3 { assert!(errs[k + 1] < errs[k]); }
    assert!(errs[3] < 5e-3);
    Ok(())
}

// flow

#[test]
fn flow_variable_density_steady() -> Result<(), VFVChemError> {
    // two domains with a density jump; inlet on the right (flow in -x); rho u = -G everywhere
    let g = 0.7;
    let mut vars = Variables::new_segments(&[0.0, 0.4, 1.0], 0.02, 0.08, 0.2)?;
    let dom_0 = vars.add_dom("dom_0")?;
    let dom_1 = vars.add_dom("dom_1")?;
    let bnd_r = vars.add_bnd(dom_1, "bnd_r")?;
    let itf = vars.add_itf(dom_0, dom_1, "itf_0")?;
    let vel_0 = vars.add_scldom_unkn(dom_0, 0.0, "")?;
    let vel_1 = vars.add_scldom_unkn(dom_1, 0.0, "")?;
    let den_0 = vars.add_scldom_func(dom_0, Box::new(|_t, x, _u| 2.0 + x * x), vec![], "")?;
    let den_1 = vars.add_scldom_func(dom_1, Box::new(|_t, x, _u| 0.5 + (3.0 * x).sin().abs()), vec![], "")?;
    let mfr_r = vars.add_sclbnd_cnst(bnd_r, g, "")?;
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom_0, vel_0, den_0)?;
    phys.add_flow_dom(dom_1, vel_1, den_1)?;
    phys.add_mfr_bnd(bnd_r, mfr_r)?;
    phys.add_fcnt_itf(itf)?;
    phys.solve(&mut vars, &mut lin(), 10, 1e-12, 1.0)?;
    let mut err: f64 = 0.0;
    for (vel, den) in [(vel_0, den_0), (vel_1, den_1)] {
        for v in 0..vars.scldom[vel].value.len() { err = err.max((vars.scldom[den].value[v] * vars.scldom[vel].value[v] + g).abs()); }
    }
    // outward face flux is +G toward -x (face_norm = -1) and -G toward +x
    for flw in &vars.sclflw {
        let seg = &vars.seg[flw.dom_id];
        for (vid, flux_n) in flw.flux.iter().enumerate() {
            for (loc, &f) in flux_n.iter().enumerate() { err = err.max((f - seg.face_norm[vid][loc] * (-g)).abs()); }
        }
    }
    println!("flow_variable_density_steady: max |rho u + G| = {err:.3e}");
    assert!(err < 1e-12);
    Ok(())
}

#[test]
fn flow_transient_density() -> Result<(), VFVChemError> {
    // rho = 1 + 0.5 t (uniform), u(0) = 1: rho u = rho - rho_t x, u = 1 - 0.5 x / rho (exact for the trapezoid rule)
    let mut vars = Variables::new_bounds(0.0, 1.0, 0.02, 0.2, 0.3)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let vel = vars.add_scldom_unkn(dom, 1.0, "")?;
    let den = vars.add_scldom_func(dom, Box::new(|t, _x, _u| 1.0 + 0.5 * t), vec![], "")?;
    let vel_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "")?;
    let mut phys = TransientFlow::new()?;
    phys.add_flow_dom(dom, vel, den)?;
    phys.add_vel_bnd(bnd_l, vel_l)?;
    phys.solve(&mut vars, &mut lin(), 0.1, 10, 10, 10, 1e-12, 1.0)?;
    let rho = 1.5;
    let err = max_err(&vars, vel, &|x| 1.0 - 0.5 * x / rho);

    // flux leaving the outlet control volume; the outlet node has one face, toward -x
    let vid_r = bnd_vid(&vars, bnd_r);
    let flux_out = -vars.sclflw[0].flux[vid_r][0] - vars.seg[dom].cvol_len[vid_r] * 0.5;
    let flux_err = (flux_out - rho * vars.scldom[vel].value[vid_r]).abs();
    println!("flow_transient_density: max velocity error = {err:.3e}; outlet flux mismatch = {flux_err:.3e}");
    assert!(err < 1e-12 && flux_err < 1e-12);
    Ok(())
}

#[test]
fn massflow_transient_conservation() -> Result<(), VFVChemError> {
    // density varying in x and t, inlet w = 1 (Dirichlet), outlet with zero diffusive flux, constant source
    // balance over all control volumes but the inlet: d/dt sum(rho w V) = J_in - F_out w_out + sum(msrc V)
    let (d, q) = (1e-2, 0.3);
    let mut vars = Variables::new_bounds(0.0, 1.0, 0.025, 0.025, 0.025)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let vel = vars.add_scldom_unkn(dom, 0.5, "")?;
    let den = vars.add_scldom_func(dom, Box::new(|t, x, _u| 1.0 + 0.5 * x * (2.0 * t).sin()), vec![], "")?;
    let conc = vars.add_scldom_unkn(dom, 0.0, "")?;
    let diff = vars.add_scldom_cnst(dom, d, "")?;
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;
    let msrc = vars.add_scldom_cnst(dom, q, "")?;
    let vel_l = vars.add_sclbnd_cnst(bnd_l, 0.5, "")?;
    let conc_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "")?;
    let mut phys = TransientMassFlow::new(1)?;
    phys.add_mass_dom(0, dom, conc, HashMap::from([(0, diff)]), zero, msrc)?;
    phys.add_conc_bnd(0, bnd_l, conc_l)?;
    phys.add_mout_bnd(0, bnd_r)?;
    phys.add_flow_dom(dom, vel, den)?;
    phys.add_vel_bnd(bnd_l, vel_l)?;
    let mut opt = AdaptiveOptions::new(0.02, 2.0, &tmp_dir("massflow_cons"))?;
    opt.fixed = true;
    opt.tol = 1e-12;
    let mut worst: f64 = 0.0;
    phys.solve_adaptive_with(&mut vars, &mut lin(), &opt, &mut |v: &Variables, _t, dt| {
        let seg = &v.seg[dom];
        let vid_ord = order(v, dom);
        let (w, wp) = (&v.scldom[conc].value, &v.scldom[conc].value_prev);
        let (r, rp) = (&v.scldom[den].value, &v.scldom[den].value_prev);
        let flux = &v.sclflw[0].flux;
        let (vid_in, vid_1, vid_out) = (vid_ord[0], vid_ord[1], *vid_ord.last().unwrap());
        let dmass: f64 = vid_ord[1..].iter().map(|&i| (r[i] * w[i] - rp[i] * wp[i]) * seg.cvol_len[i] / dt).sum();

        // conservative flux into the first interior CV through the inlet face (power law), consistent with the operators
        let f = flux[vid_in][0];
        let cond = d / seg.ngh_dist[vid_in][0];
        let pl = (1.0 - 0.1 * (f / cond).abs()).max(0.0).powi(5);
        let j_in = (cond * pl + f.max(0.0)) * (w[vid_in] - w[vid_1]) + f * w[vid_1];

        // flux leaving the outlet control volume
        let f_out = -flux[vid_out][0] - seg.cvol_len[vid_out] * (r[vid_out] - rp[vid_out]) / dt;
        let src: f64 = vid_ord[1..].iter().map(|&i| q * seg.cvol_len[i]).sum();
        let bal = dmass - (j_in - f_out * w[vid_out] + src);
        worst = worst.max(bal.abs());
    })?;
    println!("massflow_transient_conservation: worst balance = {worst:.3e}");
    assert!(worst < 1e-10);
    Ok(())
}

#[test]
fn heatflow_steady_exponential() -> Result<(), VFVChemError> {
    // heated gas: rho cp u T' = k T'', T(0) = 300, T(L) = 400; constant G; power law on a fine mesh
    let (l, g, cp, k) = (0.5, 0.01, 1000.0, 0.5);
    let pe = g * cp * l / k;
    let exact = |x: f64| 300.0 + 100.0 * ((pe * x / l).exp() - 1.0) / (pe.exp() - 1.0);
    let mut vars = Variables::new_bounds(0.0, l, 0.0025, 0.0025, 0.0025)?;
    let dom = vars.add_dom("dom")?;
    let bnd_l = vars.add_bnd(dom, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom, "bnd_r")?;
    let temp = vars.add_scldom_unkn(dom, 300.0, "")?;
    let vel = vars.add_scldom_unkn(dom, 0.0, "")?;
    let den = vars.add_scldom_func(dom, Box::new(|_t, _x, u| 101325.0 * 0.029 / (8.314 * u[0])), vec![temp], "")?;  // ideal gas
    let vlcp = vars.add_scldom_func(dom, Box::new(move |_t, _x, u| 101325.0 * 0.029 / (8.314 * u[0]) * cp), vec![temp], "")?;
    let cond = vars.add_scldom_cnst(dom, k, "")?;
    let zero = vars.add_scldom_cnst(dom, 0.0, "")?;
    let temp_l = vars.add_sclbnd_cnst(bnd_l, 300.0, "")?;
    let temp_r = vars.add_sclbnd_cnst(bnd_r, 400.0, "")?;
    let mfr_l = vars.add_sclbnd_cnst(bnd_l, g, "")?;
    let mut phys = SteadyHeatFlow::new()?;
    phys.add_heat_dom(dom, temp, vlcp, cond, zero, zero)?;
    phys.add_temp_bnd(bnd_l, temp_l)?;
    phys.add_temp_bnd(bnd_r, temp_r)?;
    phys.add_flow_dom(dom, vel, den)?;
    phys.add_mfr_bnd(bnd_l, mfr_l)?;
    phys.solve(&mut vars, &mut lin(), 50, 1e-10, 1.0)?;
    let err = max_err(&vars, temp, &exact);
    let num_tot = vars.dom[dom].num_tot();
    let g_err = (0..num_tot).map(|v| (vars.scldom[den].value[v] * vars.scldom[vel].value[v] - g).abs()).fold(0.0, f64::max);
    println!("heatflow_steady_exponential: max T error = {err:.3e} K; max |rho u - G| = {g_err:.3e}");
    assert!(err < 0.5 && g_err < 1e-12);
    Ok(())
}

// errors

#[test]
fn error_checks() -> Result<(), VFVChemError> {
    // flow domain without an inlet
    let mut vars = Variables::new_bounds(0.0, 1.0, 0.25, 0.25, 0.25)?;
    let dom = vars.add_dom("dom")?;
    let vel = vars.add_scldom_unkn(dom, 0.0, "")?;
    let den = vars.add_scldom_cnst(dom, 1.0, "")?;
    let mut phys = SteadyFlow::new()?;
    phys.add_flow_dom(dom, vel, den)?;
    let res = phys.solve(&mut vars, &mut lin(), 5, 1e-8, 1.0);
    assert!(matches!(res, Err(VFVChemError::FlowNoInlet { .. })));

    // component index out of range
    let mut phys = SteadyMass::new(1)?;
    let res = phys.add_conc_bnd(1, 0, 0);
    assert!(matches!(res, Err(VFVChemError::MassInvalidComponent { .. })));

    // mass interface between a domain with flow and one without
    let mut vars = Variables::new_segments(&[0.0, 0.5, 1.0], 0.1, 0.1, 0.1)?;
    let dom_0 = vars.add_dom("dom_0")?;
    let dom_1 = vars.add_dom("dom_1")?;
    let bnd_l = vars.add_bnd(dom_0, "bnd_l")?;
    let bnd_r = vars.add_bnd(dom_1, "bnd_r")?;
    let itf = vars.add_itf(dom_0, dom_1, "itf_0")?;
    let vel_0 = vars.add_scldom_unkn(dom_0, 1.0, "")?;
    let den_0 = vars.add_scldom_cnst(dom_0, 1.0, "")?;
    let conc_0 = vars.add_scldom_unkn(dom_0, 0.0, "")?;
    let conc_1 = vars.add_scldom_unkn(dom_1, 0.0, "")?;
    let diff_0 = vars.add_scldom_cnst(dom_0, 1.0, "")?;
    let diff_1 = vars.add_scldom_cnst(dom_1, 1.0, "")?;
    let zero_0 = vars.add_scldom_cnst(dom_0, 0.0, "")?;
    let zero_1 = vars.add_scldom_cnst(dom_1, 0.0, "")?;
    let vel_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "")?;
    let conc_l = vars.add_sclbnd_cnst(bnd_l, 1.0, "")?;
    let conc_r = vars.add_sclbnd_cnst(bnd_r, 0.0, "")?;
    let mut phys = SteadyMassFlow::new(1)?;
    phys.add_mass_dom(0, dom_0, conc_0, HashMap::from([(0, diff_0)]), zero_0, zero_0)?;
    phys.add_mass_dom(0, dom_1, conc_1, HashMap::from([(0, diff_1)]), zero_1, zero_1)?;
    phys.add_conc_bnd(0, bnd_l, conc_l)?;
    phys.add_conc_bnd(0, bnd_r, conc_r)?;
    phys.add_mcnt_itf(0, itf)?;
    phys.add_flow_dom(dom_0, vel_0, den_0)?;
    phys.add_vel_bnd(bnd_l, vel_l)?;
    let res = phys.solve(&mut vars, &mut lin(), 5, 1e-8, 1.0);
    assert!(matches!(res, Err(VFVChemError::MassInterfaceMixedFlow { .. })));
    Ok(())
}
