use crate::base::prelude::*;
use faer::Col;

// Per-field scaling of convergence measures (Picard residual and adaptive error norm); as in 2D, without vectors and
// pressure (1D solves no momentum equation and no pressure).
//
// Every unknown field (scalar domain) is measured against its own reference scale, and the measure is the maximum
// over fields, so a field of large magnitude (e.g., T ~ 300 K) does not hide the others.
//   Picard residual: max_f rms_f(x_new - x_old) / s_f
// Reference scale s_f: max |v| over both states; a value set with set_scldom_scale is a lower bound.
// A tiny absolute floor (SCALE_FLOOR) avoids division by zero for fields that are identically zero.

pub const SCALE_FLOOR: f64 = 1e-30;

pub struct FieldNorm {
    pub name: String,   // e.g. "scldom 0", for output
    pub id: usize,      // scldom id
    pub start: usize,   // start in solution vector
    pub len: usize,     // number of vertices
}

pub fn field_norms(vars: &Variables, exclude: &[usize]) -> Vec<FieldNorm> {
    // unknown fields in the solution vector; scldom ids in exclude are left out
    let mut fields = Vec::new();
    for (id, scldom) in vars.scldom.iter().enumerate() {
        if exclude.contains(&id) { continue; }
        if let SclDomType::Unknown { start } = scldom.scldom_type {
            fields.push(FieldNorm { name: format!("scldom {id}"), id, start, len: scldom.value.len() });
        }
    }
    fields
}

pub fn field_scale(vars: &Variables, field: &FieldNorm, x_a: &Col<f64>, x_b: &Col<f64>) -> f64 {
    // reference scale of a field from two states
    let (s, n) = (field.start, field.len);
    let mut vmax: f64 = 0.0;
    for x in [x_a, x_b] { for i in s..s + n { vmax = vmax.max(x[i].abs()); } }
    vmax.max(vars.scale_scldom.get(&field.id).copied().unwrap_or(0.0)).max(SCALE_FLOOR)
}

pub fn scaled_residual(vars: &Variables, fields: &[FieldNorm], x_old: &Col<f64>, x_new: &Col<f64>) -> (f64, usize) {
    // Picard residual: max over fields of rms(x_new - x_old) / s_f; returns the residual and the index of the largest field
    // with scale_off, the unscaled residual |x_new - x_old| / |x_old| (field index 0)
    // NaN if x_new is not finite
    if (0..x_new.nrows()).any(|i| !x_new[i].is_finite()) {
        return (f64::NAN, 0);
    }
    if vars.scale_off || fields.is_empty() {
        return ((x_new - x_old).norm_l2() / (x_old.norm_l2() + 1e-10), 0);
    }
    let mut res: f64 = 0.0;
    let mut fid = 0;
    for (k, field) in fields.iter().enumerate() {
        let (s, n) = (field.start, field.len);
        let mut sum = 0.0;
        for i in s..s + n { sum += (x_new[i] - x_old[i]).powi(2); }
        let rms = (sum / (n.max(1) as f64)).sqrt();
        let r = rms / field_scale(vars, field, x_old, x_new);
        if r > res || r.is_nan() { res = r; fid = k; }
    }
    (res, fid)
}
