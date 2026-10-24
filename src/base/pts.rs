use crate::base::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct Points {
    // spacing parameters
    // spacing grows linearly from h_min at the nodes to h_max at a distance ref_len
    pub h_min: f64,
    pub h_max: f64,
    pub ref_len: f64,

    // points
    // points are located between the nodes of a domain
    pub pts_pos: Vec<f64>,  // gpt (global point id) -> x

    // 1d groups
    // uses the same grp1d_id in nodes
    pub grp1d_gpt: HashMap<String, Vec<usize>>,  // grp1d_id -> [gpt, ...] in order of x
}

impl Points {
    pub fn new(nds: &Nodes, h_min: f64, h_max: f64, ref_len: f64) -> Result<Points, VFVChemError> {
        // input checking
        if !(h_min > 0.0 && h_min.is_finite()) { return Err(VFVChemError::PointsInvalidMinSpacing { h_min }) }
        if !(h_max >= h_min && h_max.is_finite()) { return Err(VFVChemError::PointsInvalidMaxSpacing { h_min, h_max }) }
        if !(ref_len > 0.0 && ref_len.is_finite()) { return Err(VFVChemError::PointsInvalidRefinementLength { ref_len }) }

        // default struct
        let mut pts = Points::default();
        pts.h_min = h_min;
        pts.h_max = h_max;
        pts.ref_len = ref_len;

        // iterate through 1d groups
        for (grp1d_id, nds_gnd) in nds.grp1d_gnd.iter() {
            // create internal points
            let x_l = nds.nds_pos[nds_gnd[0]];
            let x_r = nds.nds_pos[nds_gnd[1]];
            let pts_pos = _point_place(x_l, x_r, h_min, h_max, ref_len);
            let mut pts_id = Vec::new();

            // register internal points
            for &pos in pts_pos.iter() {
                let gpt = pts.pts_pos.len();
                pts.pts_pos.push(pos);
                pts_id.push(gpt);
            }
            pts.grp1d_gpt.insert(grp1d_id.clone(), pts_id);
        }

        Ok(pts)
    }

    #[inline]
    pub fn num_pts(&self) -> usize {
        self.pts_pos.len()
    }
}

fn _point_place(x_l: f64, x_r: f64, h_min: f64, h_max: f64, ref_len: f64) -> Vec<f64> {
    // places points between x_l and x_r with the local spacing
    // h(d) = h_min + (h_max - h_min) * min(d / ref_len, 1)
    // d is the distance to the nearer end, so interfaces are also refined
    //
    // the stretched coordinate s = int_0^d 1 / h(d') dd' counts the spacings from an end
    // the number of segments is the stretched length rounded to an integer, and the points are equally spaced in s

    // get stretched length of the domain
    let half = 0.5 * (x_r - x_l);
    let s_half = _stretch(half, h_min, h_max, ref_len);
    let num_seg = ((2.0 * s_half).round() as usize).max(1);

    // place points symmetrically from both ends
    let mut pts_pos = Vec::with_capacity(num_seg - 1);
    for k in 1..num_seg {
        let s = 2.0 * s_half * (k as f64) / (num_seg as f64);
        let x = if s <= s_half {
            x_l + _unstretch(s, h_min, h_max, ref_len)
        } else {
            x_r - _unstretch(2.0 * s_half - s, h_min, h_max, ref_len)
        };
        pts_pos.push(x);
    }

    pts_pos
}

fn _stretch(d: f64, h_min: f64, h_max: f64, ref_len: f64) -> f64 {
    // s(d) = int_0^d 1 / h(d') dd'

    // uniform spacing
    if h_max == h_min { return d / h_min; }

    // graded layer, then uniform spacing h_max
    let g = (h_max - h_min) / ref_len;
    let d_lay = d.min(ref_len);
    (1.0 + g * d_lay / h_min).ln() / g + (d - ref_len).max(0.0) / h_max
}

fn _unstretch(s: f64, h_min: f64, h_max: f64, ref_len: f64) -> f64 {
    // inverse of _stretch

    // uniform spacing
    if h_max == h_min { return s * h_min; }

    // graded layer, then uniform spacing h_max
    let g = (h_max - h_min) / ref_len;
    let s_lay = (h_max / h_min).ln() / g;
    if s <= s_lay {
        h_min * ((g * s).exp() - 1.0) / g
    } else {
        ref_len + (s - s_lay) * h_max
    }
}
