use crate::base::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct Variables {
    // geometry
    pub pts: Points,
    pub nds: Nodes,
    pub seg: Vec<Segment>,

    // grouping
    pub dom: Vec<Domain>,
    pub bnd: Vec<Boundary>,
    pub itf: Vec<Interface>,

    // scalar
    pub scldom: Vec<ScalarDomain>,
    pub sclbnd: Vec<ScalarBoundary>,
    pub sclitf: Vec<ScalarInterface>,
    pub sclflw: Vec<ScalarFlow>,

    // convergence scaling
    // reference scales of unknown fields used in the Picard residual and the adaptive error norm (see scaling.rs)
    pub scale_scldom: HashMap<usize, f64>,  // scldom id -> reference scale
    pub scale_off: bool,                    // if true, the unscaled residual |x_new - x_old| / |x_old| is used
}

impl Variables {
    // geometry
    // point spacing grows linearly from h_min at the nodes to h_max at a distance ref_len
    pub fn new_bounds(x_min: f64, x_max: f64, h_min: f64, h_max: f64, ref_len: f64) -> Result<Variables, VFVChemError> {
        // one domain "dom" with boundaries "bnd_l" and "bnd_r"
        let mut vars = Variables::default();
        vars.nds = Nodes::new_from_bounds(x_min, x_max)?;
        vars.pts = Points::new(&vars.nds, h_min, h_max, ref_len)?;
        Ok(vars)
    }

    pub fn new_segments(nds_pos: &[f64], h_min: f64, h_max: f64, ref_len: f64) -> Result<Variables, VFVChemError> {
        // domains "dom_0", "dom_1", ... between consecutive positions
        // boundaries "bnd_l" and "bnd_r" at the ends; interfaces "itf_0", "itf_1", ... between consecutive domains
        let mut vars = Variables::default();
        vars.nds = Nodes::new_from_list(nds_pos)?;
        vars.pts = Points::new(&vars.nds, h_min, h_max, ref_len)?;
        Ok(vars)
    }

    pub fn new_nodes(pos: Vec<f64>) -> Result<Variables, VFVChemError> {
        // one domain "dom" with boundaries "bnd_l" and "bnd_r"
        // pos lists every vertex in increasing order; the first and last are the nodes and the rest are points
        // check that positions are increasing
        if pos.len() < 2 { return Err(VFVChemError::NodesNotIncreasing { idx: 0 }) }
        for (idx, &x) in pos.iter().enumerate() {
            if !x.is_finite() || (idx > 0 && !(x > pos[idx - 1])) { return Err(VFVChemError::NodesNotIncreasing { idx }) }
        }

        // store nodes
        let mut vars = Variables::default();
        vars.nds = Nodes::new_from_bounds(pos[0], pos[pos.len() - 1])?;

        // store internal points
        let len: Vec<f64> = pos.windows(2).map(|p| p[1] - p[0]).collect();
        vars.pts.h_min = len.iter().cloned().fold(f64::INFINITY, f64::min);
        vars.pts.h_max = len.iter().cloned().fold(0.0, f64::max);
        vars.pts.ref_len = pos[pos.len() - 1] - pos[0];
        vars.pts.pts_pos = pos[1..pos.len() - 1].to_vec();
        vars.pts.grp1d_gpt.insert("dom".to_string(), (0..pos.len() - 2).collect());
        Ok(vars)
    }

    // groupings
    pub fn add_dom(&mut self, grp1d_id: &str) -> Result<usize, VFVChemError> {
        let dom_id = self.dom.len();
        let dom = Domain::new(dom_id, &self.pts, &self.nds, grp1d_id)?;
        let seg = Segment::new(dom_id, &dom)?;  // seg_id = dom_id;
        self.dom.push(dom);
        self.seg.push(seg);
        Ok(dom_id)
    }

    pub fn add_bnd(&mut self, dom_id: usize, grp0d_id: &str) -> Result<usize, VFVChemError> {
        let bnd_id = self.bnd.len();
        let dom = &self.dom[dom_id];
        let bnd = Boundary::new(bnd_id, &self.nds, dom, grp0d_id)?;
        self.bnd.push(bnd);
        Ok(bnd_id)
    }

    pub fn add_itf(&mut self, dom_a_id: usize, dom_b_id: usize, grp0d_id: &str) -> Result<usize, VFVChemError> {
        let itf_id = self.itf.len();
        let dom_a = &self.dom[dom_a_id];
        let dom_b = &self.dom[dom_b_id];
        let itf = Interface::new(itf_id, &self.nds, dom_a, dom_b, grp0d_id)?;
        self.itf.push(itf);
        Ok(itf_id)
    }

    // scalar
    pub fn add_scldom_cnst(&mut self, dom_id: usize, value: f64, file_path: &str) -> Result<usize, VFVChemError> {
        let scldom_id = self.scldom.len();
        let dom = &self.dom[dom_id];
        let scldom = ScalarDomain::new_cnst(scldom_id, dom, value, file_path)?;
        self.scldom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_scldom_func(&mut self, dom_id: usize, func: SclDomFunc, unkn_id: Vec<usize>, file_path: &str) -> Result<usize, VFVChemError> {
        let scldom_id = self.scldom.len();
        let dom = &self.dom[dom_id];
        let scldom = ScalarDomain::new_func(scldom_id, dom, &self.scldom, func, unkn_id, file_path)?;
        self.scldom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_scldom_unkn(&mut self, dom_id: usize, value_init: f64, file_path: &str) -> Result<usize, VFVChemError> {
        let scldom_id = self.scldom.len();
        let dom = &self.dom[dom_id];
        let scldom = ScalarDomain::new_unkn(scldom_id, dom, value_init, file_path)?;
        self.scldom.push(scldom);
        Ok(scldom_id)
    }

    pub fn add_sclbnd_cnst(&mut self, bnd_id: usize, value: f64, file_path: &str) -> Result<usize, VFVChemError> {
        let sclbnd_id = self.sclbnd.len();
        let bnd = &self.bnd[bnd_id];
        let sclbnd = ScalarBoundary::new_cnst(sclbnd_id, bnd, value, file_path)?;
        self.sclbnd.push(sclbnd);
        Ok(sclbnd_id)
    }

    pub fn add_sclbnd_func(&mut self, bnd_id: usize, func: SclBndFunc, unkn_id: Vec<usize>, file_path: &str) -> Result<usize, VFVChemError> {
        let sclbnd_id = self.sclbnd.len();
        let bnd = &self.bnd[bnd_id];
        let sclbnd = ScalarBoundary::new_func(sclbnd_id, bnd, &self.scldom, func, unkn_id, file_path)?;
        self.sclbnd.push(sclbnd);
        Ok(sclbnd_id)
    }

    pub fn add_sclitf_cnst(&mut self, itf_id: usize, value: f64, file_path: &str) -> Result<usize, VFVChemError> {
        let sclitf_id = self.sclitf.len();
        let itf = &self.itf[itf_id];
        let sclitf = ScalarInterface::new_cnst(sclitf_id, itf, value, file_path)?;
        self.sclitf.push(sclitf);
        Ok(sclitf_id)
    }

    pub fn add_sclitf_func(&mut self, itf_id: usize, func: SclItfFunc, unkn_a_id: Vec<usize>, unkn_b_id: Vec<usize>, file_path: &str) -> Result<usize, VFVChemError> {
        let sclitf_id = self.sclitf.len();
        let itf = &self.itf[itf_id];
        let sclitf = ScalarInterface::new_func(sclitf_id, itf, &self.scldom, func, unkn_a_id, unkn_b_id, file_path)?;
        self.sclitf.push(sclitf);
        Ok(sclitf_id)
    }

    pub fn add_sclflw(&mut self, dom_id: usize, vel_id: usize, den_id: usize) -> Result<usize, VFVChemError> {
        let sclflw_id = self.sclflw.len();
        let dom = &self.dom[dom_id];
        let sclflw = ScalarFlow::new(sclflw_id, dom, vel_id, den_id)?;
        self.sclflw.push(sclflw);
        Ok(sclflw_id)
    }

    // convergence scaling
    pub fn set_scldom_scale(&mut self, scldom_id: usize, scale: f64) -> Result<(), VFVChemError> {
        if scale.is_nan() || scale <= 0.0 { return Err(VFVChemError::VariablesNonPositiveScale { scale }); }
        self.scale_scldom.insert(scldom_id, scale);
        Ok(())
    }

}
