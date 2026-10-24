use crate::base::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct Nodes {
    // nodes
    // nodes are on boundaries and interfaces
    pub nds_pos: Vec<f64>,  // gnd (global node id) -> x

    // 0d groups
    // defines boundaries and interfaces
    pub grp0d_gnd: HashMap<String, Vec<usize>>, // grp0d_id -> [gnd, ...]

    // 1d groups
    // defines domain and outward normals
    // left end first, then right end
    pub grp1d_gnd: HashMap<String, Vec<usize>>, // grp1d_id -> [gnd_l, gnd_r]
    pub grp1d_norm: HashMap<String, Vec<f64>>, // grp1d_id -> [-1.0, 1.0]
}

impl Nodes {
    pub fn new_from_bounds(x_min: f64, x_max: f64) -> Result<Nodes, VFVChemError> {
        // creates nodes at the ends of one domain "dom"
        // boundaries "bnd_l" and "bnd_r"

        // input checking
        if !(x_max > x_min && x_min.is_finite() && x_max.is_finite()) { return Err(VFVChemError::NodesInvalidBoundsX { x_min, x_max }) }

        // default struct
        let mut nds = Nodes::default();

        // store nodes
        nds.nds_pos = vec![x_min, x_max];

        // nodes in 0d groups
        nds.grp0d_gnd.insert("bnd_l".to_string(), vec![0]);
        nds.grp0d_gnd.insert("bnd_r".to_string(), vec![1]);

        // nodes and outward normals in 1d groups
        nds.grp1d_gnd.insert("dom".to_string(), vec![0, 1]);
        nds.grp1d_norm.insert("dom".to_string(), vec![-1.0, 1.0]);

        Ok(nds)
    }

    pub fn new_from_list(nds_pos: &[f64]) -> Result<Nodes, VFVChemError> {
        // creates nodes at the ends of domains "dom_0", "dom_1", ... between consecutive positions
        // boundaries "bnd_l" and "bnd_r" at the ends; interfaces "itf_0", "itf_1", ... between consecutive domains

        // input checking
        if nds_pos.len() < 2 { return Err(VFVChemError::NodesNotIncreasing { idx: 0 }) }
        for (idx, &x) in nds_pos.iter().enumerate() {
            if !x.is_finite() || (idx > 0 && !(x > nds_pos[idx - 1])) { return Err(VFVChemError::NodesNotIncreasing { idx }) }
        }

        // default struct
        let mut nds = Nodes::default();
        let num_dom = nds_pos.len() - 1;

        // store nodes
        nds.nds_pos = nds_pos.to_vec();

        // nodes in 0d groups
        nds.grp0d_gnd.insert("bnd_l".to_string(), vec![0]);
        nds.grp0d_gnd.insert("bnd_r".to_string(), vec![num_dom]);
        for k in 1..num_dom {
            nds.grp0d_gnd.insert(format!("itf_{}", k - 1), vec![k]);
        }

        // nodes and outward normals in 1d groups
        for k in 0..num_dom {
            nds.grp1d_gnd.insert(format!("dom_{k}"), vec![k, k + 1]);
            nds.grp1d_norm.insert(format!("dom_{k}"), vec![-1.0, 1.0]);
        }

        Ok(nds)
    }

    #[inline]
    pub fn num_nds(&self) -> usize {
        self.nds_pos.len()
    }
}
