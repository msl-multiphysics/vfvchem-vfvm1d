use crate::base::prelude::*;
use std::collections::HashMap;

#[derive(Default)]
pub struct Boundary {
    // ids
    pub bnd_id: usize,
    pub dom_id: usize,  // domain this boundary is attached to
    pub grp0d_id: String,

    // nodes on boundary
    pub nds_gnd: Vec<usize>,  // idx (index) -> gnd; local to global
    pub nds_nid: HashMap<usize, usize>,  // gnd -> idx; global to local
    pub nds_idx: HashMap<usize, usize>,  // gnd -> idx; only nodes on this boundary
    pub nds_norm: Vec<f64>,  // idx -> outward normal (-1 at the left end and +1 at the right end)
}

impl Boundary {
    pub fn new(bnd_id: usize, nds: &Nodes, dom: &Domain, grp0d_id: &str) -> Result<Boundary, VFVChemError> {
        // input checking
        if !nds.grp0d_gnd.contains_key(grp0d_id) { return Err(VFVChemError::BoundaryGroupIDNotFound { grp0d_id: grp0d_id.to_string() })}

        // check that grp0d_id nodes are a subset of dom
        let grp0d_nds = &nds.grp0d_gnd[grp0d_id];
        let not_subset = !grp0d_nds.iter().all(|gnd| dom.nds_nid.contains_key(gnd));
        if not_subset { return Err(VFVChemError::BoundaryNoCommonBoundary { dom_id: dom.dom_id, grp0d_id: grp0d_id.to_string() })}

        // default struct
        let mut bnd = Boundary::default();
        bnd.bnd_id = bnd_id;
        bnd.dom_id = dom.dom_id;
        bnd.grp0d_id = grp0d_id.to_string();

        // get nodes
        bnd.nds_gnd = nds.grp0d_gnd[grp0d_id].clone();
        bnd.nds_nid = dom.nds_gnd.iter().enumerate().map(|(idx, &gnd)| (gnd, idx)).collect();
        bnd.nds_idx = bnd.nds_gnd.iter().enumerate().map(|(idx, &gnd)| (gnd, idx)).collect();

        // get outward normals
        let nds_norm = &nds.grp1d_norm[&dom.grp1d_id];
        bnd.nds_norm = bnd.nds_gnd.iter().map(|&gnd| nds_norm[bnd.nds_nid[&gnd]]).collect();

        Ok(bnd)
    }

    #[inline]
    pub fn num_nds(&self) -> usize {
        self.nds_gnd.len()
    }
}
