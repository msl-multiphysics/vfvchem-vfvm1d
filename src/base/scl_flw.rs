use crate::base::prelude::*;

#[derive(Default)]
pub struct ScalarFlow {
    // ids
    pub sclflw_id: usize,
    pub dom_id: usize,  // domain this is attached to

    // scalars
    pub vel_id: usize,  // velocity (u)
    pub den_id: usize,  // density (rho)

    // flow data
    // shared by the continuity and transport operators so that both use the same face mass flux
    pub flux: Vec<Vec<f64>>,  // vid, loc -> outward face mass flux
}

impl ScalarFlow {
    pub fn new(sclflw_id: usize, dom: &Domain, vel_id: usize, den_id: usize) -> Result<ScalarFlow, VFVChemError> {
        // default class
        let mut sclflw = ScalarFlow::default();
        sclflw.sclflw_id = sclflw_id;
        sclflw.dom_id = dom.dom_id;

        // field ids
        sclflw.vel_id = vel_id;
        sclflw.den_id = den_id;

        // flow data is computed by the solver after the segments; see flow_helper
        sclflw.flux = Vec::new();

        Ok(sclflw)
    }

}
