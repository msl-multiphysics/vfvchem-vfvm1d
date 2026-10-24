use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndMassFlowRate {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub vel_id: usize,  // velocity (u); unknown
    pub den_id: usize,  // density (rho)
    pub mfr_id: usize,  // inlet mass flux (G_set); positive into the domain
}

impl OpSclBndMassFlowRate {
    pub fn new(vars: &Variables, bnd_id: usize, vel_id: usize, den_id: usize, mfr_id: usize) -> Result<OpSclBndMassFlowRate, VFVChemError> {
        // sets the inlet mass flux
        // rho * u * n = -G_set

        // input checking
        check_bnd_scldom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, den_id)?;
        check_bnd_sclbnd(vars, bnd_id, mfr_id)?;

        // create struct
        let mut oper = OpSclBndMassFlowRate::default();
        oper.bnd_id = bnd_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.mfr_id = mfr_id;

        Ok(oper)
    }

    pub fn apply_initial(&self, vars: &mut Variables) {
        // loads the inlet velocity into the unknown
        // u = -G_set * n / rho

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get vertex properties
            let norm = bnd.nds_norm[idx_v];
            let den_v = vars.scldom[self.den_id].value[vid_v];
            let mfr_v = vars.sclbnd[self.mfr_id].value[idx_v];

            // load into unknown
            if den_v != 0.0 { vars.scldom[self.vel_id].value[vid_v] = -norm * mfr_v / den_v; }
        }
    }
}

impl OperatorBase for OpSclBndMassFlowRate {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let den_dom = &vars.scldom[self.den_id];
        let mfr_bnd = &vars.sclbnd[self.mfr_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get vertex properties
            let norm = bnd.nds_norm[idx_v];
            let den_v = den_dom.value[vid_v];
            let mfr_v = mfr_bnd.value[idx_v];

            // load into matrix
            self.add_sclscl_a(vars, a_triplet, self.vel_id, vid_v, self.vel_id, vid_v, den_v)?;
            self.add_scl_b(vars, b_vec, self.vel_id, vid_v, -norm * mfr_v)?;
        }

        Ok(())
    }
}
