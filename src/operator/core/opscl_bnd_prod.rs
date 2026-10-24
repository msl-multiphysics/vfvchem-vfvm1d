use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndProduction {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub prd_id: usize,  // production coefficient (s_i)
}

impl OpSclBndProduction {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, prd_id: usize) -> Result<OpSclBndProduction, VFVChemError> {
        // adds the production term to the control volumes of boundary nodes
        // d(m_i * c_i)/dt + d(m_i * c_i * u)/dx = -d/dx( -G_ij * dc_j/dx ) + R_i + s_i * c_i
        //                                                                          ^^^^^^^^^
        // s_i is positive if production and negative if consumption
        // s_i is lagged (evaluated with the previous iterate); s_i * c_i is implicit on the diagonal

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, prd_id)?;

        // create struct
        let mut oper = OpSclBndProduction::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.prd_id = prd_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclBndProduction {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_seg s * c dx
        // s_v * c_v * dx

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let seg = &vars.seg[bnd.dom_id];
        let prd_dom = &vars.scldom[self.prd_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let len = seg.cvol_len[vid_v];

            // get vertex properties
            let prd_v = prd_dom.value[vid_v];

            // load into matrix
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, prd_v * len)?;
        }

        Ok(())
    }
}
