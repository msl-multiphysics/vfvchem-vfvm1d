use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfProduction {
    // interface
    pub itf_id: usize,

    // scalars - domain a
    pub unk_a_id: usize,  // unknown scalar (c_i^a)
    pub prd_a_id: usize,  // production coefficient (s_i^a)

    // scalars - domain b
    pub unk_b_id: usize,  // unknown scalar (c_i^b)
    pub prd_b_id: usize,  // production coefficient (s_i^b)
}

impl OpSclItfProduction {
    pub fn new(vars: &Variables, itf_id: usize, unk_a_id: usize, prd_a_id: usize, unk_b_id: usize, prd_b_id: usize) -> Result<OpSclItfProduction, VFVChemError> {
        // adds the production term to the combined control volumes of interface nodes
        // d(m_i * c_i)/dt + d(m_i * c_i * u)/dx = -d/dx( -G_ij * dc_j/dx ) + R_i + s_i * c_i
        //                                                                          ^^^^^^^^^
        // s_i is positive if production and negative if consumption
        // s_i is lagged (evaluated with the previous iterate); s_i * c_i is implicit on the diagonal
        // each domain contributes its half of the control volume; loads into the rows of domain b
        // c_i^a = c_i^b by continuity, so both halves act on c_i^b

        // input checking
        check_itf_scldom(vars, itf_id, unk_a_id, unk_b_id)?;
        check_itf_scldom(vars, itf_id, prd_a_id, prd_b_id)?;

        // create struct
        let mut oper = OpSclItfProduction::default();
        oper.itf_id = itf_id;
        oper.unk_a_id = unk_a_id;
        oper.prd_a_id = prd_a_id;
        oper.unk_b_id = unk_b_id;
        oper.prd_b_id = prd_b_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclItfProduction {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_seg s * c dx
        // s_v * c_v * dx

        // get objects
        let itf = &vars.itf[self.itf_id];
        let dom_b = &vars.dom[vars.scldom[self.unk_b_id].dom_id];

        // iterate through both sides of the interface
        for prd_id in [self.prd_a_id, self.prd_b_id] {
            // get objects on this side
            let dom_id = vars.scldom[prd_id].dom_id;
            let dom = &vars.dom[dom_id];
            let seg = &vars.seg[dom_id];
            let prd_dom = &vars.scldom[prd_id];

            // iterate through nodes
            for idx_v in 0..itf.num_nds() {
                // get indices
                let gnd_v = itf.nds_gnd[idx_v];
                let vid_v = dom.num_pts() + dom.nds_nid[&gnd_v];
                let vid_b_v = dom_b.num_pts() + dom_b.nds_nid[&gnd_v];

                // get geometry
                let len = seg.cvol_len[vid_v];

                // get vertex properties
                let prd_v = prd_dom.value[vid_v];

                // load into matrix
                self.add_sclscl_a(vars, a_triplet, self.unk_b_id, vid_b_v, self.unk_b_id, vid_b_v, prd_v * len)?;
            }
        }

        Ok(())
    }
}
