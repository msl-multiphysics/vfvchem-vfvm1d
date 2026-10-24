use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndSource {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub src_id: usize,  // source (R_i)
}

impl OpSclBndSource {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, src_id: usize) -> Result<OpSclBndSource, VFVChemError> {
        // adds the source term to the control volumes of boundary nodes
        // d(m_i * c_i)/dt + d(m_i * c_i * u)/dx = -d/dx( -G_ij * dc_j/dx ) + R_i
        //                                                                    ^^^

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, src_id)?;

        // create struct
        let mut oper = OpSclBndSource::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.src_id = src_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclBndSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_seg R dx
        // R_v * dx

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let seg = &vars.seg[bnd.dom_id];
        let src_dom = &vars.scldom[self.src_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let len = seg.cvol_len[vid_v];

            // get vertex properties
            let src_v = src_dom.value[vid_v];

            // load into vector
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -src_v * len)?;
        }

        Ok(())
    }
}
