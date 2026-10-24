use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomSource {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
    pub src_id: usize,  // source (R_i)
}

impl OpSclDomSource {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, src_id: usize) -> Result<OpSclDomSource, VFVChemError> {
        // adds the source term to scalar transport equations
        // d(m_i * c_i)/dt + d(m_i * c_i * u)/dx = -d/dx( -G_ij * dc_j/dx ) + R_i
        //                                                                    ^^^

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, src_id)?;

        // create struct
        let mut oper = OpSclDomSource::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.src_id = src_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclDomSource {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_seg R dx
        // R_v * dx

        // get objects
        let dom = &vars.dom[self.dom_id];
        let seg = &vars.seg[self.dom_id];
        let src_dom = &vars.scldom[self.src_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
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
