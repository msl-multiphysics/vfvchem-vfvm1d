use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomTimeUnity {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,  // unknown scalar (c_i)
}

impl OpSclDomTimeUnity {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize) -> Result<OpSclDomTimeUnity, VFVChemError> {
        // adds the time derivative term to scalar transport equations with unit weight
        // dc_i/dt + d(c_i * u)/dx = -d/dx( -G_ij * dc_j/dx ) + R_i
        // ^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;

        // create struct
        let mut oper = OpSclDomTimeUnity::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_seg -dc/dt dx
        // -(c_v - c_v_prev) * dx / dt
        //
        // lumped at the vertex

        // get objects
        let dom = &vars.dom[self.dom_id];
        let seg = &vars.seg[self.dom_id];
        let unk_dom = &vars.scldom[self.unk_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let len = seg.cvol_len[vid_v];

            // get vertex properties
            let unk_v_prev = unk_dom.value_prev[vid_v];

            // load into matrix and vector
            self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -len / dt)?;
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, -unk_v_prev * len / dt)?;
        }

        Ok(())
    }
}

impl OperatorBase for OpSclDomTimeUnity {
    fn apply(&self, _vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        Err(VFVChemError::OperatorUseApplyTime {})
    }
}
