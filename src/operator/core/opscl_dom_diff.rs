use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomDiffusion {
    // domain
    pub dom_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (c_i)
    pub diff_id: usize,  // diffusion coefficient (G_ij)
    pub drv_id: usize,   // driving scalar (c_j)
}

impl OpSclDomDiffusion {
    pub fn new(vars: &Variables, dom_id: usize, unk_id: usize, diff_id: usize, drv_id: usize) -> Result<OpSclDomDiffusion, VFVChemError> {
        // adds the diffusion term to scalar transport equations
        // d(m_i * c_i)/dt + d(m_i * c_i * u)/dx = -d/dx( -G_ij * dc_j/dx ) + R_i
        //                                         ^^^^^^^^^^^^^^^^^^^^^^^^

        // input checking
        check_dom_scldom(vars, dom_id, unk_id)?;
        check_dom_scldom(vars, dom_id, diff_id)?;
        check_dom_scldom(vars, dom_id, drv_id)?;

        // create struct
        let mut oper = OpSclDomDiffusion::default();
        oper.dom_id = dom_id;
        oper.unk_id = unk_id;
        oper.diff_id = diff_id;
        oper.drv_id = drv_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclDomDiffusion {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_seg -d/dx( -G * dc/dx ) dx
        // [G * dc/dx * n]_faces
        // sum_ngh G * dc/dx * n
        //
        // apply normal derivative
        // sum_ngh G_f * (c_n - c_v) / dn

        // get objects
        let dom = &vars.dom[self.dom_id];
        let seg = &vars.seg[self.dom_id];
        let diff_dom = &vars.scldom[self.diff_id];

        // iterate through points
        // do not iterate through nodes to skip boundaries
        for vid_v in 0..dom.num_pts() {
            // get geometry
            let list_n = &seg.ngh_vid[vid_v];
            let dist_n = &seg.ngh_dist[vid_v];

            // get vertex properties
            let diff_v = diff_dom.value[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let dist = dist_n[loc];

                // get neighbor properties
                let diff_n = diff_dom.value[vid_n];

                // get face diffusion coefficient
                let diff_f = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean

                // load into matrix
                let coeff = diff_f / dist;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.drv_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.drv_id, vid_v, -coeff)?;
            }

        }

        Ok(())
    }
}
