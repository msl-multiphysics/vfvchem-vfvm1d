use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndAdvection {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (c_i); unknown
    pub wgt_id: usize,   // weight coefficient (m_i)
    pub diff_id: usize,  // diffusion coefficient (G_ii); for the Peclet number
    pub vel_id: usize,   // velocity (u)
    pub den_id: usize,   // density (rho); converts the mass flux

    // flow
    pub flw_id: usize,  // flow data; for the face mass flux
}

impl OpSclBndAdvection {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, wgt_id: usize, diff_id: usize, vel_id: usize, den_id: usize, flw_id: usize) -> Result<OpSclBndAdvection, VFVChemError> {
        // adds the advection term to the control volumes of boundary nodes
        // m_i * dc_i/dt + m_i * u * dc_i/dx = -d/dx( -G_ii * dc_i/dx ) + R_i
        //                 ^^^^^^^^^^^^^^^^^
        // must be used together with OpSclBndDiffusion; the power law scheme corrects its central difference
        // the boundary face is upwinded from the node and contributes nothing

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_scldom(vars, bnd_id, wgt_id)?;
        check_bnd_scldom(vars, bnd_id, diff_id)?;
        check_bnd_scldom(vars, bnd_id, vel_id)?;
        check_bnd_scldom(vars, bnd_id, den_id)?;
        check_bnd_sclflw(vars, bnd_id, flw_id)?;

        // create struct
        let mut oper = OpSclBndAdvection::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.wgt_id = wgt_id;
        oper.diff_id = diff_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.flw_id = flw_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclBndAdvection {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_seg -m * u * dc/dx dx
        // [-m * (c - c_v) * u * n]_faces
        // sum_ngh -m * (c - c_v) * u * n + sum_bnd -m * (c - c_v) * u * n
        //
        // introduce weight flux from the mass flux
        // F = m * u * n = (m / rho) * rho * u * n
        // sum_ngh -F_f * (c_f - c_v) + sum_bnd -F_b * (c_v - c_v)
        // sum_ngh -F_f * (c_f - c_v)
        //
        // apply power law scheme
        // A(|Pe|) = max(0, (1 - 0.1 * |Pe|)^5)
        // sum_ngh (G_f * (A(|Pe_f|) - 1) + max(-F_f, 0)) * (c_n - c_v)
        //
        // non-conservative form; uniform c is preserved even if continuity is not satisfied

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let seg = &vars.seg[bnd.dom_id];
        let wgt_dom = &vars.scldom[self.wgt_id];
        let diff_dom = &vars.scldom[self.diff_id];
        let den_dom = &vars.scldom[self.den_id];
        let flow_dom = &vars.sclflw[self.flw_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get geometry
            let list_n = &seg.ngh_vid[vid_v];
            let dist_n = &seg.ngh_dist[vid_v];
            let flux_n = &flow_dom.flux[vid_v];

            // get vertex properties
            let wgt_v = wgt_dom.value[vid_v];
            let diff_v = diff_dom.value[vid_v];
            let den_v = den_dom.value[vid_v];

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let dist = dist_n[loc];

                // get neighbor properties
                let wgt_n = wgt_dom.value[vid_n];
                let diff_n = diff_dom.value[vid_n];
                let den_n = den_dom.value[vid_n];

                // get face weight, density, and diffusion coefficient
                let wgt_f = 0.5 * (wgt_v + wgt_n);  // arithmetic mean
                let den_f = 0.5 * (den_v + den_n);  // arithmetic mean
                let diff_f = if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 };  // harmonic mean

                // get power law term
                let cond = diff_f / dist;
                let flux = if den_f != 0.0 { wgt_f / den_f * flux_n[loc] } else { 0.0 };
                let peclet = if cond != 0.0 { flux / cond } else { 0.0 };
                let pow_law = (1.0 - 0.1 * peclet.abs()).max(0.0).powi(5);

                // load into matrix
                let coeff = cond * (pow_law - 1.0) + (-flux).max(0.0);
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_n,  coeff)?;
                self.add_sclscl_a(vars, a_triplet, self.unk_id, vid_v, self.unk_id, vid_v, -coeff)?;
            }

        }

        Ok(())
    }
}
