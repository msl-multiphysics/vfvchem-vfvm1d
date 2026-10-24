use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndAdvectionUnity {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,           // unknown scalar (c_i); unknown
    pub diff_id: Option<usize>,  // diffusion coefficient (G_ii); for the Peclet number
    pub vel_id: usize,           // velocity (u)

    // flow
    pub flw_id: usize,  // flow data; for the face mass flux
}

impl OpSclBndAdvectionUnity {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, diff_id: Option<usize>, vel_id: usize, flw_id: usize) -> Result<OpSclBndAdvectionUnity, VFVChemError> {
        // adds the advection term in terms of the mass flux to the control volumes of boundary nodes
        // rho * dc_i/dt + rho * u * dc_i/dx = -d/dx( -G_ij * dc_j/dx ) + R_i
        //                 ^^^^^^^^^^^^^^^^^
        // must be used together with OpSclBndDiffusion; the power law scheme corrects its central difference
        // a missing diffusion coefficient sets G_f to zero and falls back to pure upwinding
        // the boundary face is upwinded from the node and contributes nothing

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        if let Some(diff_id) = diff_id {
            check_bnd_scldom(vars, bnd_id, diff_id)?;
        }
        check_bnd_scldom(vars, bnd_id, vel_id)?;
        check_bnd_sclflw(vars, bnd_id, flw_id)?;

        // create struct
        let mut oper = OpSclBndAdvectionUnity::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.diff_id = diff_id;
        oper.vel_id = vel_id;
        oper.flw_id = flw_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclBndAdvectionUnity {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_seg -rho * u * dc/dx dx
        // [-rho * (c - c_v) * u * n]_faces
        // sum_ngh -rho * (c - c_v) * u * n + sum_bnd -rho * (c - c_v) * u * n
        //
        // introduce the mass flux
        // F = rho * u * n
        // sum_ngh -F_f * (c_f - c_v) + sum_bnd -F_b * (c_v - c_v)
        // sum_ngh -F_f * (c_f - c_v)
        //
        // apply power law scheme
        // A(|Pe|) = max(0, (1 - 0.1 * |Pe|)^5)
        // sum_ngh (G_f * (A(|Pe_f|) - 1) + max(-F_f, 0)) * (c_n - c_v)

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let seg = &vars.seg[bnd.dom_id];
        let diff_dom = self.diff_id.map(|diff_id| &vars.scldom[diff_id]);
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

            // iterate through faces
            for (loc, &vid_n) in list_n.iter().enumerate() {
                // get geometry
                let dist = dist_n[loc];

                // get face diffusion coefficient; zero if absent
                let diff_f = match diff_dom {
                    Some(diff_dom) => {
                        let diff_v = diff_dom.value[vid_v];
                        let diff_n = diff_dom.value[vid_n];
                        if diff_v + diff_n != 0.0 { 2.0 * diff_v * diff_n / (diff_v + diff_n) } else { 0.0 }  // harmonic mean
                    }
                    None => 0.0,
                };

                // get power law term
                let cond = diff_f / dist;
                let flux = flux_n[loc];
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
