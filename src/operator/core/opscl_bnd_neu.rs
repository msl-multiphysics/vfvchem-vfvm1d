use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclBndNeumann {
    // boundary
    pub bnd_id: usize,

    // scalars
    pub unk_id: usize,   // unknown scalar (c_i)
    pub flux_id: usize,  // flux (N_set)
}

impl OpSclBndNeumann {
    pub fn new(vars: &Variables, bnd_id: usize, unk_id: usize, flux_id: usize) -> Result<OpSclBndNeumann, VFVChemError> {
        // sets the outward boundary flux
        // N_i * n = N_set

        // input checking
        check_bnd_scldom(vars, bnd_id, unk_id)?;
        check_bnd_sclbnd(vars, bnd_id, flux_id)?;

        // create struct
        let mut oper = OpSclBndNeumann::default();
        oper.bnd_id = bnd_id;
        oper.unk_id = unk_id;
        oper.flux_id = flux_id;

        Ok(oper)
    }
}

impl OperatorBase for OpSclBndNeumann {
    fn apply(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // -N_v
        // the boundary face has unit area

        // get objects
        let bnd = &vars.bnd[self.bnd_id];
        let dom = &vars.dom[bnd.dom_id];
        let flux_bnd = &vars.sclbnd[self.flux_id];

        // iterate through nodes
        for idx_v in 0..bnd.num_nds() {
            // get indices
            let gnd_v = bnd.nds_gnd[idx_v];
            let nid_v = dom.nds_nid[&gnd_v];
            let vid_v = dom.num_pts() + nid_v;

            // get vertex properties
            let flux_v = flux_bnd.value[idx_v];

            // load into vector
            self.add_scl_b(vars, b_vec, self.unk_id, vid_v, flux_v)?;
        }

        Ok(())
    }
}
