use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclItfDivergence {
    // interface
    pub itf_id: usize,

    // domain a
    pub vel_a_id: usize,  // velocity (u^a); unknown
    pub den_a_id: usize,  // density (rho^a)

    // domain b
    pub vel_b_id: usize,  // velocity (u^b); unknown
    pub den_b_id: usize,  // density (rho^b)

    // flow direction
    pub up_is_a: bool,  // true if the flow goes from domain a to domain b
}

impl OpSclItfDivergence {
    pub fn new(vars: &Variables, itf_id: usize, vel_a_id: usize, den_a_id: usize, vel_b_id: usize, den_b_id: usize, up_is_a: bool) -> Result<OpSclItfDivergence, VFVChemError> {
        // adds the continuity of the mass flux through the interface
        // rho_dn * u_dn = rho_up * u_up
        // loaded into the row of the downstream side; the velocity jumps if the density does
        // the row of the upstream side is the outlet of its domain; see OpSclDomDivergence

        // input checking
        check_itf_scldom(vars, itf_id, vel_a_id, vel_b_id)?;
        check_itf_scldom(vars, itf_id, den_a_id, den_b_id)?;

        // create struct
        let mut oper = OpSclItfDivergence::default();
        oper.itf_id = itf_id;
        oper.vel_a_id = vel_a_id;
        oper.den_a_id = den_a_id;
        oper.vel_b_id = vel_b_id;
        oper.den_b_id = den_b_id;
        oper.up_is_a = up_is_a;

        Ok(oper)
    }

    fn _get_sides(&self) -> ((usize, usize), (usize, usize)) {
        // ((vel_up, den_up), (vel_dn, den_dn))
        if self.up_is_a {
            ((self.vel_a_id, self.den_a_id), (self.vel_b_id, self.den_b_id))
        } else {
            ((self.vel_b_id, self.den_b_id), (self.vel_a_id, self.den_a_id))
        }
    }
}

impl OperatorBase for OpSclItfDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // rho_dn * u_dn - rho_up * u_up

        // get objects
        let itf = &vars.itf[self.itf_id];
        let ((vel_up_id, den_up_id), (vel_dn_id, den_dn_id)) = self._get_sides();
        let dom_up = &vars.dom[vars.scldom[vel_up_id].dom_id];
        let dom_dn = &vars.dom[vars.scldom[vel_dn_id].dom_id];
        let den_up_dom = &vars.scldom[den_up_id];
        let den_dn_dom = &vars.scldom[den_dn_id];

        // iterate through nodes
        for idx_v in 0..itf.num_nds() {
            // get indices
            let gnd_v = itf.nds_gnd[idx_v];
            let vid_up_v = dom_up.num_pts() + dom_up.nds_nid[&gnd_v];
            let vid_dn_v = dom_dn.num_pts() + dom_dn.nds_nid[&gnd_v];

            // get vertex properties
            let den_up_v = den_up_dom.value[vid_up_v];
            let den_dn_v = den_dn_dom.value[vid_dn_v];

            // load into matrix
            self.add_sclscl_a(vars, a_triplet, vel_dn_id, vid_dn_v, vel_dn_id, vid_dn_v,  den_dn_v)?;
            self.add_sclscl_a(vars, a_triplet, vel_dn_id, vid_dn_v, vel_up_id, vid_up_v, -den_up_v)?;
        }

        Ok(())
    }
}
