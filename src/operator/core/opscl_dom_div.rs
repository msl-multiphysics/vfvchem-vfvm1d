use crate::base::prelude::*;
use crate::operator::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

#[derive(Default)]
pub struct OpSclDomDivergence {
    // domain
    pub dom_id: usize,

    // scalars
    pub vel_id: usize,  // velocity (u); unknown
    pub den_id: usize,  // density (rho)

    // flow direction
    pub dir: f64,  // +1 if the flow goes in +x; -1 otherwise
}

impl OpSclDomDivergence {
    pub fn new(vars: &Variables, dom_id: usize, vel_id: usize, den_id: usize, dir: f64) -> Result<OpSclDomDivergence, VFVChemError> {
        // adds the continuity equation of every segment of the domain
        // d(rho)/dt = -d(rho * u)/dx
        //             ^^^^^^^^^^^^^^
        // each segment is integrated with the trapezoid rule and loaded into the row of its downstream vertex
        // this covers the points and the outlet node; the inlet node row is left for the inlet condition
        // the system is lower triangular in marching order, so there is no checkerboarding

        // input checking
        check_dom_scldom(vars, dom_id, vel_id)?;
        check_dom_scldom(vars, dom_id, den_id)?;

        // create struct
        let mut oper = OpSclDomDivergence::default();
        oper.dom_id = dom_id;
        oper.vel_id = vel_id;
        oper.den_id = den_id;
        oper.dir = dir;

        Ok(oper)
    }

    pub fn apply_time(&self, vars: &Variables, _a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // int_seg d(rho)/dt dx from the upstream vertex n to the downstream vertex v
        // ((rho_v - rho_v_prev) + (rho_n - rho_n_prev)) * dist / (2 * dt)
        //
        // multiply the equation by dir as in apply and move to the right-hand side

        // get objects
        let dom = &vars.dom[self.dom_id];
        let seg = &vars.seg[self.dom_id];
        let den_dom = &vars.scldom[self.den_id];

        // iterate through vertices
        for vid_v in 0..dom.num_tot() {
            // get upstream neighbor
            let Some(loc) = seg.face_norm[vid_v].iter().position(|&norm| norm == -self.dir) else { continue };
            let vid_n = seg.ngh_vid[vid_v][loc];
            let dist = seg.ngh_dist[vid_v][loc];

            // get vertex properties
            let den_t_v = (den_dom.value[vid_v] - den_dom.value_prev[vid_v]) / dt;
            let den_t_n = (den_dom.value[vid_n] - den_dom.value_prev[vid_n]) / dt;

            // load into vector
            self.add_scl_b(vars, b_vec, self.vel_id, vid_v, -self.dir * 0.5 * dist * (den_t_v + den_t_n))?;
        }

        Ok(())
    }
}

impl OperatorBase for OpSclDomDivergence {
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, _b_vec: &mut Col<f64>) -> Result<(), VFVChemError> {
        // int_seg d(rho * u)/dx dx from the upstream vertex n to the downstream vertex v
        // dir * (rho_v * u_v - rho_n * u_n)
        //
        // multiply the equation by dir so that the diagonal is positive
        // rho_v * u_v - rho_n * u_n

        // get objects
        let dom = &vars.dom[self.dom_id];
        let seg = &vars.seg[self.dom_id];
        let den_dom = &vars.scldom[self.den_id];

        // iterate through vertices
        for vid_v in 0..dom.num_tot() {
            // get upstream neighbor
            let Some(loc) = seg.face_norm[vid_v].iter().position(|&norm| norm == -self.dir) else { continue };
            let vid_n = seg.ngh_vid[vid_v][loc];

            // get vertex properties
            let den_v = den_dom.value[vid_v];
            let den_n = den_dom.value[vid_n];

            // load into matrix
            self.add_sclscl_a(vars, a_triplet, self.vel_id, vid_v, self.vel_id, vid_v,  den_v)?;
            self.add_sclscl_a(vars, a_triplet, self.vel_id, vid_v, self.vel_id, vid_n, -den_n)?;
        }

        Ok(())
    }
}
