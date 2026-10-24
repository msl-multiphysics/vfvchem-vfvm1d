use crate::base::prelude::*;
use crate::operator::core::flow_helper;
use crate::operator::core::flow_helper::FlowChain;
use crate::operator::prelude::*;
use crate::transient::prelude::*;
use faer::Col;
use faer::sparse::Triplet;
use std::collections::HashMap;

#[derive(Default)]
pub struct TransientFlow {
    // internal data
    pub itr_dom: Vec<usize>,             // domain
    pub itr_vel: HashMap<usize, usize>,  // velocity
    pub itr_den: HashMap<usize, usize>,  // density
    pub itr_flw: HashMap<usize, usize>,  // flow data; created in initial_start

    // boundary data
    pub vel_bnd: Vec<usize>,             // inlet boundary with velocity BC
    pub vel_vel: HashMap<usize, usize>,  // velocity
    pub mfr_bnd: Vec<usize>,             // inlet boundary with mass flow rate BC
    pub mfr_mfr: HashMap<usize, usize>,  // mass flow rate per unit area; positive into the domain

    // interface data
    pub fcnt_itf: Vec<usize>,            // interface with continuity IC

    // flow chains; built in initial_start
    pub flw_chain: Vec<FlowChain>,

    // operators
    pub oper_itr: Vec<OpSclDomDivergence>,
    pub oper_bnd_vel: Vec<OpSclBndDirichlet>,
    pub oper_bnd_mfr: Vec<OpSclBndMassFlowRate>,
    pub oper_itf_fcnt: Vec<OpSclItfDivergence>,
}

impl TransientFlow {
    pub fn new() -> Result<TransientFlow, VFVChemError> {
        let phys = TransientFlow::default();
        Ok(phys)
    }

    pub fn add_flow_dom(&mut self, dom_id: usize, vel_id: usize, den_id: usize) -> Result<(), VFVChemError> {
        self.itr_dom.push(dom_id);
        self.itr_vel.insert(dom_id, vel_id);
        self.itr_den.insert(dom_id, den_id);
        Ok(())
    }

    pub fn add_vel_bnd(&mut self, bnd_id: usize, vel_id: usize) -> Result<(), VFVChemError> {
        self.vel_bnd.push(bnd_id);
        self.vel_vel.insert(bnd_id, vel_id);
        Ok(())
    }

    pub fn add_mfr_bnd(&mut self, bnd_id: usize, mfr_id: usize) -> Result<(), VFVChemError> {
        self.mfr_bnd.push(bnd_id);
        self.mfr_mfr.insert(bnd_id, mfr_id);
        Ok(())
    }

    pub fn add_fcnt_itf(&mut self, itf_id: usize) -> Result<(), VFVChemError> {
        self.fcnt_itf.push(itf_id);
        Ok(())
    }

}

impl TransientBase for TransientFlow {
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError> {
        // iterate over unknown scalars
        let mut xid = 0;  // index in solution vector
        for (&dom_id, &vel_id) in &self.itr_vel {
            vars.scldom[vel_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_tot();
        }

        // create flow data
        // replace in place on a re-solve to avoid orphaned flow data
        for &dom_id in &self.itr_dom {
            let vel_id = self.itr_vel[&dom_id];
            let den_id = self.itr_den[&dom_id];
            match self.itr_flw.get(&dom_id) {
                Some(&flw_id) => { vars.sclflw[flw_id] = ScalarFlow::new(flw_id, &vars.dom[dom_id], vel_id, den_id)?; }
                None => { self.itr_flw.insert(dom_id, vars.add_sclflw(dom_id, vel_id, den_id)?); }
            }
        }

        // create flow chains from the inlets
        let flow_dom: Vec<(usize, usize, usize, usize)> = self.itr_dom.iter().map(|&dom_id| (dom_id, self.itr_vel[&dom_id], self.itr_den[&dom_id], self.itr_flw[&dom_id])).collect();
        let inlet_bnd: Vec<usize> = self.vel_bnd.iter().chain(self.mfr_bnd.iter()).copied().collect();
        self.flw_chain = flow_helper::build_flow_chains(vars, &flow_dom, &inlet_bnd, &self.fcnt_itf)?;

        // last xid is the matrix size
        Ok(xid)
    }

    fn initial_operator(&mut self, vars: &Variables) -> Result<(), VFVChemError> {
        // purge operators from a previous solve
        self.oper_itr.clear();
        self.oper_bnd_vel.clear();
        self.oper_bnd_mfr.clear();
        self.oper_itf_fcnt.clear();

        // internal and interface
        // continuity of each segment goes into the row of its downstream vertex, so the flow direction comes from the flow chains
        for chain in &self.flw_chain {
            for link in &chain.links {
                let oper_div = OpSclDomDivergence::new(vars, link.dom_id, link.vel_id, link.den_id, link.dir)?;
                self.oper_itr.push(oper_div);
            }
            for &(itf_id, up_is_a) in &chain.itfs {
                let dom_a_id = vars.itf[itf_id].dom_a_id;
                let dom_b_id = vars.itf[itf_id].dom_b_id;
                let (vel_a_id, vel_b_id) = (self.itr_vel[&dom_a_id], self.itr_vel[&dom_b_id]);
                let (den_a_id, den_b_id) = (self.itr_den[&dom_a_id], self.itr_den[&dom_b_id]);
                let oper_div = OpSclItfDivergence::new(vars, itf_id, vel_a_id, den_a_id, vel_b_id, den_b_id, up_is_a)?;
                self.oper_itf_fcnt.push(oper_div);
            }
        }

        // boundary (inlet)
        for &bnd_id in &self.vel_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.itr_vel[&dom_id];
            let bnd_vel_id = self.vel_vel[&bnd_id];
            let oper_dir = OpSclBndDirichlet::new(vars, bnd_id, dom_vel_id, bnd_vel_id)?;
            self.oper_bnd_vel.push(oper_dir);
        }
        for &bnd_id in &self.mfr_bnd {
            let dom_id = vars.bnd[bnd_id].dom_id;
            let dom_vel_id = self.itr_vel[&dom_id];
            let dom_den_id = self.itr_den[&dom_id];
            let bnd_mfr_id = self.mfr_mfr[&bnd_id];
            let oper_mfr = OpSclBndMassFlowRate::new(vars, bnd_id, dom_vel_id, dom_den_id, bnd_mfr_id)?;
            self.oper_bnd_mfr.push(oper_mfr);
        }

        Ok(())
    }

    fn assemble_matrix(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>, dt: f64) -> Result<(), VFVChemError> {
        // internal
        for oper_div in &self.oper_itr {
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_div.apply_time(vars, a_triplet, b_vec, dt)?;
        }

        // boundary
        for oper_dir in &self.oper_bnd_vel {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }
        for oper_mfr in &self.oper_bnd_mfr {
            oper_mfr.apply(vars, a_triplet, b_vec)?;
        }

        // interface
        for oper_div in &self.oper_itf_fcnt {
            oper_div.apply(vars, a_triplet, b_vec)?;
        }

        Ok(())
    }

    fn initial_dirichlet(&self, vars: &mut Variables) {
        // boundary
        for oper_dir in &self.oper_bnd_vel {
            oper_dir.apply_initial(vars);
        }
        for oper_mfr in &self.oper_bnd_mfr {
            oper_mfr.apply_initial(vars);
        }
    }

    fn update_flow(&self, vars: &mut Variables, dt: f64) {
        // face mass flux marched from the inlets; balances the mass change of every control volume
        flow_helper::update_flow(vars, &self.flw_chain, Some(dt));
    }
}
