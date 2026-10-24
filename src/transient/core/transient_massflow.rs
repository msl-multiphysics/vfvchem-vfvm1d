use crate::base::prelude::*;
use crate::operator::core::flow_helper;
use crate::operator::core::flow_helper::FlowChain;
use crate::operator::prelude::*;
use crate::transient::prelude::*;
use faer::Col;
use faer::sparse::Triplet;
use std::collections::HashMap;

#[derive(Default)]
pub struct TransientMassFlow {
    // mass internal data
    // indexed by component, then by domain
    pub num_comp: usize,                                       // number of components
    pub itr_conc: Vec<HashMap<usize, usize>>,                  // concentration
    pub itr_diff: Vec<HashMap<usize, HashMap<usize, usize>>>,  // diffusion coefficient; component j -> D_ij
    pub itr_mprd: Vec<HashMap<usize, usize>>,                  // mass production coefficient
    pub itr_msrc: Vec<HashMap<usize, usize>>,                  // mass source

    // mass boundary data
    pub conc_bnd: Vec<Vec<usize>>,              // boundary with concentration BC
    pub conc_conc: Vec<HashMap<usize, usize>>,  // concentration
    pub mflx_bnd: Vec<Vec<usize>>,              // boundary with diffusive mass flux BC
    pub mflx_mflx: Vec<HashMap<usize, usize>>,  // diffusive mass flux
    pub mout_bnd: Vec<Vec<usize>>,              // boundary with mass outflow BC

    // mass interface data
    pub mcnt_itf: Vec<Vec<usize>>,              // interface with continuity IC

    // flow internal data
    pub itr_fdom: Vec<usize>,            // domain
    pub itr_vel: HashMap<usize, usize>,  // velocity
    pub itr_den: HashMap<usize, usize>,  // density
    pub itr_flw: HashMap<usize, usize>,  // flow data; created in initial_start

    // flow boundary data
    pub vel_bnd: Vec<usize>,             // inlet boundary with velocity BC
    pub vel_vel: HashMap<usize, usize>,  // velocity
    pub mfr_bnd: Vec<usize>,             // inlet boundary with mass flow rate BC
    pub mfr_mfr: HashMap<usize, usize>,  // mass flow rate per unit area; positive into the domain

    // flow interface data
    pub fcnt_itf: Vec<usize>,            // interface with continuity IC

    // flow chains; built in initial_start
    pub flw_chain: Vec<FlowChain>,

    // mass operators
    pub oper_itr_mass: Vec<(Vec<OpSclDomDiffusion>, OpSclDomProduction, OpSclDomSource)>,
    pub oper_itr_mtim: Vec<OpSclDomTimeUnity>,
    pub oper_itr_mtpr: Vec<OpSclDomTimePrevious>,
    pub oper_itr_madv: Vec<OpSclDomAdvectionUnity>,
    pub oper_bnd_conc: Vec<OpSclBndDirichlet>,
    pub oper_bnd_mflx: Vec<(Vec<OpSclBndDiffusion>, OpSclBndProduction, OpSclBndSource, OpSclBndNeumann)>,
    pub oper_bnd_mtim: Vec<OpSclBndTimeUnity>,
    pub oper_bnd_mtpr: Vec<OpSclBndTimePrevious>,
    pub oper_bnd_madv: Vec<OpSclBndAdvectionUnity>,
    pub oper_bnd_mout: Vec<(OpSclBndTimePrevious, Vec<OpSclBndDiffusion>, OpSclBndAdvectionUnity, OpSclBndProduction, OpSclBndSource)>,
    pub oper_itf_mcnt: Vec<(OpSclItfContinuity, Vec<OpSclItfDiffusion>, OpSclItfProduction, OpSclItfSource)>,
    pub oper_itf_mtim: Vec<OpSclItfTimeUnity>,
    pub oper_itf_mtpr: Vec<OpSclItfTimePrevious>,
    pub oper_itf_madv: Vec<OpSclItfAdvectionUnity>,

    // flow operators
    pub oper_itr_flow: Vec<OpSclDomDivergence>,
    pub oper_bnd_vel: Vec<OpSclBndDirichlet>,
    pub oper_bnd_mfr: Vec<OpSclBndMassFlowRate>,
    pub oper_itf_fcnt: Vec<OpSclItfDivergence>,
}

impl TransientMassFlow {
    pub fn new(num_comp: usize) -> Result<TransientMassFlow, VFVChemError> {
        // mass transfer is solved in conservative form minus c_i times continuity
        // rho * dc_i/dt + rho * u * dc_i/dx = -d/dx( -D_ij * dc_j/dx ) + R_i
        // on domains with flow, c_i is per unit mass and D_ij includes the density (rho * D_ij)
        // uniform c_i is preserved even if discrete continuity is not satisfied; c_i is conserved once it is
        // on domains without flow, c_i is per unit volume and the time derivative has unit weight

        // input checking
        if num_comp == 0 { return Err(VFVChemError::MassNonPositiveNumComponents { num_comp }) }

        // create struct with data for each component
        let mut phys = TransientMassFlow::default();
        phys.num_comp = num_comp;
        phys.itr_conc = vec![HashMap::new(); num_comp];
        phys.itr_diff = vec![HashMap::new(); num_comp];
        phys.itr_mprd = vec![HashMap::new(); num_comp];
        phys.itr_msrc = vec![HashMap::new(); num_comp];
        phys.conc_bnd = vec![Vec::new(); num_comp];
        phys.conc_conc = vec![HashMap::new(); num_comp];
        phys.mflx_bnd = vec![Vec::new(); num_comp];
        phys.mflx_mflx = vec![HashMap::new(); num_comp];
        phys.mout_bnd = vec![Vec::new(); num_comp];
        phys.mcnt_itf = vec![Vec::new(); num_comp];
        Ok(phys)
    }

    pub fn add_mass_dom(&mut self, comp_id: usize, dom_id: usize, conc_id: usize, diff_ids: HashMap<usize, usize>, mprd_id: usize, msrc_id: usize) -> Result<(), VFVChemError> {
        for &comp in diff_ids.keys().chain([&comp_id]) {
            if comp >= self.num_comp { return Err(VFVChemError::MassInvalidComponent { comp, num_comp: self.num_comp }) }
        }
        self.itr_conc[comp_id].insert(dom_id, conc_id);
        self.itr_diff[comp_id].insert(dom_id, diff_ids);
        self.itr_mprd[comp_id].insert(dom_id, mprd_id);
        self.itr_msrc[comp_id].insert(dom_id, msrc_id);
        Ok(())
    }

    pub fn add_conc_bnd(&mut self, comp_id: usize, bnd_id: usize, conc_id: usize) -> Result<(), VFVChemError> {
        if comp_id >= self.num_comp { return Err(VFVChemError::MassInvalidComponent { comp: comp_id, num_comp: self.num_comp }) }
        self.conc_bnd[comp_id].push(bnd_id);
        self.conc_conc[comp_id].insert(bnd_id, conc_id);
        Ok(())
    }

    pub fn add_mflx_bnd(&mut self, comp_id: usize, bnd_id: usize, mflx_id: usize) -> Result<(), VFVChemError> {
        if comp_id >= self.num_comp { return Err(VFVChemError::MassInvalidComponent { comp: comp_id, num_comp: self.num_comp }) }
        self.mflx_bnd[comp_id].push(bnd_id);
        self.mflx_mflx[comp_id].insert(bnd_id, mflx_id);
        Ok(())
    }

    pub fn add_mout_bnd(&mut self, comp_id: usize, bnd_id: usize) -> Result<(), VFVChemError> {
        if comp_id >= self.num_comp { return Err(VFVChemError::MassInvalidComponent { comp: comp_id, num_comp: self.num_comp }) }
        self.mout_bnd[comp_id].push(bnd_id);
        Ok(())
    }

    pub fn add_mcnt_itf(&mut self, comp_id: usize, itf_id: usize) -> Result<(), VFVChemError> {
        if comp_id >= self.num_comp { return Err(VFVChemError::MassInvalidComponent { comp: comp_id, num_comp: self.num_comp }) }
        self.mcnt_itf[comp_id].push(itf_id);
        Ok(())
    }

    pub fn add_flow_dom(&mut self, dom_id: usize, vel_id: usize, den_id: usize) -> Result<(), VFVChemError> {
        self.itr_fdom.push(dom_id);
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

impl TransientBase for TransientMassFlow {
    fn initial_start(&mut self, vars: &mut Variables) -> Result<usize, VFVChemError> {
        // iterate over unknown scalars
        let mut xid = 0;  // index in solution vector
        for (&dom_id, &vel_id) in &self.itr_vel {
            vars.scldom[vel_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_tot();
        }
        for itr_conc in &self.itr_conc {
        for (&dom_id, &conc_id) in itr_conc {
            vars.scldom[conc_id].scldom_type = SclDomType::Unknown { start: xid };
            xid += vars.dom[dom_id].num_tot();
        }}

        // create flow data
        // replace in place on a re-solve to avoid orphaned flow data
        for &dom_id in &self.itr_fdom {
            let vel_id = self.itr_vel[&dom_id];
            let den_id = self.itr_den[&dom_id];
            match self.itr_flw.get(&dom_id) {
                Some(&flw_id) => { vars.sclflw[flw_id] = ScalarFlow::new(flw_id, &vars.dom[dom_id], vel_id, den_id)?; }
                None => { self.itr_flw.insert(dom_id, vars.add_sclflw(dom_id, vel_id, den_id)?); }
            }
        }

        // create flow chains from the inlets
        let flow_dom: Vec<(usize, usize, usize, usize)> = self.itr_fdom.iter().map(|&dom_id| (dom_id, self.itr_vel[&dom_id], self.itr_den[&dom_id], self.itr_flw[&dom_id])).collect();
        let inlet_bnd: Vec<usize> = self.vel_bnd.iter().chain(self.mfr_bnd.iter()).copied().collect();
        self.flw_chain = flow_helper::build_flow_chains(vars, &flow_dom, &inlet_bnd, &self.fcnt_itf)?;

        // last xid is the matrix size
        Ok(xid)
    }

    fn initial_operator(&mut self, vars: &Variables) -> Result<(), VFVChemError> {
        // purge operators from a previous solve
        self.oper_itr_mass.clear();
        self.oper_itr_mtim.clear();
        self.oper_itr_mtpr.clear();
        self.oper_itr_madv.clear();
        self.oper_bnd_conc.clear();
        self.oper_bnd_mflx.clear();
        self.oper_bnd_mtim.clear();
        self.oper_bnd_mtpr.clear();
        self.oper_bnd_madv.clear();
        self.oper_bnd_mout.clear();
        self.oper_itf_mcnt.clear();
        self.oper_itf_mtim.clear();
        self.oper_itf_mtpr.clear();
        self.oper_itf_madv.clear();
        self.oper_itr_flow.clear();
        self.oper_bnd_vel.clear();
        self.oper_bnd_mfr.clear();
        self.oper_itf_fcnt.clear();

        for comp_i in 0..self.num_comp {
            // mass internal
            for (&dom_id, &conc_id) in &self.itr_conc[comp_i] {
                let mprd_id = self.itr_mprd[comp_i][&dom_id];
                let msrc_id = self.itr_msrc[comp_i][&dom_id];
                let mut oper_diff = Vec::new();
                for (&comp_j, &diff_id) in &self.itr_diff[comp_i][&dom_id] {
                    let conc_j_id = self.itr_conc[comp_j][&dom_id];
                    oper_diff.push(OpSclDomDiffusion::new(vars, dom_id, conc_id, diff_id, conc_j_id)?);
                }
                let oper_prd = OpSclDomProduction::new(vars, dom_id, conc_id, mprd_id)?;
                let oper_src = OpSclDomSource::new(vars, dom_id, conc_id, msrc_id)?;
                self.oper_itr_mass.push((oper_diff, oper_prd, oper_src));

                // advection only on domains with flow; the Peclet number uses the self-diffusivity
                // the time derivative is weighted by the previous density on domains with flow
                match self.itr_flw.get(&dom_id) {
                    Some(&flw_id) => {
                        let diff_id = self.itr_diff[comp_i][&dom_id].get(&comp_i).copied();
                        let den_id = self.itr_den[&dom_id];
                        let oper_time = OpSclDomTimePrevious::new(vars, dom_id, conc_id, den_id)?;
                        let oper_adv = OpSclDomAdvectionUnity::new(vars, dom_id, conc_id, diff_id, flw_id)?;
                        self.oper_itr_mtpr.push(oper_time);
                        self.oper_itr_madv.push(oper_adv);
                    }
                    None => {
                        let oper_time = OpSclDomTimeUnity::new(vars, dom_id, conc_id)?;
                        self.oper_itr_mtim.push(oper_time);
                    }
                }
            }

            // mass boundary
            for &bnd_id in &self.conc_bnd[comp_i] {
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let bnd_conc_id = self.conc_conc[comp_i][&bnd_id];
                let oper_dir = OpSclBndDirichlet::new(vars, bnd_id, dom_conc_id, bnd_conc_id)?;
                self.oper_bnd_conc.push(oper_dir);
            }
            for &bnd_id in &self.mflx_bnd[comp_i] {
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let dom_mprd_id = self.itr_mprd[comp_i][&dom_id];
                let dom_msrc_id = self.itr_msrc[comp_i][&dom_id];
                let bnd_mflx_id = self.mflx_mflx[comp_i][&bnd_id];
                let mut oper_diff = Vec::new();
                for (&comp_j, &dom_diff_id) in &self.itr_diff[comp_i][&dom_id] {
                    let dom_conc_j_id = self.itr_conc[comp_j][&dom_id];
                    oper_diff.push(OpSclBndDiffusion::new(vars, bnd_id, dom_conc_id, dom_diff_id, dom_conc_j_id)?);
                }
                let oper_prd = OpSclBndProduction::new(vars, bnd_id, dom_conc_id, dom_mprd_id)?;
                let oper_src = OpSclBndSource::new(vars, bnd_id, dom_conc_id, dom_msrc_id)?;
                let oper_neu = OpSclBndNeumann::new(vars, bnd_id, dom_conc_id, bnd_mflx_id)?;
                self.oper_bnd_mflx.push((oper_diff, oper_prd, oper_src, oper_neu));

                // advection only on domains with flow; the Neumann value is the diffusive flux
                // the time derivative is weighted by the previous density on domains with flow
                match self.itr_flw.get(&dom_id) {
                    Some(&dom_flw_id) => {
                        let dom_diff_id = self.itr_diff[comp_i][&dom_id].get(&comp_i).copied();
                        let dom_vel_id = self.itr_vel[&dom_id];
                        let dom_den_id = self.itr_den[&dom_id];
                        let oper_time = OpSclBndTimePrevious::new(vars, bnd_id, dom_conc_id, dom_den_id)?;
                        let oper_adv = OpSclBndAdvectionUnity::new(vars, bnd_id, dom_conc_id, dom_diff_id, dom_vel_id, dom_flw_id)?;
                        self.oper_bnd_mtpr.push(oper_time);
                        self.oper_bnd_madv.push(oper_adv);
                    }
                    None => {
                        let oper_time = OpSclBndTimeUnity::new(vars, bnd_id, dom_conc_id)?;
                        self.oper_bnd_mtim.push(oper_time);
                    }
                }
            }
            for &bnd_id in &self.mout_bnd[comp_i] {
                // zero diffusive flux; mass leaves by advection only
                let dom_id = vars.bnd[bnd_id].dom_id;
                let dom_conc_id = self.itr_conc[comp_i][&dom_id];
                let dom_mprd_id = self.itr_mprd[comp_i][&dom_id];
                let dom_msrc_id = self.itr_msrc[comp_i][&dom_id];
                let dom_vel_id = self.itr_vel[&dom_id];
                let dom_den_id = self.itr_den[&dom_id];
                let dom_flw_id = self.itr_flw[&dom_id];
                let oper_time = OpSclBndTimePrevious::new(vars, bnd_id, dom_conc_id, dom_den_id)?;
                let mut oper_diff = Vec::new();
                for (&comp_j, &dom_diff_id) in &self.itr_diff[comp_i][&dom_id] {
                    let dom_conc_j_id = self.itr_conc[comp_j][&dom_id];
                    oper_diff.push(OpSclBndDiffusion::new(vars, bnd_id, dom_conc_id, dom_diff_id, dom_conc_j_id)?);
                }
                let dom_diff_id = self.itr_diff[comp_i][&dom_id].get(&comp_i).copied();
                let oper_adv = OpSclBndAdvectionUnity::new(vars, bnd_id, dom_conc_id, dom_diff_id, dom_vel_id, dom_flw_id)?;
                let oper_prd = OpSclBndProduction::new(vars, bnd_id, dom_conc_id, dom_mprd_id)?;
                let oper_src = OpSclBndSource::new(vars, bnd_id, dom_conc_id, dom_msrc_id)?;
                self.oper_bnd_mout.push((oper_time, oper_diff, oper_adv, oper_prd, oper_src));
            }

            // mass interface
            for &itf_id in &self.mcnt_itf[comp_i] {
                let dom_a_id = vars.itf[itf_id].dom_a_id;
                let dom_b_id = vars.itf[itf_id].dom_b_id;

                let dom_conc_a_id = self.itr_conc[comp_i][&dom_a_id];
                let dom_conc_b_id = self.itr_conc[comp_i][&dom_b_id];
                let dom_mprd_a_id = self.itr_mprd[comp_i][&dom_a_id];
                let dom_mprd_b_id = self.itr_mprd[comp_i][&dom_b_id];
                let dom_msrc_a_id = self.itr_msrc[comp_i][&dom_a_id];
                let dom_msrc_b_id = self.itr_msrc[comp_i][&dom_b_id];
                let oper_cont = OpSclItfContinuity::new(vars, itf_id, dom_conc_a_id, dom_conc_b_id)?;
                let mut oper_diff = Vec::new();
                for (&comp_j, &dom_diff_a_id) in &self.itr_diff[comp_i][&dom_a_id] {
                    let dom_diff_b_id = self.itr_diff[comp_i][&dom_b_id][&comp_j];
                    let dom_conc_j_a_id = self.itr_conc[comp_j][&dom_a_id];
                    let dom_conc_j_b_id = self.itr_conc[comp_j][&dom_b_id];
                    oper_diff.push(OpSclItfDiffusion::new(vars, itf_id, dom_conc_a_id, dom_diff_a_id, dom_conc_j_a_id, dom_conc_b_id, dom_diff_b_id, dom_conc_j_b_id)?);
                }
                let oper_prd = OpSclItfProduction::new(vars, itf_id, dom_conc_a_id, dom_mprd_a_id, dom_conc_b_id, dom_mprd_b_id)?;
                let oper_src = OpSclItfSource::new(vars, itf_id, dom_conc_a_id, dom_msrc_a_id, dom_conc_b_id, dom_msrc_b_id)?;
                self.oper_itf_mcnt.push((oper_cont, oper_diff, oper_prd, oper_src));

                // the time derivative is weighted by the previous density if both sides have flow
                // c is per unit mass with flow and per unit volume without, so both sides must match
                match (self.itr_den.get(&dom_a_id), self.itr_den.get(&dom_b_id)) {
                    (Some(&den_a_id), Some(&den_b_id)) => {
                        let oper_time = OpSclItfTimePrevious::new(vars, itf_id, dom_conc_a_id, den_a_id, dom_conc_b_id, den_b_id)?;
                        self.oper_itf_mtpr.push(oper_time);
                    }
                    (None, None) => {
                        let oper_time = OpSclItfTimeUnity::new(vars, itf_id, dom_conc_a_id, dom_conc_b_id)?;
                        self.oper_itf_mtim.push(oper_time);
                    }
                    _ => return Err(VFVChemError::MassInterfaceMixedFlow { itf_id }),
                }

                // advection on each side with flow
                for dom_id in [dom_a_id, dom_b_id] {
                    if let Some(&flw_id) = self.itr_flw.get(&dom_id) {
                        let diff_id = self.itr_diff[comp_i][&dom_id].get(&comp_i).copied();
                        let vel_id = self.itr_vel[&dom_id];
                        let oper_adv = OpSclItfAdvectionUnity::new(vars, itf_id, dom_conc_a_id, dom_conc_b_id, diff_id, vel_id, flw_id)?;
                        self.oper_itf_madv.push(oper_adv);
                    }
                }
            }
        }

        // flow internal and interface
        // continuity of each segment goes into the row of its downstream vertex, so the flow direction comes from the flow chains
        for chain in &self.flw_chain {
            for link in &chain.links {
                let oper_div = OpSclDomDivergence::new(vars, link.dom_id, link.vel_id, link.den_id, link.dir)?;
                self.oper_itr_flow.push(oper_div);
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

        // flow boundary (inlet)
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
        // flow internal
        for oper_div in &self.oper_itr_flow {
            oper_div.apply(vars, a_triplet, b_vec)?;
            oper_div.apply_time(vars, a_triplet, b_vec, dt)?;
        }

        // flow boundary
        for oper_dir in &self.oper_bnd_vel {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }
        for oper_mfr in &self.oper_bnd_mfr {
            oper_mfr.apply(vars, a_triplet, b_vec)?;
        }

        // flow interface
        for oper_div in &self.oper_itf_fcnt {
            oper_div.apply(vars, a_triplet, b_vec)?;
        }

        // mass internal
        for (oper_diff, oper_prd, oper_src) in &self.oper_itr_mass {
            for oper in oper_diff { oper.apply(vars, a_triplet, b_vec)?; }
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }
        for oper_time in &self.oper_itr_mtim {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_time in &self.oper_itr_mtpr {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_adv in &self.oper_itr_madv {
            oper_adv.apply(vars, a_triplet, b_vec)?;
        }

        // mass boundary
        for oper_dir in &self.oper_bnd_conc {
            oper_dir.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_diff, oper_prd, oper_src, oper_neu) in &self.oper_bnd_mflx {
            for oper in oper_diff { oper.apply(vars, a_triplet, b_vec)?; }
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
            oper_neu.apply(vars, a_triplet, b_vec)?;
        }
        for oper_time in &self.oper_bnd_mtim {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_time in &self.oper_bnd_mtpr {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_adv in &self.oper_bnd_madv {
            oper_adv.apply(vars, a_triplet, b_vec)?;
        }
        for (oper_time, oper_diff, oper_adv, oper_prd, oper_src) in &self.oper_bnd_mout {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
            for oper in oper_diff { oper.apply(vars, a_triplet, b_vec)?; }
            oper_adv.apply(vars, a_triplet, b_vec)?;
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }

        // mass interface
        for (oper_cont, oper_diff, oper_prd, oper_src) in &self.oper_itf_mcnt {
            oper_cont.apply(vars, a_triplet, b_vec)?;
            for oper in oper_diff { oper.apply(vars, a_triplet, b_vec)?; }
            oper_prd.apply(vars, a_triplet, b_vec)?;
            oper_src.apply(vars, a_triplet, b_vec)?;
        }
        for oper_time in &self.oper_itf_mtim {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_time in &self.oper_itf_mtpr {
            oper_time.apply_time(vars, a_triplet, b_vec, dt)?;
        }
        for oper_adv in &self.oper_itf_madv {
            oper_adv.apply(vars, a_triplet, b_vec)?;
        }

        Ok(())
    }

    fn initial_dirichlet(&self, vars: &mut Variables) {
        // boundary
        for oper_dir in &self.oper_bnd_conc {
            oper_dir.apply_initial(vars);
        }
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
