use crate::base::prelude::*;

// A flow chain is a sequence of flow domains joined end to end by flow continuity interfaces, with one inlet
// boundary at one end of the chain. Flow goes from the inlet to the other end, which is the outlet.
//
// The face mass flux used by the transport operators is not interpolated from the vertex velocities; it is marched
// from the inlet so that every control volume satisfies continuity exactly
// F_out = F_in - V_v * d(rho_v)/dt
// The inlet flux is rho * u of the inlet node. At steady state F = rho_in * u_in on every face. In transient, the
// last flux equals rho * u of the outlet node up to the difference between the trapezoid (segment) and the
// control volume integration of d(rho)/dt, which vanishes for linear d(rho)/dt.

#[derive(Default, Clone)]
pub struct FlowLink {
    pub dom_id: usize,
    pub vel_id: usize,
    pub den_id: usize,
    pub flw_id: usize,
    pub dir: f64,  // +1 if the flow goes in +x; -1 otherwise
}

#[derive(Default, Clone)]
pub struct FlowChain {
    pub inlet_bnd_id: usize,
    pub links: Vec<FlowLink>,
    pub itfs: Vec<(usize, bool)>,  // (interface, upstream side is a) between consecutive links
}

pub fn build_flow_chains(vars: &Variables, flow_dom: &[(usize, usize, usize, usize)], inlet_bnd: &[usize], flow_itf: &[usize]) -> Result<Vec<FlowChain>, VFVChemError> {
    // flow_dom: (dom, vel, den, flw)
    // every flow domain must be in exactly one chain, reached from exactly one inlet

    // get the outward normal of a domain at a node
    let norm_at = |dom_id: usize, gnd: usize| -> f64 {
        let dom = &vars.dom[dom_id];
        vars.nds.grp1d_norm[&dom.grp1d_id][dom.nds_nid[&gnd]]
    };

    // follow every inlet downstream
    let mut visited = vec![false; flow_dom.len()];
    let mut chains = Vec::new();
    for &bnd_id in inlet_bnd {
        // get inlet node
        let bnd = &vars.bnd[bnd_id];
        let gnd_in = bnd.nds_gnd[0];

        // the inlet must not be on a flow interface
        let on_itf = flow_itf.iter().any(|&itf_id| {
            let itf = &vars.itf[itf_id];
            (itf.dom_a_id == bnd.dom_id || itf.dom_b_id == bnd.dom_id) && itf.nds_gnd.contains(&gnd_in)
        });
        if on_itf { return Err(VFVChemError::FlowInvalidInlet { dom_id: bnd.dom_id }); }

        // march through the domains of the chain
        let mut chain = FlowChain { inlet_bnd_id: bnd_id, ..Default::default() };
        let mut dom_id = bnd.dom_id;
        let mut dir = -bnd.nds_norm[0];  // the flow enters through the inlet
        loop {
            // register domain
            let Some(k) = flow_dom.iter().position(|f| f.0 == dom_id) else { return Err(VFVChemError::FlowMissingDomain { dom_id }); };
            if visited[k] { return Err(VFVChemError::FlowInvalidInlet { dom_id: bnd.dom_id }); }
            visited[k] = true;
            let (_, vel_id, den_id, flw_id) = flow_dom[k];
            chain.links.push(FlowLink { dom_id, vel_id, den_id, flw_id, dir });

            // get outlet node of this domain
            let dom = &vars.dom[dom_id];
            let Some(&gnd_out) = dom.nds_gnd.iter().find(|&&gnd| norm_at(dom_id, gnd) == dir) else { break };

            // find the next domain through a flow interface
            let next = flow_itf.iter().find_map(|&itf_id| {
                let itf = &vars.itf[itf_id];
                if !itf.nds_gnd.contains(&gnd_out) { return None; }
                if itf.dom_a_id == dom_id { Some((itf_id, true, itf.dom_b_id)) }
                else if itf.dom_b_id == dom_id { Some((itf_id, false, itf.dom_a_id)) }
                else { None }
            });
            match next {
                Some((itf_id, up_is_a, dom_next_id)) => {
                    chain.itfs.push((itf_id, up_is_a));
                    dir = -norm_at(dom_next_id, gnd_out);
                    dom_id = dom_next_id;
                }
                None => break,
            }
        }
        chains.push(chain);
    }

    // every flow domain must have an inlet
    if let Some(k) = visited.iter().position(|&v| !v) {
        return Err(VFVChemError::FlowNoInlet { dom_id: flow_dom[k].0 });
    }

    Ok(chains)
}

pub fn update_flow(vars: &mut Variables, chains: &[FlowChain], dt: Option<f64>) {
    // updates the face mass flux of every flow domain; marched from the inlet of each chain
    // dt is None at steady state

    // iterate through chains
    for chain in chains {
        // get inlet flux in +x
        let link_in = &chain.links[0];
        let bnd = &vars.bnd[chain.inlet_bnd_id];
        let dom_in = &vars.dom[bnd.dom_id];
        let vid_in = dom_in.num_pts() + dom_in.nds_nid[&bnd.nds_gnd[0]];
        let mut flux = vars.scldom[link_in.den_id].value[vid_in] * vars.scldom[link_in.vel_id].value[vid_in];

        // iterate through domains
        for link in &chain.links {
            // get objects
            let dom = &vars.dom[link.dom_id];
            let seg = &vars.seg[link.dom_id];
            let den_dom = &vars.scldom[link.den_id];

            // get vertices in order of x
            let mut vid_ord = Vec::with_capacity(dom.num_tot());
            vid_ord.push(seg.seg_vid[0][0]);
            vid_ord.extend(seg.seg_vid.iter().map(|&[_, vid_1]| vid_1));

            // get mass change rate of every control volume
            let mass_t: Vec<f64> = vid_ord.iter().map(|&vid| match dt {
                Some(dt) => seg.cvol_len[vid] * (den_dom.value[vid] - den_dom.value_prev[vid]) / dt,
                None => 0.0,
            }).collect();

            // march face flux in +x from the inlet end
            // face k is between vid_ord[k] and vid_ord[k + 1]
            let num_face = seg.seg_vid.len();
            let mut face_flux = vec![0.0; num_face];
            if link.dir > 0.0 {
                for k in 0..=num_face {
                    flux -= mass_t[k];
                    if k < num_face { face_flux[k] = flux; }
                }
            } else {
                for k in (0..=num_face).rev() {
                    flux += mass_t[k];
                    if k > 0 { face_flux[k - 1] = flux; }
                }
            }

            // get outward face mass flux of every control volume
            let mut flux_vec: Vec<Vec<f64>> = seg.ngh_vid.iter().map(|list_n| vec![0.0; list_n.len()]).collect();  // vid, loc -> F_f
            for (k, &[vid_0, vid_1]) in seg.seg_vid.iter().enumerate() {
                for (vid_v, vid_n) in [(vid_0, vid_1), (vid_1, vid_0)] {
                    let Some(loc) = seg.ngh_vid[vid_v].iter().position(|&vid| vid == vid_n) else { continue };
                    flux_vec[vid_v][loc] = seg.face_norm[vid_v][loc] * face_flux[k];
                }
            }

            // modify sclflw values
            vars.sclflw[link.flw_id].flux = flux_vec;
        }
    }
}
