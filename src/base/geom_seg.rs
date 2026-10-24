use crate::base::prelude::*;

#[derive(Default)]
pub struct Segment {
    // ids
    pub seg_id: usize,
    pub dom_id: usize,  // domain this seg refers to

    // segment data
    // vid is the index of combined point and node vectors
    // i.e., vid = if vid < dom.num_pts() {pid} else {dom.num_pts() + nid}
    pub seg_vid: Vec<[usize; 2]>,  // [vid, vid] in order of x

    // control volume data
    pub cvol_len: Vec<f64>,  // vid (vertex id) -> CV (control volume) length

    // face data
    // faces are at the midpoints of the segments and have unit area
    pub face_norm: Vec<Vec<f64>>,  // vid, loc (local neighbor index) -> outward normal on face

    // neighbor vertex data
    pub ngh_vid: Vec<Vec<usize>>,  // vid -> neighbor vid
    pub ngh_dist: Vec<Vec<f64>>,  // vid -> segment length to neighbor
}

impl Segment {
    pub fn new(seg_id: usize, dom: &Domain) -> Result<Segment, VFVChemError> {
        // default struct
        let mut seg = Segment::default();
        seg.seg_id = seg_id;
        seg.dom_id = dom.dom_id;

        Ok(seg)
    }

    pub fn update_seg(&mut self, nds: &Nodes, pts: &Points, dom: &Domain) -> Result<(), VFVChemError> {
        // clear seg data
        self.seg_vid.clear();
        self.cvol_len = vec![0.0; dom.num_tot()];
        self.face_norm = vec![Vec::new(); dom.num_tot()];
        self.ngh_vid = vec![Vec::new(); dom.num_tot()];
        self.ngh_dist = vec![Vec::new(); dom.num_tot()];

        // collect points in domain
        let mut pts_pos = Vec::with_capacity(dom.num_tot());  // vid -> x
        for &gpt in &dom.pts_gpt {
            pts_pos.push(pts.pts_pos[gpt]);
        }
        for &gnd in &dom.nds_gnd {
            pts_pos.push(nds.nds_pos[gnd]);
        }

        // sort vertices in order of x
        // the nodes must be the ends of the domain
        let mut vid_ord: Vec<usize> = (0..dom.num_tot()).collect();
        vid_ord.sort_by(|&vid_0, &vid_1| pts_pos[vid_0].total_cmp(&pts_pos[vid_1]));
        let vid_l = dom.num_pts();
        let vid_r = dom.num_pts() + 1;
        if vid_ord[0] != vid_l { return Err(VFVChemError::SegmentCoincidentVertices { dom_id: dom.dom_id, x: pts_pos[vid_ord[0]] }) }
        if vid_ord[dom.num_tot() - 1] != vid_r { return Err(VFVChemError::SegmentCoincidentVertices { dom_id: dom.dom_id, x: pts_pos[vid_r] }) }

        // build segments between consecutive vertices
        // each vertex gets half of the segments on either side of it
        for idx in 0..dom.num_tot() - 1 {
            let vid_0 = vid_ord[idx];
            let vid_1 = vid_ord[idx + 1];
            let len = pts_pos[vid_1] - pts_pos[vid_0];
            if !(len > 0.0) { return Err(VFVChemError::SegmentCoincidentVertices { dom_id: dom.dom_id, x: pts_pos[vid_0] }) }

            // add segment
            self.seg_vid.push([vid_0, vid_1]);
            self.cvol_len[vid_0] += 0.5 * len;
            self.cvol_len[vid_1] += 0.5 * len;

            // add neighbors and faces
            self.ngh_vid[vid_0].push(vid_1);
            self.ngh_vid[vid_1].push(vid_0);
            self.ngh_dist[vid_0].push(len);
            self.ngh_dist[vid_1].push(len);
            self.face_norm[vid_0].push(1.0);
            self.face_norm[vid_1].push(-1.0);
        }

        Ok(())
    }
}
