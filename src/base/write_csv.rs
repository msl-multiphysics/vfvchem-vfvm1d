use crate::base::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn write_scldom_csv(scldom: &ScalarDomain, pts: &Points, nds: &Nodes, dom: &Domain, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", scldom.file_name, ts, scldom.file_type);
    let file = File::create(&file_path).map_err(|e| _csv_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // header
    writeln!(w, "x,value").map_err(|e| _csv_err(&file_path, e))?;

    // collect positions
    let mut pos = Vec::with_capacity(dom.num_tot());  // vid -> x
    for &gpt in &dom.pts_gpt {
        pos.push(pts.pts_pos[gpt]);  // vid = pid for points
    }
    for &gnd in &dom.nds_gnd {
        pos.push(nds.nds_pos[gnd]);  // vid = num_pts + nid for nodes
    }

    // write lines in order of x
    let mut vid_ord: Vec<usize> = (0..dom.num_tot()).collect();
    vid_ord.sort_by(|&vid_0, &vid_1| pos[vid_0].total_cmp(&pos[vid_1]));
    for vid in vid_ord {
        writeln!(w, "{:.9},{:.9}", pos[vid], scldom.value[vid]).map_err(|e| _csv_err(&file_path, e))?;
    }

    Ok(())
}

pub fn write_sclbnd_csv(sclbnd: &ScalarBoundary, nds: &Nodes, bnd: &Boundary, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", sclbnd.file_name, ts, sclbnd.file_type);
    let file = File::create(&file_path).map_err(|e| _csv_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // header
    writeln!(w, "x,value").map_err(|e| _csv_err(&file_path, e))?;

    // write lines
    for (idx, &gnd) in bnd.nds_gnd.iter().enumerate() {
        let x = nds.nds_pos[gnd];
        let value = sclbnd.value[idx];
        writeln!(w, "{:.9},{:.9}", x, value).map_err(|e| _csv_err(&file_path, e))?;
    }

    Ok(())
}

pub fn write_sclitf_csv(sclitf: &ScalarInterface, nds: &Nodes, itf: &Interface, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", sclitf.file_name, ts, sclitf.file_type);
    let file = File::create(&file_path).map_err(|e| _csv_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // header
    writeln!(w, "x,value").map_err(|e| _csv_err(&file_path, e))?;

    // write lines
    for (idx, &gnd) in itf.nds_gnd.iter().enumerate() {
        let x = nds.nds_pos[gnd];
        let value = sclitf.value[idx];
        writeln!(w, "{:.9},{:.9}", x, value).map_err(|e| _csv_err(&file_path, e))?;
    }

    Ok(())
}

fn _csv_err(file_path: &str, source: std::io::Error) -> VFVChemError {
    VFVChemError::WriteCsv { file_path: file_path.to_string(), source }
}
