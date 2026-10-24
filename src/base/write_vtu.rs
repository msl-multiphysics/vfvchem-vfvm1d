use crate::base::prelude::*;
use std::fs::File;
use std::io::{BufWriter, Write};

pub fn write_scldom_vtu(scldom: &ScalarDomain, pts: &Points, nds: &Nodes, dom: &Domain, dual: &Segment, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", scldom.file_name, ts, scldom.file_type);
    let file = File::create(&file_path).map_err(|e| _vtu_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // start of file
    writeln!(w, r#"<?xml version="1.0"?>"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", dom.num_tot(), dual.seg_vid.len()).map_err(|e| _vtu_err(&file_path, e))?;

    // write points
    writeln!(w, "<Points>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for &gpt in &dom.pts_gpt {
        let x = pts.pts_pos[gpt];
        writeln!(w, "{:.9} 0.0 0.0", x).map_err(|e| _vtu_err(&file_path, e))?;
    }
    for &gnd in &dom.nds_gnd {
        let x = nds.nds_pos[gnd];
        writeln!(w, "{:.9} 0.0 0.0", x).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Points>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cells (segments)
    writeln!(w, "<Cells>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="connectivity" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for &[vid_0, vid_1] in &dual.seg_vid {
        writeln!(w, "{} {}", vid_0, vid_1).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="offsets" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for seg in 0..dual.seg_vid.len() {
        writeln!(w, "{}", 2 * (seg + 1)).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="types" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for _ in 0..dual.seg_vid.len() {
        writeln!(w, "3").map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</Cells>").map_err(|e| _vtu_err(&file_path, e))?;

    // write point data
    writeln!(w, r#"<PointData Scalars="value">"#).map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" Name="value" format="ascii">"#).map_err(|e| _vtu_err(&file_path, e))?;
    for vid in 0..dom.num_tot() {
        writeln!(w, "{:.9}", scldom.value[vid]).map_err(|e| _vtu_err(&file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</PointData>").map_err(|e| _vtu_err(&file_path, e))?;

    // write cell data (none)
    writeln!(w, "<CellData>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</CellData>").map_err(|e| _vtu_err(&file_path, e))?;

    // end of file
    writeln!(w, "</Piece>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</UnstructuredGrid>").map_err(|e| _vtu_err(&file_path, e))?;
    writeln!(w, "</VTKFile>").map_err(|e| _vtu_err(&file_path, e))?;

    Ok(())
}

pub fn write_sclbnd_vtu(sclbnd: &ScalarBoundary, nds: &Nodes, bnd: &Boundary, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", sclbnd.file_name, ts, sclbnd.file_type);
    let file = File::create(&file_path).map_err(|e| _vtu_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // write nodes as vertex cells
    _write_vertex_vtu(&mut w, &file_path, nds, &bnd.nds_gnd, &sclbnd.value)
}

pub fn write_sclitf_vtu(sclitf: &ScalarInterface, nds: &Nodes, itf: &Interface, ts: usize) -> Result<(), VFVChemError> {
    // initialize file writer
    let file_path = format!("{}_{:04}.{}", sclitf.file_name, ts, sclitf.file_type);
    let file = File::create(&file_path).map_err(|e| _vtu_err(&file_path, e))?;
    let mut w = BufWriter::new(file);

    // write nodes as vertex cells
    _write_vertex_vtu(&mut w, &file_path, nds, &itf.nds_gnd, &sclitf.value)
}

fn _write_vertex_vtu(w: &mut BufWriter<File>, file_path: &str, nds: &Nodes, nds_gnd: &[usize], value: &[f64]) -> Result<(), VFVChemError> {
    // start of file
    writeln!(w, r#"<?xml version="1.0"?>"#).map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#).map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, "<UnstructuredGrid>").map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, "<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">", nds_gnd.len(), nds_gnd.len()).map_err(|e| _vtu_err(file_path, e))?;

    // write points
    writeln!(w, "<Points>").map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" NumberOfComponents="3" format="ascii">"#).map_err(|e| _vtu_err(file_path, e))?;
    for &gnd in nds_gnd {
        let x = nds.nds_pos[gnd];
        writeln!(w, "{:.9} 0.0 0.0", x).map_err(|e| _vtu_err(file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, "</Points>").map_err(|e| _vtu_err(file_path, e))?;

    // write cells (degenerate to 1 point per cell)
    writeln!(w, "<Cells>").map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="connectivity" format="ascii">"#).map_err(|e| _vtu_err(file_path, e))?;
    for idx in 0..nds_gnd.len() {
        writeln!(w, "{}", idx).map_err(|e| _vtu_err(file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="offsets" format="ascii">"#).map_err(|e| _vtu_err(file_path, e))?;
    for idx in 0..nds_gnd.len() {
        writeln!(w, "{}", idx + 1).map_err(|e| _vtu_err(file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, r#"<DataArray type="Int32" Name="types" format="ascii">"#).map_err(|e| _vtu_err(file_path, e))?;
    for _ in 0..nds_gnd.len() {
        writeln!(w, "1").map_err(|e| _vtu_err(file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, "</Cells>").map_err(|e| _vtu_err(file_path, e))?;

    // write point data
    writeln!(w, r#"<PointData Scalars="value">"#).map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, r#"<DataArray type="Float64" Name="value" format="ascii">"#).map_err(|e| _vtu_err(file_path, e))?;
    for idx in 0..nds_gnd.len() {
        writeln!(w, "{:.9}", value[idx]).map_err(|e| _vtu_err(file_path, e))?;
    }
    writeln!(w, "</DataArray>").map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, "</PointData>").map_err(|e| _vtu_err(file_path, e))?;

    // write cell data (none)
    writeln!(w, "<CellData>").map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, "</CellData>").map_err(|e| _vtu_err(file_path, e))?;

    // end of file
    writeln!(w, "</Piece>").map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, "</UnstructuredGrid>").map_err(|e| _vtu_err(file_path, e))?;
    writeln!(w, "</VTKFile>").map_err(|e| _vtu_err(file_path, e))?;

    Ok(())
}

fn _vtu_err(file_path: &str, source: std::io::Error) -> VFVChemError {
    VFVChemError::WriteVtu { file_path: file_path.to_string(), source }
}
