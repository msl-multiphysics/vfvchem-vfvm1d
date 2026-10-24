use crate::base::prelude::*;
use faer::Col;
use faer::sparse::Triplet;

pub trait OperatorBase {
    // to implement in specific solvers
    fn apply(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, b_vec: &mut Col<f64>) -> Result<(), VFVChemError>;

    // matrix insertion convenience functions

    // scalar-scalar blocks
    fn add_sclscl_a(&self, vars: &Variables, a_triplet: &mut Vec<Triplet<usize, usize, f64>>, row_scl: usize, row_idx: usize, col_scl: usize, col_idx: usize, value: f64) -> Result<(), VFVChemError> {
        let row_start = match vars.scldom[row_scl].scldom_type {
            SclDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let col_start = match vars.scldom[col_scl].scldom_type {
            SclDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let row_xid = row_start + row_idx;
        let col_xid = col_start + col_idx;
        a_triplet.push(Triplet::new(row_xid, col_xid, value));
        Ok(())
    }
    fn add_scl_b(&self, vars: &Variables, b_vec: &mut Col<f64>, row_scl: usize, row_idx: usize, value: f64) -> Result<(), VFVChemError> {
        let row_start = match vars.scldom[row_scl].scldom_type {
            SclDomType::Unknown { start } => start,
            _ => return Err(VFVChemError::OperatorNonUnknownInLinearSystem{}),
        };
        let row_xid = row_start + row_idx;
        b_vec[row_xid] += value;
        Ok(())
    }

}

// input checking convenience functions

// domain
pub fn check_dom_scldom(vars: &Variables, dom_id: usize, scldom_id: usize) -> Result<(), VFVChemError> {
    let scldom = &vars.scldom[scldom_id];
    if scldom.dom_id != dom_id { return Err(VFVChemError::OperatorScalarDomainNotOnDomain { scldom_id, dom_id }) }
    Ok(())
}

// boundary
pub fn check_bnd_scldom(vars: &Variables, bnd_id: usize, scldom_id: usize) -> Result<(), VFVChemError> {
    let scldom = &vars.scldom[scldom_id];
    let bnd = &vars.bnd[bnd_id];
    if scldom.dom_id != bnd.dom_id { return Err(VFVChemError::OperatorScalarDomainNotOnBoundary { scldom_id, bnd_id }) }
    Ok(())
}
pub fn check_bnd_sclbnd(vars: &Variables, bnd_id: usize, sclbnd_id: usize) -> Result<(), VFVChemError> {
    let sclbnd = &vars.sclbnd[sclbnd_id];
    if sclbnd.bnd_id != bnd_id { return Err(VFVChemError::OperatorScalarBoundaryNotOnBoundary { sclbnd_id, bnd_id }) }
    Ok(())
}

// interface
pub fn check_itf_scldom(vars: &Variables, itf_id: usize, scldom_a_id: usize, scldom_b_id: usize) -> Result<(), VFVChemError> {
    let scldom_a = &vars.scldom[scldom_a_id];
    let scldom_b = &vars.scldom[scldom_b_id];
    let itf = &vars.itf[itf_id];
    let order_same = (itf.dom_a_id == scldom_a.dom_id) && (itf.dom_b_id == scldom_b.dom_id);
    let order_swap = (itf.dom_a_id == scldom_b.dom_id) && (itf.dom_b_id == scldom_a.dom_id);
    if !(order_same || order_swap) { return Err(VFVChemError::OperatorScalarDomainNotOnInterface { scldom_a_id, scldom_b_id, itf_id }) }
    Ok(())
}
pub fn check_itf_sclitf(vars: &Variables, itf_id: usize, sclitf_id: usize) -> Result<(), VFVChemError> {
    let sclitf = &vars.sclitf[sclitf_id];
    if sclitf.itf_id != itf_id { return Err(VFVChemError::OperatorScalarInterfaceNotOnInterface { sclitf_id, itf_id }) }
    Ok(())
}

// domain
pub fn check_dom_sclflw(vars: &Variables, dom_id: usize, sclflw_id: usize) -> Result<(), VFVChemError> {
    let sclflw = &vars.sclflw[sclflw_id];
    if sclflw.dom_id != dom_id { return Err(VFVChemError::OperatorScalarFlowNotOnDomain { sclflw_id, dom_id }) }
    Ok(())
}

// boundary
pub fn check_bnd_sclflw(vars: &Variables, bnd_id: usize, sclflw_id: usize) -> Result<(), VFVChemError> {
    let sclflw = &vars.sclflw[sclflw_id];
    let bnd = &vars.bnd[bnd_id];
    if sclflw.dom_id != bnd.dom_id { return Err(VFVChemError::OperatorScalarFlowNotOnBoundary { sclflw_id, bnd_id }) }
    Ok(())
}
