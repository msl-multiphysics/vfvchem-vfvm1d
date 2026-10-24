use thiserror::Error;

#[derive(Debug, Error)]
#[non_exhaustive]
pub enum VFVChemError {
    // placeholder
    #[error("Placeholder")]
    Placeholder {},

    // vertices
    #[error("Points: Minimum spacing `h_min` must be positive and finite. Got h_min = {h_min}.")]
    PointsInvalidMinSpacing {
        h_min: f64,
    },
    #[error("Points: Maximum spacing `h_max` must be finite and at least h_min. Got h_min = {h_min} and h_max = {h_max}.")]
    PointsInvalidMaxSpacing {
        h_min: f64,
        h_max: f64,
    },
    #[error("Points: Refinement layer thickness `ref_len` must be positive and finite. Got ref_len = {ref_len}.")]
    PointsInvalidRefinementLength {
        ref_len: f64,
    },
    #[error("Nodes: x_max must be greater than x_min. Got x_min = {x_min} and x_max = {x_max}.")]
    NodesInvalidBoundsX {
        x_min: f64,
        x_max: f64,
    },
    #[error("Nodes: Positions must be finite and strictly increasing; at least two are needed. Failed at index {idx}.")]
    NodesNotIncreasing {
        idx: usize,
    },

    // geometry
    #[error("Domain: Group ID {grp1d_id} not found.")]
    DomainGroupIDNotFound {
        grp1d_id: String,
    },
    #[error("Boundary: Group ID {grp0d_id} not found in Nodes.")]
    BoundaryGroupIDNotFound {
        grp0d_id: String,
    },
    #[error("Boundary: Domain {dom_id} does not end at group ID {grp0d_id}.")]
    BoundaryNoCommonBoundary {
        dom_id: usize,
        grp0d_id: String,
    },
    #[error("Interface: Group ID {grp0d_id} not found in Nodes.")]
    InterfaceGroupIDNotFound {
        grp0d_id: String,
    },
    #[error("Interface: Domains {dom_a} and {dom_b} do not share an interface at group ID {grp0d_id}.")]
    InterfaceNoCommonInterface {
        dom_a: usize,
        dom_b: usize,
        grp0d_id: String,
    },
    #[error("Segment: Vertices of domain {dom_id} coincide or are out of order at x = {x}.")]
    SegmentCoincidentVertices {
        dom_id: usize,
        x: f64,
    },

    // scalars
    #[error("Scalar{geom}: Function must depend on unknown-type ScalarDomains only. ScalarDomain ID {scldom_id} is not unknown-type.")]
    ScalarFunctionDependsOnNonUnknown {
        geom: &'static str,
        scldom_id: usize,
    },
    #[error("Scalar{geom}: Output file must be CSV or VTU. Found {file_type}.")]
    ScalarUnsupportedOutputType {
        geom: &'static str,
        file_type: String,
    },
    #[error("Scalar{geom}: Cannot deduce output file type from path: {file_path}.")]
    ScalarUnknownOutputType {
        geom: &'static str,
        file_path: String,
    },

    // read and write
    #[error("WriteCsv: Cannot write file '{file_path}': {source}")]
    WriteCsv {
        file_path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("WriteVtu: Cannot write file '{file_path}': {source}")]
    WriteVtu {
        file_path: String,
        #[source]
        source: std::io::Error,
    },

    // solvers
    #[error("Solver: Time step must be greater than zero. Got dt = {dt}.")]
    SolverNonPositiveTimeStep {
        dt: f64,
    },
    #[error("Solver: Number of time steps must be greater than zero. Got num_ts = {num_ts}.")]
    SolverNonPositiveNumTimeSteps {
        num_ts: usize,
    },
    #[error("Solver: Write frequency must be greater than zero. Got num_ts_write = {num_ts_write}.")]
    SolverNonPositiveWriteFrequency {
        num_ts_write: usize,
    },
    #[error("Solver: Maximum iterations must be greater than zero. Got max_iter = {max_iter}.")]
    SolverNonPositiveMaxIter {
        max_iter: usize,
    },
    #[error("Solver: Tolerance must be greater than zero. Got tol = {tol}.")]
    SolverNonPositiveTolerance {
        tol: f64,
    },
    #[error("Variables: Reference scale must be greater than zero. Got scale = {scale}.")]
    VariablesNonPositiveScale {
        scale: f64,
    },
    #[error("Solver: Damping must be in (0, 1]. Got damp = {damp}.")]
    SolverInvalidDamping {
        damp: f64,
    },
    #[error("SteadyBase: Nonlinear iteration did not converge within {max_iter} iterations. Residual = {res:.4e}.")]
    SteadyDidNotConverge {
        max_iter: usize,
        res: f64,
    },
    #[error("TransientBase: Nonlinear iteration did not converge at time step {ts} within {max_iter} iterations. Residual = {res:.4e}.")]
    TransientDidNotConverge {
        ts: usize,
        max_iter: usize,
        res: f64,
    },
    #[error("SteadyBase: Nonlinear iteration diverged; the solution is not finite (NaN or infinity) at iteration {iter}.")]
    SteadyNonFiniteResidual {
        iter: usize,
    },
    #[error("TransientBase: Nonlinear iteration diverged; the solution is not finite (NaN or infinity) at time step {ts}, iteration {iter}.")]
    TransientNonFiniteResidual {
        ts: usize,
        iter: usize,
    },
    #[error("Mass: Number of components must be greater than zero. Got num_comp = {num_comp}.")]
    MassNonPositiveNumComponents {
        num_comp: usize,
    },
    #[error("Mass: Component index must be less than the number of components. Got comp = {comp} and num_comp = {num_comp}.")]
    MassInvalidComponent {
        comp: usize,
        num_comp: usize,
    },
    #[error("Mass: Interface ID {itf_id} joins a domain with flow and a domain without flow. Concentrations are per unit mass on domains with flow and per unit volume otherwise.")]
    MassInterfaceMixedFlow {
        itf_id: usize,
    },
    #[error("Flow: Domain chain containing Domain ID {dom_id} has no inlet. Call add_vel_bnd or add_mfr_bnd on one end of every chain of flow domains.")]
    FlowNoInlet {
        dom_id: usize,
    },
    #[error("Flow: Domain chain containing Domain ID {dom_id} has more than one inlet, or the inlet is not at an end of the chain. Continuity only allows one inlet per chain.")]
    FlowInvalidInlet {
        dom_id: usize,
    },
    #[error("Flow: Domain ID {dom_id} has no flow data. Call add_flow_dom first.")]
    FlowMissingDomain {
        dom_id: usize,
    },

    // linear solvers
    #[error("LinearLu: Could not construct symbolic sparse matrix.")]
    LinearLuSymbolicError {
        #[source]
        source: faer::sparse::CreationError,
    },
    #[error("LinearLu: Could not assemble sparse matrix.")]
    LinearLuMatrixError {
        #[source]
        source: faer::sparse::FaerError,
    },
    #[error("LinearLu: Could not compute symbolic LU factorization.")]
    LinearLuSymbolicLuError {
        #[source]
        source: faer::sparse::FaerError,
    },
    #[error("LinearLu: Could not compute LU factorization.")]
    LinearLuLuError {
        #[source]
        source: faer::sparse::linalg::LuError,
    },
    #[error("LinearGmres: Number of threads must be greater than zero. Got num_thread = {num_thread}.")]
    LinearGmresNonPositiveThreads {
        num_thread: usize,
    },
    #[error("LinearGmres: Tolerance must be greater than zero. Got tol = {tol}.")]
    LinearGmresNonPositiveTolerance {
        tol: f64,
    },
    #[error("LinearGmres: Restart size must be greater than zero. Got restart = {restart}.")]
    LinearGmresNonPositiveRestart {
        restart: usize,
    },
    #[error("LinearGmres: Maximum outer iterations must be greater than zero. Got max_outer = {max_outer}.")]
    LinearGmresNonPositiveMaxOuter {
        max_outer: usize,
    },
    #[error("LinearGmres: Could not construct symbolic sparse matrix.")]
    LinearGmresSymbolicError {
        #[source]
        source: faer::sparse::CreationError,
    },
    #[error("LinearGmres: Could not assemble sparse matrix.")]
    LinearGmresMatrixError {
        #[source]
        source: faer::sparse::FaerError,
    },
    #[error("LinearGmres: Solver did not converge. {message}")]
    LinearGmresDidNotConverge {
        message: String,
    },

    // operators
    #[error("Operator: Use apply_time instead of apply.")]
    OperatorUseApplyTime {},
    #[error("Operator: Only unknown-type scalars can be solved for in the global linear system.")]
    OperatorNonUnknownInLinearSystem {},
    #[error("Operator: ScalarDomain ID {scldom_id} is not defined on Domain ID {dom_id}.")]
    OperatorScalarDomainNotOnDomain {
        scldom_id: usize,
        dom_id: usize,
    },
    #[error("Operator: ScalarDomain ID {scldom_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorScalarDomainNotOnBoundary {
        scldom_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: ScalarBoundary ID {sclbnd_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorScalarBoundaryNotOnBoundary {
        sclbnd_id: usize,
        bnd_id: usize,
    },
    #[error("Operator: ScalarDomain IDs {scldom_a_id} and {scldom_b_id} is not defined on Interface ID {itf_id}.")]
    OperatorScalarDomainNotOnInterface {
        scldom_a_id: usize,
        scldom_b_id: usize,
        itf_id: usize,
    },
    #[error("Operator: ScalarInterface ID {sclitf_id} is not defined on Interface ID {itf_id}.")]
    OperatorScalarInterfaceNotOnInterface {
        sclitf_id: usize,
        itf_id: usize,
    },
    #[error("Operator: ScalarFlow ID {sclflw_id} is not defined on Domain ID {dom_id}.")]
    OperatorScalarFlowNotOnDomain {
        sclflw_id: usize,
        dom_id: usize,
    },
    #[error("Operator: ScalarFlow ID {sclflw_id} is not defined on Boundary ID {bnd_id}.")]
    OperatorScalarFlowNotOnBoundary {
        sclflw_id: usize,
        bnd_id: usize,
    },
}
