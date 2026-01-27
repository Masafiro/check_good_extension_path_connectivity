use crate::graph::DirectedGraph;
use fs2::FileExt;
use serde::Serialize;
use std::error::Error;
use std::fmt;
use std::fs::OpenOptions;
use sysinfo::{Pid, System};

pub fn get_physical_mem_bytes() -> u64 {
    let mut sys = sysinfo::System::new();
    let pid = Pid::from_u32(std::process::id());
    sys.refresh_process(pid);
    sys.process(pid).map(|p| p.memory()).unwrap_or(0)
}

pub fn get_physical_memory_kb() -> u64 {
    get_physical_mem_bytes() / 1024
}

pub fn get_physical_memory_mb() -> u64 {
    get_physical_mem_bytes() / 1024 / 1024
}

pub fn get_virtual_memory_bytes() -> u64 {
    let mut sys = System::new();
    let pid = Pid::from_u32(std::process::id());
    sys.refresh_process(pid);
    sys.process(pid).map(|p| p.virtual_memory()).unwrap_or(0)
}

pub fn get_virtual_memory_kb() -> u64 {
    get_virtual_memory_bytes() / 1024
}

pub fn get_virtual_memory_mb() -> u64 {
    get_virtual_memory_kb() / 1024
}

/// Intersects two sorted vectors in-place.
///
/// The first vector `existing` is modified to contain only elements that are also present in `other`.
/// Both vectors must be sorted and contain no duplicates for this to work correctly in O(N+M).
pub fn intersect_sorted_vectors_inplace<T: Ord + Copy>(v1: &mut Vec<T>, v2: &[T]) {
    let mut i = 0; // index for existing
    let mut j = 0; // index for other
    let mut write_idx = 0;

    while i < v1.len() && j < v2.len() {
        if v1[i] == v2[j] {
            if write_idx != i {
                v1[write_idx] = v1[i];
            }
            write_idx += 1;
            i += 1;
            j += 1;
        } else if v1[i] < v2[j] {
            i += 1;
        } else {
            j += 1;
        }
    }
    v1.truncate(write_idx);
}

pub struct AnalysisResult {
    pub jacobian_count: usize,
    pub vertex_count: usize,
    pub edge_count: usize,
    pub is_strongly_connected: bool,
    pub scc_count: usize,
    pub are_all_jacobians_connected_via_main_scc: bool,
    pub are_all_jacobian_pairs_connected_via_good_paths: bool,
    pub are_all_minor_sccs_size_1: bool,
    pub are_all_minor_nodes_self_isogenies: bool,
    pub all_minor_sccs_have_preds_from_main: bool,
    pub any_minor_sccs_have_succs_to_main: bool,
    pub exists_minor_to_minor_edge: bool,
    pub condensed_graph: DirectedGraph,
}

impl fmt::Display for AnalysisResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let check = |b: bool| if b { "✅ Yes" } else { "❌ No " };

        writeln!(f, "\n--- Graph Connectivity Analysis Report ---")?;
        writeln!(f, "[Basic Statistics]")?;
        writeln!(f, "Number of Jacobians:        {}", self.jacobian_count)?;
        writeln!(f, "Number of Vertices:         {}", self.vertex_count)?;
        writeln!(f, "Number of Edges:            {}", self.edge_count)?;
        writeln!(
            f,
            "Strongly Connected:          {}",
            check(self.is_strongly_connected)
        )?;
        writeln!(f, "Number of SCCs:              {}", self.scc_count)?;

        writeln!(f, "\n[Main SCC Properties]")?;
        writeln!(
            f,
            "All Jacobians in Main SCC:   {}",
            check(self.are_all_jacobians_connected_via_main_scc)
        )?;

        writeln!(f, "\n[Jacobian Pair Connectivity]")?;
        writeln!(
            f,
            "All Jacobian Pairs Connected via Good Paths:    {}",
            check(self.are_all_jacobian_pairs_connected_via_good_paths)
        )?;

        if !self.is_strongly_connected {
            writeln!(f, "\n[Minor SCC Properties]")?;

            writeln!(
                f,
                "All Minor SCCs Size 1:       {}",
                check(self.are_all_minor_sccs_size_1)
            )?;
            writeln!(
                f,
                "All Minor Nodes Self-Isog:   {}",
                check(self.are_all_minor_nodes_self_isogenies)
            )?;
            writeln!(f, "\n[Minor SCC Detailed Connectivity]")?;
            writeln!(
                f,
                "All Minor SCCs Have Predecessors from Main SCC:  {}",
                check(self.all_minor_sccs_have_preds_from_main)
            )?;
            writeln!(
                f,
                "Any Minor SCCs Have Successors to Main SCC:      {}",
                check(self.any_minor_sccs_have_succs_to_main)
            )?;
            writeln!(
                f,
                "Exists Minor-to-Minor SCC Edge:                  {}",
                check(self.exists_minor_to_minor_edge)
            )?;
            writeln!(f, "\n[Condensed Graph]")?;
            writeln!(f, "Number of Vertices:         {}", self.condensed_graph.n)?;
            writeln!(f, "adjacency list:")?;
            let adj_str = self
                .condensed_graph
                .adj
                .iter()
                .enumerate()
                .map(|(i, adj)| format!("{}:{:?}", i, adj))
                .collect::<Vec<_>>()
                .join(", ");

            writeln!(f, "Adjacency: {}", adj_str)?;
        }

        writeln!(f, "------------------------------------------")
    }
}

#[derive(Serialize)]
pub struct RowDataSsp {
    pub p: u32,
    pub d: u32,
    pub construction_time_sec: String,
    pub jacobian_count: usize,
    pub vertex_count: usize,
    pub edge_count: usize,
    pub connectivity_check_time_sec: String,
    pub is_strongly_connected: bool,
    pub scc_count: usize,
    pub are_all_jacobians_connected_via_main_scc: bool,
    pub are_all_jacobian_pairs_connected_via_good_paths: bool,
    pub are_all_minor_sccs_size_1: bool,
    pub are_all_minor_nodes_self_isogenies: bool,
    pub all_minor_sccs_have_preds_from_main: bool,
    pub any_minor_sccs_have_succs_to_main: bool,
    pub exists_minor_to_minor_edge: bool,
    pub condensed_graph: String,
    pub construction_peak_bytes: usize,
    pub connectivity_check_peak_bytes: usize,
}

#[derive(Serialize)]
pub struct RowDataSsg {
    pub p: u32,
    pub d1: u32,
    pub d2c0: u32,
    pub d2c1: u32,
    pub construction_time_sec: String,
    pub jacobian_count: usize,
    pub vertex_count: usize,
    pub edge_count: usize,
    pub connectivity_check_time_sec: String,
    pub is_strongly_connected: bool,
    pub scc_count: usize,
    pub are_all_jacobians_connected_via_main_scc: bool,
    pub are_all_jacobian_pairs_connected_via_good_paths: bool,
    pub are_all_minor_sccs_size_1: bool,
    pub are_all_minor_nodes_self_isogenies: bool,
    pub all_minor_sccs_have_preds_from_main: bool,
    pub any_minor_sccs_have_succs_to_main: bool,
    pub exists_minor_to_minor_edge: bool,
    pub condensed_graph: String,
    pub construction_peak_bytes: usize,
    pub connectivity_check_peak_bytes: usize,
}

pub fn save_to_csv<T: serde::Serialize>(path: &str, data: T) -> Result<(), Box<dyn Error>> {
    // If the file exists, append; otherwise, create a new file
    let file = OpenOptions::new()
        .write(true)
        .create(true)
        .append(true)
        .open(path)?;

    // Lock the file exclusively to prevent concurrent writes
    file.lock_exclusive()?;

    let mut wtr = csv::WriterBuilder::new()
        .has_headers(file.metadata()?.len() == 0) // Write headers only if the file is empty
        .from_writer(&file);

    wtr.serialize(data)?;
    wtr.flush()?;

    // Unlock the file after writing
    file.unlock()?;
    Ok(())
}
