//! Exhaustive good-extension reachability test for superspecial (2,2)-isogenies.
//!
//! This auxiliary verifier is intended for cases in which the coarse
//! SCC-based analysis in `ssp_ge_graph` is inconclusive for Conjecture 8.
//! In particular, it is currently used for the `ssp_all` cases p = 17 and
//! p = 23.
//!
//! Unlike `ssp_ge_graph`, this program retains the normalized theta null-point
//! induced by the preceding isogeny and performs reachability searches on the
//! resulting marked-state graph. Therefore, compatibility of consecutive
//! good extensions is preserved throughout the search without explicitly
//! comparing isogeny kernels.
//!
//! This exhaustive method is substantially more expensive than
//! `ssp_ge_graph` and is intended mainly for exceptional small-prime cases.
//!
//! Example:
//!   cargo run --release --bin ssp_edge_reachability -- \
//!       17 data/ssp_all_edge_reachability.csv false
//!   cargo run --release --bin ssp_edge_reachability -- \
//!       23 data/ssp_all_edge_reachability.csv false

use deepsize::DeepSizeOf;
use g2_isogeny::counting::{number_of_ssp_jacobian_nodes, number_of_ssp_nodes};
use g2_isogeny::curves::{Curve, Invariants};
use g2_isogeny::fq::{Fp, Fp2};
use g2_isogeny::graph::DirectedGraph;
use g2_isogeny::theta::Theta;
use g2_isogeny::utils::save_to_csv;
use serde::Serialize;
use std::collections::{HashMap, HashSet, VecDeque};
use std::env;
use std::time::{Duration, Instant};

/// A marked directed isogeny A -> B.
///
/// `null` is the projectively normalized theta null-point on B carrying the
/// marking induced by A -> B.  This marking determines which outgoing
/// isogenies from B are good extensions of the preceding isogeny.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash, DeepSizeOf)]
struct MarkedEdge<const P: u32, const D: u32> {
    source_id: usize,
    target_id: usize,
    null: [Fp2<P, D>; 10],
}

impl<const P: u32, const D: u32> MarkedEdge<P, D> {
    fn theta(&self) -> Theta<Fp2<P, D>> {
        Theta { null: self.null }
    }

    fn coarse_pair(&self) -> (usize, usize) {
        (self.source_id, self.target_id)
    }
}

#[derive(Debug)]
pub struct EdgeReachabilityAnalysis {
    pub variety_count: usize,
    pub coarse_pair_count: usize,
    pub marked_edge_state_count: usize,
    pub marked_transition_count: usize,
    pub marked_scc_count: usize,
    pub all_marked_edges_mutually_reachable: bool,
    pub all_coarse_pairs_connected_via_good_paths: bool,
    pub first_unreachable_source_pair: Option<(usize, usize)>,
    pub first_unreachable_target_pair: Option<(usize, usize)>,
}

#[derive(Serialize)]
struct RowData {
    p: u32,
    d: u32,
    construction_time_sec: String,
    reachability_check_time_sec: String,
    variety_count: usize,
    coarse_pair_count: usize,
    marked_edge_state_count: usize,
    marked_transition_count: usize,
    marked_scc_count: usize,
    all_marked_edges_mutually_reachable: bool,
    all_coarse_pairs_connected_via_good_paths: bool,
    first_unreachable_source_pair: String,
    first_unreachable_target_pair: String,
    construction_peak_bytes: usize,
    reachability_check_peak_bytes: usize,
}

/// Remove the common projective scalar from a theta null-point.
///
/// Without this normalization, the same marked theta structure can reappear
/// with different common scalar factors, so a finite-state search may keep
/// producing duplicate states.
fn normalize_theta<const P: u32, const D: u32>(
    theta: Theta<Fp2<P, D>>,
) -> Theta<Fp2<P, D>> {
    let pivot = theta
        .null
        .iter()
        .copied()
        .find(|x| !x.is_zero())
        .expect("theta null-point is identically zero");
    let pivot_inv = pivot.inv();
    let mut null = theta.null;
    for coordinate in &mut null {
        *coordinate *= pivot_inv;
    }
    Theta { null }
}

fn keep_surface<const P: u32, const D: u32>(
    theta: &Theta<Fp2<P, D>>,
    only_jacobians: bool,
) -> bool {
    !only_jacobians || theta.find_split_index().is_none()
}

/// Enumerate the unmarked superspecial principally polarized surfaces.
fn enumerate_surfaces<const P: u32, const D: u32>(
    only_jacobians: bool,
) -> (
    HashMap<Invariants<Fp2<P, D>>, usize>,
    Vec<Curve<Fp2<P, D>>>,
) {
    let expected_count = if only_jacobians {
        number_of_ssp_jacobian_nodes(P) as usize
    } else {
        number_of_ssp_nodes(P) as usize
    };

    let initial_theta = Theta::<Fp2<P, D>>::generate_ssp_theta();
    let initial_curve = initial_theta.to_curve();
    let initial_invariants = initial_curve.invariants();

    let mut invariant_to_id =
        HashMap::<Invariants<Fp2<P, D>>, usize>::with_capacity(expected_count);
    let mut curves = Vec::<Curve<Fp2<P, D>>>::with_capacity(expected_count);
    let mut queue = VecDeque::<Theta<Fp2<P, D>>>::new();

    invariant_to_id.insert(initial_invariants, 0);
    curves.push(initial_curve);
    queue.push_back(initial_theta);

    while let Some(theta) = queue.pop_front() {
        for image in theta.compute_all_twoisogenies(only_jacobians, false) {
            let curve = image.to_curve();
            let invariants = curve.invariants();
            if let std::collections::hash_map::Entry::Vacant(entry) =
                invariant_to_id.entry(invariants)
            {
                let id = curves.len();
                entry.insert(id);
                curves.push(curve);
                queue.push_back(image);
            }
        }
    }

    assert_eq!(
        invariant_to_id.len(),
        expected_count,
        "Number of surfaces mismatch: computed {}, expected {}",
        invariant_to_id.len(),
        expected_count
    );

    (invariant_to_id, curves)
}

fn surface_id<const P: u32, const D: u32>(
    theta: &Theta<Fp2<P, D>>,
    invariant_to_id: &HashMap<Invariants<Fp2<P, D>>, usize>,
) -> usize {
    let invariants = theta.to_curve().invariants();
    *invariant_to_id
        .get(&invariants)
        .expect("isogeny image was not found in the enumerated surface list")
}

fn intern_state<const P: u32, const D: u32>(
    state: MarkedEdge<P, D>,
    state_to_id: &mut HashMap<MarkedEdge<P, D>, usize>,
    states: &mut Vec<MarkedEdge<P, D>>,
    adjacency: &mut Vec<Vec<usize>>,
    queue: &mut VecDeque<usize>,
) -> usize {
    if let Some(&id) = state_to_id.get(&state) {
        return id;
    }

    let id = states.len();
    state_to_id.insert(state, id);
    states.push(state);
    adjacency.push(Vec::new());
    queue.push_back(id);
    id
}

/// Build the refined graph whose vertices retain the preceding isogeny marking.
fn build_marked_edge_graph<const P: u32, const D: u32>(
    only_jacobians: bool,
) -> (
    Vec<Curve<Fp2<P, D>>>,
    Vec<MarkedEdge<P, D>>,
    DirectedGraph,
    Duration,
    usize,
) {
    assert!(P > 5, "P must be greater than 5");
    let started = Instant::now();
    let (invariant_to_id, curves) = enumerate_surfaces::<P, D>(only_jacobians);

    let mut state_to_id = HashMap::<MarkedEdge<P, D>, usize>::new();
    let mut states = Vec::<MarkedEdge<P, D>>::new();
    let mut adjacency = Vec::<Vec<usize>>::new();
    let mut queue = VecDeque::<usize>::new();

    // Seed all 15 outgoing (2,2)-isogenies from a canonical theta structure
    // for every surface.  These are all possible starting edges.
    for (source_id, curve) in curves.iter().enumerate() {
        let source_theta = normalize_theta(Theta::<Fp2<P, D>>::from_curve(curve));
        for kernel_index in 0..15 {
            let image = source_theta.compute_twoisogeny(kernel_index);
            if !keep_surface(&image, only_jacobians) {
                continue;
            }
            let target_id = surface_id(&image, &invariant_to_id);
            let marked_image = normalize_theta(image);
            let state = MarkedEdge {
                source_id,
                target_id,
                null: marked_image.null,
            };
            intern_state(
                state,
                &mut state_to_id,
                &mut states,
                &mut adjacency,
                &mut queue,
            );
        }
    }

    // Exhaust the closure under good extensions.  In the current theta
    // implementation, indices 0..8 are precisely the good extensions.
    while let Some(state_id) = queue.pop_front() {
        let state = states[state_id];
        let theta = state.theta();
        let mut successors = Vec::with_capacity(8);

        for good_index in 0..8 {
            let image = theta.compute_twoisogeny(good_index);
            if !keep_surface(&image, only_jacobians) {
                continue;
            }
            let next_target_id = surface_id(&image, &invariant_to_id);
            let marked_image = normalize_theta(image);
            let next_state = MarkedEdge {
                source_id: state.target_id,
                target_id: next_target_id,
                null: marked_image.null,
            };
            let next_state_id = intern_state(
                next_state,
                &mut state_to_id,
                &mut states,
                &mut adjacency,
                &mut queue,
            );
            successors.push(next_state_id);
        }

        successors.sort_unstable();
        successors.dedup();
        adjacency[state_id] = successors;
    }

    let graph = DirectedGraph::from_adj(adjacency);
    let peak_memory_bytes = invariant_to_id.deep_size_of()
        + curves.deep_size_of()
        + state_to_id.deep_size_of()
        + states.deep_size_of()
        + graph.deep_size_of()
        + queue.deep_size_of();

    (curves, states, graph, started.elapsed(), peak_memory_bytes)
}

/// Multi-source BFS from every marked state representing one coarse pair.
///
/// A coarse pair (A,B) reaches (C,D) if at least one marked isogeny A -> B
/// has a good-extension path ending at at least one marked isogeny C -> D.
/// This is the surface-pair quantification used in Conjecture 8.
fn check_coarse_pair_reachability<const P: u32, const D: u32>(
    states: &[MarkedEdge<P, D>],
    graph: &DirectedGraph,
) -> (
    bool,
    Option<((usize, usize), (usize, usize))>,
    Vec<(usize, usize)>,
    usize,
) {
    let mut fibers = HashMap::<(usize, usize), Vec<usize>>::new();
    for (state_id, state) in states.iter().enumerate() {
        fibers.entry(state.coarse_pair()).or_default().push(state_id);
    }

    let mut coarse_pairs: Vec<(usize, usize)> = fibers.keys().copied().collect();
    coarse_pairs.sort_unstable();

    let mut visited = vec![false; states.len()];
    let mut bfs_queue = VecDeque::<usize>::new();
    let mut reached_pairs = HashSet::<(usize, usize)>::with_capacity(coarse_pairs.len());

    for &source_pair in &coarse_pairs {
        visited.fill(false);
        bfs_queue.clear();
        reached_pairs.clear();

        for &state_id in &fibers[&source_pair] {
            if !visited[state_id] {
                visited[state_id] = true;
                bfs_queue.push_back(state_id);
            }
        }

        while let Some(u) = bfs_queue.pop_front() {
            reached_pairs.insert(states[u].coarse_pair());
            for &v in &graph.adj[u] {
                if !visited[v] {
                    visited[v] = true;
                    bfs_queue.push_back(v);
                }
            }
        }

        if let Some(&target_pair) = coarse_pairs
            .iter()
            .find(|&&pair| !reached_pairs.contains(&pair))
        {
            let memory_bytes = fibers.deep_size_of()
                + coarse_pairs.deep_size_of()
                + visited.deep_size_of()
                + bfs_queue.deep_size_of()
                + reached_pairs.deep_size_of();
            return (
                false,
                Some((source_pair, target_pair)),
                coarse_pairs,
                memory_bytes,
            );
        }
    }

    let memory_bytes = fibers.deep_size_of()
        + coarse_pairs.deep_size_of()
        + visited.deep_size_of()
        + bfs_queue.deep_size_of()
        + reached_pairs.deep_size_of();
    (true, None, coarse_pairs, memory_bytes)
}

pub fn analyze_edge_reachability<const P: u32, const D: u32>(
    only_jacobians: bool,
) -> (
    EdgeReachabilityAnalysis,
    Duration,
    Duration,
    usize,
    usize,
) {
    let (curves, states, graph, construction_time, construction_peak_bytes) =
        build_marked_edge_graph::<P, D>(only_jacobians);

    let check_started = Instant::now();
    let (sccs, scc_memory_bytes) = graph.scc();
    let all_marked_edges_mutually_reachable = sccs.len() == 1;
    let (all_coarse_pairs_connected, first_missing, coarse_pairs, bfs_memory_bytes) =
        check_coarse_pair_reachability(&states, &graph);
    let check_time = check_started.elapsed();

    let first_unreachable_source_pair = first_missing.map(|(source, _)| source);
    let first_unreachable_target_pair = first_missing.map(|(_, target)| target);
    let marked_transition_count = graph.adj.iter().map(Vec::len).sum();
    let check_peak_bytes = scc_memory_bytes
        + bfs_memory_bytes
        + states.deep_size_of()
        + graph.deep_size_of();

    (
        EdgeReachabilityAnalysis {
            variety_count: curves.len(),
            coarse_pair_count: coarse_pairs.len(),
            marked_edge_state_count: states.len(),
            marked_transition_count,
            marked_scc_count: sccs.len(),
            all_marked_edges_mutually_reachable,
            all_coarse_pairs_connected_via_good_paths: all_coarse_pairs_connected,
            first_unreachable_source_pair,
            first_unreachable_target_pair,
        },
        construction_time,
        check_time,
        construction_peak_bytes,
        check_peak_bytes,
    )
}

fn pair_to_string(pair: Option<(usize, usize)>) -> String {
    match pair {
        Some((u, v)) => format!("({u},{v})"),
        None => String::new(),
    }
}

macro_rules! analyze_connectivity {
    ($p:literal, $filename:expr, $only_jacobians:expr) => {{
        const D: u32 = Fp::<$p>::find_non_residue();
        let (analysis, construction_time, check_time, construction_memory, check_memory) =
            analyze_edge_reachability::<$p, D>($only_jacobians);

        println!("\n--- Marked-edge good-path reachability report ---");
        println!("p = {}, d = {}", $p, D);
        println!("Number of surfaces:              {}", analysis.variety_count);
        println!("Number of coarse surface pairs:  {}", analysis.coarse_pair_count);
        println!("Number of marked edge states:    {}", analysis.marked_edge_state_count);
        println!("Number of marked transitions:    {}", analysis.marked_transition_count);
        println!("Number of marked SCCs:           {}", analysis.marked_scc_count);
        println!(
            "All marked edges mutually reachable:      {}",
            analysis.all_marked_edges_mutually_reachable
        );
        println!(
            "All coarse pairs connected by good paths: {}",
            analysis.all_coarse_pairs_connected_via_good_paths
        );
        if let (Some(source), Some(target)) = (
            analysis.first_unreachable_source_pair,
            analysis.first_unreachable_target_pair,
        ) {
            println!("First unreachable pair: {:?} -> {:?}", source, target);
        }

        let row = RowData {
            p: $p,
            d: D,
            construction_time_sec: format!("{:.4}", construction_time.as_secs_f64()),
            reachability_check_time_sec: format!("{:.4}", check_time.as_secs_f64()),
            variety_count: analysis.variety_count,
            coarse_pair_count: analysis.coarse_pair_count,
            marked_edge_state_count: analysis.marked_edge_state_count,
            marked_transition_count: analysis.marked_transition_count,
            marked_scc_count: analysis.marked_scc_count,
            all_marked_edges_mutually_reachable: analysis.all_marked_edges_mutually_reachable,
            all_coarse_pairs_connected_via_good_paths: analysis
                .all_coarse_pairs_connected_via_good_paths,
            first_unreachable_source_pair: pair_to_string(
                analysis.first_unreachable_source_pair,
            ),
            first_unreachable_target_pair: pair_to_string(
                analysis.first_unreachable_target_pair,
            ),
            construction_peak_bytes: construction_memory,
            reachability_check_peak_bytes: check_memory,
        };

        save_to_csv($filename, row).expect("failed to save CSV result");
    }};
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 4 {
        eprintln!(
            "Usage: {} <prime> <output.csv> <only_jacobians:true|false>",
            args[0]
        );
        std::process::exit(2);
    }

    let p_arg: u32 = args[1].parse().expect("prime must be a u32");
    let filename_arg = &args[2];
    let only_jacobians: bool = args[3]
        .parse()
        .expect("only_jacobians must be true or false");

    include!(concat!(env!("OUT_DIR"), "/generated_dispatch.rs"));
}
