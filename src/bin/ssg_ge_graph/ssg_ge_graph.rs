//! Compute the supersingular non-superspecial (2,2)-good extension graph and analyze its connectivity

use deepsize::DeepSizeOf;
use g2_isogeny::counting::number_of_ssg_nodes;
use g2_isogeny::curves::{IgusaInvariants, Rosenhain};
use g2_isogeny::fq::Fp4;
use g2_isogeny::graph::DirectedGraph;
use g2_isogeny::theta::Theta;
use g2_isogeny::utils::{AnalysisResult, intersect_sorted_vectors_inplace};
use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet, VecDeque};
use std::time;

pub fn compute_ssg_good_extension_graph<
    const P: u32,
    const D1: u32,
    const D2C0: u32,
    const D2C1: u32,
>() -> (
    HashMap<IgusaInvariants<Fp4<P, D1, D2C0, D2C1>>, usize>,
    Vec<Rosenhain<Fp4<P, D1, D2C0, D2C1>>>,
    Vec<(usize, usize)>,
    DirectedGraph,
    time::Duration,
    usize,
) {
    assert!(P > 5, "P must be greater than 5");

    let num_j = number_of_ssg_nodes(P) as usize;

    let initial_ram = Rosenhain::<Fp4<P, D1, D2C0, D2C1>>::generate_ssg_rosenhain();
    let initial_theta = Theta::<Fp4<P, D1, D2C0, D2C1>>::from_rosenhain(&initial_ram);
    let initial_inv = initial_ram.igusa_invariants();

    let time_start = time::Instant::now();

    let mut inv_to_id =
        HashMap::<IgusaInvariants<Fp4<P, D1, D2C0, D2C1>>, usize>::with_capacity(num_j);
    let mut ram_list = Vec::<Rosenhain<Fp4<P, D1, D2C0, D2C1>>>::with_capacity(num_j);
    let mut edge_info = HashMap::<(usize, usize), Vec<usize>>::with_capacity(15 * num_j);
    let mut dq = VecDeque::<(Theta<Fp4<P, D1, D2C0, D2C1>>, usize)>::new();

    inv_to_id.insert(initial_inv, 0);
    ram_list.push(initial_ram);
    dq.push_back((initial_theta, 0));

    while let Some((prev, prev_id)) = dq.pop_front() {
        // Not necessarily good extension
        for curr in prev.compute_all_twoisogenies(false) {
            let curr_ram = curr.to_rosenhain();
            let curr_inv = curr_ram.igusa_invariants();
            let curr_id = *inv_to_id.entry(curr_inv).or_insert_with(|| {
                let id = ram_list.len();
                ram_list.push(curr_ram);
                dq.push_back((curr, id));
                id
            });

            let key = (prev_id, curr_id);
            let mut next_set = Vec::<usize>::with_capacity(8);

            // Good extension
            for next in curr.compute_all_twoisogenies(true) {
                let next_ram = next.to_rosenhain();
                let next_inv = next_ram.igusa_invariants();
                let next_id = *inv_to_id.entry(next_inv).or_insert_with(|| {
                    let id = ram_list.len();
                    ram_list.push(next_ram);
                    dq.push_back((next, id));
                    id
                });
                next_set.push(next_id);
            }

            next_set.sort_unstable();
            next_set.dedup();

            match edge_info.entry(key) {
                Entry::Vacant(entry) => {
                    entry.insert(next_set);
                }
                Entry::Occupied(mut entry) => {
                    // Intersect two sorted lists in-place
                    let existing_set = entry.get_mut();
                    intersect_sorted_vectors_inplace(existing_set, &next_set);
                }
            }
        }
    }

    let mut node_list: Vec<(usize, usize)> = edge_info.keys().copied().collect();
    node_list.sort_unstable(); // sort for binary search

    let peak_memory_bytes = dq.deep_size_of()
        + inv_to_id.deep_size_of()
        + ram_list.deep_size_of()
        + node_list.deep_size_of()
        + edge_info.deep_size_of();

    // Build adjacency list
    let mut adjacency_list = vec![vec![]; node_list.len()];
    for ((u, v), next_set) in edge_info {
        let from_idx = node_list.binary_search(&(u, v)).unwrap();
        for w in next_set {
            let to_idx = node_list.binary_search(&(v, w)).unwrap();
            adjacency_list[from_idx].push(to_idx);
        }
    }

    let ge_graph = DirectedGraph::from_adj(adjacency_list);

    let duration = time_start.elapsed();

    assert!(
        inv_to_id.len() == num_j,
        "Number of nodes mismatch: computed {}, expected {}",
        inv_to_id.len(),
        num_j
    );

    (
        inv_to_id,
        ram_list,
        node_list,
        ge_graph,
        duration,
        peak_memory_bytes,
    )
}

pub fn check_ge_path_connectivity<const P: u32, const D1: u32, const D2C0: u32, const D2C1: u32>(
    inv_to_id: &HashMap<IgusaInvariants<Fp4<P, D1, D2C0, D2C1>>, usize>,
    ram_list: &Vec<Rosenhain<Fp4<P, D1, D2C0, D2C1>>>,
    node_list: &Vec<(usize, usize)>,
    ge_graph: &DirectedGraph,
) -> (AnalysisResult, time::Duration, usize) {
    let time_start = time::Instant::now();

    let (sccs, condensed_graph, scc_memory_bytes) = ge_graph.scc_and_condense();

    let peak_memory_bytes = scc_memory_bytes
        + inv_to_id.deep_size_of()
        + ram_list.deep_size_of()
        + node_list.deep_size_of()
        + ge_graph.deep_size_of();

    let scc_count = sccs.len();
    let is_strongly_connected = scc_count == 1;

    let main_scc: &Vec<usize> = &sccs[0];
    let mut jacobians_in_main_scc = HashSet::with_capacity(inv_to_id.len());
    for &node_idx in main_scc {
        let (j0, j1) = node_list[node_idx];
        jacobians_in_main_scc.insert(j0);
        jacobians_in_main_scc.insert(j1);
    }

    let are_all_jacobians_connected_via_main_scc = jacobians_in_main_scc.len() == inv_to_id.len();

    // Minor SCCs
    let minor_sccs = &sccs[1..];

    let mut is_in_main = vec![false; ge_graph.n];
    for &idx in main_scc {
        is_in_main[idx] = true;
    }

    let has_good_path_from_main_scc = |i: usize, comp: &Vec<usize>| -> bool {
        // (A) グラフ内に既に辺がある場合
        if condensed_graph.adj[0].contains(&i) {
            return true;
        }

        for &node_idx in comp {
            // (B) グラフ外の探索: prev -> j0 -> j1 (node(prev, j0) が main_scc にあるか)
            let (j0_id, j1_id) = node_list[node_idx];
            let j0 = Theta::from_rosenhain(&ram_list[j0_id]);

            // j0 から戻る 2-isogeny (not necessarily good)
            for prev in j0.compute_all_twoisogenies(false) {
                let prev_inv = prev.to_rosenhain().igusa_invariants();
                let prev_id = inv_to_id[&prev_inv];
                // node(prev, j0) がメインSCCに存在するかチェック
                if !is_in_main[node_list.binary_search(&(prev_id, j0_id)).unwrap()] {
                    continue;
                }

                for curr in prev.compute_all_twoisogenies(false) {
                    let curr_inv = curr.to_rosenhain().igusa_invariants();
                    let curr_id = inv_to_id[&curr_inv];
                    if curr_id != j0_id {
                        continue;
                    }

                    for next in curr.compute_all_twoisogenies(true) {
                        let next_inv = next.to_rosenhain().igusa_invariants();
                        let next_id = inv_to_id[&next_inv];
                        if next_id == j1_id {
                            return true;
                        }
                    }
                }
            }
        }
        false
    };

    let has_good_path_to_main_scc = |i: usize, comp: &Vec<usize>| -> bool {
        // (A) グラフ内に既に辺がある場合
        if condensed_graph.adj[i].contains(&0) {
            return true;
        }

        for &node_idx in comp {
            // (B) グラフ外の探索: j0 -> j1 -> next
            let (j0_id, j1_id) = node_list[node_idx];
            let j0 = Theta::from_rosenhain(&ram_list[j0_id]);

            for curr in j0.compute_all_twoisogenies(false) {
                let curr_inv = curr.to_rosenhain().igusa_invariants();
                let curr_id = inv_to_id[&curr_inv];
                if curr_id != j1_id {
                    continue;
                }

                for next in curr.compute_all_twoisogenies(true) {
                    let next_inv = next.to_rosenhain().igusa_invariants();
                    let next_id = inv_to_id[&next_inv];
                    if is_in_main[node_list.binary_search(&(j1_id, next_id)).unwrap()] {
                        return true;
                    }
                }
            }
        }
        false
    };

    let are_all_jacobian_pairs_connected_via_good_paths =
        minor_sccs.iter().enumerate().all(|(i, comp)| {
            has_good_path_from_main_scc(i, comp) && has_good_path_to_main_scc(i, comp)
        });

    let is_size_1 = |comp: &Vec<usize>| comp.len() == 1;
    let are_all_minor_sccs_size_1 = minor_sccs.iter().all(|comp| is_size_1(comp));
    if !are_all_minor_sccs_size_1 {
        let anomalous_sizes: Vec<usize> = minor_sccs
            .iter()
            .filter(|comp| !is_size_1(comp))
            .map(|comp| comp.len())
            .collect();
        println!("Anomaly minor SCC sizes: {:?}", anomalous_sizes);
    }

    let is_self_isogeny = |node_idx: usize| {
        let (j0, j1) = node_list[node_idx];
        j0 == j1
    };
    let are_all_minor_nodes_self_isogenies = minor_sccs
        .iter()
        .flatten()
        .all(|&node_idx| is_self_isogeny(node_idx));
    if !are_all_minor_nodes_self_isogenies {
        let anomalies: Vec<(usize, (usize, usize))> = minor_sccs
            .iter()
            .flatten()
            .copied()
            .filter(|&node_idx| !is_self_isogeny(node_idx))
            .map(|node_idx| (node_idx, node_list[node_idx]))
            .collect();
        println!("Anomaly minor nodes (node_idx, (j0, j1)): {:?}", anomalies);
    }

    let all_minor_sccs_have_preds_from_main = condensed_graph
        .adj
        .iter()
        .enumerate()
        .skip(1)
        .all(|(i, _)| condensed_graph.adj[0].contains(&i));

    let any_minor_sccs_have_succs_to_main = condensed_graph
        .adj
        .iter()
        .skip(1)
        .any(|minor_adj| minor_adj.contains(&0));

    let exists_minor_to_minor_edge = condensed_graph
        .adj
        .iter()
        .skip(1)
        .any(|minor_adj| minor_adj.iter().any(|&p| p != 0));

    let duration = time_start.elapsed();

    (
        AnalysisResult {
            jacobian_count: inv_to_id.len(),
            vertex_count: node_list.len(),
            edge_count: ge_graph.adj.iter().map(|vec| vec.len()).sum(),
            is_strongly_connected,
            scc_count,
            are_all_jacobians_connected_via_main_scc,
            are_all_jacobian_pairs_connected_via_good_paths,
            are_all_minor_sccs_size_1,
            are_all_minor_nodes_self_isogenies,
            all_minor_sccs_have_preds_from_main,
            any_minor_sccs_have_succs_to_main,
            exists_minor_to_minor_edge,
            condensed_graph,
        },
        duration,
        peak_memory_bytes,
    )
}
