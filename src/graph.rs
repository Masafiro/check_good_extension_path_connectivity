//! Directed graph structure and algorithms (SCC decomposition)
use deepsize::DeepSizeOf;
use std::cmp::min;
use std::collections::HashSet;

/// Directed graph structure
#[derive(Debug, DeepSizeOf)]
pub struct DirectedGraph {
    pub n: usize,
    pub adj: Vec<Vec<usize>>, // Successors: u -> v
}

impl DirectedGraph {
    /// Construct graph from adjacency list
    pub fn from_adj(adj: Vec<Vec<usize>>) -> Self {
        let n = adj.len();
        Self { n, adj }
    }

    /// Strong connected component decomposition (SCC) using Tarjan's algorithm
    /// Returns: List of vertex sets for each SCC (topological order)
    pub fn scc(&self) -> (Vec<Vec<usize>>, usize) {
        let n = self.n;
        let mut ord = vec![None; n];
        let mut low = vec![0; n];
        let mut now_ord = 0;
        let mut scc_stack = Vec::new();
        let mut on_scc_stack = vec![false; n];
        let mut components = Vec::new();
        let mut dfs_stack = Vec::new(); // To avoid stack overflow, replace recursion with an stack

        for start_node in 0..n {
            if ord[start_node].is_some() {
                continue;
            }
            dfs_stack.push((start_node, 0));

            while let Some((u, child_idx)) = dfs_stack.pop() {
                // Pre-order processing
                if child_idx == 0 {
                    // First child
                    ord[u] = Some(now_ord);
                    low[u] = now_ord;
                    now_ord += 1;
                    scc_stack.push(u);
                    on_scc_stack[u] = true;
                }

                let mut found_new_child = false;
                let children = &self.adj[u];

                for i in child_idx..children.len() {
                    let v = children[i];
                    if ord[v].is_none() {
                        dfs_stack.push((u, i + 1));
                        dfs_stack.push((v, 0));
                        found_new_child = true;
                        break;
                    } else if on_scc_stack[v] {
                        low[u] = min(low[u], ord[v].unwrap());
                    }
                }

                if found_new_child {
                    continue;
                }

                // Post-order processing
                if ord[u] == Some(low[u]) {
                    let mut component = Vec::new();
                    loop {
                        let node = scc_stack.pop().unwrap();
                        on_scc_stack[node] = false;
                        component.push(node);
                        if node == u {
                            break;
                        }
                    }
                    components.push(component);
                }

                if let Some((parent, _)) = dfs_stack.last() {
                    low[*parent] = min(low[*parent], low[u]);
                }
            }
        }

        components.reverse(); // Reverse to get correct topological order

        let peak_memory_bytes = ord.deep_size_of()
            + low.deep_size_of()
            + scc_stack.deep_size_of()
            + on_scc_stack.deep_size_of()
            + components.deep_size_of()
            + dfs_stack.deep_size_of();

        (components, peak_memory_bytes)
    }

    /// Compute SCC and condensed graph
    pub fn scc_and_condense(&self) -> (Vec<Vec<usize>>, Self, usize) {
        let (sccs, peak_memory_bytes) = self.scc();
        let mut comp_id = vec![0; self.n];
        for (i, comp) in sccs.iter().enumerate() {
            for &v in comp {
                comp_id[v] = i;
            }
        }

        let m = sccs.len();

        let mut adj_condensed_sets = vec![HashSet::new(); m];

        for u in 0..self.n {
            let cu = comp_id[u];
            for &v in &self.adj[u] {
                let cv = comp_id[v];
                adj_condensed_sets[cu].insert(cv);
            }
        }

        let mut adj_condensed = vec![vec![]; m];
        for (i, set) in adj_condensed_sets.into_iter().enumerate() {
            let mut neighbors: Vec<usize> = set.into_iter().collect();
            neighbors.sort_unstable();
            adj_condensed[i] = neighbors;
        }

        let condensed_graph = DirectedGraph {
            n: m,
            adj: adj_condensed,
        };

        (sccs, condensed_graph, peak_memory_bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper function for tests
    /// Normalize result order (within components and between components) for easier comparison
    fn normalize(mut sccs: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
        // 1. Sort contents of each component (e.g., [1, 0] -> [0, 1])
        for comp in &mut sccs {
            comp.sort();
        }
        // 2. Sort components by their first element (e.g., [[2], [0, 1]] -> [[0, 1], [2]])
        sccs.sort_by(|a, b| a[0].cmp(&b[0]));
        sccs
    }

    #[test]
    fn test_simple_cycle() {
        // Graph: 0 -> 1 -> 2 -> 0 (one strongly connected component)
        let adj = vec![
            vec![1], // node 0
            vec![2], // node 1
            vec![0], // node 2
        ];

        // Note: method name is from_edges but implementation takes adjacency list
        let graph = DirectedGraph::from_adj(adj);
        let (sccs, _) = graph.scc();

        let result = normalize(sccs);
        assert_eq!(result, vec![vec![0, 1, 2]]);
    }

    #[test]
    fn test_disconnected_cycles() {
        // Graph: 0<->1 and 2<->3 exist disconnected
        let adj = vec![
            vec![1], // 0 -> 1
            vec![0], // 1 -> 0
            vec![3], // 2 -> 3
            vec![2], // 3 -> 2
        ];

        let graph = DirectedGraph::from_adj(adj);
        let (sccs, _) = graph.scc();

        let result = normalize(sccs);
        // Expected: two components {0, 1} and {2, 3}
        assert_eq!(result, vec![vec![0, 1], vec![2, 3]]);
    }

    #[test]
    fn test_standalone_function() {
        // Test standalone function
        // Graph: 0 -> 1 -> 2, and 2 -> 1 (1 and 2 form loop, 0 leads into it)
        // Components: {0}, {1, 2}
        let adj = vec![vec![1], vec![2], vec![1]];

        let graph = DirectedGraph::from_adj(adj);
        let (sccs, _) = graph.scc();

        let result = normalize(sccs);
        assert_eq!(result, vec![vec![0], vec![1, 2]]);
    }
}
