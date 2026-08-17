# Connectivity via Good Extensions in Supersingular $(2,2)$-Isogeny Graphs

This repository contains the research code accompanying the paper

> **An Algorithm for Verifying Connectivity via Good Extensions in Supersingular $(2,2)$-Isogeny Graphs**  
> Masahiro Inoue, Ryo Ohashi, and Tsuyoshi Takagi

The code constructs and analyzes modified good-extension graphs for principally
polarized abelian surfaces. The main implementation merges parallel
$(2,2)$-isogenies with the same domain and codomain and retains only those
transitions that are valid for every isogeny represented by a merged edge. This
allows connectivity via good extensions to be certified without explicitly
storing isogeny kernels.

The repository supports the following settings:

- the full superspecial graph, including both genus-2 Jacobians and products of
  elliptic curves;
- the superspecial subgraph induced by genus-2 Jacobians; and
- the supersingular non-superspecial graph.

The implementation is intended for research and reproducibility. It has not
been audited for use in production cryptographic software.

## Method

A vertex of the modified good-extension graph is an ordered pair of isomorphism
classes of surfaces,

$(A_{\mathrm{prev}},A_{\mathrm{curr}})$


obtained by merging parallel isogenies from $A_{\mathrm{prev}}$ to
$A_{\mathrm{curr}}$. A transition

$
(A_0,A_1)\longrightarrow(A_1,A_2)
$

is retained only when every represented isogeny from $A_0$ to $A_1$ admits
a good extension from $A_1$ to $A_2$. The implementation computes this
condition by intersecting the successor sets arising from parallel isogenies.

Connectivity is then analyzed using strongly connected components (SCCs).
Tarjan's algorithm is used for the SCC decomposition. In particular, the code
checks whether a nontrivial SCC covers every surface of the original isogeny
graph, which is the sufficient condition proved in the accompanying paper.

## Repository Layout

The main files are expected to have the following layout:

```text
.
├── Cargo.toml
├── build.rs
├── primes.txt
├── run_experiments.py
├── sort_csv.py
├── LICENSE
├── data/
└── src/
    ├── lib.rs
    ├── counting.rs
    ├── curves.rs
    ├── fq.rs
    ├── graph.rs
    ├── poly.rs
    ├── theta.rs
    ├── utils.rs
    └── bin/
        ├── ssp_ge_graph.rs
        ├── ssg_ge_graph.rs
        └── ssp_edge_reachability.rs
```

The two primary executables are:

- `ssp_ge_graph`: superspecial experiments, with or without products of
  elliptic curves;
- `ssg_ge_graph`: supersingular non-superspecial experiments.

`ssp_edge_reachability` is an auxiliary exhaustive verifier for small
superspecial cases in which the coarse SCC-based criterion is inconclusive.

## Requirements

- a Rust toolchain supporting Rust edition 2024;
- Cargo;
- Python 3;
- sufficient memory for the selected prime range.

The Rust dependencies are installed automatically by Cargo. The main experiment
driver uses only the Python standard library. The optional `sort_csv.py` utility
requires `pandas`.

## Quick Start

### 1. Configure the experiment

Edit the configuration block near the top of `run_experiments.py`:

```python
# MODE = "ssp_only_jacobians"
MODE = "ssp_all"
# MODE = "ssg"

PRIME_MIN = 7
PRIME_MAX = 100
MAX_WORKERS = 1
```

The available settings are:

| Parameter | Example | Description |
| :--- | :--- | :--- |
| `MODE` | `"ssp_all"` | Full superspecial graph: Jacobians and products of elliptic curves. |
| `MODE` | `"ssp_only_jacobians"` | Superspecial graph restricted to genus-2 Jacobians. |
| `MODE` | `"ssg"` | Supersingular non-superspecial graph. |
| `PRIME_MIN` | `7` | Smallest prime to test, inclusive. The implementation requires $p>5$. |
| `PRIME_MAX` | `100` | Largest prime to test, inclusive. |
| `MAX_WORKERS` | `1` | Maximum number of experiments run concurrently. |

### 2. Run the experiments

```bash
python3 run_experiments.py
```

For the configuration above, the script tests every prime
$p\in[7,100]$ in the selected mode, with at most one worker running at a
time.

The script performs the following steps:

1. generates `primes.txt` using the sieve of Eratosthenes;
2. builds the Rust project in release mode;
3. runs the appropriate executable for every generated prime; and
4. appends one result row per prime to a CSV file under `data/`.

The output file is selected automatically:

| Mode | Executable | Output file |
| :--- | :--- | :--- |
| `ssp_all` | `ssp_ge_graph` | `data/ssp_all_results.csv` |
| `ssp_only_jacobians` | `ssp_ge_graph` | `data/ssp_only_jacobians_results.csv` |
| `ssg` | `ssg_ge_graph` | `data/ssg_results.csv` |

CSV writes are protected by an exclusive file lock, so rows produced by
concurrent workers are not interleaved. Their order is not guaranteed when
`MAX_WORKERS > 1`; use `sort_csv.py`, after setting its input filename, if a
prime-sorted file is required.

### Compile-time prime dispatch

The finite-field parameters are Rust const generics. Consequently, `build.rs`
generates a dispatch table from `primes.txt` at compile time. If a prime is not
listed in `primes.txt`, add it and rebuild before invoking a binary directly:

```bash
cargo build --release
```

The Python driver handles this automatically by regenerating `primes.txt`
before building.

## Running a Single Prime

After ensuring that the desired prime is present in `primes.txt` and rebuilding,
the primary binaries can be called directly.

Full superspecial graph:

```bash
./target/release/ssp_ge_graph 17 data/ssp_all_results.csv false
```

Superspecial Jacobian subgraph:

```bash
./target/release/ssp_ge_graph 17 data/ssp_only_jacobians_results.csv true
```

Supersingular non-superspecial graph:

```bash
./target/release/ssg_ge_graph 7 data/ssg_results.csv false
```

The final Boolean argument is `only_jacobians`. It is `true` only for
`ssp_only_jacobians`.

## Primary CSV Output

The superspecial modes use the field parameter

\[
\mathbb F_{p^2}=\mathbb F_p[x]/(x^2-d),
\]

where $d\in\mathbb F_p$ is a quadratic nonresidue. The non-superspecial mode
also uses

\[
\mathbb F_{p^4}=\mathbb F_{p^2}[y]/(y^2-D_2),
\qquad D_2=d_{2,0}+d_{2,1}x\in\mathbb F_{p^2}.
\]

The primary result files contain the following columns. `d` occurs in the
superspecial files, whereas `d1`, `d2c0`, and `d2c1` occur in the
supersingular non-superspecial file.

| Column | Meaning |
| :--- | :--- |
| `p` | Characteristic of the base field. |
| `d` | Quadratic nonresidue defining $\mathbb F_{p^2}$ in a superspecial experiment. |
| `d1` | Quadratic nonresidue defining $\mathbb F_{p^2}$ in a non-superspecial experiment. |
| `d2c0`, `d2c1` | Coefficients of $D_2=d_{2,0}+d_{2,1}x$, used to define $\mathbb F_{p^4}$ over $\mathbb F_{p^2}$. |
| `construction_time_sec` | Wall-clock time used to enumerate the surfaces and construct the modified good-extension graph. |
| `variety_count` | Number of isomorphism classes of principally polarized abelian surfaces enumerated in the selected mode. |
| `vertex_count` | Number of ordered surface pairs in the modified good-extension graph. Parallel isogenies are merged. |
| `edge_count` | Number of directed transitions in the modified good-extension graph after intersecting the successor sets of parallel isogenies. |
| `connectivity_check_time_sec` | Wall-clock time used for SCC decomposition and the subsequent connectivity checks. |
| `is_strongly_connected` | Whether the entire modified good-extension graph consists of one SCC. |
| `scc_count` | Number of SCCs in the modified good-extension graph. |
| `are_all_varieties_covered_by_main_scc` | Whether every enumerated variety occurs as an endpoint of at least one vertex in the largest SCC of the modified good-extension graph. |
| `are_all_variety_pairs_connected_via_good_paths` | Whether the implemented main/minor-SCC compatibility checks certify good-extension reachability from every ordered variety pair represented as a modified-graph vertex to every other such pair. |
| `are_all_minor_sccs_size_1` | Whether every SCC other than the largest SCC contains exactly one modified-graph vertex. |
| `are_all_minor_nodes_self_isogenies` | Whether every vertex $(A,B)$ in a minor SCC satisfies $A=B$. |
| `all_minor_sccs_have_preds_from_main` | Whether every minor SCC has a direct predecessor edge from the largest SCC in the condensed graph. |
| `any_minor_sccs_have_succs_to_main` | Whether at least one minor SCC has a direct successor edge to the largest SCC in the condensed graph. |
| `exists_minor_to_minor_edge` | Whether a minor SCC has an outgoing condensed-graph edge to a non-main SCC, including a possible self-loop. |
| `condensed_graph` | Adjacency list of the SCC condensation graph. SCC identifiers are implementation-generated indices. |
| `construction_peak_bytes` | Estimated heap footprint, in bytes, of the principal data structures retained during graph construction. |
| `connectivity_check_peak_bytes` | Estimated heap footprint, in bytes, of the principal data structures retained during connectivity analysis. |

The two memory columns are estimates obtained from the recursively measured
sizes of the relevant Rust data structures; they are not operating-system peak
resident-set-size measurements.

### Interpreting the connectivity fields

`is_strongly_connected = true` is a strong and easy-to-interpret outcome: every
modified-graph vertex reaches every other modified-graph vertex.

The accompanying paper uses the weaker SCC-coverage criterion represented by
`are_all_jacobians_connected_via_main_scc`. When this field is `true` and the
largest SCC is nontrivial, the theorem in the paper certifies vertex-to-vertex
connectivity of the underlying isogeny graph via good extensions. Thus, the
modified graph itself need not be strongly connected for the sufficient
condition to hold.

## Exhaustive Marked-Edge Verification for Small Superspecial Cases

`ssp_edge_reachability` is an auxiliary verifier for cases in which the coarse
SCC-based analysis of `ssp_ge_graph` is inconclusive. It is currently intended
primarily for the small `ssp_all` cases $p=17$ and $p=23$.

Unlike the primary modified graph, the refined state graph retains the
projectively normalized theta null-point induced by the preceding isogeny.
Compatibility is therefore preserved throughout the search without an explicit
comparison of isogeny kernels. Because this state space is much larger, this
verifier is not intended to replace `ssp_ge_graph` for large parameter ranges.

Ensure that `17` and `23` occur in `primes.txt`, rebuild, and run:

```bash
./target/release/ssp_edge_reachability \
  17 data/ssp_all_edge_reachability_results.csv false

./target/release/ssp_edge_reachability \
  23 data/ssp_all_edge_reachability_results.csv false
```

The exhaustive output contains the following columns:

| Column | Meaning |
| :--- | :--- |
| `p`, `d` | Characteristic and quadratic nonresidue defining $\mathbb F_{p^2}$. |
| `construction_time_sec` | Time used to enumerate surfaces and close the marked state space under good extensions. |
| `reachability_check_time_sec` | Time used for SCC decomposition and coarse-pair reachability searches. |
| `variety_count` | Number of enumerated superspecial surfaces. |
| `coarse_pair_count` | Number of distinct ordered surface pairs after forgetting the theta marking. |
| `marked_edge_state_count` | Number of states $(A,B,\theta)$ retaining the normalized theta marking. |
| `marked_transition_count` | Number of good-extension transitions between marked states. |
| `marked_scc_count` | Number of SCCs in the marked-state graph. |
| `all_marked_edges_mutually_reachable` | Whether the marked-state graph is strongly connected. This is a stronger condition than coarse-pair reachability. |
| `all_coarse_pairs_connected_via_good_paths` | Whether, for every ordered pair of coarse surface pairs, some compatible marked state over the source reaches some marked state over the target. |
| `first_unreachable_source_pair`, `first_unreachable_target_pair` | The first counterexample found, left empty when all coarse pairs are connected. |
| `construction_peak_bytes` | Estimated heap footprint during refined-graph construction. |
| `reachability_check_peak_bytes` | Estimated heap footprint during the reachability checks. |

It is possible to obtain

```text
all_marked_edges_mutually_reachable = false
all_coarse_pairs_connected_via_good_paths = true
```

without a contradiction. The first field quantifies over every pair of marked
states, whereas the second permits an appropriate marked representative to be
chosen over each source and target surface pair.

If $M$ is the number of marked states, $E$ the number of marked transitions,
and $Q$ the number of coarse surface pairs, the SCC decomposition costs
$O(M+E)$. The current multi-source breadth-first searches cost
$O(Q(M+E))$ in total and use $O(M+E+Q)$ memory. Each marked state has at
most eight outgoing good-extension transitions, so $E=O(M)$.

## Reproducibility Notes

- Existing CSV files are opened in append mode. Move or remove an old result
  file before starting a clean run if duplicate rows are undesirable.
- Runtime and memory values depend on the compiler, machine, operating system,
  worker count, and background load.
- Release builds should be used for reported experiments.
- Increasing `MAX_WORKERS` can reduce total elapsed time but increases aggregate
  memory consumption approximately in proportion to the number of simultaneous
  processes.
- The non-superspecial graph grows substantially faster than the superspecial
  graph; choose its prime range conservatively.

For basic validation of a checkout, run:

```bash
cargo test
cargo fmt --check
cargo clippy --all-targets
```

## Citation

If this software contributes to published work, please cite the accompanying
paper. Until final bibliographic information is available, the following entry
may be used and updated when the venue, volume, pages, or DOI are assigned:

```bibtex
@unpublished{inoue2026goodextensions,
  author = {Masahiro Inoue and Ryo Ohashi and Tsuyoshi Takagi},
  title  = {An Algorithm for Verifying Connectivity via Good Extensions in
            Supersingular (2,2)-Isogeny Graphs},
  year   = {2026},
  note   = {Manuscript}
}
```

## License

This software is released under the [MIT License](LICENSE). The copyright
notices are:

```text
Copyright (c) 2025 Ryo Ohashi and Hiroshi Onuki
Copyright (c) 2026 Masahiro Inoue
```

The MIT License permits use, copying, modification, merging, publication,
distribution, sublicensing, and sale of copies, subject to preservation of the
copyright and permission notices. See `LICENSE` for the complete license terms,
including the warranty and liability disclaimer.

