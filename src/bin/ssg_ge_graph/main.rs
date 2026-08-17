//! Analyze the connectivity of a supersingular non-superspecial (2,2)-isogeny graph through good extension paths
//! cargo run --release --bin ssg_ge_graph -- <p> [output_file]

mod ssg_ge_graph;
use g2_isogeny::fq::{Fp, Fp2};
use g2_isogeny::utils::{RowDataSsg, save_to_csv};
use ssg_ge_graph::{check_ge_path_connectivity, compute_ssg_good_extension_graph};
use std::{env, process};

macro_rules! analyze_connectivity {
    ($P:expr, $FNAME:expr, $ONLY_JACOBIANS:expr) => {{
        const D1: u32 = Fp::<$P>::find_non_residue();
        const D2: (u32, u32) = Fp2::<$P, D1>::find_non_residue();
        const D2C0: u32 = D2.0;
        const D2C1: u32 = D2.1;

        let start_time = chrono::Local::now();
        println!(
            "[{}] Starting p = {}...",
            start_time.format("%Y-%m-%d %H:%M:%S"),
            $P
        );

        let (inv_to_id, ram_list, node_list, ge_graph, construction_time, construction_peak_bytes) =
            compute_ssg_good_extension_graph::<$P, D1, D2C0, D2C1>();

        let (result, check_time, check_peak_bytes) = check_ge_path_connectivity::<$P, D1, D2C0, D2C1>(
            &inv_to_id, &ram_list, &node_list, &ge_graph,
        );

        let row = RowDataSsg {
            p: $P,
            d1: D1,
            d2c0: D2C0,
            d2c1: D2C1,
            construction_time_sec: format!("{:.4}", construction_time.as_secs_f64()),
            variety_count: result.variety_count,
            vertex_count: result.vertex_count,
            edge_count: result.edge_count,
            connectivity_check_time_sec: format!("{:.4}", check_time.as_secs_f64()),
            is_strongly_connected: result.is_strongly_connected,
            scc_count: result.scc_count,
            are_all_varieties_covered_by_main_scc: result
                .are_all_varieties_covered_by_main_scc,
            are_all_variety_pairs_connected_via_good_paths: result
                .are_all_variety_pairs_connected_via_good_paths,
            are_all_minor_sccs_size_1: result.are_all_minor_sccs_size_1,
            are_all_minor_nodes_self_isogenies: result.are_all_minor_nodes_self_isogenies,
            all_minor_sccs_have_preds_from_main: result.all_minor_sccs_have_preds_from_main,
            any_minor_sccs_have_succs_to_main: result.any_minor_sccs_have_succs_to_main,
            exists_minor_to_minor_edge: result.exists_minor_to_minor_edge,
            condensed_graph: result
                .condensed_graph
                .adj
                .iter()
                .enumerate()
                .map(|(i, succs)| format!("{}:{:?}", i, succs))
                .collect::<Vec<_>>()
                .join(", "),
            construction_peak_bytes: construction_peak_bytes,
            connectivity_check_peak_bytes: check_peak_bytes,
        };

        let end_time = chrono::Local::now();
        let duration = end_time - start_time;
        println!(
            "[{}] Finished p = {} (Time: {}h {}m {}s)",
            end_time.format("%Y-%m-%d %H:%M:%S"),
            $P,
            duration.num_hours(),
            duration.num_minutes() % 60,
            duration.num_seconds() % 60
        );

        if let Err(e) = save_to_csv($FNAME, row) {
            eprintln!("❌ Failed to save CSV for p={}: {}", $P, e);
        } else {
            println!("✅ Result for p={} saved to {}", $P, $FNAME);
        }
    }};
}

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: cargo run --release --bin ssg_ge_graph -- <p> [output_file]");
        eprintln!("Example: cargo run --release --bin ssg_ge_graph -- 11 ssg_results.csv");
        process::exit(1);
    }

    let p_arg: u32 = match args[1].parse() {
        Ok(num) => num,
        Err(_) => {
            eprintln!("Error: Argument <p> must be an integer. Got: {}", args[1]);
            process::exit(1);
        }
    };

    let filename_arg = if args.len() > 2 {
        &args[2]
    } else {
        "ssg_results.csv"
    };

    include!(concat!(env!("OUT_DIR"), "/generated_dispatch.rs"));
}
