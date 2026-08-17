use std::env;
use std::fs;
use std::path::Path;

fn main() {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let primes_path = Path::new(&manifest_dir).join("primes.txt");

    if !primes_path.exists() {
        println!("cargo:warning=primes.txt not found, skipping generation.");
        return;
    }

    let primes_content = fs::read_to_string(&primes_path).expect("Failed to read primes.txt");

    let mut match_expr = String::from("match p_arg {\n");

    for num_str in primes_content.split(',') {
        let p = num_str.trim();
        if !p.is_empty() {
            match_expr.push_str(&format!(
                "    {} => analyze_connectivity!({}, filename_arg, only_jacobians),\n",
                p, p
            ));
        }
    }
    match_expr.push_str(
        "    _ => { 
        eprintln!(\"❌ Error: Prime {} is not in the compiled list.\", p_arg);
        eprintln!(\"👉 Action required: Add {} to 'primes.txt' and rebuild.\", p_arg);
        std::process::exit(1); 
    }\n",
    );
    match_expr.push_str("}");

    let out_dir = env::var("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("generated_dispatch.rs");
    fs::write(&dest_path, match_expr).unwrap();

    println!("cargo:rerun-if-changed=primes.txt");
}