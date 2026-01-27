import subprocess
from concurrent.futures import ThreadPoolExecutor
import time
import os
import sys
from datetime import datetime

# --- Configuration --------------------------------------
# Select mode: "ssp" for superspecial, "ssg" for supersingular non-superspecial
MODE = "ssp"  # Super Special
# MODE = "ssg"  # Super Singular

PRIME_MIN = 7  # Starting prime (inclusive)
PRIME_MAX = 100  # Ending prime (inclusive)

MAX_WORKERS = 1  # Number of parallel workers


# --- Do not modify below this line unless necessary ----
PRIMES_FILE = "primes.txt"  # File containing the list of primes
BINARY_NAME = f"{MODE}_ge_graph"  # Name of the compiled binary
BINARY_PATH = f"./target/release/{BINARY_NAME}"  # Path to the compiled binary
OUTPUT_FILE = f"data/{MODE}_results.csv"


def generate_primes_file(min_val, max_val, filename):
    """Generate primes using the Sieve of Eratosthenes and save to a file"""
    print(f"🔄 Generating primes between {min_val} and {max_val}...")

    if max_val < 2:
        print("⚠️ Max value must be >= 2.")
        return

    # Sieve of Eratosthenes
    is_prime = [True] * (max_val + 1)
    is_prime[0] = is_prime[1] = False

    for i in range(2, int(max_val**0.5) + 1):
        if is_prime[i]:
            for j in range(i * i, max_val + 1, i):
                is_prime[j] = False

    primes = [str(i) for i in range(max_val + 1) if is_prime[i] and i >= min_val]

    # Write to file (comma-separated)
    with open(filename, "w") as f:
        f.write(",".join(primes))

    print(f"💾 Saved {len(primes)} primes to {filename}.\n")


def load_primes():
    """Load the list of primes from primes.txt"""
    if not os.path.exists(PRIMES_FILE):
        print(f"Error: {PRIMES_FILE} not found.")
        return []

    with open(PRIMES_FILE, "r") as f:
        content = f.read()

    primes = [p.strip() for p in content.replace("\n", ",").split(",") if p.strip().isdigit()]
    return primes


def build_project():
    """Build the project once at the beginning"""
    print("🔨 Building project...")
    try:
        # Building is essential to reflect changes in primes.txt
        subprocess.run(["cargo", "build", "--release"], check=True)
        print("✅ Build successful.\n")
    except subprocess.CalledProcessError:
        print("❌ Build failed. Aborting experiments.")
        sys.exit(1)


def run_experiment(p):
    """Run the compiled binary directly with the given prime p"""
    start_time = time.time()

    try:
        cmd = [BINARY_PATH, str(p), OUTPUT_FILE]

        subprocess.run(cmd, check=True)

        elapsed = time.time() - start_time
        m, s = divmod(elapsed, 60)
        h, m = divmod(m, 60)

    except subprocess.CalledProcessError as e:
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] ❌ Error occurred for p={p}: {e}")
    except Exception as e:
        print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] ❌ Unexpected error for p={p}: {e}")


def main():
    output_dir = os.path.dirname(OUTPUT_FILE)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"📂 Created directory: {output_dir}")

    # 1. Generate primes file
    generate_primes_file(PRIME_MIN, PRIME_MAX, PRIMES_FILE)

    # 2. Build the project once at the beginning (this also updates generated_dispatch.rs)
    build_project()

    # 3. Load the list of primes
    primes = load_primes()
    if not primes:
        print("No primes found in primes.txt.")
        return

    print(f"📋 Loaded {len(primes)} primes: {primes}")
    print(f"🔥 Starting execution with {MAX_WORKERS} parallel workers...\n")

    # 4. Run the binary in parallel
    try:
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            executor.map(run_experiment, primes)

    except KeyboardInterrupt:
        print("\n❌ Interrupted by user! Stopping all processes...")
        os._exit(1)

    print("\n🎉 All experiments completed.")


if __name__ == "__main__":
    main()
