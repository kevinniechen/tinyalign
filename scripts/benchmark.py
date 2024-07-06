import time
import subprocess
import statistics

def run_script(script_name, reference_file, reads_file, num_runs=5):
    times = []
    for _ in range(num_runs):
        start_time = time.time()
        subprocess.run(['python', f'aligners/{script_name}', reference_file, reads_file], 
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        end_time = time.time()
        times.append(end_time - start_time)
    
    avg_time = statistics.mean(times)
    std_dev = statistics.stdev(times) if len(times) > 1 else 0
    
    print(f"{script_name} average execution time: {avg_time:.4f} seconds")
    print(f"{script_name} standard deviation: {std_dev:.4f} seconds")
    print(f"Individual run times: {', '.join([f'{t:.4f}' for t in times])}")
    print()

def main():
    reference_file = 'data/reference.fasta'
    reads_file = 'data/reads.fasta'
    
    print("Running benchmarks...\n")
    run_script('naive.py', reference_file, reads_file)
    run_script('star.py', reference_file, reads_file)

if __name__ == "__main__":
    main()
