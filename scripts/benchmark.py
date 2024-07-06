import time
import subprocess

def run_script(script_name, reference_file, reads_file):
    start_time = time.time()
    result = subprocess.run(['python', f'aligners/{script_name}', reference_file, reads_file], 
                            capture_output=True, text=True)
    end_time = time.time()
    execution_time = end_time - start_time
    
    print(f"{script_name} execution time: {execution_time:.4f} seconds")
    
    # Analyze the output
    output_lines = result.stdout.split('\n')
    spliced_reads = sum(1 for line in output_lines if 'Splice:' in line)
    total_reads = sum(1 for line in output_lines if line.startswith('Read:'))
    print(f"Spliced reads: {spliced_reads}/{total_reads}")
    print(f"Output:\n{result.stdout}")
    print()

def main():
    reference_file = 'data/reference.fasta'
    reads_file = 'data/reads.fasta'
    
    print("Running benchmarks...\n")
    run_script('naive.py', reference_file, reads_file)
    run_script('star.py', reference_file, reads_file)

if __name__ == "__main__":
    main()
