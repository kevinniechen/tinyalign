import time
import subprocess

def run_script(script_name, aligner_type, **kwargs):
    start_time = time.time()
    if aligner_type == "genomic":
        result = subprocess.run(['python', f'aligners/{script_name}', kwargs['reference_file'], kwargs['reads_file']], 
                                capture_output=True, text=True)
    elif aligner_type == "transcriptome":
        result = subprocess.run(['python', f'aligners/{script_name}', kwargs['transcriptome_file'], kwargs['reads_file']], 
                                capture_output=True, text=True)
    elif aligner_type == "sample_comparison":
        result = subprocess.run(['python', f'aligners/{script_name}', kwargs['reference_file'], kwargs['sample1_file'], kwargs['sample2_file']], 
                                capture_output=True, text=True)
    end_time = time.time()
    execution_time = end_time - start_time
    
    print(f"{script_name} execution time: {execution_time:.4f} seconds")
    
    # Analyze the output
    output_lines = result.stdout.split('\n')
    if aligner_type == "genomic":
        spliced_reads = sum(1 for line in output_lines if 'Splice:' in line)
        total_reads = sum(1 for line in output_lines if line.startswith('Read:'))
        print(f"Spliced reads: {spliced_reads}/{total_reads}")
    print(f"Output:\n{result.stdout}")
    print()

def main():
    reference_file = 'data/reference.fasta'
    reads_file = 'data/reads.fasta'
    transcriptome_file = 'data/transcriptome.fasta'
    sample1_file = 'data/reads_sample1.fasta'
    sample2_file = 'data/reads_sample2.fasta'
    
    print("Running benchmarks...\n")
    
    print("Genomic Aligners:")
    # run_script('naive.py', "genomic", reference_file=reference_file, reads_file=reads_file)
    run_script('star.py', "genomic", reference_file=reference_file, reads_file=reads_file)
    
    print("Transcriptome Aligner:")
    run_script('kallisto.py', "transcriptome", transcriptome_file=transcriptome_file, reads_file=reads_file)
    
    print("Sample Comparison:")
    run_script('splash.py', "sample_comparison", reference_file=reference_file, sample1_file=sample1_file, sample2_file=sample2_file)

if __name__ == "__main__":
    main()
