import random
import sys
import re
from collections import defaultdict, Counter
from Bio import SeqIO
import numpy as np

def hash_kmer(kmer, seed=0):
    """Simple rolling hash function for k-mers."""
    h = seed
    for c in kmer:
        h = (h * 31 + ord(c)) & 0xFFFFFFFF
    return h

def build_index(reference_file, k=15):
    """Build a hash table index of k-mers from the transcriptome."""
    index = defaultdict(list)
    transcript_lengths = {}
    for record in SeqIO.parse(reference_file, "fasta"):
        seq = str(record.seq)
        transcript_lengths[record.id] = len(seq)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            kmer_hash = hash_kmer(kmer)
            index[kmer_hash].append((record.id, i))
    print(f"Index size: {len(index)}")
    print(f"First 5 index entries: {list(index.items())[:5]}")
    return index, transcript_lengths

def parse_read_header(header):
    """Parse the read header to extract relevant information."""
    parts = header.split()
    read_id = parts[0]
    read_type = parts[1].split(':')[0] if len(parts) > 1 else "Unknown"
    return read_id, read_type

def pseudo_align(read, index, k):
    """Perform pseudo-alignment of a read to the transcriptome."""
    compatible_transcripts = set()
    for i in range(len(read) - k + 1):
        kmer = read[i:i+k]
        kmer_hash = hash_kmer(kmer)
        if kmer_hash in index:
            compatible_transcripts.update(tx for tx, _ in index[kmer_hash])
            print(f"Found match for k-mer: {kmer}")
    print(f"Read: {read[:20]}... Compatible transcripts: {compatible_transcripts}")
    return compatible_transcripts

def expectation_maximization(compatibility_matrix, transcript_lengths, max_iter=100, epsilon=1e-5):
    """Perform Expectation-Maximization to estimate transcript abundances."""
    num_transcripts = len(transcript_lengths)
    num_reads = len(compatibility_matrix)
    
    # Create a mapping from transcript IDs to indices
    transcript_to_index = {tx: i for i, tx in enumerate(transcript_lengths.keys())}
    
    # Initialize abundances
    abundances = np.ones(num_transcripts) / num_transcripts
    
    print(f"Initial abundances: {abundances}")
    
    for iteration in range(max_iter):
        old_abundances = abundances.copy()
        
        # Expectation step
        probabilities = np.zeros((num_reads, num_transcripts))
        for i, compatible_transcripts in enumerate(compatibility_matrix):
            if compatible_transcripts:
                indices = [transcript_to_index[tx] for tx in compatible_transcripts]
                prob_sum = np.sum(abundances[indices])
                if prob_sum > 0:
                    prob = abundances[indices] / prob_sum
                    probabilities[i, indices] = prob
        
        # Maximization step
        abundances = np.sum(probabilities, axis=0)
        transcript_lengths_array = np.array(list(transcript_lengths.values()))
        abundances = np.divide(abundances, transcript_lengths_array, out=np.zeros_like(abundances), where=transcript_lengths_array!=0)
        abundance_sum = np.sum(abundances)
        if abundance_sum > 0:
            abundances /= abundance_sum  # Ensure abundances sum to 1
        else:
            abundances = np.ones(num_transcripts) / num_transcripts  # Reset if all abundances are zero
        
        print(f"Iteration {iteration + 1} abundances: {abundances}")
        
        # Check for convergence
        if np.max(np.abs(abundances - old_abundances)) < epsilon:
            print(f"Converged after {iteration + 1} iterations")
            break
    
    return abundances

def quantify(index, transcript_lengths, reads_file, k=15):
    """Quantify transcript abundance using pseudo-alignment and EM."""
    compatibility_matrix = []
    total_reads = 0
    
    for record in SeqIO.parse(reads_file, "fasta"):
        read_id, read_type = parse_read_header(record.description)
        read = str(record.seq)
        compatible_transcripts = pseudo_align(read, index, k)
        compatibility_matrix.append(compatible_transcripts)
        total_reads += 1
        print(f"Read ID: {read_id}, Type: {read_type}")
    
    print(f"Compatibility matrix: {compatibility_matrix}")
    
    # Perform EM
    abundances = expectation_maximization(compatibility_matrix, transcript_lengths)
    
    print(f"Final abundances: {abundances}")
    
    # Calculate TPM
    tpm = {}
    rpk = {tx: abundance / (transcript_lengths[tx] / 1000) for tx, abundance in zip(transcript_lengths.keys(), abundances)}
    scaling_factor = 1e6 / sum(rpk.values())
    for tx, rpk_value in rpk.items():
        tpm[tx] = rpk_value * scaling_factor
    
    return tpm, total_reads

def main(reference_file, reads_file):
    print("Building index...")
    index, transcript_lengths = build_index(reference_file)
    print("Index built.")
    
    print("Quantifying...")
    tpm, total_reads = quantify(index, transcript_lengths, reads_file)
    print("Quantification complete.")
    
    print(f"Processed {total_reads} reads.")
    print("\nTranscript abundances (TPM):")
    for tx, abundance in tpm.items():
        print(f"{tx}: {abundance:.2f}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python kallisto.py <reference_file> <reads_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
