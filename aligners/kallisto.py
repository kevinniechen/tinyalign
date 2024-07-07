import random
from collections import defaultdict, Counter
from Bio import SeqIO
import numpy as np

def hash_kmer(kmer, seed=0):
    """Simple rolling hash function for k-mers."""
    h = seed
    for c in kmer:
        h = (h * 31 + ord(c)) & 0xFFFFFFFF
    return h

def build_index(reference_file, k=31):
    """Build a hash table index of k-mers from the transcriptome."""
    index = defaultdict(list)
    transcript_lengths = {}
    for record in SeqIO.parse(reference_file, "fasta"):
        seq = str(record.seq)
        transcript_lengths[record.id] = len(seq)
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            index[hash_kmer(kmer)].append((record.id, i))
    return index, transcript_lengths

def pseudo_align(read, index, k):
    """Perform pseudo-alignment of a read to the transcriptome."""
    compatible_transcripts = set()
    for i in range(min(len(read) - k + 1, 5)):
        kmer = read[i:i+k]
        kmer_hash = hash_kmer(kmer)
        if kmer_hash in index:
            compatible_transcripts.update(tx for tx, _ in index[kmer_hash])
    return compatible_transcripts

def expectation_maximization(compatibility_matrix, transcript_lengths, max_iter=100, epsilon=1e-5):
    """Perform Expectation-Maximization to estimate transcript abundances."""
    num_transcripts = len(transcript_lengths)
    num_reads = len(compatibility_matrix)
    
    # Create a mapping from transcript IDs to indices
    transcript_to_index = {tx: i for i, tx in enumerate(transcript_lengths.keys())}
    
    # Initialize abundances
    abundances = np.ones(num_transcripts) / num_transcripts
    
    for _ in range(max_iter):
        old_abundances = abundances.copy()
        
        # Expectation step
        probabilities = np.zeros((num_reads, num_transcripts))
        for i, compatible_transcripts in enumerate(compatibility_matrix):
            if compatible_transcripts:
                indices = [transcript_to_index[tx] for tx in compatible_transcripts]
                prob = abundances[indices] / np.sum(abundances[indices])
                probabilities[i, indices] = prob
        
        # Maximization step
        abundances = np.sum(probabilities, axis=0)
        abundances /= np.array(list(transcript_lengths.values()))  # Normalize by transcript length
        abundances /= np.sum(abundances)  # Ensure abundances sum to 1
        
        # Check for convergence
        if np.max(np.abs(abundances - old_abundances)) < epsilon:
            break
    
    return abundances

def quantify(index, transcript_lengths, reads_file, k=31):
    """Quantify transcript abundance using pseudo-alignment and EM."""
    compatibility_matrix = []
    total_reads = 0
    
    for record in SeqIO.parse(reads_file, "fasta"):
        read = str(record.seq)
        compatible_transcripts = pseudo_align(read, index, k)
        compatibility_matrix.append(compatible_transcripts)
        total_reads += 1
    
    # Perform EM
    abundances = expectation_maximization(compatibility_matrix, transcript_lengths)
    
    # Calculate TPM
    tpm = {}
    scaling_factor = 1e6 / np.sum(abundances)
    for tx, abundance in zip(transcript_lengths.keys(), abundances):
        tpm[tx] = abundance * scaling_factor
    
    return tpm, total_reads

def run(reference_file, reads_file):
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

    return tpm

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python kallisto.py <reference_file> <reads_file>")
        sys.exit(1)
    run(sys.argv[1], sys.argv[2])