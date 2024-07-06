import sys
from Bio import SeqIO

def create_suffix_array(reference):
    """Create a suffix array for the given reference."""
    suffixes = sorted(range(len(reference)), key=lambda i: reference[i:])
    return suffixes

def align_with_suffix_array(read, reference, suffix_array):
    """Align the read using the suffix array."""
    best_pos = -1
    best_score = -1
    
    for i in suffix_array:
        if i + len(read) > len(reference):
            continue
        score = sum(1 for j in range(len(read)) if read[j] == reference[i+j])
        if score > best_score or (score == best_score and i < best_pos):
            best_score = score
            best_pos = i
    
    return best_pos, best_score

def main(reference_file, reads_file):
    # Read the reference genome
    reference = str(next(SeqIO.parse(reference_file, "fasta")).seq)

    # Create the suffix array
    print("Building suffix array...")
    suffix_array = create_suffix_array(reference)
    print("Suffix array built.")

    # Process each read
    with open(reads_file, "r") as reads:
        for record in SeqIO.parse(reads, "fasta"):
            read = str(record.seq)
            pos, score = align_with_suffix_array(read, reference, suffix_array)
            print(f"Read: {record.id}, Position: {pos}, Score: {score}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python star.py data/reference.fasta data/reads.fasta")
        sys.exit(1)
    
    reference_file = sys.argv[1]
    reads_file = sys.argv[2]
    main(reference_file, reads_file)
