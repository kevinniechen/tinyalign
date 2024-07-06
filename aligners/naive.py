import sys
from Bio import SeqIO

def simple_align(read, reference):
    best_pos = -1
    best_score = -1
    
    for i in range(len(reference) - len(read) + 1):
        score = sum(1 for j in range(len(read)) if read[j] == reference[i+j])
        if score > best_score:
            best_score = score
            best_pos = i
    
    return best_pos, best_score

def main(reference_file, reads_file):
    # Read the reference genome
    reference = next(SeqIO.parse(reference_file, "fasta")).seq

    # Process each read
    with open(reads_file, "r") as reads:
        for record in SeqIO.parse(reads, "fasta"):
            read = record.seq
            pos, score = simple_align(read, reference)
            print(f"Read: {record.id}, Position: {pos}, Score: {score}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python naive.py data/reference.fasta data/reads.fasta")
        sys.exit(1)
    
    reference_file = sys.argv[1]
    reads_file = sys.argv[2]
    main(reference_file, reads_file)
