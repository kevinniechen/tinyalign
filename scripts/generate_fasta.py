import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_random_sequence(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def main():
    # Generate reference sequence
    reference_length = 10000
    reference_seq = generate_random_sequence(reference_length)
    
    # Create reference.fasta
    with open("reference.fasta", "w") as f:
        SeqIO.write(SeqRecord(Seq(reference_seq), id="ref", description=""), f, "fasta")
    
    # Generate reads
    read_length = 20
    num_reads = 10
    reads = []
    
    for i in range(num_reads):
        start = random.randint(0, reference_length - read_length)
        read_seq = reference_seq[start:start+read_length]
        reads.append(SeqRecord(Seq(read_seq), id=f"read{i+1}", description=""))
    
    # Create reads.fasta
    with open("reads.fasta", "w") as f:
        SeqIO.write(reads, f, "fasta")

if __name__ == "__main__":
    main()
    print("Generated new reference.fasta and reads.fasta files.")
