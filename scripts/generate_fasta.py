import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_random_sequence(length, case='upper'):
    seq = ''.join(random.choice('ACGT') for _ in range(length))
    return seq.upper() if case == 'upper' else seq.lower()

def generate_genome_with_introns_exons(num_exons, exon_length, intron_length):
    genome = ""
    exon_positions = []
    for i in range(num_exons):
        exon = generate_random_sequence(exon_length, 'upper')
        genome += exon
        exon_start = len(genome) - exon_length
        exon_end = len(genome)
        exon_positions.append((exon_start, exon_end))
        if i < num_exons - 1:
            intron = generate_random_sequence(intron_length, 'lower')
            genome += intron
    return genome, exon_positions

def introduce_mismatches(seq, num_mismatches):
    seq_list = list(seq)
    positions = random.sample(range(len(seq)), num_mismatches)
    mismatches = []
    for pos in positions:
        original = seq_list[pos]
        new = random.choice([base for base in 'ACGT' if base != original])
        seq_list[pos] = new
        mismatches.append(f"{original}{pos+1}{new}")
    return ''.join(seq_list), mismatches

def generate_reads(genome, exon_positions, num_reads, read_length, max_mismatches=1):
    reads = []
    spliced_count = 0
    unspliced_count = 0
    for i in range(num_reads):
        if random.random() < 0.5:  # 50% chance of spliced read
            exon1, exon2 = sorted(random.sample(exon_positions, 2))
            start1 = random.randint(exon1[0], exon1[1] - read_length // 2)
            start2 = random.randint(exon2[0], exon2[1] - read_length // 2)
            read_seq = genome[start1:start1 + read_length // 2] + genome[start2:start2 + read_length // 2]
            description = f"exon1:{start1}-{start1+read_length//2},exon2:{start2}-{start2+read_length//2}"
            read_type = "S"  # S for spliced
            spliced_count += 1
        else:  # unspliced read
            start = random.randint(0, len(genome) - read_length)
            read_seq = genome[start:start + read_length]
            description = f"pos:{start}-{start+read_length}"
            read_type = "U"  # U for unspliced
            unspliced_count += 1
        
        # Introduce mismatches
        num_mismatches = random.randint(0, max_mismatches)
        read_seq_with_mismatches, mismatches = introduce_mismatches(read_seq, num_mismatches)
        
        # Create read ID with simplified mismatch information
        mismatch_info = f"M{num_mismatches}" if num_mismatches > 0 else "M0"
        read_id = f"R{i+1}-{read_type}-{mismatch_info}"
        
        reads.append(SeqRecord(Seq(read_seq_with_mismatches), id=read_id, description=description))
    return reads, spliced_count, unspliced_count

def main():
    # Generate genome
    num_exons = 5
    exon_length = 200
    intron_length = 800
    genome, exon_positions = generate_genome_with_introns_exons(num_exons, exon_length, intron_length)
    
    # Create reference.fasta
    with open("data/reference.fasta", "w") as f:
        SeqIO.write(SeqRecord(Seq(genome), id="ref", description=""), f, "fasta")
    
    # Generate reads
    num_reads = 100
    read_length = 100
    max_mismatches = 1
    reads, spliced_count, unspliced_count = generate_reads(genome, exon_positions, num_reads, read_length, max_mismatches)
    
    # Create reads.fasta
    with open("data/reads.fasta", "w") as f:
        SeqIO.write(reads, f, "fasta")

    print(f"Generated genome with {num_exons} exons and {num_exons - 1} introns")
    print(f"Exon length: {exon_length}, Intron length: {intron_length}")
    print(f"Generated {num_reads} reads of length {read_length}")
    print(f"Spliced reads: {spliced_count}, Unspliced reads: {unspliced_count}")
    print(f"Maximum mismatches per read: {max_mismatches}")
    print("Files 'reference.fasta' and 'reads.fasta' have been created in the 'data' directory")

    # Verify spliced read
    spliced_read = next(read for read in reads if "-S-" in read.id)
    parts = spliced_read.description.split(",")
    exon1_start, exon1_end = map(int, parts[0].split(":")[1].split("-"))
    exon2_start, exon2_end = map(int, parts[1].split(":")[1].split("-"))
    
    print("\nVerifying spliced read:")
    print(f"Read ID: {spliced_read.id}")
    print(f"Read sequence: {spliced_read.seq}")
    print(f"Exon 1 from reference: {genome[exon1_start:exon1_end]}")
    print(f"Exon 2 from reference: {genome[exon2_start:exon2_end]}")
    print(f"Combined exons: {genome[exon1_start:exon1_end] + genome[exon2_start:exon2_end]}")
    
    mismatches = int(spliced_read.id.split("-")[-1][1])
    print(f"Number of mismatches: {mismatches}")

    # Print a section of the reference genome to show exon/intron distinction
    print("\nSection of reference genome (showing exon/intron distinction):")
    print(genome[:500])  # Print the first 500 characters to show the pattern

if __name__ == "__main__":
    main()
