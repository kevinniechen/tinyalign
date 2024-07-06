import sys
from Bio import SeqIO

def create_suffix_array(reference):
    """Create a suffix array for the given reference."""
    suffixes = sorted(range(len(reference)), key=lambda i: reference[i:])
    return suffixes

def binary_search(pattern, reference, suffix_array):
    """Perform binary search on the suffix array."""
    left, right = 0, len(suffix_array) - 1
    while left <= right:
        mid = (left + right) // 2
        suffix_start = suffix_array[mid]
        suffix = reference[suffix_start:suffix_start+len(pattern)]
        if pattern == suffix[:len(pattern)]:
            return mid
        elif pattern < suffix[:len(pattern)]:
            right = mid - 1
        else:
            left = mid + 1
    return -1

def align_with_suffix_array(read, reference, suffix_array):
    """Align the read using binary search on the suffix array."""
    best_pos = -1
    best_score = -1
    read_length = len(read)
    
    index = binary_search(read, reference, suffix_array)
    if index != -1:
        best_pos = suffix_array[index]
        best_score = read_length
    else:
        # If exact match not found, search for best partial match
        left, right = 0, len(suffix_array) - 1
        while left <= right:
            mid = (left + right) // 2
            suffix_start = suffix_array[mid]
            suffix = reference[suffix_start:suffix_start+read_length]
            score = sum(1 for i in range(min(read_length, len(suffix))) if read[i] == suffix[i])
            if score > best_score or (score == best_score and suffix_start < best_pos):
                best_score = score
                best_pos = suffix_start
            if read < suffix[:read_length]:
                right = mid - 1
            else:
                left = mid + 1
    
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
            print(f"Read: {record.id}, Length: {len(read)}, Position: {pos}, Score: {score}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python star.py data/reference.fasta data/reads.fasta")
        sys.exit(1)
    
    reference_file = sys.argv[1]
    reads_file = sys.argv[2]
    main(reference_file, reads_file)
