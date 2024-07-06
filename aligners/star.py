import sys
from Bio import SeqIO

def create_suffix_array(reference):
    """Create a suffix array for the given reference."""
    suffixes = sorted(range(len(reference)), key=lambda i: reference[i:])
    return suffixes

def binary_search_prefix(pattern, reference, suffix_array):
    """Perform binary search on the suffix array to find the longest prefix match."""
    left, right = 0, len(suffix_array) - 1
    best_match_length = 0
    best_match_position = -1

    while left <= right:
        mid = (left + right) // 2
        suffix_start = suffix_array[mid]
        suffix = reference[suffix_start:]
        
        # Find the length of the matching prefix
        match_length = 0
        for a, b in zip(pattern, suffix):
            if a != b:
                break
            match_length += 1
        
        # Update best match if necessary
        if match_length > best_match_length:
            best_match_length = match_length
            best_match_position = suffix_start
        
        # Decide which half to search next
        if pattern[:match_length+1] <= suffix[:match_length+1]:
            right = mid - 1
        else:
            left = mid + 1

    return best_match_position, best_match_length

def align_with_suffix_array(read, reference, suffix_array, min_seed_length=10):
    """Align the read using binary search on the suffix array for longest prefix matches."""
    # Find the first match
    first_pos, first_length = binary_search_prefix(read, reference, suffix_array)
    
    # If there's a perfect match or no match at all, return
    if first_length == len(read):
        return [(first_pos, first_length)] if first_length >= min_seed_length else []
    elif first_length == 0:
        return []
    
    # Look for a second match starting from the mismatch
    results = []
    if first_length >= min_seed_length:
        results.append((first_pos, first_length))
    
    remaining_length = len(read) - first_length
    for start in range(first_length, len(read) - min_seed_length + 1):
        second_pos, second_length = binary_search_prefix(read[start:], reference, suffix_array)
        if second_length >= min_seed_length:
            results.append((second_pos, second_length))
            break
        remaining_length -= 1
        if remaining_length < min_seed_length:
            break
    
    return results

def format_alignment(pos, length):
    """Format a single alignment as (pos:length)."""
    return f"(pos:{pos:6d}, len:{length:3d})"

def calculate_coverage(alignments, read_length):
    """Calculate the percentage of the read covered by alignments."""
    covered_bases = sum(length for _, length in alignments)
    coverage_percentage = (covered_bases / read_length) * 100
    return coverage_percentage

def main(reference_file, reads_file):
    # Read the reference genome
    reference = str(next(SeqIO.parse(reference_file, "fasta")).seq)

    # Create the suffix array
    print("Building suffix array...")
    suffix_array = create_suffix_array(reference)
    print("Suffix array built.")

    # Process each read
    for record in SeqIO.parse(reads_file, "fasta"):
        read = str(record.seq)
        alignments = align_with_suffix_array(read, reference, suffix_array)
        
        # Calculate coverage
        coverage = calculate_coverage(alignments, len(read))
        
        # Format the alignments
        alignment_str = " ".join([format_alignment(pos, length) for pos, length in alignments])
        alignment_str = alignment_str.ljust(40)  # Ensure space for up to 2 alignments
        
        print(f"{record.id:<10} len={len(read):<3} {alignment_str} coverage={coverage:.2f}%")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python star.py data/reference.fasta data/reads.fasta")
        sys.exit(1)
    
    reference_file = sys.argv[1]
    reads_file = sys.argv[2]
    main(reference_file, reads_file)
