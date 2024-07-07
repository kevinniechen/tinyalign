import sys
from Bio import SeqIO

GAP_THRESHOLD = 4000
MAX_DIFFERENCES = 5

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
        if pattern[: match_length + 1] <= suffix[: match_length + 1]:
            right = mid - 1
        else:
            left = mid + 1

    return best_match_position, best_match_length


def align_with_suffix_array(read, reference, suffix_array, min_seed_length=10):
    """Align the read using binary search on the suffix array for longest prefix matches."""
    results = []
    start = 0
    while start < len(read):
        pos, length = binary_search_prefix(read[start:], reference, suffix_array)
        if length >= min_seed_length:
            results.append((pos, length, start))
            start += length
        else:
            start += 1
        if start + min_seed_length > len(read):
            break
    return results


def stitch_seeds(read, reference, seeds):
    if len(seeds) < 2:
        return seeds

    # Step 1: Sort seeds by their read positions
    seeds = sorted(seeds, key=lambda x: x[1])

    # Function to check if two seeds can be stitched
    def can_stitch(seed1, seed2):
        pos1, len1, start1 = seed1
        pos2, len2, start2 = seed2
        alignment_gap = abs((start2-pos2) - (start1-pos1))
        if alignment_gap > GAP_THRESHOLD:
            return False

        differences = 0
        min_length = min(pos2 - pos1, len2 - len1, len(read) - pos1, len(reference) - len1)
        for i in range(min_length):
            if read[pos1 + i] != reference[len1 + i]:
                differences += 1
            if differences > MAX_DIFFERENCES:
                return False
        return True

    # Step 2: Stitch seeds
    stitched_seeds = [seeds[0]]
    for seed in seeds[1:]:
        if can_stitch(stitched_seeds[-1], seed):
            stitched_seeds.append(seed)

    return stitched_seeds


def calculate_coverage(seeds, read_length):
    """Calculate the percentage of the read covered by seeds."""
    covered_bases = sum(length for _, length, _ in seeds)
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
        seeds = align_with_suffix_array(read, reference, suffix_array)
        stitched_seeds = stitch_seeds(read, reference, seeds)
        coverage = calculate_coverage(stitched_seeds, len(read))
        print(f"{record.id:<10} {coverage:>6.2f}% {len(stitched_seeds):<12} {stitched_seeds} @@@ {seeds}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python star.py <reference_file> <reads_file>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
