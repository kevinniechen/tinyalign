from collections import defaultdict

def find_targets(sequence, anchor_length, target_length, offset):
    pairs = []
    for i in range(len(sequence) - anchor_length - offset - target_length + 1):
        anchor = sequence[i:i+anchor_length]
        target = sequence[i+anchor_length+offset:i+anchor_length+offset+target_length]
        pairs.append((anchor, target))
    return pairs

def count_targets(samples, anchor_length, target_length, offset):
    count_tables = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for sample_name, sequences in samples.items():
        for sequence in sequences:
            pairs = find_targets(sequence, anchor_length, target_length, offset)
            for anchor, target in pairs:
                count_tables[anchor][target][sample_name] += 1
    return count_tables

def print_count_table(anchor, count_table):
    print(f"Count table for anchor: {anchor}")
    samples = sorted(set().union(*[set(counts.keys()) for counts in count_table.values()]))
    print("Target\t" + "\t".join(samples))
    for target, counts in count_table.items():
        print(f"{target}\t" + "\t".join(str(counts.get(sample, 0)) for sample in samples))

def run_splash(samples, anchor_length, target_length, offset):
    print("Counting anchor-target pairs...")
    count_tables = count_targets(samples, anchor_length, target_length, offset)
    
    print(f"\nTotal number of unique anchors: {len(count_tables)}")
    
    print("\nPrinting count tables for each anchor:")
    for anchor, count_table in count_tables.items():
        print_count_table(anchor, count_table)
        print()

    return count_tables
