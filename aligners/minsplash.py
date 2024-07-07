import numpy as np
from collections import defaultdict

class MinSPLASH:
    def __init__(self):
        self.anchors = defaultdict(lambda: defaultdict(int))
        self.targets = defaultdict(lambda: defaultdict(int))
        
    def process_sample(self, sample, sample_id):
        """Process a single RNA-seq sample."""
        for read in sample:
            anchor, target = self._extract_anchor_target(read)
            self.anchors[anchor][sample_id] += 1
            self.targets[target][sample_id] += 1
    
    def _extract_anchor_target(self, read):
        """Extract anchor and target from a read."""
        # Simplified version: assume first 20 bp is anchor, rest is target
        return read[:20], read[20:]
    
    def generate_count_table(self):
        """Generate count table for anchors and targets across samples."""
        anchor_table = self._generate_table(self.anchors)
        target_table = self._generate_table(self.targets)
        return anchor_table, target_table
    
    def _generate_table(self, data):
        """Generate count table from anchor or target data."""
        items = list(data.keys())
        samples = set()
        for counts in data.values():
            samples.update(counts.keys())
        samples = sorted(samples)
        
        table = np.zeros((len(items), len(samples)), dtype=int)
        for i, item in enumerate(items):
            for j, sample in enumerate(samples):
                table[i, j] = data[item][sample]
        
        return table, items, samples
    
    def perform_pvalue_test(self):
        """Perform p-value test on targets."""
        # TODO: Implement p-value test
        pass

# Test function
def test_minsplash():
    splash = MinSPLASH()
    
    # Simulate some sample data
    samples = [
        ["ATCGATCGATCGATCGATCGGGTTACGTACGTA", "ATCGATCGATCGATCGATCGTTACGTACGTAC"],
        ["TGCATGCATGCATGCATGCAAAAAAAAAAAAA", "TGCATGCATGCATGCATGCATTTTTTTTTTTT"]
    ]
    
    # Process samples
    for i, sample in enumerate(samples):
        splash.process_sample(sample, f"sample_{i}")
    
    # Generate count tables
    anchor_table, target_table = splash.generate_count_table()
    
    print("Anchor table shape:", anchor_table[0].shape)
    print("Target table shape:", target_table[0].shape)
    print("Anchors:", anchor_table[1])
    print("Targets:", target_table[1])
    print("Samples:", anchor_table[2])

if __name__ == "__main__":
    test_minsplash()
