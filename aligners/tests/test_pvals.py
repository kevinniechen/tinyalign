import numpy as np
from pvals import calculate_pvalue

def splash_test(contingency_table, num_random_trials=1000, random_seed=None):
    rng = np.random.default_rng(random_seed)
    min_p_value = float('inf')

    for _ in range(num_random_trials):
        random_sample_coeffs = rng.choice([-1, 1], size=contingency_table.shape[1])
        random_target_coeffs = rng.choice([0, 1], size=contingency_table.shape[0])
        
        p_value = calculate_pvalue(contingency_table, random_sample_coeffs, random_target_coeffs)
        if p_value < min_p_value:
            min_p_value = p_value

    baseline_p_value = min(min_p_value * num_random_trials, 1.0)
    return baseline_p_value

def traditional_chi_square_test(contingency_table):
    from scipy.stats import chi2_contingency
    chi2, p, dof, expected = chi2_contingency(contingency_table)
    return p

def test_splash_vs_chi_square():
    # Test case 1: Simple case with small counts
    contingency_table1 = np.array([
        [10, 10],
        [10, 10]
    ])
    print("Test Case 1")
    print("SPLASH p-value:", splash_test(contingency_table1))
    print("Chi-square p-value:", traditional_chi_square_test(contingency_table1))
    print()
    
    # Test case 2: Case with one sample having an outlier count
    contingency_table2 = np.array([
        [10, 10],
        [10, 100]
    ])
    print("Test Case 2")
    print("SPLASH p-value:", splash_test(contingency_table2))
    print("Chi-square p-value:", traditional_chi_square_test(contingency_table2))
    print()
    
    # Test case 3: Case with more variability across samples
    contingency_table3 = np.array([
        [5, 15, 10],
        [10, 20, 30],
        [15, 5, 25]
    ])
    print("Test Case 3")
    print("SPLASH p-value:", splash_test(contingency_table3))
    print("Chi-square p-value:", traditional_chi_square_test(contingency_table3))
    print()
    
    # Test case 4: Large counts with sample-dependent variations
    contingency_table4 = np.array([
        [100, 200, 300],
        [150, 250, 350],
        [200, 300, 400]
    ])
    print("Test Case 4")
    print("SPLASH p-value:", splash_test(contingency_table4))
    print("Chi-square p-value:", traditional_chi_square_test(contingency_table4))
    print()
    
    # Test case 5: Extreme outliers in one sample
    contingency_table5 = np.array([
        [10, 10, 10],
        [10, 10, 1000],
        [10, 10, 10]
    ])
    print("Test Case 5")
    print("SPLASH p-value:", splash_test(contingency_table5))
    print("Chi-square p-value:", traditional_chi_square_test(contingency_table5))
    print()

if __name__ == "__main__":
    test_splash_vs_chi_square()
