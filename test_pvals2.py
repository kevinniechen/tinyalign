import numpy as np
from scipy.stats import chi2_contingency

def calculate_pvalue(contingency_table, sample_coeffs, target_coeffs):
    if np.all(sample_coeffs == sample_coeffs[0]) or np.all(target_coeffs == target_coeffs[0]):
        return 1.0

    sample_totals = contingency_table.sum(axis=0).flatten()
    total_count = sample_totals.sum()

    expected_counts = (target_coeffs[:, np.newaxis] @ sample_totals[np.newaxis, :]) / total_count
    count_difference = contingency_table - expected_counts

    test_statistic = target_coeffs @ count_difference @ sample_coeffs

    denominator = np.sum(sample_coeffs**2) - (np.sum(sample_coeffs * np.sqrt(sample_totals))**2) / total_count
    
    if denominator <= 0 or np.isclose(denominator, 0):
        return 1.0

    S2 = test_statistic**2
    sqrt_sample_totals = np.sqrt(sample_totals)
    sum_cj_nj_sqrt = np.sum(sample_coeffs * sqrt_sample_totals)
    
    if sum_cj_nj_sqrt == 0:
        return 1.0

    a = (1 + total_count * np.sum(sample_coeffs**2) / (sum_cj_nj_sqrt**2))**-1
    term1 = 2 * np.exp(-2 * (1 - a)**2 * S2 / np.sum(sample_coeffs**2))
    term2 = 2 * np.exp(-2 * a**2 * total_count * S2 / (sum_cj_nj_sqrt**2))
    p_value = term1 + term2

    return min(p_value, 1.0)

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
    chi2, p, dof, expected = chi2_contingency(contingency_table)
    return p

# Test cases
def test_cases():
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

# Run test cases
test_cases()
