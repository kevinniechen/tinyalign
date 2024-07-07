import numpy as np

def calculate_pvalue(contingency_table, sample_coeffs, target_coeffs):
    # Ensure variation in coefficients
    if np.all(sample_coeffs == sample_coeffs[0]) or np.all(target_coeffs == target_coeffs[0]):
        return 1.0

    # Compute total sample counts
    sample_totals = contingency_table.sum(axis=0).flatten()
    total_count = sample_totals.sum()

    # Calculate expected counts under the null hypothesis
    expected_counts = (target_coeffs[:, np.newaxis] @ sample_totals[np.newaxis, :]) / total_count
    count_difference = contingency_table - expected_counts

    # Calculate test statistic
    test_statistic = target_coeffs @ count_difference @ sample_coeffs

    # Compute denominator for scaling
    denominator = np.sum(sample_coeffs**2) - (np.sum(sample_coeffs * np.sqrt(sample_totals))**2) / total_count
    
    if denominator <= 0 or np.isclose(denominator, 0):
        return 1.0  # High p-value if denominator is invalid

    # Calculate Hoeffding's bound-based p-value
    S2 = test_statistic**2
    sqrt_sample_totals = np.sqrt(sample_totals)
    sum_cj_nj_sqrt = np.sum(sample_coeffs * sqrt_sample_totals)
    
    if sum_cj_nj_sqrt == 0:
        return 1.0  # High p-value if no variation in summed coefficients

    a = (1 + total_count * np.sum(sample_coeffs**2) / (sum_cj_nj_sqrt**2))**-1
    term1 = 2 * np.exp(-2 * (1 - a)**2 * S2 / np.sum(sample_coeffs**2))
    term2 = 2 * np.exp(-2 * a**2 * total_count * S2 / (sum_cj_nj_sqrt**2))
    p_value = term1 + term2

    return min(p_value, 1.0)


def calculate_effect_size(contingency_table, sample_coeffs, target_coeffs):
    positive_samples = sample_coeffs > 0
    negative_samples = sample_coeffs < 0

    if not np.any(positive_samples) or not np.any(negative_samples):
        return 0

    # Compute the sum of target coefficients for positive and negative sample groups
    positive_sum = target_coeffs @ contingency_table[:, positive_samples]
    negative_sum = target_coeffs @ contingency_table[:, negative_samples]
    
    # Calculate mean values for positive and negative sample groups
    positive_mean = np.sum(positive_sum * sample_coeffs[positive_samples]) / positive_sum.sum()
    negative_mean = np.sum(negative_sum * sample_coeffs[negative_samples]) / negative_sum.sum()

    return abs(positive_mean - negative_mean)

def compute_statistics(contingency_table, num_random_trials=1000, random_seed=None):
    rng = np.random.default_rng(random_seed)
    
    min_p_value = float('inf')
    optimal_sample_coeffs = None
    optimal_target_coeffs = None

    for _ in range(num_random_trials):
        # Randomly choose sample and target coefficients
        random_sample_coeffs = rng.choice([-1, 1], size=contingency_table.shape[1])
        random_target_coeffs = rng.choice([0, 1], size=contingency_table.shape[0])
        
        # Compute p-value for the random coefficients
        p_value = calculate_pvalue(contingency_table, random_sample_coeffs, random_target_coeffs)
        if p_value < min_p_value:
            min_p_value = p_value
            optimal_sample_coeffs = random_sample_coeffs
            optimal_target_coeffs = random_target_coeffs

    # Apply Bonferroni correction for multiple trials
    baseline_p_value = min(min_p_value * num_random_trials, 1.0)
    optimal_p_value = min_p_value
    
    # Calculate effect size for the optimal coefficients
    effect_size = calculate_effect_size(contingency_table, optimal_sample_coeffs, optimal_target_coeffs)

    return {
        'baseline_p_value': baseline_p_value,
        'optimal_p_value': optimal_p_value,
        'effect_size': effect_size
    }
