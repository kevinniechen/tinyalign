import numpy as np
from scipy import sparse
from scipy.stats import norm

def calculate_pvalue(contingency_table, sample_coeffs, target_coeffs):
    # Check for no variation in coefficients
    if np.all(sample_coeffs == sample_coeffs[0]) or np.all(target_coeffs == target_coeffs[0]):
        return 1.0

    sample_totals = contingency_table.sum(axis=0).A1
    total_count = sample_totals.sum()

    # Calculate expected counts and difference
    expected_counts = (target_coeffs[:, np.newaxis] @ sample_totals[np.newaxis, :]) / total_count
    count_difference = contingency_table - expected_counts

    # Calculate test statistic
    test_statistic = target_coeffs @ count_difference @ sample_coeffs

    # Calculate denominator
    denominator = np.sum(sample_coeffs**2) - (np.sum(sample_coeffs * np.sqrt(sample_totals))**2) / total_count
    
    # Handle perfect correlation
    if denominator <= 0:
        return 0.0
    
    # Calculate z-score and p-value
    z_score = test_statistic / np.sqrt(denominator / 2)
    p_value = 2 * (1 - norm.cdf(abs(z_score)))

    return p_value

def calculate_effect_size(contingency_table, sample_coeffs, target_coeffs):
    positive_samples = sample_coeffs > 0
    negative_samples = sample_coeffs < 0

    if not np.any(positive_samples) or not np.any(negative_samples):
        return 0

    target_sample_product = target_coeffs @ contingency_table
    positive_sum = contingency_table[:, positive_samples] @ sample_coeffs[positive_samples]
    negative_sum = contingency_table[:, negative_samples] @ sample_coeffs[negative_samples]
    
    positive_effect = target_sample_product[positive_samples].dot(sample_coeffs[positive_samples]) / positive_sum.sum()
    negative_effect = target_sample_product[negative_samples].dot(sample_coeffs[negative_samples]) / negative_sum.sum()

    return abs(positive_effect - negative_effect)

def compute_statistics(contingency_table, num_random_trials=1000, random_seed=None):
    rng = np.random.default_rng(random_seed)
    
    # Calculate baseline p-value
    min_p_value = float('inf')
    for _ in range(num_random_trials):
        random_sample_coeffs = rng.choice([-1, 1], size=contingency_table.shape[1])
        random_target_coeffs = rng.choice([0, 1], size=contingency_table.shape[0])
        p_value = calculate_pvalue(contingency_table, random_sample_coeffs, random_target_coeffs)
        if p_value < min_p_value:
            min_p_value = p_value

    baseline_p_value = min(min_p_value * num_random_trials, 1.0)
    
    # Calculate optimal p-value and effect size
    optimal_sample_coeffs = rng.choice([-1, 1], size=contingency_table.shape[1])
    optimal_target_coeffs = rng.choice([0, 1], size=contingency_table.shape[0])
    optimal_p_value = calculate_pvalue(contingency_table, optimal_sample_coeffs, optimal_target_coeffs)
    effect_size = calculate_effect_size(contingency_table, optimal_sample_coeffs, optimal_target_coeffs)

    return {
        'baseline_p_value': baseline_p_value,
        'optimal_p_value': optimal_p_value,
        'effect_size': effect_size
    }
