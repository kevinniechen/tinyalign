import pytest
import sys
import os
import numpy as np
from scipy.stats import chi2_contingency

# Add the parent directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from splash import count_targets, calculate_pvalue

@pytest.fixture
def toy_data():
    return {
        "sample1": [
            "AAAAABBBBBCCCCCDDDDDEEEEE",
            "AAAAABBBBBDDDDDEEEEEFFFF",
            "AAAAABBBBBCCCCCDDDDDEEEEE"
        ],
        "sample2": [
            "AAAAABBBBBDDDDDEEEEEFFFF",
            "AAAAABBBBBCCCCCDDDDDEEEEE",
            "AAAAABBBBBDDDDDEEEEEFFFF"
        ]
    }

@pytest.fixture
def parameters():
    return {
        'anchor_length': 5,
        'target_length': 5,
        'offset': 5
    }

def test_count_targets(toy_data, parameters):
    count_tables = count_targets(toy_data, parameters['anchor_length'], parameters['target_length'], parameters['offset'])
    
    assert len(count_tables) == 15, "Should have 15 unique anchors"

    assert count_tables["AAAAA"] == \
           {"CCCCC": {"sample1": 2, "sample2": 1},
            "DDDDD": {"sample1": 1, "sample2": 2}}, "Incorrect count table for AAAAA"

    assert count_tables["BBBBB"] == \
           {"DDDDD": {"sample1": 2, "sample2": 1},
            "EEEEE": {"sample1": 1, "sample2": 2}}, "Incorrect count table for BBBBB"

    assert count_tables["CCCCC"] == \
           {"EEEEE": {"sample1": 2, "sample2": 1}}, "Incorrect count table for CCCCC"

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

@pytest.mark.parametrize("contingency_table, expected_relation", [
    (np.array([[10, 10], [10, 10]]), 'similar'),
    (np.array([[10, 10], [10, 100]]), 'splash_lower'),
    (np.array([[5, 15, 10], [10, 20, 30], [15, 5, 25]]), 'similar'),
    (np.array([[100, 200, 300], [150, 250, 350], [200, 300, 400]]), 'similar'),
    (np.array([[10, 10, 10], [10, 10, 1000], [10, 10, 10]]), 'splash_lower'),
])
def test_splash_vs_chi_square(contingency_table, expected_relation):
    splash_p_value = splash_test(contingency_table)
    chi_square_p_value = traditional_chi_square_test(contingency_table)
    
    assert 0 <= splash_p_value <= 1, f"SPLASH p-value {splash_p_value} is not between 0 and 1"
    assert 0 <= chi_square_p_value <= 1, f"Chi-square p-value {chi_square_p_value} is not between 0 and 1"
    
    if expected_relation == 'similar':
        assert np.isclose(splash_p_value, chi_square_p_value, rtol=0.1, atol=0.1), \
            f"Expected similar p-values, but got SPLASH: {splash_p_value}, Chi-square: {chi_square_p_value}"
    elif expected_relation == 'splash_lower':
        assert splash_p_value < chi_square_p_value, \
            f"Expected SPLASH p-value to be lower, but got SPLASH: {splash_p_value}, Chi-square: {chi_square_p_value}"

def test_calculate_pvalue():
    contingency_table = np.array([[10, 20], [30, 40]])
    sample_coeffs = np.array([1, -1])
    target_coeffs = np.array([0, 1])
    
    p_value = calculate_pvalue(contingency_table, sample_coeffs, target_coeffs)
    assert 0 <= p_value <= 1, f"Expected p-value between 0 and 1, but got {p_value}"
