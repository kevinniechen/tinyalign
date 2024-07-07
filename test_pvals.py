import numpy as np
from scipy import sparse
from pvals import calculate_pvalue, compute_statistics

def test_pvals():
    # Test case 1: Simple 2x2 matrix
    X1 = sparse.csc_matrix([[1, 2], [3, 4]])
    c1 = np.array([1, -1])
    f1 = np.array([0, 1])
    pval1 = calculate_pvalue(X1, c1, f1)
    print(f"Test case 1 p-value: {pval1}")
    assert 0 <= pval1 <= 1, "P-value should be between 0 and 1"

    # Test case 2: Larger matrix
    X2 = sparse.csc_matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
    c2 = np.array([1, -1, 1])
    f2 = np.array([0, 1, 0, 1])
    pval2 = calculate_pvalue(X2, c2, f2)
    print(f"Test case 2 p-value: {pval2}")
    assert 0 <= pval2 <= 1, "P-value should be between 0 and 1"

    # Test case 3: Edge case with all cOpt values the same
    X3 = sparse.csc_matrix([[1, 2], [3, 4]])
    c3 = np.array([1, 1])
    f3 = np.array([0, 1])
    pval3 = calculate_pvalue(X3, c3, f3)
    print(f"Test case 3 p-value: {pval3}")
    assert pval3 == 1, "P-value should be 1 when all cOpt values are the same"

    # Test case 4: Full stats computation
    X4 = sparse.csc_matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    stats = compute_statistics(X4, num_random_trials=100, random_seed=42)
    print(f"Test case 4 stats: {stats}")
    assert 0 <= stats['baseline_p_value'] <= 1, "Base p-value should be between 0 and 1"
    assert 0 <= stats['optimal_p_value'] <= 1, "Optimal p-value should be between 0 and 1"
    assert 0 <= stats['effect_size'] <= 1, "Effect size should be between 0 and 1"

    # Test case 5: Perfectly correlated data
    X5 = sparse.csc_matrix([[1, 0], [0, 1]])
    c5 = np.array([1, -1])
    f5 = np.array([1, 0])
    pval5 = calculate_pvalue(X5, c5, f5)
    print(f"Test case 5 p-value: {pval5}")
    assert pval5 < 0.05, "P-value should be small for perfectly correlated data"

    # Test case 6: Uncorrelated data
    X6 = sparse.csc_matrix([[1, 1], [1, 1]])
    c6 = np.array([1, -1])
    f6 = np.array([0.5, 0.5])
    pval6 = calculate_pvalue(X6, c6, f6)
    print(f"Test case 6 p-value: {pval6}")
    assert pval6 > 0.05, "P-value should be large for uncorrelated data"

    print("All test cases passed successfully!")

if __name__ == "__main__":
    test_pvals()
