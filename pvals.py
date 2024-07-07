import numpy as np
from scipy import sparse
from math import erfc

class Calculator_S:
    def __init__(self, X, f, c):
        self.X = X
        self.f = f
        self.c = c
        self.X_col_sum = X.sum(axis=0).A1
        self.M = self.X_col_sum.sum()

    def mult_fT_Xtild_c(self):
        Xtild = self.X - (self.f[:, np.newaxis] @ self.X_col_sum[np.newaxis, :]) / self.M
        return (self.f @ Xtild @ self.c)

    def get_X_col_sum(self):
        return self.X_col_sum

    def get_M(self):
        return self.M

    def set_f(self, f):
        self.f = f

    def set_c(self, c):
        self.c = c

def testPval(X, cOpt, fOpt):
    if np.all(cOpt == cOpt[0]):
        return 1.0

    min_f, max_f = np.min(fOpt), np.max(fOpt)

    if min_f == max_f:
        return 1.0

    if not (max_f <= 1 and min_f >= 0):
        raise ValueError("fOpt values must be between 0 and 1")

    calc_S = Calculator_S(X, fOpt, cOpt)

    S = calc_S.mult_fT_Xtild_c()

    denom = np.sum(cOpt**2)

    X_col_sum = calc_S.get_X_col_sum()
    denom_part2 = np.sum(cOpt * np.sqrt(X_col_sum))

    denom -= denom_part2**2 / calc_S.get_M()

    pval = 2 * np.exp(-2 * S**2 / denom)

    return min(pval, 1.0)

def effectSize_bin(X, cOpt, fOpt):
    pos_vec = cOpt > 0
    neg_vec = cOpt < 0

    if not np.any(pos_vec) or not np.any(neg_vec):
        return 0

    cPos = cOpt[pos_vec]
    cNeg = cOpt[neg_vec]

    fX = fOpt @ X
    X_cPos = X[:, pos_vec] @ cPos
    X_cNeg = X[:, neg_vec] @ cNeg
    
    fX_pos = fX[pos_vec]
    fX_neg = fX[neg_vec]
    
    return abs((fX_pos.dot(cPos) / X_cPos.sum()) - (fX_neg.dot(cNeg) / X_cNeg.sum()))

def calc_pval_base(sp_anch_contingency_table, num_rand_cf, random_state=None):
    rng = np.random.default_rng(random_state)
    min_pval = float('inf')

    for _ in range(num_rand_cf):
        randCs_row = rng.choice([-1, 1], size=sp_anch_contingency_table.shape[1])
        randFs_row = rng.choice([0, 1], size=sp_anch_contingency_table.shape[0])

        pval = testPval(sp_anch_contingency_table, randCs_row, randFs_row)
        if pval < min_pval:
            min_pval = pval

    pval_base = min_pval * num_rand_cf
    return min(pval_base, 1.0)

def compute_stats(anchor_contingency_table, num_rand_cf=1000, random_state=None):
    sp_anch_contingency_table = sparse.csc_matrix(anchor_contingency_table)

    pval_base = calc_pval_base(sp_anch_contingency_table, num_rand_cf, random_state)

    # For simplicity, we're not implementing the full alt_maximize function here
    # Instead, we'll use random c and f vectors for demonstration
    rng = np.random.default_rng(random_state)
    cOpt = rng.choice([-1, 1], size=sp_anch_contingency_table.shape[1])
    fOpt = rng.choice([0, 1], size=sp_anch_contingency_table.shape[0])

    pval_opt = testPval(sp_anch_contingency_table, cOpt, fOpt)
    effect_size_bin = effectSize_bin(sp_anch_contingency_table, cOpt, fOpt)

    return {
        'pval_base': pval_base,
        'pval_opt': pval_opt,
        'effect_size_bin': effect_size_bin
    }

# Test cases
def test_pvals():
    # Test case 1: Simple 2x2 matrix
    X1 = sparse.csc_matrix([[1, 2], [3, 4]])
    c1 = np.array([1, -1])
    f1 = np.array([0, 1])
    pval1 = testPval(X1, c1, f1)
    print(f"Test case 1 p-value: {pval1}")

    # Test case 2: Larger matrix
    X2 = sparse.csc_matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]])
    c2 = np.array([1, -1, 1])
    f2 = np.array([0, 1, 0, 1])
    pval2 = testPval(X2, c2, f2)
    print(f"Test case 2 p-value: {pval2}")

    # Test case 3: Edge case with all cOpt values the same
    X3 = sparse.csc_matrix([[1, 2], [3, 4]])
    c3 = np.array([1, 1])
    f3 = np.array([0, 1])
    pval3 = testPval(X3, c3, f3)
    print(f"Test case 3 p-value: {pval3}")

    # Test case 4: Full stats computation
    X4 = sparse.csc_matrix([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    stats = compute_stats(X4, num_rand_cf=100, random_state=42)
    print(f"Test case 4 stats: {stats}")

if __name__ == "__main__":
    test_pvals()
