import numpy as np

def union_(breaks, indices):
    """`breaks` is sorted in ascending order. `indices` indexes the intervals
    between consecutive points in `breaks`.
    """
    beg, end = [breaks[indices[0]]], [breaks[indices[0] + 1]]
    for j in indices[1:]:
        if end[-1] < breaks[j]:
            beg.append(breaks[j])
            end.append(breaks[j + 1])
        else:
            end[-1] = max(end[-1], breaks[j + 1])
    return np.array([beg, end]).T

def intersect_(*collections):
    intervals_ = list()
    for I in zip(*collections):
        breaks_ = np.unique(np.vstack(I))
        ints_ = [j for j, (z0, z1) in enumerate(zip(breaks_[:-1], breaks_[1:]))
                 if all(any(i0 <= z0 <= z1 <= i1 for i0, i1 in i) for i in I)]
        intervals_.append(union_(breaks_, ints_))
    return intervals_

def confidence_region(knots, frequency, levels):
## Consolidate intervals
    intervals_ = list()
    for level in levels:
        ints_ = np.flatnonzero(frequency > level)
        intervals_.append(union_(knots, ints_))
    return intervals_

def sided_CCR(A_i, B_i, A_n, B_n):
    N = A_i.shape[0]
## Compute the knots
    dA, dB = A_i - A_n, B_n - B_i
## Handle degenerate intervals
    mp_inf = sum((dB == 0) & (dA >= 0))
## Now process the rays
    dA = dA[dB != 0] ; dB = dB[dB != 0]
    ratio_ = dA / dB
    breaks = np.r_[-np.inf, np.unique(ratio_), np.inf]
    coverage = np.full(len(breaks) - 1, mp_inf, dtype=np.float)
##  Leftward rays: (-\infty, \beta]
    last_ = np.searchsorted(breaks, ratio_[dB > 0], side="right") - 1
    for r_ in last_: coverage[:r_] += 1
##  Rightward rays: [\alpha, +\infty)
    first_ = np.searchsorted(breaks, ratio_[dB < 0], side="left")
    for l_ in first_: coverage[l_:] += 1
## compute the frequency
    coverage /= N
    return breaks, coverage

def absolute_CCR(A_i, B_i, A_n, B_n):
    N = A_i.shape[0]
    if B_n == 0:
## Count the trivial cases
        A_n = np.abs(A_n)
        mp_inf = sum((B_i == 0) & (np.abs(A_i) >= A_n))
        A_i = A_i[B_i > 0] ; B_i = B_i[B_i > 0]
        R1 = -(A_n + A_i) / B_i
        R2 = (A_n - A_i) / B_i
        breaks = np.unique(np.r_[-np.inf, R1, R2, np.inf])
        coverage = np.full(len(breaks) - 1, mp_inf, dtype=np.float)
##  Leftward rays: (-\infty, \beta]
        last_ = np.searchsorted(breaks, R1, side="right") - 1
        for r_ in last_: coverage[:r_] += 1
##  Rightward rays: [\alpha, +\infty)
        first_ = np.searchsorted(breaks, R2, side="left")
        for l_ in first_: coverage[l_:] += 1
    elif B_n > 0:
## Handle the trivial cases
        trivial_ = (B_i >= B_n) & (A_i * B_n == A_n * B_i)
        mp_inf = sum(trivial_)
        A_i = A_i[~trivial_] ; B_i = B_i[~trivial_]
## Precompute
        dA, dB = A_i - A_n, B_n - B_i
        sA, sB = A_i + A_n, B_i + B_n
        AiBn_c = A_i * B_n - A_n * B_i
## Now, create breakpoints
        R1 = dA.copy()
        R1[dB != 0] /= dB[dB != 0]
        R2 = - sA / sB
        breaks = np.unique(np.r_[-np.inf, R1, R2, - A_i[B_i > 0] / B_i[B_i > 0], np.inf])
        coverage = np.full(len(breaks) - 1, mp_inf, dtype=np.float)
##  Rightward rays: [\alpha, +\infty)
        first_ = np.searchsorted(breaks, R2[(dB <= 0) & (AiBn_c > 0)], side="left")
        for l_ in first_: coverage[l_:] += 1
##  Rightward rays: [\alpha, +\infty)
        first_ = np.searchsorted(breaks, R1[(dB < 0) & (AiBn_c > 0)], side="left")
        for l_ in first_: coverage[l_:] += 1
##  Leftward rays: (-\infty, \beta]
        last_ = np.searchsorted(breaks, R2[(dB <= 0) & (AiBn_c < 0)], side="right") - 1
        for r_ in last_: coverage[:r_] += 1
##  Leftward rays: (-\infty, \beta]
        last_ = np.searchsorted(breaks, R1[(dB < 0) & (AiBn_c < 0)], side="right") - 1
        for r_ in last_: coverage[:r_] += 1
##  Type9 intervals: [\alpha, \beta]
        mask_ = (dB > 0) & (AiBn_c >= 0)
        first_ = np.searchsorted(breaks, R2[mask_], side="left")
        last_ = np.searchsorted(breaks, R1[mask_], side="right") - 1
        for l_, r_ in zip(first_, last_): coverage[l_:r_] += 1
##  Type10 intervals: [\beta, \alpha]
        mask_ = (dB > 0) & (AiBn_c < 0)
        first_ = np.searchsorted(breaks, R1[mask_], side="left")
        last_ = np.searchsorted(breaks, R2[mask_], side="right") - 1
        for l_, r_ in zip(first_, last_): coverage[l_:r_] += 1
    else:
        raise ValueError("""`B_n` cannot be negative.""")
## compute the frequency
    return breaks, coverage / N

def RRCM(A, B, levels):
    A, B = A.copy(), B.copy()
    A[B < 0] *= -1
    B[B < 0] *= -1
    N, A_n, B_n = A.shape[0], A[-1], B[-1]
    
    ## compute P and Q
    P, Q = - (A + A_n) / (B + B_n), (A - A_n) / (B_n - B)

    breaks_ = np.unique(np.r_[P, Q, -np.inf, np.inf])
    breaks_ = breaks_[~np.isnan(breaks_)]

    ## locate the associated endpoints
    Pi_, Qi_ = np.searchsorted(breaks_, P), np.searchsorted(breaks_, Q)
    delta_ = np.zeros(breaks_.shape[0] + 1, dtype=np.float)
    for p, q, a_, b_ in zip(Pi_, Qi_, A, B):
        if b_ == B_n:
            if abs(a_) >= abs(A_n) and B_n == 0:
                delta_[0] += 1
                delta_[-1] -= 1
            elif B_n > 0:
                s, t = 0, -1
                if a_ > A_n:
                    s = p
                elif a_ > A_n:
                    t = p + 1
                delta_[s] += 1
                delta_[t] -= 1
        else:
            if a_ * B_n >= b_ * A_n:
                delta_[q + 1] -= 1
                delta_[p] += 1
            else:
                delta_[p + 1] -= 1
                delta_[q] += 1
    ## 
            if b_ > B_n:
                delta_[0] += 1
                delta_[-1] -= 1
    ## compute the frequency
    coverage_ = np.cumsum(delta_) / N
    ## Consolidate intervals
    intervals_ = list()
    for level in levels:
        indices_ = np.flatnonzero(coverage_ > level)
        beg, end = [indices_[0]], [indices_[0]+1]
        for j in indices_:
            if end[-1] < j:
                end[-1] -= 1
                beg.append(j)
                end.append(j + 1)
            else:
                end[-1] = max(end[-1], j + 1)
        end[-1] -= 1
        intervals_.append(breaks_[np.array([beg, end]).T])
    return intervals_

def CCR(A, B, levels):
    hi = confidence_region(*sided_CCR(A, B, A[-1:], B[-1:]), levels=levels / 2)
    lo = confidence_region(*sided_CCR(-A, -B, -A[-1:], -B[-1:]), levels=levels / 2)
    return intersect_(lo, hi)

def CRR(A, B, levels):
    A, B = A.copy(), B.copy()
    N, A_n, B_n = A.shape[0], A[-1], B[-1]
    ## compute P and Q
    Q = (A - A_n) / (B_n - B)

    breaks_ = np.unique(np.r_[Q, -np.inf, np.inf])
    breaks_ = breaks_[~np.isnan(breaks_)]

    ## locate the associated endpoints
    Qi_ = np.searchsorted(breaks_, Q)
    d_U_ = np.zeros(breaks_.shape[0] + 1, dtype=np.float)
    d_L_ = np.zeros(breaks_.shape[0] + 1, dtype=np.float)
    for q, a_, b_ in zip(Qi_, A, B):
        if b_ == B_n:
            if a_ >= A_n:
                d_U_[0] += 1
                d_U_[-1] -= 1
            if a_ <= A_n:
                d_L_[0] += 1
                d_L_[-1] -= 1
        elif b_ > B_n:
            d_L_[0] += 1
            d_L_[q+1] -= 1
            d_U_[q] += 1
            d_U_[-1] -= 1
        else:
            d_U_[0] += 1
            d_U_[q+1] -= 1
            d_L_[q] += 1
            d_L_[-1] -= 1

    ## compute the frequency
    coverage_ = np.minimum(np.cumsum(d_U_), np.cumsum(d_L_)) / N

    ## Consolidate intervals
    intervals_ = list()
    for level in levels:
        indices_ = np.flatnonzero(coverage_ > level / 2)
        beg, end = [indices_[0]], [indices_[0]+1]
        for j in indices_:
            if end[-1] < j:
                end[-1] -= 1
                beg.append(j)
                end.append(j + 1)
            else:
                end[-1] = max(end[-1], j + 1)
        end[-1] -= 1
        intervals_.append(breaks_[np.array([beg, end]).T])
    return intervals_
