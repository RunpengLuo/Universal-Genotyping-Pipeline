"""Call germline heterozygous sites from tumor-only allele counts via a latent LOH-state chain.

Last update: 2026-08-16

Observed, per site i:
    - ``k_i`` = min(REF_COUNT, ALT_COUNT)
    - ``n_i`` = REF_COUNT + ALT_COUNT
    - ``d_i``, bp to the previous site on the same arm

Latent, per site i:
    - ``s_i`` in 0..K-1, the segment state, a Markov chain along the arm, where
      K = 1 + ``n_retained``. State 0 is clonal LOH. Returned as ``p_loh`` and ``state``.
    - ``g_i`` in {het, hom}, the germline genotype. Marginalised out inside the emission;
      returned as ``p_het``.

Transition:
    - the chain restarts at every chromosome arm, with ``s_1`` uniform over the K states
    - breakpoints are Poisson at rate ``t`` per BASE PAIR, so

          P(s_i = s_i-1 | d_i) = exp(-t * d_i)

      and the rest of the mass is uniform over the other K-1 states

Emission:
    - ``b_s``, held in ``means[s]``, is the minor fraction a germline het shows in state
      ``s``; a germline hom shows ``b_hom`` in every state
    - ``g_i ~ Bernoulli(pi)``
    - ``b_hom`` is a point mass at ``hom_laf``
    - ``b_0`` is a point mass at ``loh_laf``, roughly (1 - purity) / (2 - purity);
      ``b_1..b_K-1`` are flat on (b_0 + ``margin``, 0.5]
    - each site draws its own rate from ``Beta(tau*b, tau*(1-b))``, integrated out to give

          BB(k; n, b) = BetaBinom(k; n, tau*b, tau*(1-b))

      where ``tau`` = inf leaves the binomial

Likelihood, over one arm:
    - Joint likelihood

          P(k, s, g | n, pos) = P(s_1) prod_i P(s_i | s_i-1, d_i) prod_i P(g_i) P(k_i | n_i, s_i, g_i)

    - Conditional likelihood, summed over both phases since ``k`` is a minimum

          P(k_i | n_i, s_i, g_i) = BB(k_i; n_i, b) + BB(k_i; n_i, 1 - b)
                                   b = b_s_i for g_i = het, b_hom for g_i = hom

    - Marginal likelihood, what the chain sees

          P(k_i | n_i, s_i) = pi [BB(k_i; n_i, b_s_i) + BB(k_i; n_i, 1 - b_s_i)]
                              + (1 - pi) [BB(k_i; n_i, b_hom) + BB(k_i; n_i, 1 - b_hom)]

Learnable parameters, fitted by EM:
    - ``b_1..b_K-1``, over at most ``n_iter`` iterations, stopping below ``em_tol``
    - ``pi``, when ``learn_pi``

Returned per site:
    - ``p_loh`` = P(s_i = 0 | data)
    - ``p_het`` = P(g_i = het | data), marginalised over s_i
    - ``state`` = the Viterbi MAP s_i

Calling rule, applied after decoding:
    - below ``min_dp`` depth the site is ./.
    - a site is 0/1 unless BOTH ``p_het`` < ``p_het_min`` and ``p_loh`` < ``loh_min``

Notes:
    - derived from hetdetect's ``hmm_decode``
    - emission is a beta-binomial on counts, not a Gaussian on LAF, so depth enters it
    - transitions are per bp, not per site, so SNP density does not set segment length
    - the clonal-LOH state is explicit and pinned at ``b_0``
    - the genotype test sits inside the emission, not in a post-hoc filter
    - ``b_hom`` and ``b_0`` stay separate: sequencing error against normal admixture
"""

import logging

import numpy as np
from numba import njit
from scipy.special import betaln, gammaln, logsumexp

GT_CODES = ("0/0", "0/1", "1/1", "./.")
GT_DTYPE = "<U3"


def log_binom_coef(k, n):
    """``log C(n, k)``, the state-independent part of every ``BB``.

    Computed once per fit and passed down as *log_coef*.

    Args:
        k: Minor counts.
        n: Depths, aligned to ``k``.

    Returns:
        float64 array.
    """
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def _bb_logpmf(k, n, b, tau, log_coef=None, eps=1e-16):
    """``log BB(k; n, b)``.

    Args:
        k: Minor counts.
        n: Depths, aligned to ``k``.
        b: Scalar allele fraction.
        tau: Concentration; ``np.inf`` is the binomial.
        log_coef: Precomputed :func:`log_binom_coef`.
        eps: ``b`` is clipped to ``[eps, 1 - eps]`` on the binomial branch.

    Returns:
        float64 array of log probabilities.
    """
    if log_coef is None:
        log_coef = log_binom_coef(k, n)
    if not np.isfinite(tau):
        p = min(max(b, eps), 1.0 - eps)
        return log_coef + k * np.log(p) + (n - k) * np.log1p(-p)
    a1, a2 = b * tau, (1.0 - b) * tau
    return log_coef + betaln(k + a1, n - k + a2) - betaln(a1, a2)


def bb_logpmf(k, n, b, tau=np.inf, log_coef=None):
    """``log [BB(k; n, b) + BB(k; n, 1 - b)]``.

    Args:
        k: Minor counts.
        n: Depths, aligned to ``k``.
        b: Scalar allele fraction.
        tau: Concentration; ``np.inf`` is the binomial.
        log_coef: Precomputed :func:`log_binom_coef`.

    Returns:
        float64 array of log probabilities.
    """
    if log_coef is None:
        log_coef = log_binom_coef(k, n)
    return np.logaddexp(
        _bb_logpmf(k, n, b, tau, log_coef), _bb_logpmf(k, n, 1.0 - b, tau, log_coef)
    )


def phase_weight(k, n, b, tau=np.inf, log_coef=None):
    """Posterior that ``k`` came from the ``b`` phase rather than ``1-b``.

    Needed only by the M-step: estimating ``b`` from ``k/n`` is biased low because the min of
    a symmetric pair sits below the mean.

    Args:
        k: Minor counts.
        n: Depths.
        b: Scalar allele fraction.
        tau: Concentration; ``np.inf`` is the binomial.
        log_coef: Precomputed :func:`log_binom_coef`.

    Returns:
        float64 array in [0, 1].
    """
    if log_coef is None:
        log_coef = log_binom_coef(k, n)
    a = _bb_logpmf(k, n, b, tau, log_coef)
    c = _bb_logpmf(k, n, 1.0 - b, tau, log_coef)
    return np.exp(a - np.logaddexp(a, c))


def _pi_array(pi, shape, eps=1e-9):
    """Broadcast a scalar or per-site het prior to *shape*, clipped to ``[eps, 1 - eps]``."""
    return np.clip(
        np.broadcast_to(np.asarray(pi, dtype="float64"), shape), eps, 1 - eps
    )


def emission_logp(k, n, means, hom_laf, pi, tau=np.inf, log_coef=None):
    """``log P(k_i | n_i, s_i)``, the marginal likelihood, per site and state.

    Args:
        k: Minor counts.
        n: Depths, aligned to ``k``.
        means: ``b_s`` per state.
        hom_laf: ``b_hom``.
        pi: Het prior; scalar or per-site.
        tau: Concentration; ``np.inf`` is the binomial.
        log_coef: Precomputed :func:`log_binom_coef`.

    Returns:
        ``(n_sites, len(means))`` float64 array.
    """
    if log_coef is None:
        log_coef = log_binom_coef(k, n)
    pi_a = _pi_array(pi, np.shape(k))
    log_pi = np.log(pi_a)
    log_1mpi = np.log1p(-pi_a)
    hom = log_1mpi + bb_logpmf(k, n, hom_laf, tau, log_coef)
    out = np.empty((k.size, len(means)))
    for j, b in enumerate(means):
        np.logaddexp(log_pi + bb_logpmf(k, n, b, tau, log_coef), hom, out=out[:, j])
    return out


def het_posterior(k, n, means, hom_laf, pi, tau=np.inf, log_coef=None):
    """``P(g_i = het | k_i, n_i, s_i)`` per site and state.

    Args:
        k: Minor counts.
        n: Depths.
        means: ``b_s`` per state.
        hom_laf: ``b_hom``.
        pi: Het prior; scalar or per-site.
        tau: Concentration; ``np.inf`` is the binomial.
        log_coef: Precomputed :func:`log_binom_coef`.

    Returns:
        ``(n_sites, len(means))`` float64 array.
    """
    if log_coef is None:
        log_coef = log_binom_coef(k, n)
    pi_a = _pi_array(pi, np.shape(k))
    log_pi = np.log(pi_a)
    hom = np.log1p(-pi_a) + bb_logpmf(k, n, hom_laf, tau, log_coef)
    out = np.empty((k.size, len(means)))
    for j, b in enumerate(means):
        het = bb_logpmf(k, n, b, tau, log_coef)
        out[:, j] = pi_a * np.exp(het - np.logaddexp(log_pi + het, hom))
    return out


def arm_batches(group_id, pos):
    """Sort sites into their arms and pad the arms into one rectangular index.

    Args:
        group_id: Integer arm per site.
        pos: Site positions.

    Returns:
        Tuple of ``idx`` ``(n_arms, longest_arm)`` int64, padded with ``n_sites``, and ``mask``
        marking its real entries.
    """
    order = np.lexsort((pos, group_id))
    bounds = np.flatnonzero(np.diff(group_id[order]) != 0) + 1
    segs = [s for s in np.split(order, bounds) if s.size]
    n = pos.size
    length = max(s.size for s in segs)
    idx = np.full((len(segs), length), n, dtype="int64")
    mask = np.zeros((len(segs), length), dtype=bool)
    for a, s in enumerate(segs):
        idx[a, : s.size] = s
        mask[a, : s.size] = True
    return idx, mask


def transition_logs(pos, idx, mask, t, n_states, floor=1e-300):
    """Stay/switch log probabilities for the step into each site.

    Args:
        pos: Site positions.
        idx: Arm index from :func:`arm_batches`.
        mask: Its real entries.
        t: Breakpoint rate per bp.
        n_states: K.
        floor: Lower bound on the switch probability before taking its log.

    Returns:
        Tuple of ``log_stay`` and ``log_switch``, each length ``n_sites + 1``.
    """
    n = pos.size
    p = np.zeros(n + 1)
    prev = np.full(idx.shape, n, dtype="int64")
    prev[:, 1:] = idx[:, :-1]
    step = mask.copy()
    step[:, 0] = False
    cur_i = idx[step]
    prev_i = prev[step]
    d = np.abs(pos[cur_i].astype("int64") - pos[prev_i].astype("int64"))
    p[cur_i] = -np.expm1(-t * d)
    p = np.minimum(p, (n_states - 1) / n_states)
    log_stay = np.log1p(-p)
    log_switch = np.log(np.maximum(p / (n_states - 1), floor))
    return log_stay, log_switch


@njit(cache=True)
def _fb_kernel(e, log_diff, log_switch, idx, mask, n_sites):
    """Forward-backward recursions over the padded arm batch.

    Args:
        e: ``(n_sites + 1, K)`` log emissions; padded positions index the trailing row.
        log_diff: ``log(stay - switch)`` per site, length ``n_sites + 1``.
        log_switch: ``log(switch)`` per site, length ``n_sites + 1``.
        idx: ``(n_arms, longest_arm)`` arm index, padded with ``n_sites``.
        mask: Its real entries.
        n_sites: Number of real sites.

    Returns:
        Tuple of ``alpha``, ``beta`` (both ``(n_sites + 1, K)``) and the log likelihood.
    """
    A, L = idx.shape
    K = e.shape[1]
    alpha = np.zeros((n_sites + 1, K))
    beta = np.zeros((n_sites + 1, K))
    cur = np.empty((A, K))
    nxt = np.empty((A, K))
    log_k = np.log(K)
    loglik = 0.0

    for a in range(A):
        r = idx[a, 0]
        mx = -np.inf
        for j in range(K):
            v = e[r, j] - log_k
            cur[a, j] = v
            if v > mx:
                mx = v
        acc = 0.0
        for j in range(K):
            acc += np.exp(cur[a, j] - mx)
        norm = mx + np.log(acc)
        for j in range(K):
            cur[a, j] -= norm
            alpha[r, j] = cur[a, j]
        loglik += norm

    for i in range(1, L):
        for a in range(A):
            if not mask[a, i]:
                continue
            r = idx[a, i]
            mx = -np.inf
            for j in range(K):
                if cur[a, j] > mx:
                    mx = cur[a, j]
            acc = 0.0
            for j in range(K):
                acc += np.exp(cur[a, j] - mx)
            tot = mx + np.log(acc)
            sw = log_switch[r] + tot
            mx = -np.inf
            for j in range(K):
                v = log_diff[r] + cur[a, j]
                hi = v if v > sw else sw
                lo = sw if v > sw else v
                v = e[r, j] + (hi + np.log1p(np.exp(lo - hi)) if hi > -np.inf else hi)
                nxt[a, j] = v
                if v > mx:
                    mx = v
            acc = 0.0
            for j in range(K):
                acc += np.exp(nxt[a, j] - mx)
            norm = mx + np.log(acc)
            for j in range(K):
                cur[a, j] = nxt[a, j] - norm
                alpha[r, j] = cur[a, j]
            loglik += norm

    for a in range(A):
        for j in range(K):
            cur[a, j] = 0.0
    for i in range(L - 1, 0, -1):
        for a in range(A):
            if not mask[a, i]:
                continue
            r = idx[a, i]
            mx = -np.inf
            for j in range(K):
                v = e[r, j] + cur[a, j]
                nxt[a, j] = v
                if v > mx:
                    mx = v
            acc = 0.0
            for j in range(K):
                acc += np.exp(nxt[a, j] - mx)
            tot = mx + np.log(acc)
            sw = log_switch[r] + tot
            mx = -np.inf
            for j in range(K):
                v = log_diff[r] + nxt[a, j]
                hi = v if v > sw else sw
                lo = sw if v > sw else v
                v = hi + np.log1p(np.exp(lo - hi)) if hi > -np.inf else hi
                nxt[a, j] = v
                if v > mx:
                    mx = v
            acc = 0.0
            for j in range(K):
                acc += np.exp(nxt[a, j] - mx)
            norm = mx + np.log(acc)
            pr = idx[a, i - 1]
            for j in range(K):
                cur[a, j] = nxt[a, j] - norm
                beta[pr, j] = cur[a, j]
    return alpha, beta, loglik


def forward_backward(log_e, log_stay, log_switch, idx, mask):
    """Exact forward-backward over every arm at once.

    Uniform off-diagonal mass makes each step ``O(K)`` rather than ``O(K^2)``. Arms are
    independent, so stepping the padded batch over the longest arm is exact.

    Args:
        log_e: ``(n_sites, K)`` log-emission matrix.
        log_stay: ``log(stay)`` per site.
        log_switch: ``log(switch)`` per site.
        idx: Arm index from :func:`arm_batches`.
        mask: Its real entries.

    Returns:
        Tuple of ``gamma`` ``(n_sites, K)`` and the total log likelihood.
    """
    n, K = log_e.shape
    e = np.empty((n + 1, K))
    e[:n] = log_e
    e[n] = 0.0
    with np.errstate(divide="ignore"):
        log_diff = np.log(np.maximum(np.exp(log_stay) - np.exp(log_switch), 0.0))

    alpha, beta, loglik = _fb_kernel(e, log_diff, log_switch, idx, mask, n)
    g = alpha[:n] + beta[:n]
    g -= logsumexp(g, axis=1, keepdims=True)
    return np.exp(g), loglik


@njit(cache=True)
def _viterbi_kernel(e, log_stay, log_switch, idx, mask, n_sites):
    """MAP path over the padded arm batch.

    Uniform off-diagonal mass means the best predecessor of state ``j`` is either ``j`` itself
    or the best other state, so the top two of ``delta`` are all the step needs.

    Args:
        e: ``(n_sites + 1, K)`` log emissions; padded positions index the trailing row.
        log_stay: ``log(stay)`` per site, length ``n_sites + 1``.
        log_switch: ``log(switch)`` per site, length ``n_sites + 1``.
        idx: ``(n_arms, longest_arm)`` arm index, padded with ``n_sites``.
        mask: Its real entries.
        n_sites: Number of real sites.

    Returns:
        ``(n_sites,)`` int8 MAP state per site.
    """
    A, L = idx.shape
    K = e.shape[1]
    delta = np.empty((A, K))
    back = np.zeros((n_sites + 1, K), dtype=np.int8)
    path = np.zeros(n_sites, dtype=np.int8)
    log_k = np.log(K)

    for a in range(A):
        r = idx[a, 0]
        for j in range(K):
            delta[a, j] = e[r, j] - log_k

    for i in range(1, L):
        for a in range(A):
            if not mask[a, i]:
                continue
            r = idx[a, i]
            best = 0
            top1 = -np.inf
            for j in range(K):
                if delta[a, j] > top1:
                    top1 = delta[a, j]
                    best = j
            second = 0
            top2 = -np.inf
            for j in range(K):
                if j != best and delta[a, j] > top2:
                    top2 = delta[a, j]
                    second = j
            for j in range(K):
                if j == best:
                    other_val = top2
                    other_arg = second
                else:
                    other_val = top1
                    other_arg = best
                stay = delta[a, j] + log_stay[r]
                switch = other_val + log_switch[r]
                if stay >= switch:
                    delta[a, j] = e[r, j] + stay
                    back[r, j] = j
                else:
                    delta[a, j] = e[r, j] + switch
                    back[r, j] = other_arg

    state = np.zeros(A, dtype=np.int64)
    for a in range(A):
        best = 0
        top1 = -np.inf
        for j in range(K):
            if delta[a, j] > top1:
                top1 = delta[a, j]
                best = j
        state[a] = best
    for i in range(L - 1, -1, -1):
        for a in range(A):
            if not mask[a, i]:
                continue
            r = idx[a, i]
            path[r] = state[a]
            if i:
                state[a] = back[r, state[a]]
    return path


def viterbi(log_e, log_stay, log_switch, idx, mask):
    """MAP state path over every arm at once.

    Args:
        log_e: ``(n_sites, K)`` log-emission matrix.
        log_stay: ``log(stay)`` per site.
        log_switch: ``log(switch)`` per site.
        idx: Arm index from :func:`arm_batches`.
        mask: Its real entries.

    Returns:
        ``(n_sites,)`` int8 MAP state per site.
    """
    n, K = log_e.shape
    e = np.empty((n + 1, K))
    e[:n] = log_e
    e[n] = 0.0
    return _viterbi_kernel(e, log_stay, log_switch, idx, mask, n)


def m_step_means(
    k, n, gamma, means, hom_laf, pi, margin, n_fixed=1, tau=np.inf, log_coef=None
):
    """Weighted-moment update of ``b_s``, leaving the first ``n_fixed`` states pinned.

    Each site is weighted by ``gamma`` times ``P(g_i = het | k_i, n_i, s_i)``, then corrected
    by the phase posterior so the estimate is not biased low by ``k`` being a minimum.

    Args:
        k: Minor counts.
        n: Depths.
        gamma: ``(n_sites, K)`` state posteriors.
        means: Current ``b_s`` per state.
        hom_laf: ``b_hom``.
        pi: Het prior.
        margin: Learned states are bounded below by ``means[n_fixed - 1] + margin``.
        n_fixed: Number of leading states held fixed; 1 pins ``b_0``.
        tau: Concentration; ``np.inf`` is the binomial.
        log_coef: Precomputed :func:`log_binom_coef`.

    Note:
        ``b = 0.5`` is a fixed point reached over iterations, not in one step: below it the
        phase posterior favours the observed minor side, so the estimate climbs toward 0.5.

        Under a beta-binomial ``E[k] = n*b`` still holds, so this update stays consistent but
        is no longer the exact MLE, which is why :func:`clonal_loh_hmm` steps back whenever an
        iteration lowers the log likelihood.

    Returns:
        Updated ``b_s``: the pinned prefix unchanged, the rest sorted and clipped.
    """
    if log_coef is None:
        log_coef = log_binom_coef(k, n)
    pi_a = _pi_array(pi, np.shape(k))
    log_pi = np.log(pi_a)
    hom = np.log1p(-pi_a) + bb_logpmf(k, n, hom_laf, tau, log_coef)
    means = np.asarray(means, dtype="float64")
    lo = means[n_fixed - 1] + margin if n_fixed else hom_laf + margin
    out = np.empty(len(means) - n_fixed)
    for j in range(n_fixed, len(means)):
        b = means[j]
        het = bb_logpmf(k, n, b, tau, log_coef)
        w = pi_a * np.exp(het - np.logaddexp(log_pi + het, hom))
        w *= gamma[:, j]
        q = phase_weight(k, n, b, tau, log_coef)
        eff = q * k + (1.0 - q) * (n - k)
        den = float((w * n).sum())
        out[j - n_fixed] = float((w * eff).sum()) / den if den > 0 else b
    return np.concatenate([means[:n_fixed], np.sort(np.clip(out, lo, 0.5))])


def m_step_pi(gamma, het_post):
    """Update the global het prior from the joint responsibilities.

    Args:
        gamma: ``(n_sites, K)`` state posteriors.
        het_post: ``(n_sites, K)`` ``P(het | data, state)``.

    Returns:
        Scalar het prior.
    """
    return float((gamma * het_post).sum() / gamma.shape[0])


def init_means(k, n, n_retained, loh_laf, margin):
    """Starting ``b_1..b_K-1``, from evenly spaced quantiles of the per-site minor fraction.

    Args:
        k: Minor counts.
        n: Depths.
        n_retained: Number of states to initialise.
        loh_laf: ``b_0``.
        margin: Minimum separation from ``b_0``.

    Returns:
        Sorted array of ``n_retained`` levels.
    """
    maf = k / np.maximum(n, 1)
    pool = maf[maf > loh_laf + margin]
    qs = np.linspace(0.1, 0.9, n_retained)
    m = np.quantile(pool, qs) if pool.size else np.linspace(0.1, 0.5, n_retained)
    return np.sort(np.clip(m, loh_laf + margin, 0.5))


def clonal_loh_hmm(
    rc,
    ac,
    pos,
    group_id,
    hom_laf=0.005,
    loh_laf=0.05,
    t=1e-6,
    pi=0.5,
    tau=np.inf,
    n_retained=5,
    n_iter=10,
    margin=0.02,
    learn_pi=True,
    em_tol=1e-5,
    loh_min=0.1,
    p_het_min=0.1,
    min_dp=10,
    init_b=None,
):
    """Fit the chain by EM, decode it, and call each site's germline genotype.

    Args:
        rc: REF counts, in any order; the chain is built by sorting on ``group_id``, ``pos``.
        ac: ALT counts, aligned to ``rc``.
        pos: Site positions.
        group_id: Integer arm per site.
        hom_laf: ``b_hom``.
        loh_laf: ``b_0``.
        t: Breakpoint rate per bp.
        pi: Het prior; scalar or per-site.
        tau: Concentration; ``np.inf`` is the binomial.
        n_retained: Number of states beyond ``b_0``.
        n_iter: Maximum EM iterations; 0 decodes from the initial levels.
        margin: Minimum separation of a learned level from ``b_0``.
        learn_pi: Re-estimate ``pi``. It also absorbs model misspecification, so the fitted
            value settles below the true het rate and is not an estimate of it.
        em_tol: Relative log-likelihood improvement below which EM stops.
        loh_min: Keep ``0/1`` wherever ``p_loh`` reaches this.
        p_het_min: Keep ``0/1`` wherever ``p_het`` reaches this.
        min_dp: Depth below which the site is ``./.``.
        init_b: Starting ``b_1..b_K-1``; quantile-initialised when None.

    Returns:
        Dict with ``gt`` (from :data:`GT_CODES`), ``p_loh``, ``p_het``, ``state``, ``maf``,
        ``means``, ``pi``, ``loglik``, ``n_iter`` (iterations run) and ``converged``.
    """
    rc = np.asarray(rc)
    ac = np.asarray(ac)
    pos = np.asarray(pos)
    n = (rc + ac).astype("int64")
    k = np.minimum(rc, ac).astype("int64")
    with np.errstate(divide="ignore", invalid="ignore"):
        maf = np.where(n > 0, k / np.maximum(n, 1), np.nan)

    idx, mask = arm_batches(np.asarray(group_id), pos)
    log_coef = log_binom_coef(k, n)
    log_stay, log_switch = transition_logs(pos, idx, mask, t, 1 + n_retained)
    logging.info(
        f"clonal-LOH HMM: {k.size} sites over {idx.shape[0]} arm(s), "
        f"longest {int(mask.sum(axis=1).max())} sites"
    )
    retained = (
        init_means(k, n, n_retained, loh_laf, margin)
        if init_b is None
        else np.asarray(init_b)
    )
    means = np.concatenate([[loh_laf], retained])
    prev = -np.inf
    prev_means, prev_pi, prev_gamma = means.copy(), pi, None
    gamma = loglik = None
    converged = stale = False
    for it in range(max(n_iter, 0) + 1):
        log_e = emission_logp(k, n, means, hom_laf, pi, tau, log_coef)
        gamma, loglik = forward_backward(log_e, log_stay, log_switch, idx, mask)
        pi_mean = pi if np.isscalar(pi) else float(np.mean(pi))
        logging.info(
            f"[em {it:>2}] loglik {loglik:,.0f}  pi {pi_mean:.4f}  means {np.round(means, 4)}"
            f"  occ {np.round(gamma.mean(axis=0), 3)}"
        )
        if loglik < prev:
            # NB: the moment M-step is not the exact MLE under a beta-binomial; step back
            means, pi, gamma, loglik = prev_means, prev_pi, prev_gamma, prev
            logging.info(
                f"[em] log likelihood fell at iteration {it}, reverted and stopped"
            )
            converged = stale = True
            break
        if it == n_iter:
            break
        if prev > -np.inf and abs(loglik - prev) < em_tol * abs(prev):
            logging.info(f"[em] converged at iteration {it}")
            converged = True
            break
        prev, prev_means, prev_pi, prev_gamma = loglik, means.copy(), pi, gamma
        means = m_step_means(
            k,
            n,
            gamma,
            means,
            hom_laf,
            pi,
            margin,
            n_fixed=1,
            tau=tau,
            log_coef=log_coef,
        )
        if learn_pi:
            pi = m_step_pi(
                gamma, het_posterior(k, n, means, hom_laf, pi, tau, log_coef)
            )

    if stale:
        log_e = emission_logp(k, n, means, hom_laf, pi, tau, log_coef)
    state = viterbi(log_e, log_stay, log_switch, idx, mask)
    het_post = het_posterior(k, n, means, hom_laf, pi, tau, log_coef)
    p_loh = gamma[:, 0]
    p_het = (gamma * het_post).sum(axis=1)

    gt = np.full(rc.size, "0/1", dtype=GT_DTYPE)
    hom = (p_het < p_het_min) & (p_loh < loh_min)
    gt[hom & (ac < rc)] = "0/0"
    gt[hom & (rc <= ac)] = "1/1"
    gt[n < min_dp] = "./."

    pi_mean = pi if np.isscalar(pi) else float(np.mean(pi))
    logging.info(
        f"fitted means={np.round(means, 4)}, pi={pi_mean:.4f}, loglik={loglik:,.0f}"
    )
    logging.info(
        f"genotypes: het={int((gt == '0/1').sum())}, hom_alt={int((gt == '1/1').sum())}, "
        f"hom_ref={int((gt == '0/0').sum())}, no_call={int((gt == './.').sum())}; "
        f"sites in a Viterbi LOH segment={int((state == 0).sum())}"
    )
    return {
        "gt": gt,
        "p_loh": p_loh,
        "p_het": p_het,
        "state": state,
        "maf": maf,
        "means": means,
        "pi": pi,
        "loglik": loglik,
        "n_iter": it,
        "converged": converged,
    }
