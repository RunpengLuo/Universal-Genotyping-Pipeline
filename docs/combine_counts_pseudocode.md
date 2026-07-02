# `combine_counts` (bulk) — pseudo-code

SNP-informed adaptive binning of a fixed window grid, then allele/depth aggregation + RDR.

## 1. Parameters & inputs

- **Inputs:** SNPs $S$ (`#CHR, POS0, PS`, opt. `feature_id`); allele matrices $T,A,B \in \mathbb{R}^{N\times M}$; window grid $\mathcal W$ with corrected depth $D\in\mathbb{R}^{W\times M}$; tumor cols $\mathcal T$.
- **Params:** $\rho$=`min_snp_reads`, $\kappa$=`min_snp_per_block`, $L$=`max_blocksize`, gene-aware $g$, group keys $G$ = `region_id, PS[, phase_group]`.
- **Setup:** assign each SNP to its window; roll up per window $n_w=\#\text{SNPs}$, $\;r_w=\sum_{s\in w}T_{s,\mathcal T}$.

## 2. Adaptive binning loop

Within each group (group change ⇒ forced cut), greedily merge consecutive windows:

```
acc, n ← r_w0, n_w0                       # accumulate reads + SNP count
for next window w:
    close = ( (min_j acc[j] ≥ ρ  AND  n ≥ κ)  OR  span ≥ L )  AND  unit(w) ≠ unit(w-1)
    if close:  emit bin; acc, n ← r_w, n_w
    else:      acc += r_w;  n += n_w
# last block: keep if it meets ρ,κ (or L); else merge into previous bin
```

Gives $K$ bins and SNP→bin map $\beta$. Aggregate:

$$T^{bb}_{k}=\!\!\sum_{\beta(s)=k}\!\!T_s,\quad \text{BAF}_k=\frac{B^{bb}_k}{T^{bb}_k},\quad D^{bb}_k=\sum_{w\in k}D_w,\quad \text{RDR}_{k}=\frac{D^{bb}_k/\!\sum D^{bb}}{\;\text{(matched normal)}}$$

Drop NaN bins, re-index, assign per-bin switch probs.

## 3. Output (`bb_dir/MSR{ρ}/{stream}/`, one `MSR{ρ}/` subdir per $\rho$)

`bb.tsv.gz` (bin table) · `bb.{Tallele,Aallele,Ballele,depth,rdr}.npz` · `sample_ids.tsv` · QC PDF `combine_counts.{stream}.MSR{ρ}.pdf`.

## 4. Multiple bin sizes — pick one

`min_snp_reads` is a **list** ($\rho$ controls bin size). One `combine_counts` job loads/preprocesses **once** (SNP→window, phase-flip, gene blocks are $\rho$-independent) and loops the binning over each $\rho$, writing one `MSR{ρ}/` subdir. Default: $\rho\in\{500,1000\}$.

**Noise metric** (no ground truth — adjacent bins are mostly same CN, so jitter = noise):
$$V(\rho)=\operatorname{median}_k\,|x_{k+1}-x_k|,\qquad x=\text{RDR or }|\text{BAF}-0.5|$$

**Selection (simplest):** BAF SE $\approx 0.5/\sqrt{\rho}$, so for target noise $\sigma^*$ take $\rho\ge(0.5/\sigma^*)^2$ (e.g. $\sigma^*{=}0.05\Rightarrow\rho{\approx}100$). No sweep needed.

**Or elbow:** plot $V(\rho)$ vs #bins $K(\rho)$, normalize both to $[0,1]$, pick the knee (point farthest from the endpoint chord) — best resolution-vs-noise trade.

> **TODO:** auto-recommend a default $\rho$ at the end of `combine_counts*` (compute $V(\rho)$/$K(\rho)$ over the swept subdirs, report the elbow). Marked as a `TODO` in both scripts; for now, pick by eye from the per-$\rho$ QC PDFs.
