"""Pseudobulk SNP annotation helpers.

Runpeng Luo
Last update: 2026-08-06

Dependencies:
  pandas; script_utils (io_utils, utils).

Notes/References:
  Caller: workflow/scripts/annotate_snps_pseudobulk.py
"""

import pandas as pd  # noqa: F401


def compute_snp_statistics(
    raw_snps_list,
    modalities,
    base_snps,
    chrom_sizes,
    chroms,
):
    """LEGACY, unused: per-chromosome SNP genotyping statistics.

    No caller since annotate_snps_pseudobulk stopped writing
    pseudobulk_snp_statistics.tsv. Kept for reference; delete when nothing needs it.

    Returns a DataFrame with columns:
    #CHR, LENGTH, #SNPS-DNA, #SNPS-RNA, #SNPS-shared, #SNPS-union,
    #SNPS-het, #SNPS-hom_alt, #SNPS-hom_ref, #SNPS-dropped, #SNPS-kept.
    """
    modality_keys = {}
    for idx, modality in enumerate(modalities):
        modality_keys[modality] = set(raw_snps_list[idx]["KEY"])

    has_dna = "DNA" in modality_keys
    has_rna = "RNA" in modality_keys
    if has_dna and has_rna:
        shared_keys = modality_keys["DNA"] & modality_keys["RNA"]
    else:
        shared_keys = set()

    rows = []
    for chrom in chroms:
        length = chrom_sizes.get(chrom, 0)
        ch_mask = base_snps["#CHROM"] == chrom

        n_dna = 0
        n_rna = 0
        if has_dna:
            n_dna = int(
                (raw_snps_list[modalities.index("DNA")]["#CHROM"] == chrom).sum()
            )
        if has_rna:
            n_rna = int(
                (raw_snps_list[modalities.index("RNA")]["#CHROM"] == chrom).sum()
            )

        ch_keys = set(base_snps.loc[ch_mask, "KEY"])
        n_shared = len(ch_keys & shared_keys)
        n_union = int(ch_mask.sum())

        ch_snps = base_snps.loc[ch_mask]
        n_het = int((ch_snps["SAMPLE"] == "0/1").sum())
        n_hom_alt = int((ch_snps["SAMPLE"] == "1/1").sum())
        n_hom_ref = int((ch_snps["SAMPLE"] == "0/0").sum())
        n_dropped = n_union - n_het - n_hom_alt - n_hom_ref
        n_kept = n_het + n_hom_alt

        rows.append(
            [
                chrom,
                length,
                n_dna,
                n_rna,
                n_shared,
                n_union,
                n_het,
                n_hom_alt,
                n_hom_ref,
                n_dropped,
                n_kept,
            ]
        )

    # TOTAL row
    total_len = sum(chrom_sizes.get(c, 0) for c in chroms)
    total_dna = sum(r[2] for r in rows)
    total_rna = sum(r[3] for r in rows)
    total_shared = sum(r[4] for r in rows)
    total_union = sum(r[5] for r in rows)
    total_het = sum(r[6] for r in rows)
    total_hom_alt = sum(r[7] for r in rows)
    total_hom_ref = sum(r[8] for r in rows)
    total_dropped = sum(r[9] for r in rows)
    total_kept = sum(r[10] for r in rows)
    rows.append(
        [
            "TOTAL",
            total_len,
            total_dna,
            total_rna,
            total_shared,
            total_union,
            total_het,
            total_hom_alt,
            total_hom_ref,
            total_dropped,
            total_kept,
        ]
    )

    columns = [
        "#CHR",
        "LENGTH",
        "#SNPS-DNA",
        "#SNPS-RNA",
        "#SNPS-shared",
        "#SNPS-union",
        "#SNPS-het",
        "#SNPS-hom_alt",
        "#SNPS-hom_ref",
        "#SNPS-dropped",
        "#SNPS-kept",
    ]
    return pd.DataFrame(rows, columns=columns)
