import os
import logging

snakemake_handle = snakemake

t = int(getattr(snakemake_handle, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import pandas as pd

from utils import sort_df_chr, REFVER2SEXCHROM

##################################################
"""
Parse genetic map files from Shapeit or Eagle resources.
chr-prefix will always be added to comply with other tools.
"""

log_file = snakemake_handle.log[0]
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

# inputs
gmap_files = list(snakemake_handle.input["gmap_files"])

# parameters
chrnames = list(snakemake_handle.params["chrnames"])
phaser = snakemake_handle.params["phaser"]
reference_version = snakemake_handle.params["reference_version"]

# outputs
gmap_tsv = snakemake_handle.output["gmap_tsv"]

logging.info(
    f"parse genetic map files, phaser={phaser}, reference_version={reference_version}"
)

required_columns = ["#CHR", "POS", "cM"]
if phaser == "eagle":
    # Eagle map(s): whitespace-delimited; one file for all chroms or one per chrom.
    rename_cols = {
        "chr": "#CHR",
        "position": "POS",
        "COMBINED_rate(cM/Mb)": "recomb_rate",
        "Genetic_Map(cM)": "cM",
    }
    files = [gmap_files[0]] if len(set(gmap_files)) == 1 else gmap_files
    genetic_map = pd.concat(
        [
            pd.read_table(f, sep=" ", index_col=None, dtype={"chr": str}).rename(
                columns=rename_cols
            )
            for f in files
        ],
        ignore_index=True,
    )
    missing = [c for c in required_columns if c not in genetic_map.columns]
    assert not missing, (
        f"eagle gmap missing expected columns {missing}; "
        f"got {list(genetic_map.columns)}"
    )

    # Strip any 'chr' prefix, filter to requested chroms (as bare strings), re-add.
    wanted = {str(c) for c in chrnames}
    genetic_map["#CHR"] = (
        genetic_map["#CHR"].astype(str).str.replace(r"^chr", "", regex=True)
    )
    # Relabel numeric sex chroms to letters: via REFVER2SEXCHROM for a known
    # reference, else the "largest non-autosome int is X" heuristic. Labels are
    # strings, so compare by int (key=int).
    labels = set(genetic_map["#CHR"])
    sexmap = REFVER2SEXCHROM.get(reference_version)
    if sexmap is not None:
        for sex_name, sex_num in sexmap.items():
            if sex_name not in wanted or sex_name in labels:
                continue
            if str(sex_num) in labels:
                genetic_map.loc[genetic_map["#CHR"] == str(sex_num), "#CHR"] = sex_name
                logging.info(
                    f"eagle: relabeled chromosome {sex_num} as {sex_name} "
                    f"(reference_version={reference_version})"
                )
            else:
                logging.warning(
                    f"eagle: {sex_name} requested but chromosome {sex_num} absent "
                    f"in gmap (reference_version={reference_version})"
                )
    elif "X" in wanted and "X" not in labels:
        autosomes = {c for c in wanted if c.isdigit()}
        cand = [c for c in labels if c.isdigit() and c not in autosomes]
        if cand:
            x_label = max(cand, key=int)
            genetic_map.loc[genetic_map["#CHR"] == x_label, "#CHR"] = "X"
            logging.info(
                f"eagle: relabeled largest int chromosome '{x_label}' as 'X' "
                f"(heuristic; reference_version={reference_version!r})"
            )
        else:
            logging.warning(
                f"eagle: X requested but no X-like chromosome found in gmap "
                f"(reference_version={reference_version!r})"
            )

    genetic_map = genetic_map[genetic_map["#CHR"].isin(wanted)].reset_index(drop=True)
    assert len(genetic_map) > 0, (
        f"no eagle gmap rows match requested chromosomes {sorted(wanted)}; "
        f"map chromosome labels were {sorted(labels)}"
    )
    genetic_map["#CHR"] = "chr" + genetic_map["#CHR"]
    genetic_map = sort_df_chr(genetic_map, ch="#CHR", pos="POS")
    logging.info(
        f"eagle: #rows={len(genetic_map)}, "
        f"chroms={sorted(genetic_map['#CHR'].unique().tolist())}, "
        f"cM range=[{genetic_map['cM'].min():.4f}, {genetic_map['cM'].max():.4f}]"
    )
    genetic_map[required_columns].to_csv(gmap_tsv, sep="\t", header=True, index=False)

if phaser == "shapeit":
    genetic_maps = []
    for chrname, gmap_file in zip(chrnames, gmap_files):
        genetic_map = pd.read_csv(
            gmap_file,
            sep="\t",
            comment="#",
        )
        assert "pos" in genetic_map.columns and "cM" in genetic_map.columns, (
            "gmap.gz file is invalid"
        )
        genetic_map["#CHR"] = f"chr{chrname}"
        genetic_map["POS"] = genetic_map["pos"]

        genetic_maps.append(genetic_map[["#CHR", "POS", "cM"]].reset_index(drop=True))

    genetic_map = pd.concat(genetic_maps, ignore_index=True)
    genetic_map = sort_df_chr(genetic_map, ch="#CHR", pos="POS")
    logging.info(
        f"shapeit: #rows={len(genetic_map)}, "
        f"chroms={sorted(genetic_map['#CHR'].unique().tolist())}, "
        f"cM range=[{genetic_map['cM'].min():.4f}, {genetic_map['cM'].max():.4f}]"
    )
    genetic_map[required_columns].to_csv(gmap_tsv, sep="\t", header=True, index=False)
