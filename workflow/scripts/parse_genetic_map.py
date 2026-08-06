"""Parse Eagle/SHAPEIT genetic map files into one chr-prefixed table.

Inputs
  gmap_files: one map for all chromosomes (Eagle) or one per chromosome (SHAPEIT).
Parameters
  chroms: this run's chromosomes, chr-prefixed
  phaser: eagle | shapeit
  species: selects the sex-chromosome numbering used to relabel numeric labels
Outputs:
  gmap_tsv: #CHR POS cM, chr-prefixed, genomically sorted
"""

import logging

snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, sort_df_chr, SPECIES2SEXCHROM

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import pandas as pd


# inputs
gmap_files = list(snakemake_handle.input["gmap_files"])

# parameters
chroms = list(snakemake_handle.params["chroms"])
phaser = snakemake_handle.params["phaser"]
species = snakemake_handle.params["species"]

# outputs
gmap_tsv = snakemake_handle.output["gmap_tsv"]

logging.info(f"parse genetic map files, phaser={phaser}, species={species}")

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
    wanted = {str(c) for c in chroms}
    genetic_map["#CHR"] = (
        genetic_map["#CHR"].astype(str).str.replace(r"^chr", "", regex=True)
    )
    labels = set(genetic_map["#CHR"])
    if not wanted <= labels:
        sexmap = SPECIES2SEXCHROM.get(species)
        if sexmap is None:
            raise ValueError(
                f"eagle gmap lacks requested chromosome(s) {sorted(wanted - labels)} "
                f"and species={species!r} has no sex-chromosome numbering to relabel "
                f"by; known species: {sorted(SPECIES2SEXCHROM)}"
            )
        for sex_name, sex_num in sexmap.items():
            if sex_name not in wanted or sex_name in labels:
                continue
            if str(sex_num) in labels:
                genetic_map.loc[genetic_map["#CHR"] == str(sex_num), "#CHR"] = sex_name
                labels = set(genetic_map["#CHR"])
                logging.info(
                    f"eagle: relabeled chromosome {sex_num} as {sex_name} "
                    f"(species={species})"
                )
        missing_chroms = sorted(wanted - labels)
        if missing_chroms:
            raise ValueError(
                f"eagle gmap has no rows for requested chromosome(s) {missing_chroms}; "
                f"map labels are {sorted(labels)} (species={species}). Supply a genetic "
                "map that covers them, or relabel it to match."
            )
    genetic_map = genetic_map[genetic_map["#CHR"].isin(wanted)].reset_index(drop=True)
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
    for chrom, gmap_file in zip(chroms, gmap_files):
        genetic_map = pd.read_csv(
            gmap_file,
            sep="\t",
            comment="#",
        )
        assert "pos" in genetic_map.columns and "cM" in genetic_map.columns, (
            "gmap.gz file is invalid"
        )
        genetic_map["#CHR"] = f"chr{chrom}"
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
