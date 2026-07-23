# Resources

External data required by the pipeline.

---

## SNP Panels

VCF format. Set via `snp_panel` or `snp_targets` in config.

| Panel | Download |
|-------|----------|
| 1kGP phase3 AF>=5e-2 (~92 MB, hg38) | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz) |
| 1kGP phase3 AF>=5e-4 (~568 MB, hg38) | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz) |
| 1kGP n=3,202 (hg38) | [FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) — use [`scripts/process_1kGP_3202_panel.sh --ref hg38`](scripts/process_1kGP_3202_panel.sh) |
| 1kGP n=3,202 (chm13v2.0, biallelic) | [S3](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/) — use [`scripts/process_1kGP_3202_panel.sh --ref chm13v2`](scripts/process_1kGP_3202_panel.sh) |
| MGP v5 strain SNPs (mm10, biallelic) | [UCSC](https://hgdownload.soe.ucsc.edu/gbdb/mm10/mouseStrains/mgpV5MergedSNPsAlldbSNP142.vcf.gz) — use [`build_mouse_b6_129_panel.sh`](../../LLM%20plans/allen_collab/build_mouse_b6_129_panel.sh) (produces sites-only `snp_panel`, multi-strain `phasing_panel/`, and `target_positions/` in one pass) |

### Building `snp_targets` from any panel VCF

The `genotype_snps_bulk` rule requires `config["snp_targets"]` — a directory of per-chromosome position files. To build these from any SNP panel VCF:

```bash
bash resources/scripts/build_snp_targets.sh /path/to/snp_panel.vcf.gz /path/to/snp_targets
```

This produces `target.chr{1..22,X}.pos.gz` + `.tbi` index files. Existing chromosomes are skipped (idempotent). Requires `bcftools`, `bgzip`, `tabix`.

---

## Phasing Panels

BCF format, one per chromosome. Set via `phasing_panel` in config.

| Panel | Download |
|-------|----------|
| 1kGP phase3 (n=2,504, hg38) | [download](http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip) |
| 1kGP phase3 (n=3,202, hg38) | produced by `process_1kGP_3202_panel.sh --ref hg38` (see SNP Panels) |
| 1kGP n=3,202 (chm13v2.0) | produced by `process_1kGP_3202_panel.sh --ref chm13v2` (see SNP Panels). Phased with SHAPEIT5 v1.1 + T2T-native maps ([phasing_T2T](https://github.com/JosephLalli/phasing_T2T)) |
| gnomAD HGDP+1KG (n=4,099, hg38) | `gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes` |
| TOPMed (n=97,256, hg38) | via [imputation server](https://imputation.biodatacatalyst.nhlbi.nih.gov) |
| MGP v5 strains (mm10) | produced by [`build_mouse_b6_129_panel.sh`](../../LLM%20plans/allen_collab/build_mouse_b6_129_panel.sh); inbred strain GTs phased trivially (`0/0`→`0|0`, `1/1`→`1|1`, strain-het→`./.`) |

### Genetic Maps

Set the full gmap path via `gmap_path` in config. Use `{chrname}` placeholder for per-chrom files (SHAPEIT5), or a literal path for the single-file case (Eagle2).

| Reference | Source | Example `gmap_path` |
|---|---|---|
| **hg38** | bundled with Eagle2 (`tables/`) and SHAPEIT5 (`resources/maps/b38/`) | SHAPEIT5: `/path/to/shapeit5/resources/maps/b38/chr{chrname}.b38.gmap.gz`<br>Eagle2: `/path/to/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz` |
| **chm13v2** | download [T2T-native scaled maps](https://github.com/JosephLalli/phasing_T2T/tree/main/resources/recombination_maps/t2t_native_scaled_maps); convert via [`scripts/convert_gmap_to_eagle.py`](scripts/convert_gmap_to_eagle.py) for Eagle2 | SHAPEIT5: `/path/to/chm13v2_maps/chr{chrname}.t2t.scaled.gmap.gz`<br>Eagle2: `/path/to/eagle_chm13v2/genetic_map_chm13v2_withX.txt.gz` |
| **mm10** | build from Karl Broman's CoxMapV3 (see [`build_mouse_gmap_mm10.sh`](../../LLM%20plans/allen_collab/build_mouse_gmap_mm10.sh)) — produces both SHAPEIT5 per-chrom files and a single Eagle2 file | SHAPEIT5: `/path/to/mm10_gmap/shapeit5/chr{chrname}.mm10.gmap.gz`<br>Eagle2: `/path/to/mm10_gmap/eagle2/genetic_map_mm10_withX.txt.gz` |

---

## Gene Annotation (GTF)

Set via `gtf_file` in config.

| Source | Download |
|--------|----------|
| GENCODE v38 (hg38) | `wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz` |
| 10x GRCh38-2024-A (hg38) | `curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz` → `genes/genes.gtf.gz` |
| 10x mm10-2020-A (mm10) | `curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz` → `genes/genes.gtf.gz` |
| UCSC ncbiRefSeq (chm13v2.0, chr-style) | `wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz` |
| NCBI RefSeq (chm13v2.0, accession-style) | [NCBI FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/) — `GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz` |

**Converting NCBI accession-style GTF to chr-style:** The NCBI GTF uses RefSeq accessions (e.g., `NC_060925.1`) instead of `chr1`. To rename, download the [assembly report](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt) and run:

```bash
# 1. Build accession→chr mapping (tab-separated) from the assembly report
grep -v '^#' GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt \
  | tr -d '\r' \
  | awk -F'\t' '$10 ~ /^chr([0-9]+|[XY])$/ {print $7 "\t" $10}' > accession_to_chr.tsv

# 2. Rename and filter chromosomes in the GTF
zcat GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz \
  | tr -d '\r' \
  | awk -F'\t' -v OFS='\t' 'NR==FNR{m[$1]=$2;next} /^#/{print;next} ($1 in m){$1=m[$1];print}' accession_to_chr.tsv - \
  | gzip > GCF_009914755.1_T2T-CHM13v2.0_genomic.chr.gtf.gz
```

Build 10x Cell Ranger ARC reference with [`scripts/build_cellranger_arc_ref_chm13v2.sh`](scripts/build_cellranger_arc_ref_chm13v2.sh) ([10x guide](https://kb.10xgenomics.com/s/article/29207065679501-Building-a-Custom-T2T-reference-for-Cell-Ranger-ARC)).

---

## Window BED (Bias Correction)

Built in-workflow per stream by `workflow/rules/build_windows.smk` (bulk only); there is no `window_bed` config key. `build_segment_bed` first derives `aux/segment.bed` (region_id arm + seg_id chunk) from `region_bed`, splitting each arm at the union of `files.breakpoint_bedpe` cuts. Then per stream: WGS tiles fixed-size windows within the segment BED (`aux/wgs_windows.bed.gz`); WES adaptively tiles the `files.wes_targets_bed` capture targets within it (`aux/wes_windows.bed.gz`). Each window gets a `GC` column, plus optional `MAP` (from `mappability_bed`) and `REPLI` (Repli-seq) columns.

### Mappability track

Optional; set `mappability_bed` in config to add the `MAP` column.

- [k100.Umap.MultiTrackMappability.bw](http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw) — bigWig format, convert to BED with `bigWigToBedGraph` (UCSC tools).

### Replication timing (Repli-seq)

Added automatically for `hg19`/`hg38` (auto-detected from `reference_version`) as the `REPLI` column. The build downloads 16 ENCODE Repli-seq WaveSignal bigWig files from [UCSC](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/) (hg19), converts via `bigWigToBedGraph`, and lifts to hg38 using [hg19ToHg38.over.chain.gz](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz). Requires `bigWigToBedGraph` and `liftOver` (UCSC tools).

### WES exon capture targets

WES mode needs a vendor exon capture BED as `files.wes_targets_bed` on a bulkWES record. Example (IDT xGen):

- [xgen-exome-research-panel-targets-hg38.bb](https://hgdownload.soe.ucsc.edu/gbdb/hg38/exomeProbesets/xgen-exome-research-panel-targets-hg38.bb) — bigBed format, convert to BED with `bigBedToBed` (UCSC tools).

---

## Blacklist BED

Set via `blacklist_bed` in config. Pre-built:

- `data/hg38-blacklist.v2.bed.gz` — [ENCODE blacklist v2](https://github.com/Boyle-Lab/Blacklist)
- `data/mm10-blacklist.v2.bed.gz` — [ENCODE blacklist v2 (mm10)](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz)

---

## Phasing Tools

One required for `bulk_genotyping` and `single_cell_genotyping` modes.

| Tool | Source |
|------|--------|
| Eagle2 | [download](https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz) |
| SHAPEIT5 | [GitHub](https://github.com/odelaneau/shapeit5) |
| LongPhase | [GitHub](https://github.com/twolinin/longphase) |
