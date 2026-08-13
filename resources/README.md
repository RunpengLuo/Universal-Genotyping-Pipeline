# Resources

## Table of Contents
- [Things You Should Know](#things-you-should-know)
- [Genome Reference](#genome-reference)
  - [Converting NCBI accession-style GTF to Chr notations](#converting-ncbi-accession-style-gtf-to-chr-notations)
- [Pre-built Window BED Files](#pre-built-window-bed-files)
- [SNP Panels](#snp-panels)
  - [Support external SNP panel](#support-external-snp-panel)
- [Population-haplotype Panels](#population-haplotype-panels)
  - [Genetic Maps](#genetic-maps)
- [ENCODE Blacklist](#encode-blacklist)
- [Sequencing bias correction](#sequencing-bias-correction)
  - [Mappability track](#mappability-track)
  - [Replication timing](#replication-timing)

---

## Things You Should Know

> [!IMPORTANT]
> Every input must use `chr`-prefixed (UCSC-style) contig names: the reference FASTA,
> alignments (BAM/CRAM), `snp_panel` / `phasing_panel`, and `region_bed` / `window_bed`.
> The pipeline addresses each chromosome as `chr{chrname}`, where config `chromosomes` stays
> bare (e.g. `[1, 2, ..., 22, X]`). GRCh37/b37-style data with bare contigs
> (`1`, `MT`) will fail - a `chr1` region query does not match a `1` contig - and note that
> many "hg19" BAMs are actually b37. Use the UCSC `chr`-prefixed build, or rename contigs
> first (`bcftools annotate --rename-chrs`, `samtools reheader`) before running.

> [!NOTE]
> Small references (`genome_size`, `region_bed`, blacklist, window BEDs) are bundled under
> [`data/`](data/). Larger resources - SNP / phasing panels and mm10 genetic maps - are not
> bundled; build them on demand with the scripts in [`scripts/`](scripts/) (see below).

---

## Genome Reference

| Species | Reference | Gene Annotation | `genome_size` | `region_bed` |
|---------|-----------|-----------------|---------------|--------------|
| Human | [hg19](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) | - | [hg19.chrom.sizes](data/hg19.chrom.sizes) | [hg19.regions.bed](data/hg19.regions.bed) |
| Human | [hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) | [GENCODE v38](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz), [10x GRCh38-2024-A](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz) | [hg38.chrom.sizes](data/hg38.chrom.sizes) | [hg38.regions.bed](data/hg38.regions.bed) |
| Human | [T2T-chm13v2.0](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz) | [UCSC hs1 ncbiRefSeq](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz), [NCBI accession-style](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz) | [T2T-CHM13v2.0.sizes](data/T2T-CHM13v2.0.sizes) | [T2T-CHM13v2.0.regions.bed](data/T2T-CHM13v2.0.regions.bed) |
| Mouse | [mm10](https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz) | [10x mm10-2020-A](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz) | [mm10.chrom.sizes](data/mm10.chrom.sizes) | [mm10.regions.bed](data/mm10.regions.bed) |

> [!NOTE]
> - 10x annotations sit inside the Cell Ranger tarball at `genes/genes.gtf.gz`.
> - The NCBI accession-style T2T GTF uses RefSeq accessions (e.g. `NC_060925.1`); rename to `chr*` first, see [Converting NCBI accession-style GTF to Chr notations](#converting-ncbi-accession-style-gtf-to-chr-notations).

### Converting NCBI accession-style GTF to `chr` notations 
The NCBI GTF uses RefSeq accessions (e.g., `NC_060925.1`) instead of `chr1`. To rename, download the [assembly report](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_assembly_report.txt) and run:

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

---

## Pre-built Window BED Files

| Species | Reference | Grid | Covariate columns | File |
|---------|-----------|------|-------------------|------|
| Human | hg19 | 1 kb | GC, MAP, REPLI | [windows.1kbp.hg19.bed.gz](data/windows.1kbp.hg19.bed.gz) |
| Human | hg38 | 1 kb | GC, MAP, REPLI | [windows.1kbp.hg38.bed.gz](data/windows.1kbp.hg38.bed.gz) |
| Human | T2T-chm13v2.0 | 1 kb | GC | [windows.1kbp.chm13v2.bed.gz](data/windows.1kbp.chm13v2.bed.gz) |
| Mouse | mm10 | 1 kb | GC | [windows.1kbp.mm10.bed.gz](data/windows.1kbp.mm10.bed.gz) |

> [!TIP]
> - The pre-built windows already accounts for `region_bed` and masked by `blacklist_bed`.
> - Set `params_build_windows.window_size` and let the pipeline to auto-build non pre-built windows.
> - A custom BED file can be provided via `window_bed`, See [reference.md](../docs/reference.md#configuration).

---

## SNP Panels

VCF format. Set via `snp_panel` or `snp_targets` in config. User needs to download and preprocess them.

| Panel | Species | Reference | Source |
|-------|---------|-----------|----------|
| 1kGP phase3 AF>=5e-2 | Human | hg38 | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz) |
| 1kGP phase3 AF>=5e-4 | Human | hg38 | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz) |
| 1kGP n=3,202 | Human | hg38 | [FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) — use [`scripts/process_1kGP_3202_panel.sh --ref hg38`](scripts/process_1kGP_3202_panel.sh) |
| 1kGP n=3,202 | Human | T2T-chm13v2.0 | [S3](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/) — use [`scripts/process_1kGP_3202_panel.sh --ref chm13v2`](scripts/process_1kGP_3202_panel.sh) |
| MGP v5 | Mouse | mm10 | [UCSC MGP v5 VCF](https://hgdownload.soe.ucsc.edu/gbdb/mm10/mouseStrains/mgpV5MergedSNPsAlldbSNP142.vcf.gz) -> [`scripts/build_mouse_mgp_panel.sh`](scripts/build_mouse_mgp_panel.sh) (all strains; `--strains` to subset; biallelic SNP panel + phasing panel + targets) |

### Support external SNP panel

The `genotype_snps_bulk` rule requires `config["snp_targets"]` — a directory of per-chromosome position files. To build these from any SNP panel VCF:

We provide a preparation script [scripts/build_snp_targets.sh](./scripts/build_snp_targets.sh) to generate targeted germline SNP positions required by the pipeline. The script requires `bcftools`, `bgzip`, and `tabix` on `$PATH`.

```bash
bash resources/scripts/build_snp_targets.sh /path/to/snp_panel.vcf.gz /path/to/snp_targets
# /path/to/snp_targets/target.chr{1..22,X}.pos.gz + .tbi
```

> [!TIP]
> We recommend user to use the full SNP panel without AF cutoff to achieve
> best germline SNP genotyping performance.

---

## Population-haplotype Panels

| Panel | Species | Reference | Source |
|-------|---------|-----------|----------|
| 1kGP phase3 (n=2,504) | Human | hg38 | [download](http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip) |
| 1kGP phase3 (n=3,202) | Human | hg38 | produced by `process_1kGP_3202_panel.sh --ref hg38` (see SNP Panels) |
| 1kGP n=3,202 | Human | T2T chm13v2.0 | produced by `process_1kGP_3202_panel.sh --ref chm13v2` (see SNP Panels). Phased with SHAPEIT5 v1.1 + T2T-native maps ([phasing_T2T](https://github.com/JosephLalli/phasing_T2T)) |
| gnomAD HGDP+1KG (n=4,099) | Human | hg38 | `gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes` |
| TOPMed (n=97,256) | Human | hg38 | via [imputation server](https://imputation.biodatacatalyst.nhlbi.nih.gov) |
| MGP v5 | Mouse | mm10 | per-chrom phased BCFs from the [UCSC MGP v5 VCF](https://hgdownload.soe.ucsc.edu/gbdb/mm10/mouseStrains/mgpV5MergedSNPsAlldbSNP142.vcf.gz), built by [`scripts/build_mouse_mgp_panel.sh`](scripts/build_mouse_mgp_panel.sh) (all 36 strains by default, `--strains` to subset; inbred GTs phased trivially, strain-het masked) |

### Genetic Maps

| Species | Reference | Source |
|---|---|---|
| Human | **hg38** | bundled with [Eagle2](https://github.com/poruloh/Eagle) (`tables/`) and [SHAPEIT5](https://github.com/odelaneau/shapeit) (`resources/maps/b38/`) |
| Human | **chm13v2** | download [T2T-native scaled maps](https://github.com/JosephLalli/phasing_T2T/tree/main/resources/recombination_maps/t2t_native_scaled_maps); convert via [`scripts/convert_gmap_to_eagle.py`](scripts/convert_gmap_to_eagle.py) for Eagle2 |
| Mouse | **mm10** | download CoxMapV3 [`OrigMaps/CoxMaps_rev_build38.csv`](https://raw.githubusercontent.com/kbroman/CoxMapV3/main/OrigMaps/CoxMaps_rev_build38.csv), then build via [`scripts/build_mouse_gmap_mm10.py`](scripts/build_mouse_gmap_mm10.py); writes `eagle2/genetic_map_mm10_withX.txt.gz` + `shapeit5/chr{N}.mm10.gmap.gz` |

Here is an example config for `hg38`:

```yaml
# Eagle2 (single file)
phaser: "eagle"
gmap_path: "/path/to/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"

# SHAPEIT5 (per-chromosome, {chrname} substituted per chromosome by pipeline)
phaser: "shapeit"
gmap_path: "/path/to/shapeit5/resources/maps/b38/chr{chrname}.b38.gmap.gz"

# mm10 (Eagle2; download CoxMaps_rev_build38.csv, then build:
#   python resources/scripts/build_mouse_gmap_mm10.py CoxMaps_rev_build38.csv /path/to/mm10_gmap)
phaser: "eagle"
gmap_path: "/path/to/mm10_gmap/eagle2/genetic_map_mm10_withX.txt.gz"
```

---

## ENCODE Blacklist

| Blacklist | Species | Reference | Source |
|-----------|---------|-----------|----------|
| ENCODE blacklist v2 | Human | hg38 | [hg38-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz) |
| ENCODE blacklist v2 | Mouse | mm10 | [mm10-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz) |

---

## Sequencing bias correction
### Mappability track

Optional; set `mappability_bed` in config to add the `MAP` column. Convert the bigWig to BED with `bigWigToBedGraph` (UCSC tools).

| Track | Species | Reference | Source |
|-------|---------|-----------|----------|
| k100 Umap multi-track | Human | hg38 | [k100.Umap.MultiTrackMappability.bw](http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw) |

### Replication timing

This is automatically handled by the pipeline via rule `repliseq_bigwig_to_bedgraph`, then `repliseq_liftover` when the run is not hg19 (`REPLI_LIFTOVER` in `config/const.py`).

| Track | Species | Reference | Source |
|-------|---------|-----------|----------|
| ENCODE UW Repli-seq WaveSignal (16 bigWig) | Human | hg19 | [UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwRepliSeq/) |
| liftOver chain hg19 -> hg38 | Human | hg38 | [hg19ToHg38.over.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz) |
| liftOver chain hg19 -> chm13v2 (`hs1`) | Human | chm13v2 | [hg19ToHs1.over.chain.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHs1.over.chain.gz) |

