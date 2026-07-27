# Resources

## Table of Contents
- [SNP Panels](#snp-panels)
  - [Support external SNP panel](#support-external-snp-panel)
- [Population-haplotype Panels](#population-haplotype-panels)
  - [Genetic Maps](#genetic-maps)
- [Gene Annotation (GTF)](#gene-annotation-gtf)
  - [Converting NCBI accession-style GTF to Chr notations](#converting-ncbi-accession-style-gtf-to-chr-notations)
- [Genome Sizes & Regions](#genome-sizes--regions)
- [ENCODE Blacklist](#encode-blacklist)
- [Sequencing bias correction (Optional)](#sequencing-bias-correction-optional)
  - [Mappability track](#mappability-track)
  - [Replication timing](#replication-timing)

---

## SNP Panels

VCF format. Set via `snp_panel` or `snp_targets` in config.

| Panel | Species | Reference | Download |
|-------|---------|-----------|----------|
| 1kGP phase3 AF>=5e-2 | Human | hg38 | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz) |
| 1kGP phase3 AF>=5e-4 | Human | hg38 | [download](https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz) |
| 1kGP n=3,202 | Human | hg38 | [FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) — use [`scripts/process_1kGP_3202_panel.sh --ref hg38`](scripts/process_1kGP_3202_panel.sh) |
| 1kGP n=3,202 | Human | T2T-chm13v2.0 | [S3](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/) — use [`scripts/process_1kGP_3202_panel.sh --ref chm13v2`](scripts/process_1kGP_3202_panel.sh) |
| MGP v5 strain | Mouse | mm10 | [UCSC](https://hgdownload.soe.ucsc.edu/gbdb/mm10/mouseStrains/mgpV5MergedSNPsAlldbSNP142.vcf.gz) — TODO |

### Support external SNP panel

The `genotype_snps_bulk` rule requires `config["snp_targets"]` — a directory of per-chromosome position files. To build these from any SNP panel VCF:

We provide a preparation script [scripts/build_snp_targets.sh](./scripts/build_snp_targets.sh) to generate targeted germline SNP positions required by the pipeline. The script requires `bcftools`, `bgzip`, and `tabix` on `$PATH`.

```bash
bash resources/scripts/build_snp_targets.sh /path/to/snp_panel.vcf.gz /path/to/snp_targets
# /path/to/snp_targets/target.chr{1..22,X}.pos.gz + .tbi
```

---

## Population-haplotype Panels

| Panel | Species | Reference | Download |
|-------|---------|-----------|----------|
| 1kGP phase3 (n=2,504) | Human | hg38 | [download](http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip) |
| 1kGP phase3 (n=3,202) | Human | hg38 | produced by `process_1kGP_3202_panel.sh --ref hg38` (see SNP Panels) |
| 1kGP n=3,202 | Human | T2T chm13v2.0 | produced by `process_1kGP_3202_panel.sh --ref chm13v2` (see SNP Panels). Phased with SHAPEIT5 v1.1 + T2T-native maps ([phasing_T2T](https://github.com/JosephLalli/phasing_T2T)) |
| gnomAD HGDP+1KG (n=4,099) | Human | hg38 | `gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes` |
| TOPMed (n=97,256) | Human | hg38 | via [imputation server](https://imputation.biodatacatalyst.nhlbi.nih.gov) |
| MGP v5 strains | Mouse | mm10 | TODO |

### Genetic Maps

| Reference | Source | Example `gmap_path` |
|---|---|---|
| **hg38** | bundled with Eagle2 (`tables/`) and SHAPEIT5 (`resources/maps/b38/`) | SHAPEIT5: `/path/to/shapeit5/resources/maps/b38/chr{chrname}.b38.gmap.gz`<br>Eagle2: `/path/to/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz` |
| **chm13v2** | download [T2T-native scaled maps](https://github.com/JosephLalli/phasing_T2T/tree/main/resources/recombination_maps/t2t_native_scaled_maps); convert via [`scripts/convert_gmap_to_eagle.py`](scripts/convert_gmap_to_eagle.py) for Eagle2 | SHAPEIT5: `/path/to/chm13v2_maps/chr{chrname}.t2t.scaled.gmap.gz`<br>Eagle2: `/path/to/eagle_chm13v2/genetic_map_chm13v2_withX.txt.gz` |
| **mm10** | build from Karl Broman's CoxMapV3 (`build_mouse_gmap_mm10.sh`, TODO) — produces both SHAPEIT5 per-chrom files and a single Eagle2 file | SHAPEIT5: `/path/to/mm10_gmap/shapeit5/chr{chrname}.mm10.gmap.gz`<br>Eagle2: `/path/to/mm10_gmap/eagle2/genetic_map_mm10_withX.txt.gz` |

---

## Gene Annotation (GTF)

Set via `gtf_file` in config.

| Source | Species | Reference | Download |
|--------|---------|-----------|----------|
| GENCODE v38 | Human | hg38 | [gencode.v38.annotation.gtf.gz](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz) |
| 10x GRCh38-2024-A | Human | hg38 | [refdata-gex-GRCh38-2024-A.tar.gz](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz): `genes/genes.gtf.gz` |
| 10x mm10-2020-A | Mouse | mm10 | [refdata-gex-mm10-2020-A.tar.gz](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz): `genes/genes.gtf.gz` |
| UCSC ncbiRefSeq (chr-style) | Human | T2T-chm13v2.0 | [hs1.ncbiRefSeq.gtf.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz) |
| NCBI RefSeq (accession-style) | Human | T2T-chm13v2.0 | [GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz) |

### Converting NCBI accession-style GTF to Chr notations 
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

## Genome Sizes & Regions

Set via `genome_size` (two-column `chrom<TAB>size`) and `region_bed` (arm-level whitelist) in config; pre-built under `data/`.

| Species | Reference | `genome_size` | `region_bed` |
|---------|-----------|---------------|--------------|
| Human | hg19 | [hg19.chrom.sizes](data/hg19.chrom.sizes) | [hg19.regions.bed](data/hg19.regions.bed) |
| Human | hg38 | [hg38.chrom.sizes](data/hg38.chrom.sizes) | [hg38.regions.bed](data/hg38.regions.bed) |
| Human | T2T-chm13v2.0 | [T2T-CHM13v2.0.sizes](data/T2T-CHM13v2.0.sizes) | [T2T-CHM13v2.0.regions.bed](data/T2T-CHM13v2.0.regions.bed) |

---

## ENCODE Blacklist

Set via `blacklist_bed` in config; both pre-built under `data/`.

| Blacklist | Species | Reference | Download |
|-----------|---------|-----------|----------|
| ENCODE blacklist v2 | Human | hg38 | [hg38-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz) |
| ENCODE blacklist v2 | Mouse | mm10 | [mm10-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz) |

---

## Sequencing bias correction (Optional)
### Mappability track

Optional; set `mappability_bed` in config to add the `MAP` column. Convert the bigWig to BED with `bigWigToBedGraph` (UCSC tools).

| Track | Species | Reference | Download |
|-------|---------|-----------|----------|
| k100 Umap multi-track | Human | hg38 | [k100.Umap.MultiTrackMappability.bw](http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw) |

### Replication timing

This is automatically handled by the pipeline via rule `repliseq_bigwig_to_bedgraph` and lift-over to `hg38` via [hg19ToHg38.over.chain.gz](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz) using `repliseq_liftover` if needed.

| Track | Species | Reference | Download |
|-------|---------|-----------|----------|
| ENCODE UW Repli-seq WaveSignal (16 bigWig) | Human | hg19 | [UCSC](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/) |

