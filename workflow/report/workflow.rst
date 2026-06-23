Universal Genotyping Pipeline
=============================

SNP genotyping, phasing, and allele/depth counting across bulk WGS/WES, scRNA,
scATAC, and VISIUM assays. Outputs (binned allele counts, depth, RDR / CNV
segments) feed HATCHet3 and CalicoST.

This report collects per-stage **QC statistics** (SNP/phasing/depth tables) and
**QC plots** (allele-frequency, read-depth correction, and segmentation figures)
for ``sample_id = {{ snakemake.config["sample_id"] }}`` run in
``workflow_mode = {{ snakemake.config["workflow_mode"] }}``.

Use the *Statistics* tab for runtime/provenance and the *Results* categories
(``QC stats``, ``QC plots``) for the generated artifacts.
