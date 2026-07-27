# Universal Genotyping Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake->=9.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](VERSION)

Universal Genotyping Pipeline is a Snakemake preprocessing pipeline for downstream allele-specific CNA inference softwares including:
- [HATCHet](https://github.com/raphael-group/hatchet) - bulk **short-read** WGS/WES, **long-read** PacBio HiFi/ONT,
- [Copy-typing](https://github.com/raphael-group/Copy-typing) - scRNA-seq, scATAC-seq,
- [CalicoST](https://github.com/raphael-group/CalicoST) - Visium ST.

---

## Documentation

| Document | Description |
|----------|-------------|
| [docs/pipeline.md](docs/pipeline.md) | Install, configure, and run the Snakemake pipeline. |
| [docs/bulk_genotyping.md](docs/bulk_genotyping.md) | Bulk genotyping workflow |
| [docs/single_cell_genotyping.md](docs/single_cell_genotyping.md) | Single-cell/spatial genotyping workflow |
| [docs/copytyping_preprocess.md](docs/copytyping_preprocess.md) | Copy-typing preprocessing workflow |
| [docs/workflow.md](docs/workflow.md) | Detailed description per workflow |
| [docs/sample_sheet.md](docs/sample_sheet.md) | Sample sheet specification |
| [docs/reference.md](docs/reference.md) | Manual for input, output, and (hyper-)parameters |
| [config/config.yaml](config/config.yaml) | Default configuration (auto-loaded by the Snakefile) |
| [resources/templates/](resources/templates/) | User config and sample-file templates |
| [resources/README.md](resources/README.md) | External data resources |
| [CHANGELOG.md](CHANGELOG.md) | Release notes and version history |
