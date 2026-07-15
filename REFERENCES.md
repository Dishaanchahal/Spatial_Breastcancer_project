# References

Methods and biological references underpinning the analyses in this repository.

## Dataset
- Wu SZ, Al-Eryani G, Roden DL, et al. A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics* 2021;53:1334–1347. doi:10.1038/s41588-021-00911-1 (GEO: GSE176078)

## Deconvolution & integration
- Kleshchevnikov V, Shmatko A, Dann E, et al. Cell2location maps fine-grained cell types in spatial transcriptomics. *Nature Biotechnology* 2022;40:661–671. doi:10.1038/s41587-021-01139-4
- Cable DM, Murray E, Zou LS, et al. Robust decomposition of cell type mixtures in spatial transcriptomics (RCTD/spacexr). *Nature Biotechnology* 2022;40:517–526. doi:10.1038/s41587-021-00830-w
- Hao Y, Hao S, Andersen-Nissen E, et al. Integrated analysis of multimodal single-cell data (Seurat v4). *Cell* 2021;184:3573–3587. doi:10.1016/j.cell.2021.04.048
- Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology* 2018;19:15. doi:10.1186/s13059-017-1382-0
- Li H, et al. A comprehensive benchmarking with practical guidelines for cellular deconvolution of spatial transcriptomics. *Nature Communications* 2023;14:1548.
- Li B, et al. Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution. *Nature Methods* 2022;19:662–670.
- Sang-aram C, et al. Spotless: a reproducible pipeline for benchmarking cell type deconvolution in spatial transcriptomics. *eLife* 2024;12:RP88431.

## Cell–cell communication
- Dimitrov D, Türei D, Garrido-Rodriguez M, et al. Comparison of methods and resources for cell–cell communication inference from single-cell RNA-Seq data (LIANA). *Nature Communications* 2022;13:3224. doi:10.1038/s41467-022-30755-0
- Jin S, Guerrero-Juarez CF, Zhang L, et al. Inference and analysis of cell–cell communication using CellChat. *Nature Communications* 2021;12:1088.
- Cang Z, Zhao Y, Almet AA, et al. Screening cell–cell communication in spatial transcriptomics via collective optimal transport (COMMOT). *Nature Methods* 2023;20:218–228. doi:10.1038/s41592-022-01728-4

## Spatial niches & tissue architecture
- Palla G, Spitzer H, Klein M, et al. Squidpy: a scalable framework for spatial omics analysis. *Nature Methods* 2022;19:171–178. doi:10.1038/s41592-021-01358-2
- Goltsev Y, Samusik N, Kennedy-Darling J, et al. Deep profiling of mouse splenic architecture with CODEX multiplexed imaging (cellular neighborhoods). *Cell* 2018;174:968–981. doi:10.1016/j.cell.2018.07.010
- Schürch CM, Bhate SS, Barlow GL, et al. Coordinated cellular neighborhoods orchestrate antitumoral immunity at the colorectal cancer invasive front. *Cell* 2020;182:1341–1359. doi:10.1016/j.cell.2020.07.005
- Janesick A, Shelansky R, Gottscho AD, et al. High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis. *Nature Communications* 2023;14:8353.
- Oliveira MF, et al. High-definition spatial transcriptomic profiling of immune cell populations in colorectal cancer. *Nature Genetics* 2025.
- Wang X, et al. Spatial transcriptomics reveals substantial heterogeneity in triple-negative breast cancer. *Nature Communications* 2024.

## Pathway / functional activity
- Schubert M, Klinger B, Klünemann M, et al. Perturbation-response genes reveal signaling footprints in cancer gene expression (PROGENy). *Nature Communications* 2018;9:20. doi:10.1038/s41467-017-02391-6
- Badia-i-Mompel P, Vélez Santiago J, Braunger J, et al. decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinformatics Advances* 2022;2:vbac016. doi:10.1093/bioadv/vbac016

## Tertiary lymphoid structures (TLS)
- Cabrita R, Lauss M, Sanna A, et al. Tertiary lymphoid structures improve immunotherapy and survival in melanoma (9-gene TLS signature). *Nature* 2020;577:561–565. doi:10.1038/s41586-019-1914-8
- Helmink BA, Reddy SM, Gao J, et al. B cells and tertiary lymphoid structures promote immunotherapy response. *Nature* 2020;577:549–555. doi:10.1038/s41586-019-1922-8

## Tumor microenvironment & CAF biology
- Costa A, Kieffer Y, Scholer-Dahirel A, et al. Fibroblast heterogeneity and immunosuppressive environment in human breast cancer (CAF-S1/CXCL12). *Cancer Cell* 2018;33:463–479. doi:10.1016/j.ccell.2018.01.011
- Jenkins L, Jungwirth U, Avgustinova A, et al. Cancer-associated fibroblasts suppress CD8+ T cell infiltration and confer resistance to immune checkpoint blockade. *Cancer Research* 2022;82:2904–2917.
- Tharp KM, et al. Tumor-associated macrophages restrict CD8+ T cell function through collagen deposition and metabolic reprogramming. *Nature Cancer* 2024.
- Sabit H, et al. The role of tumor microenvironment and immune cell crosstalk in triple-negative breast cancer. *Cancer Letters* 2025.

## Clinical context (TILs / immunotherapy)
- Wood SJ, et al. High tumor infiltrating lymphocytes are associated with pathological complete response in TNBC treated with neoadjuvant KEYNOTE-522 chemoimmunotherapy. *Breast Cancer Research and Treatment* 2024.
- Barroso-Sousa R, et al. Prediction of pathologic complete response to chemoimmunotherapy in TNBC using tumor-infiltrating lymphocytes. *Journal of Clinical Oncology* 2024.

## External validation dataset
- 10x Genomics. Human Breast Cancer (Block A, Sections 1–2), Visium Spatial Gene Expression, v1.1.0 (invasive ductal carcinoma). https://www.10xgenomics.com/resources/datasets/human-breast-cancer-block-a-section-1-1-standard-1-1-0

## Clinical cohort & survival analysis
- The Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. *Nature* 2012;490:61–70. doi:10.1038/nature11412
- Goldman MJ, Craft B, Hastie M, et al. Visualizing and interpreting cancer genomics data via the Xena platform. *Nature Biotechnology* 2020;38:675–678. doi:10.1038/s41587-020-0546-8 (data source for TCGA-BRCA HiSeqV2 + survival)
- Davidson-Pilon C. lifelines: survival analysis in Python. *Journal of Open Source Software* 2019;4(40):1317. doi:10.21105/joss.01317

## METABRIC validation cohort
- Curtis C, Shah SP, Chin SF, et al. The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature* 2012;486:346–352. doi:10.1038/nature10983
- Pereira B, Chin SF, Rueda OM, et al. The somatic mutation profiles of 2,433 breast cancers refine their genomic and transcriptomic landscapes. *Nature Communications* 2016;7:11479. doi:10.1038/ncomms11479
- METABRIC expression + clinical obtained via a public GitHub mirror of the cBioPortal `brca_metabric` study (cBioPortal datahub was unreachable at analysis time).
