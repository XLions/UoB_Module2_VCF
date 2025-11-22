# VCF Trio Challenges
Group 13, **Task E - Delly 3**, Module 2 Genomics & Next Generation Sequencing Group Task, *MSc Bioinfomatics, University of Birmingham*. The task is provided by Dr.Deena.
## Task Description
<ol>
<li>Summarize the various combination of variants, their location, type and function.</li>
<li>Identify de novo mutation in the trio you are assigned, if any.</li>
<li>Annotate the variants of interest and explain their importance.</li>
</ol>

## R Packages Used
| Package | Version | Usage |
|------|------|------|
| [stringr](https://cran.r-project.org/web/packages/stringr/index.html) | 1.5.1 | Deal with characters. |
| [progress](https://cran.r-project.org/web/packages/progress/index.html) | 1.2.3 | Visualize the process of loop. |
| [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html) | 2.0.0 | Load many packages. |
| [ggrepel](https://cran.r-project.org/web/packages/ggrepel/index.html) | 0.9.6 | Make labels in plot look better. |
| [RCircos](https://cran.r-project.org/web/packages/RCircos/index.html) | 1.2.2 | Creat the plot of chromsome locations. |
| [rtracklayer](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html) | 1.68.0 | Creat the plot of chromsome locations. |
| [dplyr](https://dplyr.tidyverse.org/) | 1.1.4 | Pipe operation. |
| [ggplot2](https://ggplot2.tidyverse.org/) | 3.5.2 | Create and save plots. |

## Repository Structure
```bash
UoB_Module2_VCF/
├── Data/                   # Input data
│   ├── Homo_sapiens.GRCh38.114.chr.gtf.gz
│   ├── gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz
│   ├── VCF_SplitedDetailedLongerDataFrame.RDS
│   ├── vcf.RDS
│   └── DellyVariation.vcf  # Example VCF file
├── Figs/
│   ├── A_PrimaryExploration/
│   │   └── ...
│   ├── B_FindCild/
│   │   ├── 1.Sample_Freq.pdf
│   │   ├── 2.DeNovoMutations.csv
│   │   └── FindChildFrame.RDS
│   └── C_InterestedGenes
│       └── ...
├── R.R                     # R Script
└── README.md               # Project overview

