# External tools and attribution

| Tool | Role | Setup |
|---|---|---|
| [CrispHunter](https://github.com/2021JohnSheng/CrispHunter) | Workflow inspiration | Reference project |
| CaSilico | Local Cas12a 20 nt crRNA generation and structural checks | Existing R installation, configurable Rscript/library paths |
| [EasyDesign](https://github.com/scRNA-Compt/EasyDesign) | Official Cas12a 21 nt designs and activity prediction | Separate Python 3.8 / TensorFlow 2.3 environment and official model files |
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | Nuclear 18S alignment | WSL Conda environment, tested 7.526 |
| [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) | Local alignment, databases and low-complexity screening | WSL Conda environment, tested 2.16 |
| [BRAKER4](https://github.com/Gaius-Augustus/BRAKER4) | Optional genome annotation | Official Snakemake workflow, Apptainer and required evidence/databases |
| [OrthoFinder](https://github.com/davidemms/OrthoFinder) | Optional protein orthogroup identification | Shared-gene mode with multiple assemblies |
| [XlsxWriter](https://xlsxwriter.readthedocs.io/) | Portable Excel export | Main Python requirements |
| [ViennaRNA/RNAfold](https://www.tbi.univie.ac.at/RNA/) | No-PAM spacer RNA self-folding | Separate structure environment, tested 2.8.0 |
| NumPy | Packed short-guide background screening | Windows and WSL Python environments |

Install dependencies from their upstream distributions and observe their licenses and citation guidance. Model weights, third-party source trees, binaries, R libraries, containers, environments, downloaded biological datasets and local analysis results are external to this repository.

The `environment.yml` supplies the primary 18S Linux tools. Optional BRAKER4, Snakemake, Apptainer and OrthoFinder must be installed separately according to the selected workflow. EasyDesign requires its own compatible legacy environment and is not installed into the main Python 3.11 controller environment.
