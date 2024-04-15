![workflow](https://github.com/naturalis/barcode-constrained-phylogeny/actions/workflows/python-package-conda.yml/badge.svg)
[![License: Apache-2.0](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10519081.svg)](https://doi.org/10.5281/zenodo.10519081)

![Logo](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/doc/logo-small.png?raw=true)

# Bactria: BarCode TRee Inference and Analysis
This repository contains code and data for building very large, topologically-constrained 
barcode phylogenies through a divide-and-conquer strategy. Such trees are useful as 
reference materials for curating barcode data by detecting rogue terminals (indicating
incorrect taxonomic annotation) and in the comparable calculation of alpha and beta 
biodiversity metrics across metabarcoding assays. 

The input data for the approach we develop here currently comes from BOLD data dumps. 
The international database [BOLD Systems](https://www.boldsystems.org/index.php) 
contains DNA barcodes for hundreds of thousands of species, with multiple barcodes per 
species. The data dumps we use here are TSV files whose columns conform to the nascent
BCDM (barcode data model) vocabulary. As such, other data sources that conform to this
vocabulary could in the future be used as well, such as [UNITE](https://unite.ut.ee/).

Theoretically, such data could be filtered and aligned per DNA marker to make 
phylogenetic trees. However, there are two limiting factors: building very large 
phylogenies is computationally intensive, and barcodes are not considered ideal for 
building big trees because they are short (providing insufficient signal to resolve large 
trees) and because they tend to saturate across large patristic distances.

![concept](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/doc/concept.png)

Both problems can be mitigated by using the 
[Open Tree of Life](https://tree.opentreeoflife.org/opentree/argus/opentree13.4@ott93302) 
as a further source of phylogenetic signal. The BOLD data can be split into chunks that 
correspond to Open Tree of Life clades. These chunks can be made into alignments and 
subtrees. The OpenTOL can be used as a constraint in the algorithms to make these. The 
chunks are then combined in a large synthesis by grafting them on a backbone made from 
exemplar taxa from the subtrees. Here too, the OpenTOL is a source of phylogenetic 
constraint.

In this repository this concept is developed for the COI-5P marker, but the aim is to 
achieve equivalent functionality for plant barcoding markers (matK, rbcL) and for some
part of the ITS region.

## Installation

The pipeline and its dependencies are managed using conda. On a linux or osx system, you 
can follow these steps to set up the `bactria` Conda environment using an `environment.yml` 
file and a `requirements.txt` file:

1. **Clone the Repository:**  
   Clone the repository containing the environment files to your local machine:
   ```bash
   git clone https://github.com/naturalis/barcode-constrained-phylogeny.git
   cd barcode-constrained-phylogeny
   ```
2. **Create the Conda Environment:**
   Create the bactria Conda/Mamba environment using the environment.yml file with the following 
   command:
   ```bash
   mamba env create -f workflow/envs/environment.yml
   ```
   This command will create a new Conda environment named bactria with the packages 
   specified in the environment.yml file. This step is largely a placeholder because
   most of the dependency management is handled at the level of individual pipeline
   steps, which each have their own environment specification.
3. **Activate the Environment:**
   After creating the environment, activate it using the conda activate command:
   ```bash
   mamba activate bactria
   ```
4. **Verify the Environment:**
   Verify that the bactria environment was set up correctly and that all packages were 
   installed using the conda list command:
   ```bash
   mamba list
   ```
   This command will list all packages installed in the active conda environment. You should 
   see all the packages specified in the environment.yml file and the requirements.txt file.

It is recommended that the mamba environment is configured to use strict channel priorities. 
This is crucial for having robust and correct environments (for details, see 
[here](https://conda-forge.org/docs/user/tipsandtricks.html)). Consider configuring strict 
priorities by executing `conda config --set channel_priority strict`.

## How to run

The pipeline is implemented using snakemake, which is available within the conda 
environment that results from the installation. Important before running the snakemake pipeline 
is to change in [config/config.yaml](config/config.yaml) the number of threads available on your 
computer. Which marker gene is used in the pipeline is also specified in the config.yaml (default 
COI-5P). Prior to execution, the BOLD data package (we used the 
[release of 18 December 2023](https://www.boldsystems.org/index.php/datapackage?id=BOLD_Public.18-Dec-2023))
must be downloaded manually and stored in the [resources/](resources/) directory.
The same goes for the OpenTOL datapackage (we used the 
[13.4 release](https://tree.opentreeoflife.org/about/synthesis-release/v13.4))
If a BOLD or OpenTOL release from another date is used the file names in config.yaml need to be updated.

How to run the entire pipeline:

```bash 
snakemake -j {number of threads} --use-conda
```

Snakemake rules can be performed separately:
```bash 
snakemake -R {Rule} -j {number of threads} --use-conda
```

Enter the same number at {number of threads} as you filled in previously in src/config.yaml.
In {Rule} insert the rule to be performed.

Here is an overview of all the rules in the Snakefile:

![graphviz (1)](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/doc/dag.svg)
(zoomed view is available [here](https://raw.githubusercontent.com/naturalis/barcode-constrained-phylogeny/main/doc/dag.svg))

More detailed documentation of the individual rules is provided [here](workflow/documentation.md).

## Repository layout

Below is the top-level layout of the repository. This layout is in line with 
[community standards](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) and must be adhered to.
All of these subfolders contains further explanatory READMEs to explain their contents in more detail.

- [config](config/) - configuration files
- [doc](doc/) - documentation and background literature
- [logs](logs/) - where log files are written during pipeline runtime
- [resources](resources/) - external data resources (from BOLD and OpenTree) are downloaded here
- [results](results/) - intermediate and final results are generated here
- [workflow](workflow/) - script source code and driver snakefile 

## License

&copy; 2023-2024 Naturalis Biodiversity Center

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except 
in compliance with the License. You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)
   
Unless required by applicable law or agreed to in writing, software distributed under the License 
is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express 
or implied. See the License for the specific language governing permissions and limitations under 
the License.
