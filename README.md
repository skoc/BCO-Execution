# BCO Execution

Repository to share data, workflows, apps, and scripts working toward the goal of execution of BioCompute Objects (BCO).

![](https://img.shields.io/badge/Status-under--dev-red.svg) ![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg) ![License: AGPL-3.0](https://img.shields.io/github/license/skoc/BCO-Execution.svg)

<img src="https://raw.githubusercontent.com/skoc/BCO-Execution/master/img/bco_execution_flow.png" align="center" alt="summary" />

This repository also includes two example workflows that can be used to generate a BCO.
One, is a simple CWL tool definition using the common NGS quality control program
[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). The second
is a part of the GATK workflows, the Broad Best Practice Data Pre-processing Workflow.
This workflow is intended to prepare whole-exome or whole-genome sequencing data
for variant calling. The included workflow is publicly available from the Broad
and Seven Bridges' public Open Workflows repository [here](https://github.com/sevenbridges-openworkflows/Broad-Best-Practice-Data-pre-processing-CWL1.0-workflow).
It has been packed into a single .cwl file for convenience. The included .yml
files describe jobs that were validated through [cwltool](https://github.com/common-workflow-language/cwltool) and [toil](https://toil.readthedocs.io/en/latest/).

# Usage Notes

This repository include all the required apps/tools to crate BCO from scratch or CWL workflow, and provides an easy-to-use script to validate and excute given BCO with provided workflow language in its 'Execution Domain'.

Begin by cloning the repo to your local machine.

## Software requirements

### BCO APP Installation

#### Getting the BCO R Shiny app

1. Navigate to the BCO App Github Page: https://github.com/sbg/bco-app
2. Clone the repo, or follow the directions for running through the Docker container

Using the Docker version is recommended unless doing development work on the BCO App itself.

### BCO Runner Script Installation [[1](https://github.com/HadleyKing/bcotool)]

#### To install:

Clone the repository and install the enviroment:

`git clone https://github.com/skoc/BCO-Execution.git`

`cd BCO-Execution/bco_runner`

Create a conda environment from it as follows:

`conda env create -f environment.yml`

`conda activate bco`

#### Supported commands:

```
python bco_runner.py --help

positional arguments:
  {functions,license,validate,run_cwl}
    functions           list all available functions
    license             Prints BCO License
    validate            Validation options. Used to test a BCO against a JSON
                        schema. If no schema is supplied the ieee-2791-schema
                        is used as the default
    run_cwl             run a CWL

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

```

#### Validate BCO
```
python bco_runner.py validate -b <BCO_FILE>
``` 

#### Run CWL of BCO
```
python bco_runner.py run_cwl -b <BCO_FILE>

```


# References

1. [bcotool by Hadley King](https://github.com/HadleyKing/bcotool)

2. [BioCompute Objects](https://biocomputeobject.org/)
