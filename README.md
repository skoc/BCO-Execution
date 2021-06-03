# BCO Execution

Repository to share data and scripts working toward the goal of execution of BCOs
as CWL workflows.

![](https://img.shields.io/badge/Status-under--dev-red.svg) ![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)![License: AGPL-3.0](https://img.shields.io/github/license/skoc/BCO-Execution.svg)

<img src="https://raw.githubusercontent.com/skoc/BCO-Execution/master/img/bco_execution_flow.png" align="center" alt="summary" />

# Usage Notes

This repository documents rough, in-progress, work toward creating an executable
BCO. To accomplish this, test workflows as well as test datasets are required.
Two workflows are (so far) included here. 

Begin by cloning the repo to your local machine.

## Software requirements

### BCO APP Installation

#### Getting the BCO R Shiny app

1. Navigate to the BCO App Github Page: https://github.com/sbg/bco-app
2. Clone the repo, or follow the directions for running through the Docker container

Using the Docker version is recommended unless doing development work on the BCO App itself.

### BCO Runner Script Installation

#### To install:

Clone the repository and install the enviroment:

`git clone https://github.com/skoc/BCO-Execution.git`

`cd BCO-Execution/bco_runner`

Create a conda environment from it as follows:

`conda env create -f environment.yml`

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

#### Run CWL of BCO via run_cwl
```
python bco_runner.py cwl_runner -b <BCO_FILE>

```

## Getting the required test data

The test dataset is whole-genome sequencing supplied by the broad for testing
GATK best practice workflows. These are accessible through either AWS or GCP.

1. Navigate to the test_data directory in the repository and either run the bash
or execute the same command contained within using `gsutil` to download gatk test
data to your local machine.
2. Look in the file test_data directory for the file which outlines the required
files you'll need to download from the SB Platform. These are available as public
data files, as well as in the PDXNet Datapool.


# References

[bcotool by Hadley King](https://github.com/HadleyKing/bcotool)

[BioCompute Objects](https://biocomputeobject.org/)
