# BCO Execution

Repository to share data and scripts working toward the goal of execution of BCOs
as CWL workflows.

## Usage Notes

This repository documents rough, in-progress, work toward creating an executable
BCO. To accomplish this, test workflows as well as test datasets are required.
Two workflows are (so far) included here. One is from the SilicoFCM project, and
is the intended final use case (bwa-gatk-annotation-bix2mus-workflow.cwl). The
other approximates many of the same steps, but is a public workflow available on
all SBG platforms (gatk-best-practice-data-pre-processing-4.1.0.0.cwl).

Begin by cloning the repo to your local machine.

### Software requirements

1. cwltool - The CWL reference runner available through pypi `pip install cwltool`
2. toil - An open-source pure-Python workflow engine available through pypi `pip install toil`
3. cwl extension for toil `pip install 'toil[cwl]'`
4. Docker - https://docs.docker.com/docker-for-mac/install/ (or linux)
5. gsutil - cli utility for interfacing with Google Cloud buckets `pip install gsutil` (to get test data)

### Getting the required test data

The test dataset is whole-genome sequencing supplied by the broad for testing
GATK best practice workflows. These are accessible through either AWS or GCP.

1. Navigate to the test_data directory in the repository and either run the bash
or execute the same command contained within using `gsutil` to download gatk test
data to your local machine.
2. Look in the file test_data directory for the file which outlines the required
files you'll need to download from the SB Platform. These are available as public
data files, as well as in the PDXNet Datapool.

### Getting the BCO R Shiny app

1. Navigate to the BCO App Github Page: https://github.com/sbg/bco-app
2. Clone the repo, or follow the directions for running through the Docker container

Using the Docker version is recommended unless doing development work on the BCO App itself.

### Goals

1. Identify test workflows
2. Identify test datasets
3. Create BCO from test workflows
4. Execute test workflow, or subset of steps from workflow, locally
5. Determine implementation details for executable BCO
6. Create executable BCO for test workflows

### Current status

#### 2021-06-02
- Test datasets located
- BCO App creates BCO from test workflows
- Test workflows not yet executing locally
