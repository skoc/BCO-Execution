# BCO Runner Script

## To install:

Clone the repository and install the enviroment:

`git clone https://github.com/skoc/BCO-Execution.git`

`cd BCO-Execution/bco_runner`

Create a conda environment from it as follows:

`conda env create -f environment.yml`

`conda activate bco`

## Supported commands:

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

### Validate BCO
```
python bco_runner.py validate -b <BCO_FILE>
``` 

### Run CWL of BCO via run_cwl
```
python bco_runner.py cwl_runner -b <BCO_FILE>

```

## Reference

[bcotool by Hadley King](https://github.com/HadleyKing/bcotool)

[BioCompute Objects](https://biocomputeobject.org/)