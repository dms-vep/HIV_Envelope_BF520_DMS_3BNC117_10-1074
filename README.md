# Deep mutational scanning using HIV Envelopes TRO11 and BF520 to identify mutations that escape neutralization by antibodies 3BNC117 and 10-1074

For documentation of the analysis, see [https://dms-vep.github.io/HIV_Envelope_BF520_DMS_3BNC117_10-1074/](https://dms-vep.github.io/HIV_Envelope_BF520_DMS_3BNC117_10-1074/) and [https://dms-vep.github.io/HIV_Envelope_TRO11_DMS_3BNC117_10-1074/](https://dms-vep.github.io/HIV_Envelope_TRO11_DMS_3BNC117_10-1074/). 

## Organization of this repo

This repository is meant to be run in the same directory as [https://dms-vep.github.io/HIV_Envelope_TRO11_DMS_3BNC117_10-1074/](https://dms-vep.github.io/HIV_Envelope_TRO11_DMS_3BNC117_10-1074/) . This analysis should be run by following the instructions in this README before running the analysis in the other repository.

### `dms-vep-pipeline-3` submodule

Most of the analysis is done by the [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline-3), which was added as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to this pipeline via:

    git submodule add https://github.com/dms-vep/dms-vep-pipeline-3

This added the file [.gitmodules](.gitmodules) and the submodule [dms-vep-pipeline-3](dms-vep-pipeline-3), which was then committed to the repo.
Note that if you want a specific commit or tag of [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) or to update to a new commit, follow the [steps here](https://stackoverflow.com/a/10916398), basically:

    cd dms-vep-pipeline-3
    git checkout <commit>

and then `cd ../` back to the top-level directory, and add and commit the updated `dms-vep-pipeline` submodule.
You can also make changes to the [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) that you commit back to that repo.

### Code and configuration
The [snakemake](https://snakemake.readthedocs.io/) pipeline itself is run by the [Snakefile](Snakefile), which includes [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) reads its configuration from [config.yaml](config.yaml). There are a few additional configuration files in the [./data/](data) folder. 
The [conda](https://docs.conda.io/) environment used by the pipeline is that specified in the `environment.yml` file in [dms-vep-pipeline-3](dms-vep-pipeline-3).

Additional scripts and notebooks that are specific to this analysis and not part of [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) are in [./scripts/](scripts) and [./notebooks/](notebooks). These analysis are included in the snakemake pipeline by adding them to [custom_rules.smk](custom_rules.smk).

### Input data
Input data for the pipeline are in [./data/](data).

### Results and documentation

The results of running the pipeline are placed in [./results/](results).
Only some of these results are tracked to save space (see [.gitignore](.gitignore)).

The pipeline builds HTML documentation for the pipeline in [./docs/](docs), which is rendered via GitHub Pages at [https://dms-vep.github.io/HIV_Envelope_BF520_DMS_3BNC117_10-1074/](https://dms-vep.github.io/HIV_Envelope_BF520_DMS_3BNC117_10-1074/).

## Running the pipeline
To run the pipeline, build the conda environment `dms-vep-pipeline-3` in the `environment.yml` file of [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3), activate it, and run [snakemake](https://snakemake.readthedocs.io/), such as:

    conda activate dms-vep-pipeline
    snakemake -c 1 -j 1 -s dms-vep-pipeline-3/Snakefile --use-conda --rerun-incomplete

To run on the Hutch cluster via [slurm](https://slurm.schedmd.com/), you can run the file [run_Hutch_cluster.bash](run_Hutch_cluster.bash):

    sbatch run_Hutch_cluster.bash 