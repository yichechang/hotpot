# hotpot
Collection of random helper functions and scripts, mostly for image analysis.

## Dependencies

### For snakemake pipelines

1. Install `snakemake`. See [Snakemake's official documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
   for details. Please install version 7.32.4 but not 8.0.0 or later,
   which should still be compatible with the pipelines in this repo but
   the instructions of how to run the pipelines are different.

2. Have snakemake create the conda environments for you on a node with
   internet access.

   ```
   snakemake -s {path/to/Snakefile} \
       --use-conda \
       --conda-prefix {path/to/automatically/created/envs} \
       -c1 all_envs
   ```

   We generally don't use `--conda-create-envs-only` when we want
   to install all dependencies ahead of time, independent of if certain
   rules have become discoverable (e.g., if they are after a
   checkpoint). But for simpler pipelines without checkpoints and if
   `all_env` is not provided, `--conda-create-envs-only` is a good
   option. Or, simply let snakemake create the conda environments on
   the fly if the pipeline is meant to be run entirely on a node with
   internet access.
