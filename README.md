# hotpot
Collection of random helper functions and scripts, mostly for image analysis. 

## Dependencies

### For snakemake pipelines

1. Install `snakemake`. See [Snakemake's official documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
for details.

2. Have snakemake create the conda environments for you on a node with 
   internet access.
   
   ```
   snakemake -s {path/to/Snakefile} \
       --use-conda \
       --conda-prefix {path/to/automatically/created/envs} \
       -c1 all_envs
   ```

   We don't use `--conda-create-envs-only` because we want to install 
   all dependencies ahead of time, independent of if certain rules 
   have become discoverable (e.g., if they are after a checkpoint).