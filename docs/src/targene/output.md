# Outputs

All outputs are generated in the `$(OUTDIR)` (default: `results`) directory. Here we succintly describe the most important ones:

- `summary.csv`: is the main output of the pipeline, it contains all summary statistics and information for each estimand of interest. Those are further described [here](https://targene.github.io/TargetedEstimation.jl/stable/tmle_estimation/#Output-file).
- `hdf5files/inf_curves`: contains influence curves in HDF5 format if those were requested (see the `SAVE_IC` Nextflow parameter).
- `tmle_inputs/final.data.csv`: contains the input dataset to all TMLE processes.
- Other sub-directories contain intermediate results that may still be of interest for debugging purposes.
