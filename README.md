# Predicting the global invasion of Drosophila suzukii to improve biosecurity preparedness

Scripts for the reserach paper in the Jounal of Applied Ecology

James Maino, Rafael Schouten

Scripts in `src` folder:

To set up the scripts here:

1. Install julia 1.5

2. run julia with 

```bash
julia --project='/path/to/this/repo'
```

3. Then in the julia REPL run:

```julia
] instantiate
```

Where `]` just gets you into Pkg mode

## Establishment

## Get the data

- Download the SMAP L4_SM_gph_v4 dataset from
  https://nsidc.org/data/smap/smap-data.html for all of 2017. It will be a few
  hundred GB
- If you can, store it on a SD hard disk. Disk reads are lot of the load time.

### Using a GPU

This greatly reduces the time for running sensitivity analysis

- If you have an Nvidia GPU, make sure it's set up and manually install
  CUDA.jl and its system dependencies following the installation 
  [documentaion](https://juliagpu.gitlab.io/CUDA.jl/installation/overview/)

### Run the script

- establishment.jl 

This builds growth rates datasets using GrowthMaps.jl, including sensitivity
analysis. You can run it as a script or step through it, and plot some of the 
outputs (plotting is commented out for use as a script). 

There are some user-specific settings in the first few lines of the file, which
you will need to edit.

## Spread Optimisation, validation and sensitivity

These scripts can be run after establishment.jl generates growth-rate datasets.

At the top of the script, set your working directory (where this file is),
set the `humanactity` type to `:gdp` or `:pop`.

Then you can run the whole file. By default it will load the pre-computed
optimisation results, as optimising all the models will take a number of days on
a desktop computer. The whole file may take an hour to run.

## Spread predictions and plots

- spread.jl generates spread outputs used in the paper
- again, there are settings at the top of the file
- the script also contains a live app interface (ElectronOutput) you can uncomment to use
  to explore the model
