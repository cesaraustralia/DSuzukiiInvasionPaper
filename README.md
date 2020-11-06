# Predicting the global invasion of Drosophila suzukii to improve biosecurity preparedness

Scripts for the research paper in the Jounal of Applied Ecology,
By James Maino, Rafael Schouten and Paul Umina

Scripts and code by Rafael Schouten

Scripts in `src` folder:

To set up the scripts here:

1. Install julia 1.5

2. run julia with 

```bash
julia --project='/path/to/this/repo'
```

3. Then to install the required pacakge, in the julia REPL run:

```julia
] instantiate
```

After that, you should be able to run the scripts in this repository 
with exactly the same dependencies used.

The scripts included here can be run from an IDE like VScode or Atom,
or run directly from the command line (instructions below).


## Establishment

### Get the data

- Download the SMAP L4_SM_gph_v4 dataset from
  https://nsidc.org/data/smap/smap-data.html for all of 2017. It will be a few hundred GB
- If you can, store it on an SD or M.2 hard drive. Disk reads are lot of the load time.

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

To run, use the terminal command line:

```bash
SCRIPTPATH="/path/to/this/repo"
julia --project=$SCRIPTPATH $SCRIPTPATH/src/establishment.jl
```

Or step through the file interactively in Atom or VScode.

## Spread Optimisation, validation and sensitivity

These scripts can be run after establishment.jl generates growth-rate datasets.
All other input data is included in this repository.

At the top of the script, set the `humanactivity` type to `:gdp` or `:pop`
for to use either GDP or population density datasets to drive human dispersal
dynamics.

Then you can run the file. By default it will load the pre-computed
optimisation results, and run only a few replicates, as optimising all the 
models will take a number of days on a desktop computer. 
It may still take a while to run.

```bash
SCRIPTPATH="/path/to/this/repo"
julia --project=$SCRIPTPATH $SCRIPTPATH/src/spreadoptimisation.jl
```

## Spread predictions and inputs for plots

- spread.jl generates spread outputs used in the paper
- again, there are settings at the top of the file
- the script also contains a live app interface (ElectronOutput) can be
  uncommented to explore the model predictions for spread in Australia

```bash
SCRIPTPATH="/path/to/this/repo"
julia --project=$SCRIPTPATH $SCRIPTPATH/src/spread.jl
```

