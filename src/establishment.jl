
### Set user-specific settings ###

# Set the repo directory. You may need to change this manually in
# some contexts
basedir = joinpath(@__DIR__, "..")
smapfolder = "/home/raf/Data/SMAP_L4_SM_gph_v4"
# Only set to true if you have an Nvidia GPU and CUDA.jl properly set up.
use_cuda = true
# Set the days of the month to use for a quality/performance tradoff
# The final results use all days, but that may take a long time.
days = (1,15) # Use three days a month for performance
# days = 1:31 # Use all days
 

# Load required packages

using GrowthMaps, Unitful, Dates, Setfield, Statistics, Flatten
using GeoData, NCDatasets, HDF5
using Unitful: °C, K, cal, mol
# Optionally load CUDA
if use_cuda 
    using CUDA
    CUDA.allowscalar(false)
    arraytype = CuArray
else
    arraytype = Array
end


### Define model components ###

# Growth model based on Schoolfield (1981).

p = 3.626804f-01
ΔH_A = 3.586625f4cal/mol
ΔH_H = 1.431237f5cal/mol
Thalf_H = 2.988454f2K
ΔH_L = -1.108988f5cal/mol
Thalf_L = 2.459949f2K
T_ref = K(25.0f0°C)
growthresponse = Layer(:surface_temp, K,
    SchoolfieldIntrinsicGrowth(p, ΔH_A, ΔH_L, Thalf_L, ΔH_H, Thalf_H, T_ref)
)
# Optionally plot the growth rate curve:
# using Plots, UnitfulRecipes
# temprange = collect(0.0:0.1:40.0)u"°C"
# plot(x -> GrowthMaps.conditionalrate(growthresponse, x), temprange)

coldthresh = -10.0f0°C |> K  # Stephens2015
coldmort = -log(1.23f0)K^-1 # Stephens2015
coldstress = Layer(:surface_temp, K, ColdStress(coldthresh, coldmort))
heatthresh = 30.0f0°C |> K # Enriquez2017
heatmort = -log(1.15f0)K^-1 # Enriquez2017
heatstress = Layer(:surface_temp, K, HeatStress(heatthresh, heatmort))
wiltthresh = 0.5f0
wiltmort = -log(1.1f0)
wiltstress = Layer(:land_fraction_wilting, WiltStress(wiltthresh, wiltmort));
wetthresh = 0.8f0
wetmort = -10.0f0
wetstress = Layer(:sm_surface_wetness, WetStress(wetthresh, wetmort));


### Set the simulation time span ###

tspan = DateTime(2017,1):Month(1):DateTime(2017,12)


### Load the SMAP dataset ###

filter = Where(t -> dayofmonth(t) in days && t >= first(tspan) && t < last(tspan) + Month(1))
smapseries = SMAPseries(smapfolder)[filter]



### Run main model and sensitivity analysis ####

basemodel = growthresponse, wiltstress, wetstress, coldstress, heatstress;

series = smapseries
outputs = Dict()
stressors = (heatstress=heatstress, coldstress=coldstress, wiltstress=wiltstress, wetstress=wetstress)
stressorkeys = (:heatstress, :coldstress, :wiltstress, :wetstress)
s80s = (Symbol(s, "_", 80) => (s, 0.8f0) for s in stressorkeys)  
s120s = (Symbol(s, "_", 120) => (s, 1.2f0) for s in stressorkeys)  
scenariokeys  = (; s80s..., s120s...)
sensitivityseries = smapseries
stressscenarios = map(scenariokeys) do (key, change)
    stressor = stressors[key]
    i = findfirst(k -> k == key, keys(stressors))
    altered_stressor = deepcopy(stressor)
    println("$key $val with newval: ", altered_stressor.model.mortalityrate * change)
    @set! altered_stressor.model.mortalityrate = altered_stressor.model.mortalityrate * change
    altered_model = (growthresponse, values(stressors)[1:i-1]..., altered_stressor, values(stressors)[i+1:end]...)
end
scenarios = (; standard=basemodel, stressscenarios...)

@time outputs = mapgrowth(
     scenarios;
     series=smapseries,
     tspan=tspan,
     arraytype=CuArray,
);

# Save NetCDF files for sensitivity analysys

sensitivitydir = joinpath(basedir, "output/sensitivity")
mkpath(sensitivitydir)
for scenario in keys(outputs)
    output = GeoArray(outputs[scenario]; name=:growthrates)
    write(joinpath(sensitivitydir, "growthrates_$(scenario).ncd"), NCDarray, output)
end

# Optionall plot the stored result
# using Plots
# NCDarray(joinpath(sensitivitydir, "growthrates_standard.ncd"))[Ti(1)] |> plot
# NCDarray(joinpath(sensitivitydir, "growthrates_heatstress_80.ncd"))[Ti(1)] |> plot
# NCDarray(joinpath(sensitivitydir, "growthrates_heatstress_120.ncd"))[Ti(1)] |> plot
