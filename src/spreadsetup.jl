
### Load the required packages ###

using DimensionalData, GeoData, ArchGDAL, NCDatasets, HDF5
using CSV, DataFrames, DataFramesMeta, TableView
using DynamicGrids, Dispersal 
using Statistics, Dates, Unitful, Setfield, Flatten, FieldMetadata
using Plots, ColorSchemes, Colors
using DimensionalData: setdims, rebuild, Between

const DG = DynamicGrids


### Override default parameter range bounds ###

# For both interfaces and optimiser

import FieldMetadata: @bounds, bounds, @flattenable, flattenable
using Flatten: flatten

@bounds HumanDispersal begin
    human_exponent  | (0.0, 3.0)
    dist_exponent   | (0.0, 3.0)
    dispersalperpop | (0.0, 1e-8)
    max_dispersers  | (1e1, 1e5)
end
# Constant growth
@bounds ExactLogisticGrowth begin
    intrinsicrate | (0.0, 10.0)
end
# Alee
@bounds AlleeExtinction begin
    minfounders | (2e0, 1e3)
end

# Dont parametrise carrycap
@flattenable ExactLogisticGrowthMap begin
    carrycap | false
end
@flattenable ExactLogisticGrowth begin
    carrycap | false
end
# Don't parametrise λ for main models, it will just be the top bound 
# so we need to estimate it. This is reverted below for local-dispersal
# only models.
@flattenable ExponentialKernel begin
    λ | false
end



### Define simulation settings ###

timestep = Month(1)

us_tspan = DateTime(2008, 5):timestep:DateTime(2013, 12)
eu_tspan = DateTime(2008, 5):timestep:DateTime(2016, 12)

aus = Lon(Between(113.3402, 153.9523)), Lat(Between(-43.62234, -10.65125))
us = Lon(Between(-130.332, -66.9)), Lat(Between(18.977, 71.59))
eu = Lon(Between(-13.7, 40.15)), Lat(Between(34.9, 70.30))

aus_incursionpoints = (
    Melbourne=(-37.805896, 144.959527),
    Mildura=(-34.219504, 142.130864),
    Coffs_Harbour=(-30.287245, 153.092991),
    Sydney=(-33.839943, 151.006101),
    Adelaide=(-34.901608, 138.601547),
    Port_Augusta=(-32.466201, 137.813850),
    Devonport=(-41.180545, 146.314887),
    Hobart=(-42.881742, 147.323879),
    Brisbane=(-27.436190, 152.990588),
    Cairns=(-16.937281, 145.747709),
    Perth=(-31.9505, 115.8605),
    Geraldton=(-28.778138, 114.615632)
)

us_incursionpoints = (
    San_Francisco=(36.9453, -122.0695),
)

obs = CSV.File(joinpath(basedir, "data/Oersted_occurrence.csv")) |> DataFrame
eu_obs = @where(obs, :Year .!== missing, :Year .== 2008, :Country .!= "USA")
eu_incursionpoints = collect(zip(eu_obs.Latitude, eu_obs.Longitude))

incursionpointstatekeys = (
    Melbourne = (:Vic),
    Mildura	 = (:Vic),
    Coffs_Harbour =	(:NSW),
    Sydney = (:NSW),
    Adelaide = (:SA),
    Port_Augusta = (:SA),
    Devonport  = (:Tas),
    Hobart = (:Tas),
    Brisbane = (:QLD),
    Cairns = (:QLD),
    Perth =	(:WA),
    Geraldton =	(:WA)
)


### Load numbered states rasters ###

loadstates(path, selectors) = GDALarray(joinpath(basedir, path);
    name=:States, usercrs=EPSG(4326))[Band(1), selectors...] |>
    x->permutedims(x, (Lat, Lon)) |> x -> replace_missing(x, 0) |>
    x->trunc.(Int, x)
aus_states = loadstates(joinpath(basedir, "data/aus_states.tif"), aus)
us_states = loadstates(joinpath(basedir, "data/US_states.tif"), us)
eu_states = loadstates(joinpath(basedir, "data/EU_states.tif"), eu)
us_states = loadstates(joinpath(basedir, "data/US_states.tif"), ())
# aus_states |> plot
# eu_states |> plot
# us_states |> plot


### Define the RuleSets ###

growthratesfilepath = joinpath(basedir, "data/sensitivity/growthrates_standard.ncd")
isfile(growthratesfilepath)
growthrates = NCDarray(growthratesfilepath; crs=crs(aus_states), dimcrs=EPSG(4326)) |> GeoArray |>
    x -> setdims(x, (@set dims(x, Ti).mode.span = Regular(Month(1)))) |>
    x -> setdims(x, convertmode(Projected, dims(x, Lon))) |>
    x -> setdims(x, convertmode(Projected, dims(x, Lat))) |>
    x -> permutedims(x, (Lat, Lon, Ti)) |>
    x -> reverse(x; dims=Lat)


# Make growth rates aux grids, replacing the `missing` values with zero.
# Dispersal.jl `init` can't contain `missing` - it will spread everywhere.

grzm = floattype.(replace_missing(growthrates, 0.0))
aus_growthrates_zeromissing = grzm[aus...]
us_growthrates_zeromissing = grzm[us...]
eu_growthrates_zeromissing = grzm[eu...]
aus_growthrates_zeromissing[Ti(1)] |> plot
us_growthrates_zeromissing[Ti(1)] |> plot
eu_growthrates_zeromissing[Ti(1)] |> plot
aus_aux = (growthrates=aus_growthrates_zeromissing,)
us_aux = (growthrates=us_growthrates_zeromissing,)
eu_aux = (growthrates=eu_growthrates_zeromissing,)



### Define masking layers ###

# These define regions that are not run, specifically, the sea.

aus_boolmask = GeoData.boolmask(aus_states)
aus_missingmask = GeoData.missingmask(aus_states)
aus_missingmask = GeoData.missingmask(aus_states)
us_boolmask = GeoData.boolmask(growthrates[us..., Ti(1)])
us_missingmask = GeoData.missingmask(growthrates[us..., Ti(1)])
eu_boolmask = GeoData.boolmask(eu_states)
eu_missingmask = GeoData.missingmask(eu_states)
# plot(aus_boolmask)
# plot(aus_missingmask)
# plot(us_boolmask)
# plot(us_missingmask)
# plot(eu_boolmask)
# plot(eu_missingmask)

# Growth masks. These are the cells where the mean growth rate is above zero.

aus_growthmask = rebuild(aus_boolmask, (mean(aus_growthrates_zeromissing; dims=Ti)[Ti(1)] .* aus_boolmask) .> 0)
us_growthmask = rebuild(us_boolmask, (mean(us_growthrates_zeromissing; dims=Ti)[Ti(1)] .* us_boolmask) .> 0)
eu_growthmask = rebuild(eu_boolmask, (mean(eu_growthrates_zeromissing; dims=Ti)[Ti(1)] .* eu_boolmask) .> 0)
# aus_growthmask |> plot
# us_growthmask |> plot
# eu_growthmask |> plot



### Load human population

human_pop_path = joinpath(basedir, "data/population_density.tif")
human_pop = GDALarray(human_pop_path; usercrs=EPSG(4326))[Band(1)] |>
            permutedims |> replace_missing
# human_pop |> plot
aus_humanpop = human_pop[aus...]
us_humanpop = human_pop[us...]
eu_humanpop = human_pop[eu...]
# aus_humanpop |> plot
# us_humanpop |> plot
# eu_humanpop |> plot

gdp_path = joinpath(basedir, "data/GDP_USD_per_m2.tif")
gdp = GDALarray(gdp_path; usercrs=EPSG(4326))[Band(1)] .* 1e6 |> permutedims
gdp = (@set gdp.missingval = minimum(gdp)) |> replace_missing
# gdp |> plot
aus_gdp = gdp[aus...]
us_gdp = gdp[us...]
eu_gdp = gdp[eu...]
# aus_gdp |> plot
# us_gdp |> plot
# eu_gdp |> plot


@assert size(eu_humanpop) == size(eu_growthmask[Ti(1)]) == size(eu_states)
@assert size(us_humanpop) == size(us_growthmask[Ti(1)]) == size(us_states)
@assert size(aus_humanpop) == size(aus_growthmask[Ti(1)]) == size(aus_states)




### Define Rules ###

carrycap = floattype(1e9)
growth = ExactLogisticGrowthMap{:population,:population}(
    layerkey=Val(:growthrates),
    carrycap=carrycap,
    timestep=Day(1),
    nsteps=floattype(1.0),
);
growth

### Constant growth ###

constant_growth = ExactLogisticGrowth{:population,:population}(
    intrinsicrate=floattype(0.1),
    timestep=Day(1),
    carrycap=carrycap
)

### Local dispersal ###

λ = floattype(0.0125)
radius = 1

@time hood = DispersalKernel{radius}(
    formulation=ExponentialKernel(λ),
    distancemethod=AreaToArea(30),
    cellsize=floattype(1.0)
)
localdisp = InwardsPopulationDispersal{:population,:population}(hood)
log.(hood.kernel) |> heatmap
savefig(joinpath(basedir, "output/dispersal_kernel.png"))

### Allee effects

allee = AlleeExtinction{:population,:population}(
    minfounders=floattype(22.0)
);


### Human Dispersal ###

scale = 10
hkwargs = (
    read=:population,
    scale=scale,
    human_exponent=floattype(1.075),
    dist_exponent=floattype(1.429),
    dispersalperpop=floattype(8.703e-9),
    max_dispersers=floattype(3.264e4),
    nshortlisted=100,
)

if humanactivity == :pop
    println("human activity from population")
    aus_humandisp = HumanDispersal(; human_pop=parent(aus_humanpop), hkwargs...)
    us_humandisp = HumanDispersal(; human_pop=parent(us_humanpop), hkwargs...)
    eu_humandisp = HumanDispersal(; human_pop=parent(eu_humanpop), hkwargs...)
elseif humanactivity == :gdp
    println("human activity from gdp")
    aus_humandisp = HumanDispersal(; human_pop=parent(aus_gdp), hkwargs...)
    us_humandisp = HumanDispersal(; human_pop=parent(us_gdp), hkwargs...)
    eu_humandisp = HumanDispersal(; human_pop=parent(eu_gdp), hkwargs...)
end

I = (DimensionalData.dims2indices(aus_states, (Lat(Contains(-36.805896)), Lon(Contains(143.959527)))) .-1) .÷ scale .+ 1
melb_dests = Dispersal.populate!(zero(aus_growthrates_zeromissing[Ti(1)]), aus_humandisp, I...)
plot(melb_dests, title="Melbourne destinations")
savefig(joinpath(basedir, "output/human_dispersal_destinations.png"))

### Define initialisation grids

assign_coord!(popgrid, (lat, lon), initsize) = begin
   coord = Lat(Contains(lat)), Lon(Contains(lon))
   popgrid[coord...] = initsize
end
init_popgrid!(popgrid, (lat, lon), initsize) = begin
   popgrid .= floattype(0.0)
    assign_coord!(popgrid, (lat, lon), initsize)
   nothing
end

aus_populationgrid = zero(aus_growthrates_zeromissing[Ti(1)])
init_popgrid!(aus_populationgrid, aus_incursionpoints[:Brisbane], carrycap)
us_populationgrid = zero(us_growthrates_zeromissing[Ti(1)])
init_popgrid!(us_populationgrid, us_incursionpoints[:San_Francisco], carrycap)
eu_populationgrid = zero(eu_growthrates_zeromissing[Ti(1)])
map(ip -> assign_coord!(eu_populationgrid, ip, carrycap), eu_incursionpoints)
# p = plot(eu_boolmask, label="")
# scatter!(p, reverse.(eu_incursionpoints); title="")
# eu_populationgrid |> plot



### Define the combined rulesets ###

# HumanDispersal has location specific precalculated
# data that the optimiser needs to update, making this a little
# messier than it otherwise would be

chain = Chain(localdisp, allee, growth)

# Test rule type stability
output = ArrayOutput((population=aus_populationgrid,);
    mask=collect(aus_boolmask),
    aux=(growthrates=aus_growthrates_zeromissing,),
    tspan=DateTime(2020, 1):timestep:DateTime(2026, 1)
)
DynamicGrids.isinferred(output, localdisp)
DynamicGrids.isinferred(output, allee)
DynamicGrids.isinferred(output, growth)
DynamicGrids.isinferred(output, chain)
DynamicGrids.isinferred(output, aus_humandisp)

ruleset_kwargs = (
    timestep=timestep,
    opt=SparseOpt(),
)
aus_ruleset = Ruleset(
    aus_humandisp, chain;
    ruleset_kwargs...
);
us_ruleset = Ruleset(
    us_humandisp, chain;
    ruleset_kwargs...
);
eu_ruleset = Ruleset(
    eu_humandisp, chain;
    ruleset_kwargs...
);
DynamicGrids.isinferred(output, us_ruleset)


### Output ###

# Define some color processors to use in live simulations.

zerocolor = ARGB32(0.6)
maskcolor = ARGB32(0.3)
textconfig = TextConfig(font="cantarel", bcolor=maskcolor, fcolor=ARGB32(0.8))

greyscale = ColorProcessor(Greyscale(), zerocolor, maskcolor, textconfig)
oranges = ColorProcessor(ColorSchemes.Oranges_3, zerocolor, maskcolor, textconfig)
jet = ColorProcessor(ColorSchemes.jet, zerocolor, maskcolor, textconfig)
viridis = ColorProcessor(ColorSchemes.viridis, zerocolor, maskcolor, textconfig)
inferno = ColorProcessor(ColorSchemes.inferno, zerocolor, maskcolor, textconfig)
magma = ColorProcessor(ColorSchemes.magma, zerocolor, maskcolor, textconfig)
