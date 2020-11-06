
basedir = realpath(joinpath(@__DIR__, ".."))
basedir = realpath("/home/raf/.julia/dev/SpottedWingPaper")
@show basedir 
floattype = Float64
# Set the input data used in human dispersal models
humanactivity = :gdp
# humanactivity = :pop

# Opimisation and validation settings
cpu_cores = :multi
# cpu_cores = :single
ngroups = 5 # should be your physical cpu cores minus 1
# groupsize = 50 # Used in the paper - for 250 replicates when ngroups == 5
groupsize = 1 # Fast - ngroups total replicates, to check the scripts




### Load data and packages ###

include(joinpath(basedir, "src", "spreadsetup.jl"))
using Optim, LossFunctions, JLD2, FileIO, DataStructures, LabelledArrays, DynamicGridsGtk

### Load common sctipt functions ###

# Main parametrisation method
function parametrise_ruleset!(optimresults, parametriser, key)
    println("\n\n==== Model: ", key, now(), " ===================================================")
    ruleset = parametriser.ruleset
    initpars = flatten(rules(ruleset))
    println(initpars)
    # Make a labelled array so we can ignore the order
    parnames = fieldnameflatten(rules(ruleset), Number)
    # Get the lower and upper bounds for params with flatten
    bnds = metaflatten(rules(ruleset), FieldMetadata.bounds, Number)
    lower = first.(bnds)
    upper = last.(bnds)
    names = fieldnameflatten(rules(ruleset), Number)
    namedpars = @LArray [initpars...] parnames
    display(DataFrame(names=[names...], pars=[initpars...], lower=[lower...], upper=[upper...]))
    # Zero params to make sure parametrizer is working
    output.running = false
    ruleset.rules = reconstruct(ruleset.rules, initpars .* 0)
    println("Running optimiser...")
    @time res = Optim.optimize(
        parametriser, [lower...], [upper...], namedpars, SAMIN(),
        Optim.Options(iterations=iterations, show_trace=true, store_trace=true)
    )
    optimresults[key] = res
    println("Saving JLD2 file")
    save("output/optimresults_$key.jld2", "optimresults", optimresults)
    println("Calculating losses...")
    optimresults
end

# Method to add results to a dataframe table
function ruleset_results!(resultdf, rulesets, parametriser, params, location, scenario=:regular)
    for key in keys(params)
        println("\n==== Model: ", key, now(), " ===================================================")
        ruleparams = params[key]
        # Set the parameters to zero, to make sure they are updated in the optimizer
        ruleset = deepcopy(rulesets[key])
        ruleset.rules = reconstruct(ruleset.rules, ruleparams .* 0)
        parametriser = @set parametriser.ruleset = ruleset
        loss = parametriser(ruleparams)
        # Standard deviation of accuracy
        std_ = std(accuracy.(Ref(Dispersal.targets(parametriser.objective)),
                             vcat(parametriser.results...)))
        row = (
            location=location,
            model=key,
            scenario=scenario,
            loss=loss,
            std=std_,
            accuracy=accuracy(Dispersal.targets(parametriser.objective), loss)
        )
        push!(resultdf, row)
    end
end
accuracy(target, loss) = one(loss) - loss/length(target)

# Plot accuracy heatmap
function lossplot(output::RegionOutput, ruleset)
    sim!(output, ruleset)
    diff = (transformfunc.(output[1]) - transformfunc.(output.objective.occurrence))
    sleep(0.5)
    println("Accuracy: $(1 - (sum(abs.(diff)) * 0.5 / length(diff)))")
    heatmap(diff; clims=(-1,1))
end

# Use ColorRegionFit in GtkOutput
function regionsim(ruleset, objective::RegionObjective,
        mask, popgrid, growthrates, tspan, carrycap; fps=20, outputtype=GtkOutput, filename="output/regionsim.gif")
    truescheme = Greyscale()
    falsescheme = reverse(ColorSchemes.Reds_3)
    truezerocolor =  (0.5, 0.5, 0.5)
    falsezerocolor = falsescheme[1]
    maskcolor =      (0.0, 1.0, 1.0)
    region_processor = ColorRegionFit(
        objective, truescheme, falsescheme,
        truezerocolor, falsezerocolor, maskcolor
    )
    gtkoutput = outputtype((population=popgrid,); filename=filename, mask=BitArray(collect(mask)), aux=(growthrates=growthrates,),
        store=true,
        tspan=tspan,
        processor=region_processor,
        minval=zero(carrycap),
        maxval=carrycap,
        fps=fps,
    )
    sim!(gtkoutput, ruleset)
    gtkoutput
end


# Sensitivity routines
function growthrate_sensitivity!(resultdf, rulesets, parametriser, params, stressors, change, selectors, location)
    output = parametriser.output
    crs_ = if output isa Output
        DynamicGrids.aux(output)[:growthrates]
    elseif output isa Vector{<:Output}
        DynamicGrids.aux(first(output))[:growthrates]
    end |> crs
    for s in stressors, c in change
        scenario = Symbol(s, :_, c)
        println(scenario)
        # Load growthrates file matching this scenario
        growthratesfilepath = joinpath(basedir, "output/sensitivity/growthrates_$(scenario).ncd")
        isfile(growthratesfilepath) || error("$growthratesfilepath not found")
        growthrates = NCDarray(growthratesfilepath; crs=crs_, dimcrs=EPSG(4326))[selectors...]  |>
            x -> setdims(x, (@set dims(x, Ti).mode.span = Regular(Month(1)))) |>
            x -> convertmode(Projected, x) |>
            x -> permutedims(x, (Lat, Lon, Ti)) |>
            x -> reverse(x; dims=Lat) |>
            x -> replace_missing(x, 0)
        # set the output aux to the new growthrates
        crs_ = if output isa Output
            @set! parametriser.output.extent.aux = (growthrates=growthrates,)
        elseif output isa Vector{<:Output}
            newoutputs = map(output) do o
                @set! o.extent.aux = (growthrates=growthrates,)
            end
            @set! parametriser.output = newoutputs
        end
        println("\n\n\n######### scenario: $scenario ########################################")
        ruleset_results!(resultdf, rulesets, parametriser, params, location, scenario)
    end
end

function param_sensitivity!(resultdf, rulesets, parametrizer, params, changekeys, change, location)
    # Rebuild all rulesets with adjusted params
    for paramkey in changekeys, c in change
        scenario = Symbol(paramkey, :_, c)
        modified_rulesets = Dict()
        for rulesetkey in keys(rulesets)
            rs = rulesets[rulesetkey]
            names = fieldnameflatten(rules(rs), Number)
            pars = NamedTuple{names}(flatten(rules(rs), Number))
            if paramkey in names
                ruleset = deepcopy(rs)
                # multiply the parameter by the amount of change
                pars = @set pars[paramkey] = pars[paramkey] * c
                # rebuild
                ruleset.rules = reconstruct(rules(ruleset), pars, Number)
                modified_rulesets[rulesetkey] = ruleset
            end
        end
        println("\n\n\n######### scenario: $scenario% ########################################")
        ruleset_results!(resultdf, modified_rulesets, parametrizer, params, location, scenario)
    end
end

function years_detected(df, years, nstates)
    occurrence = zeros(Bool, nstates, length(years))
    for y in eachindex(years), state in 1:nstates
        year = years[y]
        occurrence[state, y] = year >= df.year_detected[state]
    end
    occurrence
end


# Load US data to optimise the model:

timestep = Month(1)
us_df = CSV.read(joinpath(basedir, "data/region_year_detected_US.csv"))
us_years = Dates.year(first(us_tspan)):Dates.year(last(us_tspan))
us_nstates = round(Int, maximum(us_states))
us_occurrence = years_detected(us_df, us_years, us_nstates)

heatmap(us_occurrence; title="US occurrence", xlab="Year", ylab="State ID")


# Define model combinations for comparison

kwargs = (
    timestep=timestep,
    opt=SparseOpt(),
)
us_stochastic_rulesets = (
    full=Ruleset(us_humandisp, Chain(localdisp, allee, growth); kwargs...),
    nolocal=Ruleset(us_humandisp, Chain(allee, growth); kwargs...),
    noallee=Ruleset(us_humandisp, Chain(localdisp, growth); kwargs...),
    noclimate = Ruleset(us_humandisp, Chain(localdisp, allee, constant_growth); kwargs...),
)
null_rulesets = (
    null=Ruleset(; kwargs...),
)
null_params = (null=(),)



## Set up the optimiser

framesperstep = 12
startmonth = 5
transformfunc = x -> 2x - 1 # Transform to -1, +1 instead of 0, 1
lossfunc = ZeroOneLoss()
iterations = 1000
detectionthreshold = floattype(1e7)
threading = if cpu_cores == :multi
    ThreadedReplicates()
else
    SingleCoreReplicates()
end

us_objective = RegionObjective(
    detectionthreshold,
    parent(us_states),
    us_occurrence,
    framesperstep,
    startmonth,
)
us_aux = (growthrates=us_growthrates_zeromissing,)
us_ruleset = us_stochastic_rulesets[:full]
#us_ruleset = us_deterministic_rulesets[:nomodel]
us_output = RegionOutput((population=parent(us_populationgrid),);
    aux=us_aux,
    mask=BitArray(us_boolmask),
    tspan=us_tspan,
    objective=us_objective,
)


# rs = us_full_parametrized # If you have loaded the parametrised version
# rs = us_stochastic_rulesets[:full] # The initial version
# using DynamicGridsInteract
# InteractOutput((population=parent(us_populationgrid),);
    # ruleset=us_full_parametrized,
    # aux=us_aux,
    # mask=BitArray(us_boolmask),
    # tspan=us_tspan,
    # maxval=carrycap,
    # processor=inferno,
# ) |> display

#' ## Profile model performance
#' 
#' It's good to check the everythign is running well
#' before starting a 24 hour simulation:
#' 
#+ 

# simdata = DynamicGrids.SimData(us_output.extent, us_ruleset)
# using BenchmarkTools
# @btime sim!($us_output, $us_ruleset)
# @profiler for i = 1:10 sim!(us_output, us_ruleset; simdata=simdata) end
# isinferred(us_output, us_ruleset)

## Run the optimizer

# **WARNING** this might take a few days if there is no params file already saved

stochparamsfile = "output/stochparams_$humanactivity.jld2"
optimresultsfile = "output/optimresults_$humanactivity.jld2"
optimresults = OrderedDict()
if isfile(stochparamsfile)
    stochparams = load(stochparamsfile, "stochparams_$humanactivity")
    stochparams = load(stochparamsfile, "stochparams_$humanactivity")
    stochrulesetkeys = keys(stochparams)
    stochparamvals = values(stochparams)
else
    for key in keys(us_stochastic_rulesets)
        us_parametriser = Parametriser(
            us_stochastic_rulesets[key], us_output, us_objective, transformfunc, lossfunc, ngroups, groupsize, threading
        )
        parametrise_ruleset!(optimresults, us_parametriser, key)
        # Save as we go
        stochrulesetkeys = Tuple(keys(optimresults))
        stochparamvals = ((Optim.minimizer(optimresults[k]) for k in stochrulesetkeys)...,)
        stochparams = NamedTuple{stochrulesetkeys}(stochparamvals)
        save(optimresultsfile, "optimresults", optimresults)
        save("output/stochparams_$humanactivity.jld2", "stochparams_$humanactivity", stochparams)
    end
    optimresults
end


# Test the output:

us_full_paremetrized = reconstruct(us_stochastic_rulesets[:full], stochparams[:full])
lossplot(us_output, us_full_paremetrized)
# regionoutput = regionsim(us_full_paremetrized, us_objective, us_boolmask, us_populationgrid,
    # us_growthrates_zeromissing, us_tspan, carrycap)

# Add US results to dataframes

# Add parameters to DataFrame
allparamkeys = reduce(union, map(keys, stochparams); init=Symbol[])
colnames = [:lower, :upper, stochrulesetkeys...]
missingparams = [Any[missing for pk in allparamkeys] for rk in colnames]
paramdf = DataFrame([allparamkeys, missingparams...], [:param, colnames...])
fullbounds = metaflatten(us_stochastic_rulesets[:full].rules, FieldMetadata.bounds, Union{Real,Quantity})
fullkeys = fieldnameflatten(us_stochastic_rulesets[:full].rules, Union{Real,Quantity})
for key in stochrulesetkeys
    stochrulesetresults = stochparams[key]
    for (row, rowkey) in enumerate(allparamkeys)
        if rowkey in keys(stochparams[key])
            paramdf[key][row] = stochrulesetresults[rowkey]
        end
    end
end
# Add param bounds, using :full model
for (row, rowkey) in enumerate(allparamkeys)
    if rowkey in fullkeys
        paramdf[:lower][row], paramdf[:upper][row] = fullbounds[row]
    end
end
# Set intrinsicrate bounds from ExactLogisticGrowth model
select = paramdf.param .== :intrinsicrate
paramdf[select, :lower], paramdf[select, :upper] =
    FieldMetadata.bounds(ExactLogisticGrowth, :intrinsicrate)

# Add parameter names
stochparamnames = keys(stochparams[:full])


### Get the results in a dataframe: ###

# Loss/accuracy results dataframe
resultdf = DataFrame(location=[], model=[], scenario=[], loss=Float64[], std=Float64[], accuracy=Float64[])
us_parametriser = Parametriser(
    first(values(us_stochastic_rulesets)), us_output, us_objective, transformfunc,
    lossfunc, ngroups, groupsize, threading
)
us_null_parametriser = Parametriser(
    first(values(null_rulesets)), us_output, us_objective, transformfunc,
    lossfunc, ngroups, groupsize, threading
)
us_location = :US


# Get the results for all stochastic rulsets, and the null rulesets:

ruleset_results!(resultdf, us_stochastic_rulesets, us_parametriser, stochparams, us_location)
ruleset_results!(resultdf, null_rulesets, us_null_parametriser, null_params, us_location)


# Show the output in the plot pane:

display(paramdf)
display(resultdf)
showtable(resultdf)


stressors = (:heatstress, :coldstress, :wiltstress, :wetstress)
paramnames = (:human_exponent, :dist_exponent, :dispersalperpop, :max_dispersers)
change = 80, 120
selectors = us
growthrate_sensitivity!(resultdf, us_stochastic_rulesets, us_parametriser, stochparams, stressors, change, selectors, us_location)
param_sensitivity!(resultdf, us_stochastic_rulesets, us_parametriser, stochparams, paramnames, change, us_location)
showtable(resultdf)


# Load data for validation against the EU

eu_stochastic_rulesets = (
    full=Ruleset(eu_humandisp, Chain(localdisp, allee, growth); kwargs...),
    nolocal=Ruleset(eu_humandisp, Chain(allee, growth); kwargs...),
    noallee=Ruleset(eu_humandisp, Chain(localdisp, growth); kwargs...),
    noclimate = Ruleset(eu_humandisp, Chain(localdisp, allee, constant_growth); kwargs...)
)
eu_rulesets = map((r, p) -> reconstruct(r, p, Number), eu_stochastic_rulesets, stochparams)
eu_df = CSV.read(joinpath(basedir, "data/region_year_detected_EU.csv"))
eu_years = Dates.year(first(eu_tspan)):Dates.year(last(eu_tspan))
eu_occurrence = years_detected(eu_df, eu_years, maximum(eu_states))
eu_objective = RegionObjective(detectionthreshold, eu_states, eu_occurrence, framesperstep, startmonth)
eu_aux = (growthrates=eu_growthrates_zeromissing,)
eu_init = (population=eu_populationgrid,)
eu_output = RegionOutput(eu_init;
    aux=eu_aux,
    mask=eu_boolmask,
    tspan=eu_tspan,
    objective=eu_objective,
)
eu_parametriser = Parametriser(
    first(values(eu_rulesets)), eu_output, eu_objective,
    transformfunc, lossfunc, ngroups, groupsize, threading
)
eu_null_parametriser = Parametriser(
    first(values(null_rulesets)), eu_output, eu_objective, transformfunc,
    lossfunc, ngroups, groupsize, threading
)
eu_location = :EU



# Test output

heatmap(eu_occurrence; title="EU occurrence", xlab="Year", ylab="State ID")
eu_output.extent.init.population |> plot
lossplot(eu_output, null_rulesets[:null])
savefig("output/eu_null_loss.png")

eu_full = reconstruct(eu_stochastic_rulesets[:full], stochparams[:full])
lossplot(eu_output, eu_full)
savefig("output/eu_full_loss.png")
# regionoutput = regionsim(eu_full, eu_objective, eu_boolmask, eu_populationgrid,
          # eu_growthrates_zeromissing, eu_tspan, carrycap; fps=10)#, outputtype=GifOutput, filename="eu_gdp.gif")
 


### Validate US/EU ###

selectors = eu
ruleset_results!(resultdf, null_rulesets, eu_null_parametriser, null_params, eu_location)
ruleset_results!(resultdf, eu_stochastic_rulesets, eu_parametriser, stochparams, eu_location)
growthrate_sensitivity!(resultdf, eu_stochastic_rulesets, eu_parametriser, stochparams, stressors, change, selectors, eu_location)
param_sensitivity!(resultdf, eu_stochastic_rulesets, eu_parametriser, stochparams, paramnames, change, eu_location)
display(resultdf)
@save "output/resultdf_$(humanactivity).jld2" resultdf
CSV.write("output/sensitivity_$(humanactivity).csv", resultdf)
CSV.write("output/params_$(humanactivity).csv", paramdf)


### NoHuman model ###

# This needs to be separated out so that we can
# optimise the no-human model with a free λ parameter.

println("==== Running nohuman model ==============================================")

λ = floattype(0.0125)
radius = 6

@time hood = DispersalKernel{radius}(
    formulation=ExponentialKernel(λ),
    distancemethod=AreaToArea(30),
    cellsize=floattype(1.0)
)

free_localdisp = InwardsPopulationDispersal{:population,:population}(hood)
log.(hood.kernel) |> heatmap

nohuman = Ruleset(Chain(free_localdisp, allee, growth); kwargs...)
nohuman_rulesets = (nohuman=nohuman,)
# Parametrise λ - it's now a general spread model, not local movement
@flattenable @bounds ExponentialKernel begin
    λ | true | (0.0, 100.0)
end

# This model is deteministic so we don't need replicates
ngroups = 1
groupsize = 1
threading = SingleCoreReplicates()
us_nohuman_parametriser = Parametriser(
    nohuman, us_output, us_objective, transformfunc, lossfunc, ngroups, groupsize, threading
)
eu_nohuman_parametriser = Parametriser(
    nohuman, eu_output, eu_objective, transformfunc, lossfunc, ngroups, groupsize, threading
)

# Parametrise unless an output is already saved
nohumanfile = joinpath(basedir, "output/params_$(humanactivity)_nohuman.jld2")
if isfile(nohumanfile)
    nohuman_params = load(nohumanfile, "nohuman_params")
    nohuman_paramvals = values(nohuman_params)
else
    parametrise_ruleset!(optimresults, us_nohuman_parametriser, :nohuman)
    nohuman_rulesetkeys = (:nohuman,)
    nohuman_paramvals = Optim.minimizer(optimresults[:nohuman])
    nohuman_params = NamedTuple{(:nohuman,)}((nohuman_paramvals,))
end

parametrized_nohuman = reconstruct(nohuman, nohuman_paramvals, Number)
lossplot(us_output, parametrized_nohuman)
# regionoutput = regionsim(
    # parametrized_nohuman, us_objective, us_boolmask, us_populationgrid,
    # us_growthrates_zeromissing, us_tspan, carrycap;
    #outputtype=GifOutput, filename="output/us_nohuman_regionsim.gif", fps=10,
# )
lossplot(eu_output, nohuman)
# regionoutput = regionsim(
    # parametrized_nohuman, eu_objective, eu_boolmask, eu_populationgrid,
    # eu_growthrates_zeromissing, eu_tspan, carrycap;
    #outputtype=GifOutput, filename="output/eu_nohuman_regionsim.gif", fps=10,
# )

# Regular results
ruleset_results!(resultdf, nohuman_rulesets, us_nohuman_parametriser, nohuman_params, us_location)
ruleset_results!(resultdf, nohuman_rulesets, eu_nohuman_parametriser, nohuman_params, eu_location)
# Growth rate sensitivity
growthrate_sensitivity!(resultdf, nohuman_rulesets, us_nohuman_parametriser, nohuman_params, stressors, change, selectors, us_location)
growthrate_sensitivity!(resultdf, nohuman_rulesets, eu_nohuman_parametriser, nohuman_params, stressors, change, selectors, eu_location)
# Dispersal parameter sensitivity
nohuman_paramnames = keys(nohuman_params[:nohuman])
param_sensitivity!(resultdf, nohuman_rulesets, us_nohuman_parametriser, nohuman_params, nohuman_paramnames, change, us_location)
param_sensitivity!(resultdf, nohuman_rulesets, eu_nohuman_parametriser, nohuman_params, nohuman_paramnames, change, eu_location)

nohumancol = deepcopy(paramdf.full) .= missing
insert!(paramdf, 8, nohumancol, :nohuman)
select = paramdf.param .== :minfounders
paramdf[select, :nohuman] = nohuman_params[1].minfounders
newrow = [:λ, 0.0, 100.0, missing, missing, missing, missing, nohuman_params[1].λ]
push!(paramdf, newrow)
showtable(paramdf)
showtable(resultdf)
display(paramdf)
display(resultdf)


# Save results as a csv file

@save "output/params_$(humanactivity)_nohuman.jld2" nohuman_params
CSV.write(joinpath(basedir, "output/sensitivity_$(humanactivity)_nohuman.csv"), resultdf)
CSV.write(joinpath(basedir, "output/params_$(humanactivity)_nohuman.csv"), paramdf)
