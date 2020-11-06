
basedir = realpath(joinpath(@__DIR__, ".."))
floattype = Float64
# Set the input data used in human dispersal models
humanactivity = :gdp
# humanactivity = :pop

# Set simulation replicates and grid aggregation
# Fast but innacurate
nreps = 2
scale = 10
# Slow but accurate: for the paper
# nreps = 100
# scale = 2

include(joinpath(basedir, "src", "spreadsetup.jl"))

# locname =                 :US
# locname =                 :EU
locname =                 :Aus

### Load everything required ###

if locname == :US
    tspan =                   us_tspan
    ruleset =                 us_ruleset
    boolmask =                us_boolmask
    missingmask_ =            us_missingmask
    growthrates_zeromissing = us_growthrates_zeromissing
    populationgrid =          us_populationgrid
    incursionpoints =         us_incursionpoints
    incursion =               us_incursionpoints[:San_Francisco]
elseif locname == :EU
    tspan =                   eu_tspan
    ruleset =                 eu_ruleset
    boolmask =                eu_boolmask
    missingmask_ =            eu_missingmask
    growthrates_zeromissing = eu_growthrates_zeromissing
    populationgrid =          eu_populationgrid
    incursionpoints =         eu_incursionpoints
    #incursion = eu_incursionpoints[:Rome]
elseif locname == :Aus
    tspan =                   DateTime(2021, 1):timestep:DateTime(2026, 1)
    ruleset =                 aus_ruleset
    boolmask =                aus_boolmask
    missingmask_ =            aus_missingmask
    growthrates_zeromissing = aus_growthrates_zeromissing
    populationgrid =          aus_populationgrid
    incursionpoints =         aus_incursionpoints
    incursion = aus_incursionpoints[:Brisbane]
end
aux = (growthrates=growthrates_zeromissing,)
init = (population=populationgrid,)
# Reconstruct with fitted parameters if you have them allready
using JLD2, FileIO, LabelledArrays
stochparams = load(joinpath(basedir, "output/stochparams_$humanactivity.jld2"), "stochparams_$humanactivity")
ruleset = reconstruct(ruleset, stochparams[:full])


# Run the simulation in a live interface :
# unccomment to run in a live interface
# using DynamicGridsInteract
# output = ElectronOutput((population=parent(populationgrid),);
#     ruleset=ruleset,
#     mask=BitArray(collect(boolmask)),
#     aux=aux,
#     tspan=tspan,
#     store=true,
#     processor=inferno,
#     minval=zero(carrycap),
#     maxval=carrycap,
#     fps=50,
# )
# display(output)
# And save the simulation as a gif:
# savegif(joinpath(basedir, "output/dispersal.gif"), output, ruleset; fps=10)


## Mapping incursion-point sensitivity

function incursionreps!(output, ruleset, missingmask, nreps)
    steps_established = similar(output[1][:population], Int)
    steps_established .= 0
    @set! steps_established.refdims = ()
    for rep in 1:nreps
        sim!(output, ruleset)
        println("rep: ", rep)
        map(output) do step
            steps_established .+= step[:population] .>= 1
        end
    end
    GeoArray(steps_established ./ nreps .* missingmask; missingval=missing, name=:established)
end
function incursionplot(established, key, nreps)
    plot(established;
        color=:inferno,
        xlabel="",
        ylabel="",
        title="$key: num timesteps established with $nreps reps"
    )
end
function set_incursion(output, (lat, lon))
    DynamicGrids.init(output)[:population] .= 0
    DynamicGrids.init(output)[:population][Lat(Contains(lat)), Lon(Contains(lon))] = carrycap
end

# Plot the initial incursion:

output = ArrayOutput(init; tspan=tspan, aux=aux, mask=boolmask)
established = incursionreps!(output, ruleset, missingmask_, nreps)
incursionplot(established, locname, nreps)
filename = joinpath(basedir, "output/$(locname)_steps_established_$(nreps)reps")
savefig("$filename.png")
write("$filename.tif", GDALarray, replace_missing(established, typemin(Float64)))
# Check the save worked
GDALarray("$filename.tif"; usercrs=EPSG(4326)) |> plot


# Plot all the incursions in incursionpoints:
plots = []

for (key, loc) in zip(keys(incursionpoints), incursionpoints)
    set_incursion(output, loc)
    steps_established = incursionreps!(output, ruleset, missingmask_, nreps)
    filename = joinpath(basedir, "output/$(locname)_steps_established_$(nreps)reps")
    write("$filename.tif", GDALarray, replace_missing(steps_established, typemin(Float64)))
    push!(plots, incursionplot(steps_established, key, nreps))
end

# Plot all togeth
plot(plots...; size=(1000, 600))
savefig(joinpath(basedir, "output/months_established.png"))

#' We can also save the plots as individual figures.
for (i, key) in enumerate(keys(incursionpoints))
    plot(plots[i])
    savefig(joinpath(basedir, "output/months_established_$key.png"))
end



### Incursion point grid ###

function incursion_point_grid!(
    cellsinvaded6, cellsinvaded12, ruleset, extent, nreps, scale, season;
    writedir="$basedir/output/tiffs", write=false
)
    if write
        rm(writedir; recursive=true)
        mkdir(writedir)
    end
    init = DynamicGrids.init(extent)
    londim, latdim = dims(init[:population], (Lon, Lat))
    latindex = reverse(val(latdim))
    outputs = [ArrayOutput(init; extent=extent) for t in 1:Threads.nthreads()]
    simdata = [DynamicGrids.SimData(extent, ruleset) for t in 1:Threads.nthreads()]
    inits = [deepcopy(init) for t in 1:Threads.nthreads()]
    Threads.@threads for i = 1:scale:size(cellsinvaded12, 1) * scale
        acc = similar(init[:population], UInt8)
        acc .= 0
        @set! acc.missingval = 0x00
        threadinit = inits[Threads.threadid()]
        threaddata = simdata[Threads.threadid()]
        threadoutput = outputs[Threads.threadid()]
        for j = 1:scale:size(cellsinvaded12, 2) * scale
            println("i, j: ", (i, j))
            DynamicGrids.mask(extent)[i, j] || continue
            checkbounds(Bool, acc, i, j) || continue
            lat, lon = ArchGDAL.reproject([[londim[j], latindex[i]]], crs(init[:population]), EPSG(4326))[1]
            # Zero the accumulator array
            acc .= 0
            # Zero init array
            inits[Threads.threadid()][:population] .= 0
            # Set the current cell to maximum population
            inits[Threads.threadid()][:population][i, j] = carrycap
            invaded6 = 0
            invaded12 = 0
            for k in 1:nreps
                sim!(threadoutput, ruleset; init=inits[Threads.threadid()], simdata=threaddata)
                invaded6 += count(x -> x > one(x), threadoutput[7][:population])
                invaded12 += count(x -> x > one(x), threadoutput[13][:population])
                acc .+= threadoutput[13][:population] .> 0
            end
            cellsinvaded6[(i - 1) ÷ scale + 1, (j - 1) ÷ scale + 1] = invaded6 / nreps
            cellsinvaded12[(i - 1) ÷ scale + 1, (j - 1) ÷ scale + 1] = invaded12 / nreps
            if write && invaded12 > 0
                filepath = "$writedir/cells_invaded_from_$(lon)_$(lat)_$season.tif"
                # Remove underscores
                filepath = replace(filepath, r"(.*)\.(.*)\.(.*)\.tif" => s"\1_\2_\3.tif")
                GeoData.write(filepath, GDALarray, acc)
            end
        end
    end
end

function spread_probability!(output, I...; 
    gridname, 
    initval,
    ruleset=ruleset(output), 
    acc=zero(output.extent.init[gridname]), 
    init=output.extent.init, 
    nreps=1,
    simdata=nothing
)
    invaded = 0
    acc .= 0
    initgrids = deepcopy(init)
    initgrids[gridname] .= 0
    # Set the current cell to maximum population
    initgrids[gridname][I...] = carrycap
    @show nreps
    for k in 1:nreps
        sim!(output, ruleset; init=initgrids, simdata=simdata)
        invaded += count(x -> x > one(x), output[end][gridname])
        acc .+= output[end][gridname] .> 0
    end
    acc ./= nreps
    return acc
end

pyplot() 
oneyear = DateTime(2021):Month(1):DateTime(2022)
output = ArrayOutput(init; tspan=oneyear, aux=aux, mask=boolmask)
spreaddir = joinpath(basedir, "output/spread_probability")
mkpath(spreaddir)

for key in keys(incursionpoints)
    println("==== Generating spread probability for $key ======================")
    lat, lon = incursionpoints[key]
    spreadprob = spread_probability!(output, Lon(Contains(lon)), Lat(Contains(lat)); 
        initval=carrycap,                  
        gridname=:population, 
        ruleset=ruleset, 
        simdata=nothing,
        nreps=nreps,
    )
    filename = joinpath(spreaddir, "$(key)_spread_probability_$(nreps)reps")
    plot(spreadprob; title=key, colorbar_title="Spread probability")
    savefig(string(filename, ".png"))
    write(string(filename, ".tif"), GDALarray, spreadprob)
end


# Incursion locations from a spreadsheet

cellsinvaded = GeoData.aggregate(Center(), populationgrid, scale) .* 0
@set! cellsinvaded.name = "Cells invaded"
cellsinvaded6_summer = deepcopy(cellsinvaded)
cellsinvaded12_summer = deepcopy(cellsinvaded)
cellsinvaded6_winter = deepcopy(cellsinvaded)
cellsinvaded12_winter = deepcopy(cellsinvaded)
tspan_summer = DateTime(2021, 1):timestep:DateTime(2022, 1)
tspan_winter = DateTime(2021, 7):timestep:DateTime(2022, 7)
init = (population=zero(populationgrid),)

# Summer sims
extent = DynamicGrids.Extent(init, boolmask, aux, tspan)
@time incursion_point_grid!(
    cellsinvaded6_summer, cellsinvaded12_summer, ruleset, extent, nreps, scale, "summer"; 
    write=false
)

# Winter sims
@set! extent.tspan = tspan_winter
@time incursion_point_grid!(
    cellsinvaded6_winter, cellsinvaded12_winter, ruleset, extent, nreps, scale, "winter"; 
    write=false
)

# Now plot and save the combined cells invaded arrays:

function plot_cells_invaded(ci, label, months, nreps, mmask; kwargs...)
    toarea = 9*9*1e-6
    area = ci .* toarea .* GeoData.aggregate(Center(), mmask, scale)
    p = plot(area;
        colorbar_title="Area invaded from incursion pt. (million km²)",
        title="Area invaded after $months months",
        kwargs...
    )
    path = "output/cellsinvaded_$(label)_$(months)months_$(nreps)reps"
    savefig(path * ".png")
    write(path * ".tif", GDALarray, replace_missing(ci, -1.00e+10))
    p
end
plot_cells_invaded(cellsinvaded6_summer, "summer_$humanactivity", 6, nreps, missingmask_)
plot_cells_invaded(cellsinvaded12_summer, "summer_$humanactivity", 12, nreps, missingmask_)
plot_cells_invaded(cellsinvaded6_winter, "winter_$humanactivity", 6, nreps, missingmask_)
plot_cells_invaded(cellsinvaded12_winter, "winter_$humanactivity", 12, nreps, missingmask_)

