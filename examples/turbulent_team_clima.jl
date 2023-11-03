using Images
using Colors
using Oceananigans
using Oceananigans.Fields: regrid!
using GLMakie

image_filename = "20221104-DSC06270.jpg"
output_filename = "turbulent_team_clima.jld2"

include("image_to_initial_conditions.jl")

simulation = buoyanct_image_simulation(image_filename, output_filename, Nz=64)

run!(simulation)

bt = FieldTimeSeries(filename, "b")

n = Observable(1)
bn = @lift interior(bt[$n], :, 1, :)

fig = Figure(resolution=(600aspect, 600))
ax = GLMakie.Axis(fig[1, 1])

heatmap!(ax, bn, colormap=:grays)

Nt = length(bt.times)

record(fig, "turbulent_team_clima.mp4", 1:Nt, framerate=48) do nn
    @info string("Drawing frame $nn of $Nt...")
    n[] = nn
end

display(fig)

