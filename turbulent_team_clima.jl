using Images
using Colors
using Oceananigans
using Oceananigans.Fields: regrid!
using GLMakie

include("image_to_initial_conditions.jl")

image_filename = "20231029-_DSC6318.jpg" #"20221104-DSC06270.jpg"
output_filename = "turbulent_plume.jld2"

simulation = buoyant_image_simulation(image_filename, output_filename, Nz=1024)
simulation.stop_time = 2.0

run!(simulation)

bt = FieldTimeSeries(output_filename, "b")

n = Observable(1)
bn = @lift interior(bt[$n], :, 1, :)

Nx, Ny, Nz = size(simulation.model.grid)
aspect = Nx / Nz
fig = Figure(resolution=(600aspect, 600))
ax = GLMakie.Axis(fig[1, 1])

heatmap!(ax, bn, colormap=:grays)

Nt = length(bt.times)

record(fig, "turbulent_plume.mp4", 1:Nt, framerate=48) do nn
    @info string("Drawing frame $nn of $Nt...")
    n[] = nn
end

display(fig)

