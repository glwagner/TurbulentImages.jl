using Oceananigans
using TurbulentImages
using GLMakie

image_filename = "kitten.jpg"
output_name = "kitten"

output_filename = output_name * ".jld2"
output_moviename = output_name * ".mp4"
image_path = joinpath(@__DIR__, image_filename)

simulation = turbulent_image_simulation(image_path, output_filename,
                                        advection = WENO(order=9),
                                        z_pixels = 128)
run!(simulation)

# Make a nice movie

bt = FieldTimeSeries(output_filename, "b")

n = Observable(1)
bn = @lift interior(bt[$n], :, 1, :)

Nx, Ny, Nz = size(simulation.model.grid)
aspect = Nx / Nz

fig = Figure(resolution=(600aspect, 600))
ax = GLMakie.Axis(fig[1, 1])
heatmap!(ax, bn, colormap=:grays)
hidedecorations!(ax)

Nt = length(bt.times)

record(fig, output_moviename, 1:Nt, framerate=48) do nn
    @info string("Drawing frame $nn of $Nt...")
    n[] = nn
end

