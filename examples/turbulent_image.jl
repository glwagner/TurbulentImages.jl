using Oceananigans
using TurbulentImages
using GLMakie

image_filename = "20231114-UND02051.jpg"
output_name = "turbulent_dillon"
output_filename = output_name * ".jld2"
output_moviename = output_name * ".mp4"

simulation = turbulent_image_simulation(image_filename, output_filename, z_pixels=512)
run!(simulation)

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

display(fig)

