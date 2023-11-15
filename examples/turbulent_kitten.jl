using Oceananigans
using TurbulentImages
using GLMakie

image_filename = "kitten.jpg"
output_name = "kitten"

output_filename = output_name * ".jld2"
output_moviename = output_name * ".mp4"
image_path = joinpath(@__DIR__, image_filename)

simulation = turbulent_image_simulation(image_path, output_filename,
                                        architecture = CPU(),
                                        advection = WENO(order=9),
                                        z_pixels = 128)

simulation.stop_time = 2 # how turbulent does kitty want to get?
run!(simulation)

# Make a nice movie

bt = FieldTimeSeries(output_filename, "b")

n = Observable(1)
bn = @lift interior(bt[$n], :, 1, :)

Nx, Ny, Nz = size(simulation.model.grid)
aspect = Nx / Nz

fig = Figure(resolution=(600aspect, 600))
ax = GLMakie.Axis(fig[1, 1], title="meow!")
heatmap!(ax, bn, colormap=:grays, colorrange=(0, 1))
hidedecorations!(ax)

# Make a sweet movie that also goes in reverse
stillframes = 10
framerate = 24
movingframes = length(bt.times)

record(fig, output_moviename; framerate) do io
    [recordframe!(io) for _ = 1:stillframes]

    for nn in 1:movingframes
        n[] = nn
        recordframe!(io)
    end

    [recordframe!(io) for _ = 1:stillframes]

    ax.title[] = "!woem"

    for nn in movingframes:-1:1
        n[] = nn
        recordframe!(io)
    end

    [recordframe!(io) for _ = 1:stillframes]
end

