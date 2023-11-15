module TurbulentImages

using Images
using Colors
using Oceananigans
using Oceananigans.Fields: regrid!

export turbulent_image_simulation

function saturate!(vals, ϵ)
    @. vals = (1 + ϵ) * vals - ϵ/2
    vals .= clamp.(vals, 0, 1)
    return nothing
end

function image_to_initial_condition(filename;
                                    architecture = CPU(),
                                    topology = (Periodic, Flat, Bounded))
                                    
    color_img = load(filename)
    gray_img = Gray.(color_img)
    gray_vals = map(pixel -> Float64(pixel.val), gray_img)
    gray_vals = rotr90(gray_vals)

    Ix, Iz = size(gray_vals)
    aspect = Ix / Iz
    z = (0, 1)
    x = (0, aspect)
    
    native_grid = RectilinearGrid(architecture; topology,
                                  size = (Ix, Iz),
                                  z = (0, 1),
                                  x = (0, aspect))

    # Regrid image to new grid
    gray_field = CenterField(native_grid)
    gray_vals = reshape(gray_vals, Ix, 1, Iz)
    set!(gray_field, gray_vals)

    return gray_field
end

function regrid_xy(grid, native_field)

    native_grid = native_field.grid
    arch = Oceananigans.Architectures.architecture(native_grid)
    Ix, Iy, Iz = size(native_grid)
    aspect = Ix / Iz
    Nx, Ny, Nz = size(grid)

    intermediate_grid = RectilinearGrid(arch,
                                        size = (Ix, Nz),
                                        z = (0, 1),
                                        x = (0, aspect),
                                        topology = Oceananigans.Grids.topology(grid))

    intermediate_field = CenterField(intermediate_grid)
    regrid!(intermediate_field, native_field)

    field = CenterField(grid)
    regrid!(field, intermediate_field)

    return field
end

"""

    turbulent_image_simulation(img_filename,
                               output_filename;
                               advection = WENO(order=5),
                               progress_schedule = IterationInterval(10),
                               output_schedule = IterationInterval(5),
                               Nz=128)

"""
function turbulent_image_simulation(img_filename,
                                    output_filename;
                                    architecture = CPU(),
                                    advection = WENO(order=5),
                                    progress_schedule = IterationInterval(10),
                                    output_schedule = IterationInterval(5),
                                    wizard_schedule = IterationInterval(5),
                                    x_pixels = nothing,
                                    z_pixels = nothing)

    img = image_to_initial_condition(img_filename)
    Ix, Iy, Iz = size(img)
    aspect = Ix / Iz

    if isnothing(x_pixels) # Nx is not set
        if isnothing(z_pixels) # Nz is not set either!
            Nz = 128 # sensible default
        else
            Nz = z_pixels
        end
        Nx = floor(Int, aspect * Nz)
    elseif isnothing(z_pixels) # Nz is not set, but Nx is
        Nx = x_pixels
        Nz = floor(Int, Nx / aspect)
    end

    grid = RectilinearGrid(architecture,
                           size = (Nx, Nz),
                           halo = (5, 5),
                           z = (0, 1),
                           x = (0, aspect),
                           topology = (Periodic, Flat, Bounded))

    bi = regrid_xy(grid, img)

    model = NonhydrostaticModel(; grid, advection,
                                tracers = :b,
                                timestepper = :RungeKutta3,
                                buoyancy = BuoyancyTracer())

    set!(model, b=bi)

    simulation = Simulation(model, Δt=0.01, stop_time=2)
    wizard = TimeStepWizard(cfl=0.8, max_change=1.1)
    simulation.callbacks[:wizard] = Callback(wizard, wizard_schedule)

    wall_clock = Ref(time_ns())
    function progress(sim)
        elapsed = 1e-9 * (time_ns() - wall_clock[])

        @info string("Iter: ",
                     iteration(sim),
                     ", time: ", time(sim),
                     ", wall time: ", prettytime(elapsed),
                     ", Δt: ", sim.Δt)

        wall_clock[] = time_ns()

        return nothing
    end

    simulation.callbacks[:progress] = Callback(progress, progress_schedule)

    b = model.tracers.b
    outputs = (; b)
    simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                        schedule = output_schedule,
                                                        filename = output_filename,
                                                        overwrite_existing = true)

    return simulation
end

end
