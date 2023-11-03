module TurbulentImages

function saturate!(vals, ϵ)
    @. vals = (1 + ϵ) * vals - ϵ/2
    vals .= clamp.(vals, 0, 1)
    return nothing
end

function image_to_initial_condition(filename; topology = (Periodic, Flat, Bounded))
                                    
    color_img = load(filename)
    gray_img = Gray.(color_img)
    gray_vals = map(pixel -> Float64(pixel.val), gray_img)
    gray_vals = rotr90(gray_vals)

    Ix, Iz = size(gray_vals)
    aspect = Ix / Iz
    z = (0, 1)
    x = (0, aspect)
    
    native_grid = RectilinearGrid(size=(Ix, Iz), z=(0, 1), x=(0, aspect); topology)

    # Regrid image to new grid
    gray_field = CenterField(native_grid)
    set!(gray_field, gray_vals)

    return gray_field
end

function regrid_xy(grid, native_field)

    native_grid = native_field.grid
    Ix, Iy, Iz = size(native_grid)
    intermediate_grid = RectilinearGrid(size=(Ix, Nz), z=(0, 1), x=(0, aspect), topology)

    intermediate_field = CenterField(intermediate_grid)
    regrid!(intermediate_field, native_field)

    field = CenterField(grid)
    regrid!(field, intermediate_field)

    return field
end

function buoyant_image_simulation(img_filename, output_filename; Nz=128)
    img = image_to_initial_conditions(img_filename)

    Ix, Iy, Iz = size(img)
    aspect = Ix / Iz
    Nx = floor(Int, aspect * 256)
    grid = RectilinearGrid(size=(Nx, Nz), z=(0, 1), x=(0, aspect),
                           topology = (Periodic, Flat, Bounded))

    bi = regrid_xy(grid, img)

    advection = WENO(order=5)
    model = NonhydrostaticModel(; grid, advection,
                                tracers = :b,
                                timestepper = :RungeKutta3,
                                buoyancy = BuoyancyTracer())

    set!(model, b=bi)

    simulation = Simulation(model, Δt=0.01, stop_time=2)

    wizard = TimeStepWizard(cfl=0.8, max_change=1.1)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5))

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

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    b = model.tracers.b
    outputs = (; b)
    simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                        schedule = IterationInterval(5),
                                                        filename = output_filename,
                                                        overwrite_existing = true)

    return simulation
end



end
