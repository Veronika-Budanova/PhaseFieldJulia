include("Imports.jl")
include("Config.jl")
include("Grid.jl")
include("Initialization.jl")
include("BoundaryConditions.jl")
include("Operators.jl")
include("Thermodynamics.jl")
include("Stresses.jl")
include("Solver.jl")
include("Timestep.jl")
include("Plotting.jl")

function main(io)

    ##########
    # GRID
    ##########

    x, xc = CreateGrid()

    ##########
    # INITIAL CONDITIONS
    ##########
    u = InitVelocity()
    rho1, rho2, rho = InitDensitiesPistonImitationSmooth(xc)
    E = InitAlmansiTens()

    ##########
    # BOUNDARY CONDITIONS
    ##########
    
    rho1 = BoundCondScalar(rho1)
    rho2 = BoundCondScalar(rho2)
    rho = BoundCondScalar(rho)
    u = BoundCondVector(u)
    E = BoundCondScalar(E)
    
    Phi = -g .* xc

    ##########
    # TIME VARIABLES
    ##########

    t = 0.0
    step = 0

    next_frame_t = frame_dt
    frame_count = 0

    ##########
    # ARRAYS
    ##########

    energy_hist = []
    time_hist = []
    frames_data = []

    if !isdir(direct)
        mkdir(direct)
    end

    println("Начало:")
    init_mass1 = sum(rho1[fict + 1:N + fict - 1]) * h
    init_mass2 = sum(rho2[fict + 1:N + fict - 1]) * h
    println("   Initial mass 1 = ", init_mass1)
    println("   Initial mass 2 = ", init_mass2)
    println("   Initial sum of masses = ", init_mass1 + init_mass2)
    println("   t_max = $t_max")
    println("   step_max = $step_max")
    println("   frame_dt = $frame_dt")
    flush(stdout)

    mu_el = MuEl(muelB, rho1, rho2)
	lam_el = LamEl(lamelB, rho1, rho2)

    p = Pressure(E, mu_el, lam_el, rho1, rho2)

    ##########
    # INITIAL FRAME
    ##########

    save_state_frame!(xc, x, rho1, rho2, rho, u, p, E, Phi, fict,
                        t, step, frame_count, energy_hist, time_hist,
                        frames_data)

    frame_count += 1
    
    ##########
    # MAIN LOOP
    ##########

    while step < step_max

        if t >= t_max
            println("Достигнут t_max")
            flush(stdout)
            break
        end

        dt = delta_t

        if t + dt > t_max
            dt = t_max - t
        end

        rho1_new, rho2_new, u_new = ConservMassAndMomentumStep(rho1,rho2,u,E,Phi,dt)

        E_new = ConservAlmansiTensor(E,u,dt)

        rho1 = rho1_new
        rho2 = rho2_new
        rho = rho1 .+ rho2

        u = u_new
        E = E_new

        rho1 = BoundCondScalar(rho1)
        rho2 = BoundCondScalar(rho2)
        rho = BoundCondScalar(rho)
        u = BoundCondVector(u)
        E = BoundCondScalar(E)
        #=
        if step < 50 || step % 1000 == 0
            println("step=$step t=$(round(t,sigdigits=4)) min(rho1)=$(round(minimum(rho1),sigdigits=4)) min(rho2)=$(round(minimum(rho2),sigdigits=4))")
            flush(stdout)
        end
        =#

        if !CheckStability(rho1,rho2,u,step,t,E,mu_el,lam_el,Phi)
            break
        end

        if t >= next_frame_t - 1e-10

            mu_el = MuEl(muelB, rho1, rho2)
            lam_el = LamEl(lamelB, rho1, rho2)

            p = Pressure(E, mu_el, lam_el, rho1, rho2)

            save_state_frame!(xc, x, rho1, rho2, rho, u, p, E, Phi, fict, t,
                                step, frame_count, energy_hist, time_hist, frames_data)

            frame_count += 1
            next_frame_t += frame_dt

            progress_pct = 100.0 * (t / t_max)

            println(@sprintf("Прогресс %.1f%% | шаг %d | время %.5f ",
                            progress_pct, step, t))
            flush(io)
            cp("log.txt", joinpath(direct, "log.txt"), force=true)
            cp("Config.jl", joinpath(direct, "Config.jl"), force=true)
        end

        t += dt
        step += 1

    end

    println("Всего шагов: $step")
    println("Финальное время: $t")
    flush(io)
    cp("log.txt", joinpath(direct, "log.txt"), force=true)
    cp("Config.jl", joinpath(direct, "Config.jl"), force=true)

end

open("log.txt", "w") do io
    redirect_stdout(io) do
        main(io)  
    end
end


