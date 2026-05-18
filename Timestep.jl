function compute_dt_adaptive(h, rho1, rho2, u, M0, theta, lam11, lam22, eta, CFL, E)
    #=
    Mloc = M(M0, rho1, rho2)
    Mmax = maximum(Mloc)
    lambda_max = max(abs(lam11), abs(lam22))

    dt_diffusion = if Mmax * lambda_max > 1e-20
        CFL * h^4 / (Mmax * lambda_max / theta)
    else
        1e10
    end

    u_max = maximum(abs.([u[i][1] for i in eachindex(u)]))
    dt_convection = if u_max > 1e-10
        CFL * h / u_max
    else
        1e10
    end

    rho_min_val = max(minimum(rho1 .+ rho2), 1.0)
    dt_viscous = CFL * h^2 * rho_min_val / (eta + 1e-10)

    rho_min = max(minimum(rho1 .+ rho2), 1.0)

    # Ограничение от градиентных членов в импульсе
    dt_momentum_chem = CFL * h^2 / (lambda_max / rho_min + 1e-20)

    # Упругая волна с локальными параметрами
    mu_el = MuEl(muelB, rho1, rho2)
	lam_el = LamEl(lamelB, rho1, rho2)

    c2_max = maximum(lam_el .+ 2 .* mu_el) / rho_min
    dt_elastic = if c2_max > 1e-10
        CFL * h / sqrt(c2_max)
    else
        1e10
    end

    # Ограничение от градиента напряжений
    K_cur = K(E, mu_el, lam_el, rho1, rho2)
    Pi_el_cur = PiEl(E, K_cur)
    Pi_ns_cur = PiNS(u)
    d_Pi = dStar([Pi_ns_cur[i][1,1] + Pi_el_cur[i][1,1] for i in eachindex(E)])
    d_Pi_max = maximum(abs.(d_Pi)) + 1e-20
    dt_stress = CFL * h * rho_min / d_Pi_max

    dt = min(dt_diffusion, dt_convection, dt_viscous, dt_momentum_chem, dt_elastic, dt_stress)

    #println("  dt_diff=$(round(dt_diffusion,sigdigits=3)) | dt_conv=$(round(dt_convection,sigdigits=3)) | dt_visc=$(round(dt_viscous,sigdigits=3)) | dt_mom=$(round(dt_momentum_chem,sigdigits=3)) | dt_el=$(round(dt_elastic,sigdigits=3)) → dt=$(round(dt,sigdigits=3))")


    return max(dt, 1e-8)
    =#
    return 1e-6
end