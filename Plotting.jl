
function compute_energy(rho1, rho2, u, E, Phi)
    mu_el  = MuEl(muelB, rho1, rho2)
    lam_el = LamEl(lamelB, rho1, rho2)

    E_psi = sum(Psi0(
        E[fict+1:N-1+fict],
        mu_el[fict+1:N-1+fict],
        lam_el[fict+1:N-1+fict],
        rho1[fict+1:N-1+fict],
        rho2[fict+1:N-1+fict]
    )) * h

    u_x   = [u[i][1] for i in 1:N+2*fict]
    s_u_2 = s(u_x .* u_x)
    rho   = rho1 .+ rho2
    E_kinetic = 0.5 * sum(rho[fict+1:N-1+fict] .* s_u_2[fict+1:N-1+fict]) * h

    E_lam = 0.5 * sum(
        lam11 .* s(dStar(rho1) .* dStar(rho1))[fict+1:N-1+fict] .+
        lam12 .* s(dStar(rho2) .* dStar(rho1))[fict+1:N-1+fict] .+
        lam21 .* s(dStar(rho1) .* dStar(rho2))[fict+1:N-1+fict] .+
        lam22 .* s(dStar(rho2) .* dStar(rho2))[fict+1:N-1+fict]
    ) * h

    E_force = sum(rho[fict+1:N-1+fict] .* Phi[fict+1:N-1+fict]) * h
    E_total = E_psi + E_kinetic + E_lam - E_force

    return E_total, E_psi, E_kinetic, E_lam, E_force
end


function compute_momentum(rho1, rho2, u)
    u_x  = [u[i][1] for i in 1:N+2*fict]
    sRho = sStar(rho1 .+ rho2)
    return sum(sRho[fict+1:N+fict] .* u_x[fict+1:N+fict]) * h
end


function save_frame(xc, x, rho1, rho2, u, p, fict, t,
                    frame_count, time_hist, energy_hist,
                    frames_data, mu1_full, mu2_full)

    fig = plt.figure(figsize=(18,8))

    ############################
    # densities
    ############################

    ax1 = fig.add_subplot(2,3,1)

    ax1.plot(xc[fict+1:N+fict-1],rho1[fict+1:N+fict-1],label="rho1")
    ax1.plot(xc[fict+1:N+fict-1],rho2[fict+1:N+fict-1],label="rho2")
    ax1.plot(xc[fict+1:N+fict-1],rho1[fict+1:N+fict-1]+rho2[fict+1:N+fict-1],label="rho_total")

    ax1.set_xlabel("x")
    ax1.set_ylabel("Density")
    ax1.set_title(@sprintf("time = %.4f",t))
    ax1.grid(true)
    ax1.legend()

    ############################
    # energy
    ############################

    ax2 = fig.add_subplot(2,3,2)

    ax2.plot(time_hist,energy_hist,"k-",linewidth=2,label="Total energy")

    ax2.set_xlabel("time")
    ax2.set_ylabel("Total energy")
    ax2.set_title("Total energy (time)")
    ax2.grid(true)

    ############################
    # chemical potentials
    ############################

    ax3 = fig.add_subplot(2,3,3)

    mu1_min_hist = [frames_data[i][8] for i in 1:length(frames_data)]
    mu1_max_hist = [frames_data[i][9] for i in 1:length(frames_data)]

    mu2_min_hist = [frames_data[i][10] for i in 1:length(frames_data)]
    mu2_max_hist = [frames_data[i][11] for i in 1:length(frames_data)]

    times_hist = [frames_data[i][5] for i in 1:length(frames_data)]

    push!(mu1_min_hist,minimum(mu1_full[fict+1:N+fict-1]))
    push!(mu1_max_hist,maximum(mu1_full[fict+1:N+fict-1]))

    push!(mu2_min_hist,minimum(mu2_full[fict+1:N+fict-1]))
    push!(mu2_max_hist,maximum(mu2_full[fict+1:N+fict-1]))

    push!(times_hist,t)

    ax3.plot(times_hist,mu1_min_hist,"b--",linewidth=1.5,label="Mu1 min")
    ax3.plot(times_hist,mu1_max_hist,"b-",linewidth=1.5,label="Mu1 max")

    ax3.plot(times_hist,mu2_min_hist,"r--",linewidth=1.5,label="Mu2 min")
    ax3.plot(times_hist,mu2_max_hist,"r-",linewidth=1.5,label="Mu2 max")

    ax3.set_xlabel("time")
    ax3.set_ylabel("Chemical potential extrema")
    ax3.set_title("Chemical potential extrema (time)")
    ax3.grid(true)
    ax3.legend()

    ############################
    # velocity
    ############################

    ax4 = fig.add_subplot(2,3,4)

    u_x = [u[i][1] for i in fict+1:N+fict]
    ax4.plot(x[fict+1:N+fict], u_x, label="u")

    ax4.set_xlabel("x")
    ax4.set_ylabel("Velocity")
    ax4.set_title(@sprintf("time = %.4f",t))
    ax4.grid(true)
    ax4.legend()

    ############################
    # pressure
    ############################

    ax5 = fig.add_subplot(2,3,5)

    ax5.plot(xc[fict+1:N+fict-1],p[fict+1:N+fict-1],label="p")

    ax5.set_xlabel("x")
    ax5.set_ylabel("Pressure")
    ax5.set_title(@sprintf("time = %.4f",t))
    ax5.grid(true)
    ax5.legend()

    ############################
    # approximate compressibility
    ############################

    ax6 = fig.add_subplot(2,3,6)

    sStar_rho = sStar(rho1 .+ rho2)
    c_approx = [sqrt(max(0.0, Apsi * sStar_rho[i])) for i in fict+1:N+fict]
    u_x = [u[i][1] for i in fict+1:N+fict]
    ax6.plot(x[fict+1:N+fict],(u_x./c_approx) .* (u_x./c_approx),label="beta")

    ax6.set_xlabel("x")
    ax6.set_ylabel("Compressibility")
    ax6.set_title(@sprintf("time = %.4f",t))
    ax6.grid(true)
    ax6.legend()

    plt.tight_layout()
    plt.savefig(@sprintf("%s/rho_%06d.png", direct, frame_count))
    plt.close()
end

function save_state_frame!(xc, x, rho1, rho2, rho, u, p,
                            E, Phi, fict, t, step, frame_count, energy_hist,
                            time_hist, frames_data)

    ############################
    # energy
    ############################

    E_total, E_psi, E_kinetic, E_lam, E_force = compute_energy(rho1, rho2, u, E, Phi)

    println(@sprintf("E_psi=%.10e  E_kin=%.10e  E_lam=%.10e  E_force=%.10e  E_total=%.10e",
    E_psi, E_kinetic, E_lam, E_force, E_total))
    flush(stdout)

    push!(energy_hist, E_total)
    push!(time_hist, t)

    ############################
    # chemical potentials
    ############################

    mu_el  = MuEl(muelB, rho1, rho2)
    lam_el = LamEl(lamelB, rho1, rho2)

    mu1_full = Mu1Hat(E, mu_el, lam_el, rho1, rho2)
    mu2_full = Mu2Hat(E, mu_el, lam_el, rho1, rho2)

    ############################
    # frames data
    ############################

    push!(frames_data,(
                    copy(rho1[fict+1:N+fict-1]),
                    copy(rho2[fict+1:N+fict-1]),
                    copy(u[fict+1:N+fict]),
                    copy(p[fict+1:N+fict-1]),
                    t,
                    step,
                    E_total,
                    minimum(mu1_full[fict+1:N+fict-1]),
                    maximum(mu1_full[fict+1:N+fict-1]),
                    minimum(mu2_full[fict+1:N+fict-1]),
                    maximum(mu2_full[fict+1:N+fict-1])
    ))

    ############################
    # plotting
    ############################

    save_frame(xc, x, rho1, rho2, u, p, fict, t,
                frame_count, time_hist, energy_hist,
                frames_data, mu1_full, mu2_full)

    ############################
    # diagnostics
    ############################

    current_mass1 = sum(rho1[fict+1:end-fict]) * h
    current_mass2 = sum(rho2[fict+1:end-fict]) * h

    println(@sprintf("total energy = %.6e", E_total))

    println(@sprintf(
        "mass 1 = %.6f, mass 2 = %.6f, sum of masses = %.6f",
        current_mass1,
        current_mass2,
        current_mass1 + current_mass2
    ))
    flush(stdout)

end
