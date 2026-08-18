function compute_energy(rho1, rho2, u, E, Phi)
    mu_el  = MuEl!(zeros(length(rho1)), rho1, rho2)
    lam_el = LamEl!(zeros(length(rho1)), rho1, rho2)
 
    E_psi = sum(Psi0(
        E[fict+1:N-1+fict],
        mu_el[fict+1:N-1+fict],
        lam_el[fict+1:N-1+fict],
        rho1[fict+1:N-1+fict],
        rho2[fict+1:N-1+fict]
    )) * h
 
    u_x   = [u[i][1] for i in 1:N+2*fict]
    s_u_2 = s!(zeros(N + 2 * fict - 1), u_x .* u_x)
    rho   = rho1 .+ rho2
    E_kinetic = 0.5 * sum(rho[fict+1:N-1+fict] .* s_u_2[fict+1:N-1+fict]) * h
 
    dStar_rho1 = dStar!(zeros(N + 2 * fict), rho1)
    dStar_rho2 = dStar!(zeros(N + 2 * fict), rho2)
    E_lam = 0.5 * sum(
        lam11 .* s!(zeros(N + 2 * fict - 1), dStar_rho1 .* dStar_rho1)[fict+1:N-1+fict] .+
        lam12 .* s!(zeros(N + 2 * fict - 1), dStar_rho2 .* dStar_rho1)[fict+1:N-1+fict] .+
        lam21 .* s!(zeros(N + 2 * fict - 1), dStar_rho1 .* dStar_rho2)[fict+1:N-1+fict] .+
        lam22 .* s!(zeros(N + 2 * fict - 1), dStar_rho2 .* dStar_rho2)[fict+1:N-1+fict]
    ) * h
 
    E_force = sum(rho[fict+1:N-1+fict] .* Phi[fict+1:N-1+fict]) * h
    E_total = E_psi + E_kinetic + E_lam - E_force
 
    return E_total, E_psi, E_kinetic, E_lam, E_force
end

function fmt_sci(x)
    if x == 0
        return "0"
    end
    exp_val = floor(Int, log10(abs(x)))
    mantissa = x / 10.0^exp_val
    if isapprox(mantissa, round(mantissa); atol=1e-10)
        return @sprintf("%d \\cdot 10^{%d}", Int(round(mantissa)), exp_val)
    else
        return @sprintf("%.1f \\cdot 10^{%d}", mantissa, exp_val)
    end
end


function save_state_frame!(xc, x, rho1, rho2, rho, u, p,
                            E, Phi, fict, t, step, frame_count, energy_hist,
                            time_hist, frames_data)
 
    E_total, E_psi, E_kinetic, E_lam, E_force = compute_energy(rho1, rho2, u, E, Phi)
 
    push!(energy_hist, E_total)
    push!(time_hist, t)
 
    mu_el  = MuEl!(zeros(length(rho1)), rho1, rho2)
    lam_el = LamEl!(zeros(length(rho1)), rho1, rho2)
 
    buf = CreateBuffers()
    mu1_full = zeros(length(rho1))
    mu2_full = zeros(length(rho1))
    Mu1Hat!(mu1_full, buf, E, rho1, rho2)
    Mu2Hat!(mu2_full, buf, E, rho1, rho2)
 
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
 
end