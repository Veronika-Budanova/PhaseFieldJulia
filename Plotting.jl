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


function save_frame(xc, x, rho1, rho2, u, p, fict, t,
                    frame_count, time_hist, energy_hist,
                    frames_data, mu1_full, mu2_full)

    fig = plt.figure(figsize=(15,12))
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["axes.formatter.use_mathtext"] = true
    ticker = pyimport("matplotlib.ticker")

    down = -0.11


    ############################
    # глобальные оценки y_lim
    ############################
    rho_total_max = max(rho1A + rho2A, rho1B + rho2B)
    ylim_rho_max = rho_total_max * 1.5

    u_min_est = (abs(u_init) + g * t_max) * 0.5
    u_max_est = (abs(u_init) + g * t_max) * 0.8

    rho_max_single = max(rho1A, rho1B, rho2A, rho2B)
    mu_max_est = Apsi * rho_max_single * 1.3

    p_max_est = Apsi * rho_max_single^2 * 1.1

    E_psi_est = Apsi * rho_max_single^2 * L
    E_kin_est = 0.5 * rho_total_max * (abs(u_init) + g * t_max)^2 * L
    E_force_est = rho_total_max * g * L^2
    E_total_max_est = (E_psi_est + E_kin_est + E_force_est) * 0.12
    E_total_min_est = max(1.0, 0.03 * E_psi_est)

    ############################
    # densities
    ############################
    ax1 = fig.add_subplot(3,2,1)

    ax1.plot(xc[fict+1:N+fict-1],rho1[fict+1:N+fict-1],color="navy",label=raw"$\rho_1$")
    ax1.plot(xc[fict+1:N+fict-1],rho2[fict+1:N+fict-1],color="crimson",label=raw"$\rho_2$")
    ax1.plot(xc[fict+1:N+fict-1],rho1[fict+1:N+fict-1]+rho2[fict+1:N+fict-1],ls="--",color="gray",label=raw"$\rho_1+\rho_2$")
    #ax1.set_xlabel(raw"x, мм", loc="right")
    ax1.annotate(raw"x, $м \cdot 10^{-3}$", 
                xy=(1.01, -0.01),  # x=1.0 (правый край), y отрицательный (ниже оси)
                xycoords="axes fraction",
                ha="left",       
                va="top",         
                fontsize=12,
                rotation=0)
    #ax1.set_ylabel(raw"кг/м$^3$", loc="top")
    ax1.annotate(raw"$\rho$, кг/м$^3$", 
                xy=(-0.09, 1.075),  # x отрицательный (слева от оси), y=1.0 (верх)
                xycoords="axes fraction",
                ha="left",        
                va="top",         
                fontsize=12,
                rotation=0)      
    ax1.set_ylim(-50, ylim_rho_max*1.5)
    ax1.set_title(@sprintf("Плотности при t = %.1f мс", t), y=down, pad=0)
    ax1.grid(true)
    ax1.legend()

    ############################
    # energy  
    ############################
    ax2 = fig.add_subplot(3,2,2)

    #if length(energy_hist) >= 2
        #E_ref = energy_hist[end-1]
        #log_diff = [log(e#= - E_ref + 1e-30=#) for e in energy_hist]
        ax2.plot(time_hist, energy_hist, color="darkred", linewidth=1.5)
        ax2.set_yscale("log")
    #end
    #ax2.set_xlabel(raw"время, мс")
    ax2.annotate(raw"t, $c \cdot 10^{-3}$", 
                xy=(1.05, -0.01),  # x=1.0 (правый край), y отрицательный (ниже оси)
                xycoords="axes fraction",
                ha="left",       
                va="top",         
                fontsize=12,
                rotation=0)
    #ax2.set_ylabel(raw"$E$")
    ax2.annotate(raw"$E$", 
                xy=(-0.03, 1.075),  # x отрицательный (слева от оси), y=1.0 (верх)
                xycoords="axes fraction",
                ha="left",        
                va="top",         
                fontsize=12,
                rotation=0)  
    #ax2.set_ylabel(raw"$\ln(E - E_{\mathrm{предпослед}})$")
    ax2.set_xlim(0, t_max)
    ax2.set_ylim(E_total_min_est, E_total_max_est)
    ax2.set_title(raw"Полная энергия (логарифмическая шкала)", y=down, pad=0)
    ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax2.grid(true, which="major", linestyle="-", linewidth=0.8, alpha=0.7)
    ax2.grid(true, which="minor", linestyle="--", linewidth=0.5, alpha=0.5)


    ############################
    # velocity
    ############################
    ax3 = fig.add_subplot(3,2,3)

    u_x = [u[i][1] for i in fict+1:N+fict]
    ax3.plot(x[fict+1:N+fict], u_x, color="black")
    #ax3.set_xlabel(raw"x, мм")
    ax3.annotate(raw"x, $м \cdot 10^{-3}$", 
                xy=(1.01, -0.01),  # x=1.0 (правый край), y отрицательный (ниже оси)
                xycoords="axes fraction",
                ha="left",       
                va="top",         
                fontsize=12,
                rotation=0)
    #ax3.set_ylabel(raw"м/с")
    ax3.annotate(raw"$u$, м/с", 
                xy=(-0.09, 1.075),  # x отрицательный (слева от оси), y=1.0 (верх)
                xycoords="axes fraction",
                ha="left",        
                va="top",         
                fontsize=12,
                rotation=0)  
    ax3.set_ylim(-u_max_est, u_min_est)
    ax3.set_title(@sprintf("Скорость при t = %.1f мс", t), y=down, pad=0)
    ax3.grid(true)

    ############################
    # chemical potentials
    ############################
    ax4 = fig.add_subplot(3,2,4)

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

    ax4.plot(times_hist,mu1_min_hist,color="magenta",ls="-",linewidth=1.5,label=raw"Минимум $\hat{\mu}_1$")
    ax4.plot(times_hist,mu2_min_hist,color="royalblue",ls="-",linewidth=1.5,label=raw"Минимум $\hat{\mu}_2$")
    ax4.plot(times_hist,mu1_max_hist,color="magenta",ls="--",linewidth=1.5,label=raw"Максимум $\hat{\mu}_1$")
    ax4.plot(times_hist,mu2_max_hist,color="royalblue",ls="--",linewidth=1.5,label=raw"Максимум $\hat{\mu}_2$")
    #ax4.set_xlabel(raw"время, мс")
    ax4.annotate(raw"t, $c \cdot 10^{-3}$", 
                xy=(1.05, -0.01),  # x=1.0 (правый край), y отрицательный (ниже оси)
                xycoords="axes fraction",
                ha="left",       
                va="top",         
                fontsize=12,
                rotation=0)
    #ax4.set_ylabel(raw"Дж/кг")
    ax4.annotate(raw"$\mu$, Дж/кг", 
                xy=(-0.09, 1.075),  # x отрицательный (слева от оси), y=1.0 (верх)
                xycoords="axes fraction",
                ha="left",        
                va="top",         
                fontsize=12,
                rotation=0)  
    ax4.set_xlim(0, t_max)
    ax4.set_ylim(-mu_max_est, mu_max_est*1.2)
    ax4.set_title(raw"Минимумы и максимумы химических потенциалов", y=down, pad=0)
    ax4.grid(true)
    ax4.legend()

    ############################
    # pressure
    ############################
    ax5 = fig.add_subplot(3,2,5)

    ax5.plot(xc[fict+1:N+fict-1],p[fict+1:N+fict-1],color="indigo")
    #ax5.set_xlabel(raw"x, мм")
    ax5.annotate(raw"x, $м \cdot 10^{-3}$", 
                xy=(1.01, -0.01),  # x=1.0 (правый край), y отрицательный (ниже оси)
                xycoords="axes fraction",
                ha="left",       
                va="top",         
                fontsize=12,
                rotation=0)
    #ax5.set_ylabel(raw"Па")
    ax5.annotate(raw"$p$, Па", 
                xy=(-0.09, 1.075),  # x отрицательный (слева от оси), y=1.0 (верх)
                xycoords="axes fraction",
                ha="left",        
                va="top",         
                fontsize=12,
                rotation=0)  
    ax5.set_ylim(-p_max_est, p_max_est*1.5)
    ax5.set_title(@sprintf("Давление при t = %.1f мс", t), y=down, pad=0)
    ax5.grid(true)

    ############################
    # параметры расчёта
    ############################
    ax6 = fig.add_subplot(3,2,6)
    ax6.axis("off")

    u_arrow = u_init >= 0 ? raw"$\Rightarrow$" : raw"$\Leftarrow$"
    g_arrow = g >= 0 ? raw"$\Leftarrow$" : raw"$\Rightarrow$"

    line1 = @sprintf("\$|u_0| = %.2f\$ м/с    %s", abs(u_init), u_arrow)
    line2 = @sprintf("\$|g| = %.2f \\cdot 10^{3}\$ м/с\$^2\$  %s", g, g_arrow)
    line3 = "\$A_\\psi = B_\\psi = " * fmt_sci(Apsi) * "\$ м⁵/(кг·с²)"
    line4 = @sprintf("\$\\eta = %.1f \\cdot 10^{-3}\$ Па·с", eta)
    line5 = @sprintf("\$\\zeta = %.1f \\cdot 10^{-3}\$ Па·с", zeta)
    line6 = @sprintf("\$\\lambda_{el} = \\mu_{el} = %.1f\$ Па", lamelB)
    line7 = "\$\\lambda_{11} = \\lambda_{22} = " * fmt_sci(lam11) * " \\cdot 10^{-6}\$ м\$^7\$/(кг\$\\cdot\$с\$^2\$)"
    line8 = @sprintf("\$M_0 = %.1f \\cdot 10^{-3}\$ с·кг·К/м\$^3\$", M0)
    line9 = @sprintf("\$N = %d\$", N)
    line10 = @sprintf("\$L = %.1f\\cdot10^{-3}\$ м", L)

    params_text = "Параметры расчёта:\n" *
                  "\n" *
                  line1 * "\n" *
                  line2 * "\n" *
                  "\n" *
                  line3 * "\n" *
                  line4 * "\n" *
                  line5 * "\n" *
                  line6 * "\n" *
                  line7 * "\n" *
                  line8 * "\n" *
                  line9 * "\n" *
                  line10

    #=
    ax6.text(0.5, 0.5, params_text,
             transform=ax6.transAxes,
             fontsize=14,
             verticalalignment="center",
             horizontalalignment="center",
             fontfamily="serif",
             bbox=Dict("boxstyle"=>"round,pad=1.0", "facecolor"=>"white", "edgecolor"=>"gray", "alpha"=>1.0))
    =#
    # Для варианта с ручным созданием рамки:
    text_obj = ax6.text(0.5, 0.5, params_text,
                        transform=ax6.transAxes,
                        fontsize=14,
                        verticalalignment="center",
                        horizontalalignment="center",
                        fontfamily="serif")

    # Получаем bounding box текста
    renderer = fig.canvas.get_renderer()
    bbox = text_obj.get_window_extent(renderer=renderer)
    bbox = bbox.transformed(ax6.transData.inverted())

    # Создаем рамку с помощью patches.Rectangle
    padding_x = 0.165
    padding_y = 0.015
    rect = patches.Rectangle((bbox.x0 - padding_x, bbox.y0 - padding_y),
                            bbox.width + 2*padding_x,
                            bbox.height + 2*padding_y,
                            facecolor="white",
                            edgecolor="gray",
                            alpha=1.0,
                            transform=ax6.transData,
                            clip_on=false)
    ax6.add_patch(rect)
    text_obj.set_zorder(10)

    ############################
    # стрелки на осях
    ############################
    function add_axes_arrows(ax)
        # Скрываем правую и верхнюю spines
        ax.spines["right"].set_visible(false)
        ax.spines["top"].set_visible(false)
        
        # Получаем позиции spines
        spine_bottom = ax.spines["bottom"]
        spine_left = ax.spines["left"]
        
        # Добавляем стрелку на конце bottom spine
        ax.plot([1.0], [0], "k>", transform=ax.transAxes, markersize=8, clip_on=false)
        
        # Добавляем стрелку на конце left spine
        ax.plot([0], [1.0], "^k", transform=ax.transAxes, markersize=8, clip_on=false)
    end
    
    # Применяем ко всем подграфикам
    for ax in [ax1, ax2, ax3, ax4, ax5]
        add_axes_arrows(ax)
    end

    ############################
    # общее форматирование
    ############################
    for ax in fig.get_axes()
        if ax != ax2
            ax.xaxis.set_major_locator(ticker.AutoLocator())
            ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
            ax.yaxis.set_major_locator(ticker.AutoLocator())
            ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))

            ax.grid(true, which="major", linestyle="-", linewidth=0.8, alpha=0.7)
            ax.grid(true, which="minor", linestyle="--", linewidth=0.5, alpha=0.5)
        end
    end

    #ax2.grid(true, which="both", linestyle="-", linewidth=0.5, alpha=0.5)

    plt.tight_layout()
    plt.savefig(@sprintf("%s/rho_%05d.png", direct, frame_count))
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