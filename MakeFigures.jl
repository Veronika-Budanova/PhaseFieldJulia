#=
Построение графиков распределения любых величин.

Запуск: julia MakeFigures.jl <config>.jl
Результат: Figures/<config>

В <config>.jl:
    Построение аналитического решения: plot_analytics = true
    Выбор величин для построения графиков: requests = [:density, :energy, :mu, [:density, :energy, :velocity, :mu]]
    Регулировка масштаба графиков: rho_max_coeff = 1.75 и т.д.
=#

include("Imports.jl")
include("Instruments.jl")
DownloadConfig()
include("Initialization.jl")  

if plot_analytics
    include("Analytics.jl")
end

ticker = pyimport("matplotlib.ticker")
plt.rcParams["font.family"] = "serif"
plt.rcParams["axes.formatter.use_mathtext"] = true

####################################################################################################

# --- чтение данных численного решения из data.jls ---
io     = open(joinpath(direct, "data.jls"), "r")
grids  = deserialize(io)
frames = []
while !eof(io)
    try
        push!(frames, deserialize(io))
    catch
        break
    end
end
close(io)

ts = [fr[5] for fr in frames]
if plot_analytics
    xs = [sol(fr[5])[1] for fr in frames]
end

x  = grids.x
xc = grids.xc

xc_in = xc[fict+1:N+fict-1]
x_in  = x[fict+1:N+fict]

####################################################################################################

# --- подобранные вручную оптимальные масштабы графиков ---
rho_total_max = max(rho1A + rho2A, rho1B + rho2B)
ylim_rho_min  = rho_min
ylim_rho_max  = rho_total_max * rho_max_coeff

u_max_est  = (abs(u_init) + g * t_max)
ylim_u_min = -u_max_est * u_min_coeff
ylim_u_max = u_max_est * u_max_coeff

rho_max_single = max(rho1A, rho1B, rho2A, rho2B)
mu_max_est     = Apsi * rho_max_single * mu_max_coeff
ylim_mu_min    = -mu_max_est * mu_min_coeff
ylim_mu_max    = mu_max_est * mu_max_coeff

p_max_est  = Apsi * rho_max_single^2
ylim_p_min = -p_max_est * p_min_coeff
ylim_p_max = p_max_est * p_max_coeff

e_psi_est        = Apsi * rho_max_single^2 * L
e_kin_est        = 0.5 * rho_total_max * (abs(u_init) + g * t_max)^2 * L
e_force_est      = rho_total_max * g * L^2
e_total_max_est  = e_psi_est + e_kin_est + e_force_est
e_total_min_est  = max(1.0, e_psi_est)
ylim_e_total_min = e_total_min_est * E_min_coeff
ylim_e_total_max = e_total_max_est * E_max_coeff

####################################################################################################

# --- построение графиков ---
function DrawDensity(ax, fr, k)
    rho1 = fr[1]
    rho2 = fr[2]
    t = fr[5]
    ax.plot(xc_in, rho1, color="navy", label=raw"$\rho_1$")
    ax.plot(xc_in, rho2, color="crimson", label=raw"$\rho_2$")
    ax.plot(xc_in, rho1 .+ rho2, ls="--", color="gray", label=raw"$\rho_1+\rho_2$")

    if plot_analytics && k <= length(ts)
        xc_piston = x0_piston + xs[k]
        a1, a2 = BuildDensity(xc_in, xc_piston)
        ax.plot(xc_in, a1, color="navy",   ls="--", alpha=0.5, label=raw"$\rho_{1, analytic}$")
        ax.plot(xc_in, a2, color="crimson", ls="--", alpha=0.5, label=raw"$\rho_{2, analytic}$")
    end

    ax.annotate(raw"x, $м \cdot 10^{-3}$",
                xy=(1.01, -0.01),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.annotate(raw"$\rho$, кг/м$^3$",
                xy=(-0.09, 1.075),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.set_ylim(ylim_rho_min, ylim_rho_max)
    ax.set_title(@sprintf("Плотности при t = %.1f мс", t), y=down, pad=0)
    ax.grid(true)
    ax.legend()
end

function DrawVelocity(ax, fr)
    u = fr[3]
    t = fr[5]
    u_x = [u[i][1] for i in 1:length(u)]
    ax.plot(x_in, u_x, color="black")
    ax.annotate(raw"x, $м \cdot 10^{-3}$",
                xy=(1.01, -0.01),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.annotate(raw"$u$, м/с",
                xy=(-0.09, 1.075),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.set_ylim(ylim_u_min, ylim_u_max)
    ax.set_title(@sprintf("Скорость при t = %.1f мс", t), y=down, pad=0)
    ax.grid(true)
end

function DrawPressure(ax, fr)
    p = fr[4]
    t = fr[5]
    ax.plot(xc_in, p, color="indigo")
    ax.annotate(raw"x, $м \cdot 10^{-3}$",
                xy=(1.01, -0.01),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.annotate(raw"$p$, Па",
                xy=(-0.09, 1.075),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.set_ylim(ylim_p_min, ylim_p_max)
    ax.set_title(@sprintf("Давление при t = %.1f мс", t), y=down, pad=0)
    ax.grid(true)
end

function DrawEnergy(ax, k)
    times = [frames[i][5] for i in 1:k]
    energies = [frames[i][7] for i in 1:k]
    ax.plot(times, energies, color="darkred", linewidth=1.5)
    ax.set_yscale("log")
    ax.annotate(raw"t, $c \cdot 10^{-3}$",
                xy=(1.05, -0.01),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.annotate(raw"$E$",
                xy=(-0.03, 1.075),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.set_xlim(0, t_max)
    ax.set_ylim(ylim_e_total_min, ylim_e_total_max)
    ax.set_title(raw"Полная энергия (логарифмическая шкала)", y=down, pad=0)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax.grid(true, which="major", linestyle="-", linewidth=0.8, alpha=0.7)
    ax.grid(true, which="minor", linestyle="--", linewidth=0.5, alpha=0.5)
end

function DrawMu(ax, k)
    times = [frames[i][5] for i in 1:k]
    mu1_min = [frames[i][8] for i in 1:k]
    mu1_max = [frames[i][9] for i in 1:k]
    mu2_min = [frames[i][10] for i in 1:k]
    mu2_max = [frames[i][11] for i in 1:k]
    ax.plot(times, mu1_min, color="magenta", ls="-", linewidth=1.5, label=raw"Минимум $\hat{\mu}_1$")
    ax.plot(times, mu2_min, color="royalblue", ls="-", linewidth=1.5, label=raw"Минимум $\hat{\mu}_2$")
    ax.plot(times, mu1_max, color="magenta", ls="--", linewidth=1.5, label=raw"Максимум $\hat{\mu}_1$")
    ax.plot(times, mu2_max, color="royalblue", ls="--", linewidth=1.5, label=raw"Максимум $\hat{\mu}_2$")
    ax.annotate(raw"t, $c \cdot 10^{-3}$",
                xy=(1.05, -0.01),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.annotate(raw"$\mu$, Дж/кг",
                xy=(-0.09, 1.075),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.set_xlim(0, t_max)
    ax.set_ylim(ylim_mu_min, ylim_mu_max)
    ax.set_title(raw"Минимумы и максимумы химических потенциалов", y=down, pad=0)
    ax.grid(true)
    ax.legend()
end

function DrawQuantity(ax, sym, fr, k)
    if sym == :density
        DrawDensity(ax, fr, k)
    elseif sym == :velocity
        DrawVelocity(ax, fr)
    elseif sym == :pressure
        DrawPressure(ax, fr)
    elseif sym == :energy
        DrawEnergy(ax, k)
    elseif sym == :mu
        DrawMu(ax, k)
    else
        error("Неизвестная величина: $sym")
    end
end

function FinishAxis(ax, sym)
    ax.spines["right"].set_visible(false)
    ax.spines["top"].set_visible(false)
    ax.plot([1.0], [0], "k>", transform=ax.transAxes, markersize=8, clip_on=false)
    ax.plot([0], [1.0], "^k", transform=ax.transAxes, markersize=8, clip_on=false)
    if sym != :energy
        ax.xaxis.set_major_locator(ticker.AutoLocator())
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax.yaxis.set_major_locator(ticker.AutoLocator())
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax.grid(true, which="major", linestyle="-", linewidth=0.8, alpha=0.7)
        ax.grid(true, which="minor", linestyle="--", linewidth=0.5, alpha=0.5)
    end
end

####################################################################################################

# --- построение графика зависимости положения центра поршня от времени ---
function PistonCenterNumeric(rho1)
    w = clamp.((rho1 .- rho1A) ./ (rho1B - rho1A), 0.0, 1.0)
    sw = sum(w)
    return sw == 0.0 ? NaN : sum(w .* xc_in) / sw
end

function MakeCenterPlot()
    nfr = plot_analytics ? min(length(frames), length(xs)) : length(frames)
    tt = [frames[i][5] for i in 1:nfr]
    cn = [PistonCenterNumeric(frames[i][1]) for i in 1:nfr]
    fig = plt.figure(figsize=(7.5, 4))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(tt, cn, color="navy", label=raw"Численный $x_c$")
    if plot_analytics
        ca = [x0_piston + xs[i] for i in 1:nfr]
        ax.plot(tt, ca, color="crimson", ls="--", label=raw"Приближенный аналитический $x_c$")
    end
    ax.annotate(raw"t, $c \cdot 10^{-3}$",
                xy=(1.05, -0.01),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.annotate(raw"$x_c$, $м \cdot 10^{-3}$",
                xy=(-0.09, 1.075),
                xycoords="axes fraction",
                ha="left",
                va="top",
                fontsize=12,
                rotation=0)
    ax.set_xlim(0, t_max)
    ax.set_title(raw"Положение центра твердого тела", y=down, pad=0)
    ax.grid(true)
    ax.legend()
    FinishAxis(ax, :center)
    plt.tight_layout()
    plt.savefig(joinpath(direct, "piston_center.png"), bbox_inches="tight")
    plt.close(fig)
end

####################################################################################################

# --- отдельные графики и несколько сабплотов на одном графике ---
function MakeIndividual(sym, fr, k, idx)
    fig = plt.figure(figsize=(7.5, 4))
    ax = fig.add_subplot(1, 1, 1)
    DrawQuantity(ax, sym, fr, k)
    FinishAxis(ax, sym)
    plt.tight_layout()
    plt.savefig(joinpath(direct, suffix[sym], @sprintf("rho_%05d_%s.png", idx, suffix[sym])), bbox_inches="tight")
    plt.close(fig)
end

function MakeCombined(syms, fr, k, idx)
    n = length(syms)
    ncols = 2
    nrows = ceil(Int, n / ncols)
    fig = plt.figure(figsize=(7.5 * ncols, 4 * nrows))
    for (j, sym) in enumerate(syms)
        ax = fig.add_subplot(nrows, ncols, j)
        DrawQuantity(ax, sym, fr, k)
        FinishAxis(ax, sym)
    end
    plt.tight_layout()
    plt.savefig(joinpath(direct, "mp4", @sprintf("rho_%05d_mp4.png", idx)), bbox_inches="tight")
    plt.close(fig)
end

for sub in ("rho", "E", "mu", "mp4")
    mkpath(joinpath(direct, sub))
end

for (i, fr) in enumerate(frames)
    idx = i - 1
    for req in requests
        if isa(req, Symbol)
            MakeIndividual(req, fr, i, idx)
        else
            MakeCombined(req, fr, i, idx)
        end
    end
    flush(stdout)
end

MakeCenterPlot()