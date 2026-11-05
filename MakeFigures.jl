#=
Построение графиков по данным численного решения (одномерного или двумерного).

Запуск: julia MakeFigures.jl <config>.jl   (конфиг ищется в папке Configs)
Результат: каталог direct из конфига

В <config>.jl:
    Построение аналитического решения: plot_analytics = true
    Выбор графиков: requests = [:density_y, :energy, [:density_y, :velocity_y, :energy, :mu]]
        символ - отдельный график, вектор символов - сабплоты на одном графике
        профили вдоль y в столбце x = x_slice: :density_y, :velocity_y, :pressure_y
        профили вдоль x в строке y = y_slice:  :density_x, :velocity_x, :pressure_x
        двумерные поля:                        :rho1_2d, :rho2_2d, :rho_2d, :u_2d, :v_2d, :p_2d
        функции времени:                       :energy, :mu
    Профили берутся из сохраненных во фрейме профилей, а если их нет, то из двумерных полей.
    Аналитическое решение наносится на профили вдоль y (:density_y, :velocity_y) и на график центра поршня.
    Регулировка масштаба графиков: rho_max_coeff = 1.75 и т.д.
Всегда строится график положения центра поршня по столбцу x = x_slice: piston_center.png.

Одномерная задача (dim = 1): запросы вдоль x и двумерные поля игнорируются,
профили вдоль y строятся по всей одномерной области.

Фреймы читаются из data.jls потоково (в памяти одновременно только один фрейм),
функции времени собираются отдельным первым проходом.
=#

include("Imports.jl")
include("Instruments.jl")
DownloadConfig()

# --- по умолчанию задача одномерная ---
if !@isdefined(dim)
    const dim = 1
end
dim in (1, 2) || error("Недопустимое значение dim = $dim, допустимы 1 и 2")

if plot_analytics
    include("Analytics.jl")
end

ticker = pyimport("matplotlib.ticker")
plt.rcParams["font.family"] = "serif"
plt.rcParams["axes.formatter.use_mathtext"] = true

####################################################################################################

# --- известные графики и то, какие данные фрейма им нужны ---
const PROFILE_Y_REQ = (:density_y, :velocity_y, :pressure_y)
const PROFILE_X_REQ = (:density_x, :velocity_x, :pressure_x)
const FIELD_REQ     = (:rho1_2d, :rho2_2d, :rho_2d, :u_2d, :v_2d, :p_2d)
const TIME_REQ      = (:energy, :mu)

# --- в одномерной задаче запросы вдоль x и двумерные поля не строятся (игнорируются) ---
function FilterRequests1D(reqs)
    keep(sym) = !(sym in PROFILE_X_REQ || sym in FIELD_REQ)
    out = []
    for req in reqs
        if isa(req, Symbol)
            keep(req) && push!(out, req)
        else
            sub = filter(keep, req)
            isempty(sub) || push!(out, sub)
        end
    end
    return out
end

if dim == 1
    requests = FilterRequests1D(requests)
end

####################################################################################################

# --- чтение данных численного решения из data.jls ---
data_path = joinpath(direct, "data.jls")

function ReadHeader(path)
    io  = open(path, "r")
    hdr = deserialize(io)
    close(io)
    return hdr
end

# --- f(fr, k) вызывается для каждого фрейма по очереди, возвращает число фреймов ---
function ForEachFrame(f, path)
    io = open(path, "r")
    deserialize(io)   # заголовок
    k = 0
    while !eof(io)
        local fr
        try
            fr = deserialize(io)
        catch
            break
        end
        k += 1
        f(fr, k)
    end
    close(io)
    return k
end

hdr     = ReadHeader(data_path)
sl      = hdr.slices
content = hdr.frame_content

xc_in = sl.xc_in
x_in  = sl.x_in
yc_in = sl.yc_in
y_in  = sl.y_in

####################################################################################################

# --- проверка запросов до начала построения ---
function CheckRequest(sym)
    if sym in PROFILE_Y_REQ
        (:profile_y in content || :fields in content) ||
            error("Для $sym во фреймах нужны :profile_y или :fields (frame_content = $content)")
    elseif sym in PROFILE_X_REQ
        (:profile_x in content || :fields in content) ||
            error("Для $sym во фреймах нужны :profile_x или :fields (frame_content = $content)")
    elseif sym in FIELD_REQ
        :fields in content ||
            error("Для $sym во фреймах нужны :fields (frame_content = $content)")
    elseif !(sym in TIME_REQ)
        error("Неизвестная величина: $sym")
    end
    haskey(suffix, sym) || error("Для $sym нет записи в suffix")
end

for req in requests
    for sym in (isa(req, Symbol) ? (req,) : req)
        CheckRequest(sym)
    end
end
(:profile_y in content || :fields in content) ||
    error("Для графика центра поршня во фреймах нужны :profile_y или :fields")

####################################################################################################

# --- профиль вдоль y в столбце x = x_slice: из profile_y или из двумерных полей ---
# rho1, rho2, p, u -- в центрах ячеек (yc_in), v -- на гранях (y_in)
function ProfileY(fr, name)
    if fr.profile_y !== nothing
        return getproperty(fr.profile_y, name)
    end
    F  = fr.fields
    jj = sl.jj
    if name == :u
        return 0.5 .* (F.u[:, jj] .+ F.u[:, jj + 1])
    elseif name == :v
        return F.v[:, jj]
    else
        return getproperty(F, name)[:, jj]
    end
end

# --- профиль вдоль x в строке y = y_slice: из profile_x или из двумерных полей ---
# rho1, rho2, p, v -- в центрах ячеек (xc_in), u -- на гранях (x_in)
function ProfileX(fr, name)
    if fr.profile_x !== nothing
        return getproperty(fr.profile_x, name)
    end
    F  = fr.fields
    ii = sl.ii
    if name == :v
        return 0.5 .* (F.v[ii, :] .+ F.v[ii + 1, :])
    elseif name == :u
        return F.u[ii, :]
    else
        return getproperty(F, name)[ii, :]
    end
end

####################################################################################################
function PistonCenterNumeric(rho1, coord)
    N = length(rho1)
    if N < 3 return NaN end
    
    grad = zeros(Float64, N)
    grad[1] = (rho1[2] - rho1[1]) / (coord[2] - coord[1])
    for i in 2:(N-1)
        grad[i] = (rho1[i+1] - rho1[i-1]) / (coord[i+1] - coord[i-1])
    end
    grad[N] = (rho1[N] - rho1[N-1]) / (coord[N] - coord[N-1])
    
    idx_L = argmax(grad)
    idx_R = argmin(grad)
    
    x_L_num = coord[idx_L]
    x_R_num = coord[idx_R]
    
    return (x_R_num + x_L_num) / 2.0
end

ts      = Float64[]
Es      = Float64[]
mu1_min = Float64[]
mu1_max = Float64[]
mu2_min = Float64[]
mu2_max = Float64[]
cn      = Float64[]
ca      = Float64[]

ForEachFrame(data_path) do fr, k
    push!(ts, fr.t)
    push!(Es, fr.E_total)
    push!(mu1_min, fr.mu1_min)
    push!(mu1_max, fr.mu1_max)
    push!(mu2_min, fr.mu2_min)
    push!(mu2_max, fr.mu2_max)
    push!(cn, PistonCenterNumeric(ProfileY(fr, :rho1), yc_in))
    if plot_analytics == true
        push!(ca, analytic.sol(fr.t)[1])
    end
end

####################################################################################################

# --- подобранные вручную оптимальные масштабы графиков ---
rhoA = rho1A + rho2A
rhoB = rho1B + rho2B

rho_total_max = max(rhoA, rhoB)
ylim_rho_min  = rho_min
ylim_rho_max  = rho_total_max * rho_max_coeff

# --- в одномерной задаче скорость u отсутствует ---
u_max_est  = ((dim == 1 ? abs(v_init) : max(abs(u_init), abs(v_init))) + g * t_max)
ylim_u_min = -u_max_est * u_min_coeff
ylim_u_max = u_max_est * u_max_coeff

rho_max_single = max(rho1A, rho1B, rho2A, rho2B)
mu_max_est     = Apsi * rho_max_single * mu_max_coeff
ylim_mu_min    = -mu_max_est * mu_min_coeff
ylim_mu_max    = mu_max_est * mu_max_coeff

p_max_est  = Apsi * rho_max_single^2
ylim_p_min = -p_max_est * p_min_coeff
ylim_p_max = p_max_est * p_max_coeff

# --- в одномерной задаче энергия считается на единицу ширины, поэтому мера области - Ly ---
area             = dim == 1 ? Ly : Lx * Ly
e_psi_est        = Apsi * rho_max_single^2 * area
e_kin_est        = 0.5 * rho_total_max * u_max_est^2 * area
e_force_est      = rho_total_max * g * Ly * area
e_total_max_est  = e_psi_est + e_kin_est + e_force_est
e_total_min_est  = max(1.0, e_psi_est)
ylim_e_total_min = e_total_min_est * E_min_coeff
ylim_e_total_max = e_total_max_est * E_max_coeff

# --- пределы цветовой шкалы двумерных полей (постоянные для всех фреймов) ---
clim = Dict(:rho1 => (min(rho1A, rho1B), max(rho1A, rho1B)),
            :rho2 => (min(rho2A, rho2B), max(rho2A, rho2B)),
            :rho  => (min(rhoA, rhoB),   max(rhoA, rhoB)),
            :u    => (ylim_u_min, ylim_u_max),
            :v    => (ylim_u_min, ylim_u_max),
            :p    => (ylim_p_min, ylim_p_max))

####################################################################################################

# --- подписи осей ---
LabelX(ax) = ax.annotate(raw"x, $м \cdot 10^{-3}$", xy=(1.01, -0.01), xycoords="axes fraction",
                         ha="left", va="top", fontsize=12, rotation=0)
LabelY(ax) = ax.annotate(raw"y, $м \cdot 10^{-3}$", xy=(1.01, -0.01), xycoords="axes fraction",
                         ha="left", va="top", fontsize=12, rotation=0)
LabelT(ax) = ax.annotate(raw"t, $c \cdot 10^{-3}$", xy=(1.05, -0.01), xycoords="axes fraction",
                         ha="left", va="top", fontsize=12, rotation=0)
LabelV(ax, s; xpos=-0.09) = ax.annotate(s, xy=(xpos, 1.075), xycoords="axes fraction",
                                        ha="left", va="top", fontsize=12, rotation=0)

# --- профили вдоль y (dir = :y) или вдоль x (dir = :x) ---
function ProfileSetup(dir)
    if dir == :y
        place = dim == 1 ? "вдоль y" : @sprintf("вдоль y при x = %.2f", sl.x_slice)
        return ProfileY, yc_in, y_in, LabelY, place
    else
        return ProfileX, xc_in, x_in, LabelX, @sprintf("вдоль x при y = %.2f", sl.y_slice)
    end
end

function DrawDensityProfile(ax, fr, dir)
    Prof, coord_c, coord_f, Label, place = ProfileSetup(dir)
    r1 = Prof(fr, :rho1)
    r2 = Prof(fr, :rho2)
    ax.plot(coord_c, r1, color="navy", label=raw"$\rho_1$")
    ax.plot(coord_c, r2, color="crimson", label=raw"$\rho_2$")
    ax.plot(coord_c, r1 .+ r2, ls="--", color="gray", label=raw"$\rho_1+\rho_2$")

    if plot_analytics && dir == :y
        an = AnalyticProfilesY(yc_in, y_in, fr.t)
        ax.plot(yc_in, an.rho1, color="navy",   alpha=0.5, label=raw"$\rho_{1, analytic}$")
        ax.plot(yc_in, an.rho2, color="crimson",alpha=0.5, label=raw"$\rho_{2, analytic}$")
    end

    Label(ax)
    LabelV(ax, raw"$\rho$, кг/м$^3$")
    ax.set_ylim(ylim_rho_min, ylim_rho_max)
    ax.set_title(@sprintf("Плотности %s, t = %.1f мс", place, fr.t), y=down, pad=0)
    ax.grid(true)
    ax.legend()
end

# --- компонента скорости вдоль профиля на гранях, поперечная в центрах ячеек ---
# --- в одномерной задаче u отсутствует (во фрейме u = nothing) и не строится ---
function DrawVelocityProfile(ax, fr, dir)
    Prof, coord_c, coord_f, Label, place = ProfileSetup(dir)
    u = Prof(fr, :u)
    v = Prof(fr, :v)
    if dir == :y
        ax.plot(coord_f, v, color="black", label=raw"$v$")
        if u !== nothing
            ax.plot(coord_c, u, color="gray", ls="--", label=raw"$u$")
        end
        if plot_analytics
            an = AnalyticProfilesY(yc_in, y_in, fr.t)
            ax.plot(y_in, an.v, color="crimson", alpha=0.7, label=raw"$v_{analytic}$")
        end
    else
        ax.plot(coord_f, u, color="black", label=raw"$u$")
        ax.plot(coord_c, v, color="gray", ls="--", label=raw"$v$")
    end
    Label(ax)
    LabelV(ax, dim == 1 ? raw"$v$, м/с" : raw"$u, v$, м/с")
    ax.set_ylim(ylim_u_min, ylim_u_max)
    ax.set_title(@sprintf("Скорость %s, t = %.1f мс", place, fr.t), y=down, pad=0)
    ax.grid(true)
    ax.legend()
end

function DrawPressureProfile(ax, fr, dir)
    Prof, coord_c, coord_f, Label, place = ProfileSetup(dir)
    ax.plot(coord_c, Prof(fr, :p), color="indigo")
    Label(ax)
    LabelV(ax, raw"$p$, Па")
    ax.set_ylim(ylim_p_min, ylim_p_max)
    ax.set_title(@sprintf("Давление %s, t = %.1f мс", place, fr.t), y=down, pad=0)
    ax.grid(true)
end

# --- двумерное поле: rho1, rho2, p в центрах ячеек, u на гранях по x, v на гранях по y ---
function Draw2D(ax, fr, name)
    F = fr.fields
    if name == :rho
        C, X, Y, title, cmap = F.rho1 .+ F.rho2, xc_in, yc_in, raw"$\rho_1+\rho_2$, кг/м$^3$", "viridis"
    elseif name == :rho1
        C, X, Y, title, cmap = F.rho1, xc_in, yc_in, raw"$\rho_1$, кг/м$^3$", "viridis"
    elseif name == :rho2
        C, X, Y, title, cmap = F.rho2, xc_in, yc_in, raw"$\rho_2$, кг/м$^3$", "viridis"
    elseif name == :u
        C, X, Y, title, cmap = F.u, x_in, yc_in, raw"$u$, м/с", "RdBu_r"
    elseif name == :v
        C, X, Y, title, cmap = F.v, xc_in, y_in, raw"$v$, м/с", "RdBu_r"
    elseif name == :p
        C, X, Y, title, cmap = F.p, xc_in, yc_in, raw"$p$, Па", "magma"
    end
    vmin, vmax = clim[name]
    mesh = ax.pcolormesh(X, Y, C, shading="nearest", cmap=cmap, vmin=vmin, vmax=vmax)
    plt.colorbar(mesh, ax=ax)
    ax.set_aspect("equal")
    ax.set_xlabel(raw"x, $м \cdot 10^{-3}$")
    ax.set_ylabel(raw"y, $м \cdot 10^{-3}$")
    ax.set_title(@sprintf("%s при t = %.1f мс", title, fr.t))
end

function DrawEnergy(ax, k)
    ax.plot(ts[1:k], Es[1:k], color="darkred", linewidth=1.5)
    ax.set_yscale("log")
    LabelT(ax)
    LabelV(ax, raw"$E$"; xpos=-0.03)
    ax.set_xlim(0, t_max)
    ax.set_ylim(ylim_e_total_min, ylim_e_total_max)
    ax.set_title(raw"Полная энергия (логарифмическая шкала)", y=down, pad=0)
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax.grid(true, which="major", linestyle="-", linewidth=0.8, alpha=0.7)
    ax.grid(true, which="minor", linestyle="--", linewidth=0.5, alpha=0.5)
end

function DrawMu(ax, k)
    tt = ts[1:k]
    ax.plot(tt, mu1_min[1:k], color="magenta", ls="-", linewidth=1.5, label=raw"Минимум $\hat{\mu}_1$")
    ax.plot(tt, mu2_min[1:k], color="royalblue", ls="-", linewidth=1.5, label=raw"Минимум $\hat{\mu}_2$")
    ax.plot(tt, mu1_max[1:k], color="magenta", ls="--", linewidth=1.5, label=raw"Максимум $\hat{\mu}_1$")
    ax.plot(tt, mu2_max[1:k], color="royalblue", ls="--", linewidth=1.5, label=raw"Максимум $\hat{\mu}_2$")
    LabelT(ax)
    LabelV(ax, raw"$\mu$, Дж/кг")
    ax.set_xlim(0, t_max)
    ax.set_ylim(ylim_mu_min, ylim_mu_max)
    ax.set_title(raw"Минимумы и максимумы химических потенциалов", y=down, pad=0)
    ax.grid(true)
    ax.legend()
end

# --- соответствие запросов и функций построения ---
const DRAW = Dict(
    :density_y  => (ax, fr, k) -> DrawDensityProfile(ax, fr, :y),
    :density_x  => (ax, fr, k) -> DrawDensityProfile(ax, fr, :x),
    :velocity_y => (ax, fr, k) -> DrawVelocityProfile(ax, fr, :y),
    :velocity_x => (ax, fr, k) -> DrawVelocityProfile(ax, fr, :x),
    :pressure_y => (ax, fr, k) -> DrawPressureProfile(ax, fr, :y),
    :pressure_x => (ax, fr, k) -> DrawPressureProfile(ax, fr, :x),
    :rho1_2d    => (ax, fr, k) -> Draw2D(ax, fr, :rho1),
    :rho2_2d    => (ax, fr, k) -> Draw2D(ax, fr, :rho2),
    :rho_2d     => (ax, fr, k) -> Draw2D(ax, fr, :rho),
    :u_2d       => (ax, fr, k) -> Draw2D(ax, fr, :u),
    :v_2d       => (ax, fr, k) -> Draw2D(ax, fr, :v),
    :p_2d       => (ax, fr, k) -> Draw2D(ax, fr, :p),
    :energy     => (ax, fr, k) -> DrawEnergy(ax, k),
    :mu         => (ax, fr, k) -> DrawMu(ax, k),
)

function DrawQuantity(ax, sym, fr, k)
    haskey(DRAW, sym) || error("Неизвестная величина: $sym")
    DRAW[sym](ax, fr, k)
end

# --- стрелки осей и сетка (для двумерных полей не нужны) ---
function FinishAxis(ax, sym)
    if sym in FIELD_REQ
        return
    end
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
function MakeCenterPlot()
    fig = plt.figure(figsize=(9, 5))
    ax = fig.add_subplot(1, 1, 1)
    label_num = raw"Численное решение" #dim == 1 ? $y_c$" : @sprintf("Численный \$y_c\$ при x = %.2f", sl.x_slice)
    ax.plot(ts, cn, color="navy", label=label_num)
    if plot_analytics
        ca = [AnalyticState(t).y_p for t in ts]
        ax.plot(ts, ca, label=raw"Аналитическое решение")#Приближенный аналитический $y_c$, " 
    end
    LabelT(ax)
    LabelV(ax, raw"$y_c$, $м \cdot 10^{-3}$")
    ax.set_xlim(0, t_max)
    #ax.set_title(raw"Положение центра твердого тела", y=down, pad=0)
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

for req in requests
    mkpath(joinpath(direct, isa(req, Symbol) ? suffix[req] : "mp4"))
end

# --- второй проход: графики по каждому фрейму ---
# --- фреймы, дописанные в data.jls после первого прохода (расчет еще идет), пропускаются ---

ForEachFrame(data_path) do fr, k
    k <= length(ts) || return
    idx = k - 1
    for req in requests
        if isa(req, Symbol)
            MakeIndividual(req, fr, k, idx)
        else
            MakeCombined(req, fr, k, idx)
        end
    end
    flush(stdout)
end

MakeCenterPlot()

# --- сообщение о завершении построения графиков ---
println("Построение графиков завершено.")
println("Графики сохранены в директорию: ", direct)
