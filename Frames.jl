#=
Файл, содержащий функции записи фреймов численного решения в data.jls.
Включен в Main.jl через include.

Структура data.jls: заголовок FrameHeader, затем фреймы (именованные кортежи).

Двумерная задача. Содержимое фреймов задается в конфиге списком frame_content:
    :fields    - двумерные поля без фиктивных ячеек:
                 rho1, rho2, p в центрах ячеек ((Ny-1) x (Nx-1)),
                 u на гранях по x ((Ny-1) x Nx), v на гранях по y (Ny x (Nx-1))
    :profile_x - профили вдоль x в строке ячеек, ближайшей к y = y_slice:
                 rho1, rho2, p, v в центрах ячеек (координаты xc_in), u на гранях (координаты x_in)
    :profile_y - профили вдоль y в столбце ячеек, ближайшем к x = x_slice:
                 rho1, rho2, p, u в центрах ячеек (координаты yc_in), v на гранях (координаты y_in)

Одномерная задача. frame_content из конфига игнорируется, всегда пишется только :profile_y:
                 rho1, rho2, p в центрах ячеек (координаты yc_in), v на гранях (координаты y_in), u = nothing.
    Поля, не имеющие смысла в 1D (x, xc, срезы, x_in, xc_in), хранятся как nothing.

Отсутствующая часть фрейма хранится как nothing.
Скаляры (t, step, E_total, экстремумы mu_hat) пишутся всегда.
=#

const FRAME_CONTENT_KNOWN = (:fields, :profile_x, :profile_y)

####################################################################################################
# --- двумерная задача ---

# --- индексы и координаты срезов, внутренние координаты сетки ---
# i_s, j_s - индексы в полных массивах (с фиктивными ячейками), ii, jj - во внутренних
function FrameSlices(x, xc, y, yc)
    for c in frame_content
        c in FRAME_CONTENT_KNOWN || error("Неизвестный элемент frame_content: $c")
    end

    I = ficty + 1:Ny + ficty - 1
    J = fictx + 1:Nx + fictx - 1
    j_s = J[argmin(abs.(xc[J] .- x_slice))]
    i_s = I[argmin(abs.(yc[I] .- y_slice))]

    return (i_s = i_s, j_s = j_s,
            ii = i_s - ficty, jj = j_s - fictx,
            x_slice = xc[j_s], y_slice = yc[i_s],
            x_in = x[fictx + 1:Nx + fictx], xc_in = xc[J],
            y_in = y[ficty + 1:Ny + ficty], yc_in = yc[I])
end

FrameHeader(x, xc, y, yc, slices) = (x = x, xc = xc, y = y, yc = yc,
                                     slices = slices, frame_content = copy(frame_content))

# --- в frames_data хранится только последний фрейм (все фреймы уже лежат в data.jls) ---
function SaveStateFrame!(slices, rho1::AbstractMatrix, rho2, U, p, E, Phi,
                         t, step, energy_hist, time_hist, frames_data)

    E_total, E_psi, E_kinetic, E_lam, E_force = ComputeEnergy(rho1, rho2, U, E, Phi)

    push!(energy_hist, E_total)
    push!(time_hist, t)

    mu1_full = Mu1Hat!(zeros(size(rho1)), E, rho1, rho2)
    mu2_full = Mu2Hat!(zeros(size(rho1)), E, rho1, rho2)

    u  = U[1]
    v  = U[2]
    I  = ficty + 1:Ny + ficty - 1   # внутренние ячейки по y
    J  = fictx + 1:Nx + fictx - 1   # внутренние ячейки по x
    If = ficty + 1:Ny + ficty       # грани по y (узлы y_in)
    Jf = fictx + 1:Nx + fictx       # грани по x (узлы x_in)
    i_s = slices.i_s
    j_s = slices.j_s

    fields = nothing
    if :fields in frame_content
        fields = (rho1 = rho1[I, J], rho2 = rho2[I, J], p = p[I, J],
                  u = u[I, Jf], v = v[If, J])
    end

    profile_x = nothing
    if :profile_x in frame_content
        profile_x = (rho1 = rho1[i_s, J], rho2 = rho2[i_s, J], p = p[i_s, J],
                     u = u[i_s, Jf],
                     v = 0.5 .* (v[i_s, J] .+ v[i_s + 1, J]))
    end

    profile_y = nothing
    if :profile_y in frame_content
        profile_y = (rho1 = rho1[I, j_s], rho2 = rho2[I, j_s], p = p[I, j_s],
                     u = 0.5 .* (u[I, j_s] .+ u[I, j_s + 1]),
                     v = v[If, j_s])
    end

    empty!(frames_data)
    push!(frames_data, (
                    t         = t,
                    step      = step,
                    E_total   = E_total,
                    mu1_min   = minimum(mu1_full[I, J]),
                    mu1_max   = maximum(mu1_full[I, J]),
                    mu2_min   = minimum(mu2_full[I, J]),
                    mu2_max   = maximum(mu2_full[I, J]),
                    fields    = fields,
                    profile_x = profile_x,
                    profile_y = profile_y
    ))

end

####################################################################################################
# --- одномерная задача (вдоль y) ---

function FrameSlices(y, yc)
    for c in frame_content
        c in FRAME_CONTENT_KNOWN || error("Неизвестный элемент frame_content: $c")
    end

    I = ficty + 1:Ny + ficty - 1

    return (i_s = nothing, j_s = nothing,
            ii = nothing, jj = nothing,
            x_slice = nothing, y_slice = nothing,
            x_in = nothing, xc_in = nothing,
            y_in = y[ficty + 1:Ny + ficty], yc_in = yc[I])
end

FrameHeader(y, yc, slices) = (x = nothing, xc = nothing, y = y, yc = yc,
                              slices = slices, frame_content = [:profile_y])

function SaveStateFrame!(slices, rho1::AbstractVector, rho2, U, p, E, Phi,
                         t, step, energy_hist, time_hist, frames_data)

    E_total, E_psi, E_kinetic, E_lam, E_force = ComputeEnergy(rho1, rho2, U, E, Phi)

    push!(energy_hist, E_total)
    push!(time_hist, t)

    mu1_full = Mu1Hat!(zeros(size(rho1)), E, rho1, rho2)
    mu2_full = Mu2Hat!(zeros(size(rho1)), E, rho1, rho2)

    v  = U
    I  = ficty + 1:Ny + ficty - 1   # внутренние ячейки по y
    If = ficty + 1:Ny + ficty       # грани по y (узлы y_in)

    profile_y = (rho1 = rho1[I], rho2 = rho2[I], p = p[I], u = nothing, v = v[If])

    empty!(frames_data)
    push!(frames_data, (
                    t         = t,
                    step      = step,
                    E_total   = E_total,
                    mu1_min   = minimum(mu1_full[I]),
                    mu1_max   = maximum(mu1_full[I]),
                    mu2_min   = minimum(mu2_full[I]),
                    mu2_max   = maximum(mu2_full[I]),
                    fields    = nothing,
                    profile_x = nothing,
                    profile_y = profile_y
    ))

end

####################################################################################################
# --- общие для обеих размерностей ---

function MakeFrame(data_io, slices, rho1, rho2, U, E, phi,
                   t, step, frame_count, energy_hist, time_hist, frames_data)
    mu_el  = MuEl!(zeros(size(rho1)), rho1, rho2)
    lam_el = LamEl!(zeros(size(rho1)), rho1, rho2)
    p      = Pressure(E, mu_el, lam_el, rho1, rho2)

    SaveStateFrame!(slices, rho1, rho2, U, p, E, phi,
                    t, step, energy_hist, time_hist, frames_data)
    serialize(data_io, frames_data[end])
    flush(data_io)

    return frame_count + 1, mu_el, lam_el
end
