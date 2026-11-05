#=
Файл, содержащий функции задания начальных условий.
Плотности можно задать с плавными и резкими границами раздела фаз, а также в виде смеси для спинодального распада.
Только в двумерной задаче: квадратная капля ("square_drop") и шарик ("ball") фазы B в фазе A.
Включен в Main.jl через include.

Функции с суффиксом 2D/1D - для соответствующей размерности задачи (одномерная задача - вдоль y).

Тензорные поля (E, K, Pi_el, Pi_ns) - массивы матриц Tensor3 (SMatrix 3x3 из StaticArrays).
Матрицы лежат в массиве подряд, без указателей; они неизменяемые, поэтому в ячейку записывается целая матрица.
=#

# --- тип тензора в ячейке: статическая матрица 3x3 ---
const Tensor3 = SMatrix{3, 3, Float64, 9}

# --- тензорное поле из чекпоинта в виде массива Tensor3 (чекпоинты старого формата хранят массив Matrix{Float64}) ---
TensorField(A) = map(Tensor3, A)

####################################################################################################
# --- двумерная задача ---

function InitVelocity2D()
    u = fill(u_init, Ny + 2*ficty - 1, Nx + 2*fictx)
    v = fill(v_init, Ny + 2*ficty, Nx + 2*fictx - 1)
    return [u, v]
end

function InitDensityChoise(xc::AbstractVector, yc::AbstractVector)
    if rho_init == "sharp" || rho_init == "smooth"
        return InitDensitiesPistonImitation2D(yc)
    elseif rho_init == "spinodal"
        return InitDensitiesSpinodal(xc, yc)
    elseif rho_init == "square_drop"
        return InitDensitiesSquareDrop2D(xc, yc)
    elseif rho_init == "ball"
        return InitDensitiesBall2D(xc, yc)
    end
end

function InitDensitiesPistonImitation2D(yc)
    rho1, rho2 = BuildDensity2D(yc, y0_piston)
    rho = rho1 .+ rho2
    return rho1, rho2, rho
end

function BuildDensity2D(yc, yc_piston)
    rho1 = fill(0.0, Ny + 2*ficty - 1, Nx + 2*fictx - 1)
    rho2 = fill(0.0, Ny + 2*ficty - 1, Nx + 2*fictx - 1)

    for j in 1:Nx + 2 * fictx - 1
        for i in 1:Ny + 2 * ficty - 1
            f = PistonFraction(yc[i], yc_piston)
            rho1[i, j] = rho1A + (rho1B - rho1A) * f
            rho2[i, j] = rho2A + (rho2B - rho2A) * f
        end
    end
    return rho1, rho2
end

function InitDensitiesSpinodal(xc::AbstractVector, yc::AbstractVector)
    rho1 = fill(0.0, Ny + 2*ficty - 1, Nx + 2*fictx - 1)
    rho2 = fill(0.0, Ny + 2*ficty - 1, Nx + 2*fictx - 1)

    rho1_mean = 0.5 * (rho1A + rho1B)
    rho2_mean = 0.5 * (rho2A + rho2B)
    amp = 0.02
    Random.seed!(123)

    for j in 1:Nx + 2 * fictx - 1
        for i in 1:Ny + 2 * ficty - 1
            rho1[i, j] = rho1_mean * (1.0 + amp * (2.0 * rand() - 1.0))
            rho2[i, j] = 1000 - rho1[i, j]
        end
    end
    rho = rho1 .+ rho2
    return rho1, rho2, rho

end

# --- квадратная капля фазы B в фазе A: сторона piston_width, центр (x0_piston, y0_piston) ---
# --- под действием поверхностного натяжения капля должна стать круглой ---
function InitDensitiesSquareDrop2D(xc, yc)
    rho1 = fill(0.0, Ny + 2*ficty - 1, Nx + 2*fictx - 1)
    rho2 = fill(0.0, Ny + 2*ficty - 1, Nx + 2*fictx - 1)

    for j in 1:Nx + 2 * fictx - 1
        for i in 1:Ny + 2 * ficty - 1
            f = SquareDropFraction(xc[j], yc[i])
            rho1[i, j] = rho1A + (rho1B - rho1A) * f
            rho2[i, j] = rho2A + (rho2B - rho2A) * f
        end
    end
    rho = rho1 .+ rho2
    return rho1, rho2, rho
end

# --- шарик фазы B в фазе A: диаметр piston_width, центр (x0_piston, y0_piston) ---
# --- E = 0, поэтому начальная круглая форма - ненапряженное состояние, упругие силы удерживают ее при падении ---
function InitDensitiesBall2D(xc, yc)
    rho1 = fill(0.0, Ny + 2*ficty - 1, Nx + 2*fictx - 1)
    rho2 = fill(0.0, Ny + 2*ficty - 1, Nx + 2*fictx - 1)

    for j in 1:Nx + 2 * fictx - 1
        for i in 1:Ny + 2 * ficty - 1
            f = BallFraction(xc[j], yc[i])
            rho1[i, j] = rho1A + (rho1B - rho1A) * f
            rho2[i, j] = rho2A + (rho2B - rho2A) * f
        end
    end
    rho = rho1 .+ rho2
    return rho1, rho2, rho
end

# --- доля фазы B в квадрате: произведение плавных ступенек по x и по y ---
function SquareDropFraction(x_c, y_c)
    if bound_x == "periodic_x"
        distx = mod(x_c - x0_piston + 0.5 * Lx, Lx) - 0.5 * Lx
    elseif bound_x == "wall_x"
        distx = x_c - x0_piston
    end

    if bound_y == "periodic_y"
        disty = mod(y_c - y0_piston + 0.5 * Ly, Ly) - 0.5 * Ly
    elseif bound_y == "wall_y"
        disty = y_c - y0_piston
    end

    return 0.5 * (1.0 + tanh(b * (0.5 * piston_width + distx))) *
           0.5 * (1.0 + tanh(b * (0.5 * piston_width - distx))) *
           0.5 * (1.0 + tanh(b * (0.5 * piston_width + disty))) *
           0.5 * (1.0 + tanh(b * (0.5 * piston_width - disty)))
end

# --- доля фазы B в круге: плавная ступенька по расстоянию до центра ---
function BallFraction(x_c, y_c)
    if bound_x == "periodic_x"
        distx = mod(x_c - x0_piston + 0.5 * Lx, Lx) - 0.5 * Lx
    elseif bound_x == "wall_x"
        distx = x_c - x0_piston
    end

    if bound_y == "periodic_y"
        disty = mod(y_c - y0_piston + 0.5 * Ly, Ly) - 0.5 * Ly
    elseif bound_y == "wall_y"
        disty = y_c - y0_piston
    end

    r = sqrt(distx^2 + disty^2)
    return 0.5 * (1.0 + tanh(b * (0.5 * piston_width - r)))
end

function InitAlmansiTens2D()
    return fill(zero(Tensor3), Ny + 2*ficty - 1, Nx + 2*fictx - 1)
end

####################################################################################################
# --- одномерная задача (вдоль y) ---

# --- скорость v (вдоль y) на гранях ячеек ---
function InitVelocity1D()
    v = fill(v_init, Ny + 2*ficty)
    return v
end

function InitDensityChoise(yc::AbstractVector)
    if rho_init == "sharp" || rho_init == "smooth"
        return InitDensitiesPistonImitation1D(yc)
    elseif rho_init == "spinodal"
        return InitDensitiesSpinodal(yc)
    elseif rho_init == "square_drop" || rho_init == "ball"
        error("Начальная плотность rho_init = \"$rho_init\" задается только в двумерной задаче (dim = 2)")
    end
end

function InitDensitiesPistonImitation1D(yc)
    rho1, rho2 = BuildDensity1D(yc, y0_piston)
    rho = rho1 .+ rho2
    return rho1, rho2, rho
end

function BuildDensity1D(yc, yc_piston)
    rho1 = fill(0.0, Ny + 2*ficty - 1)
    rho2 = fill(0.0, Ny + 2*ficty - 1)

    for i in 1:Ny + 2 * ficty - 1
        f = PistonFraction(yc[i], yc_piston)
        rho1[i] = rho1A + (rho1B - rho1A) * f
        rho2[i] = rho2A + (rho2B - rho2A) * f
    end
    return rho1, rho2
end

function InitDensitiesSpinodal(yc::AbstractVector)
    rho1 = fill(0.0, Ny + 2*ficty - 1)
    rho2 = fill(0.0, Ny + 2*ficty - 1)

    rho1_mean = 0.5 * (rho1A + rho1B)
    rho2_mean = 0.5 * (rho2A + rho2B)
    amp = 0.02
    Random.seed!(123)

    for i in 1:Ny + 2 * ficty - 1
        rho1[i] = rho1_mean * (1.0 + amp * (2.0 * rand() - 1.0))
        rho2[i] = 1000 - rho1[i]
    end
    rho = rho1 .+ rho2
    return rho1, rho2, rho

end

function InitAlmansiTens1D()
    return fill(zero(Tensor3), Ny + 2*ficty - 1)
end

####################################################################################################
# --- общие для обеих размерностей ---

function PistonFraction(y_c, yc_piston)
    if bound_y == "periodic_y"
        disty = mod(y_c - yc_piston + 0.5 * Ly, Ly) - 0.5 * Ly
    elseif bound_y == "wall_y"
        disty = y_c - yc_piston
    end

    if rho_init == "sharp"
        return abs(disty) <= 0.5 * piston_width ? 1.0 : 0.0
    elseif rho_init == "smooth"
        return 0.5 * (1.0 + tanh(b * (0.5 * piston_width + disty))) *
               0.5 * (1.0 + tanh(b * (0.5 * piston_width - disty)))
    end
end