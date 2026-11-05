#=
Файл, содержащий функции применения разностных и интерполирующих операторов.
Включен в Main.jl через include.

Операторы с ! записывают результат в заранее выделенный массив out (используются на каждом шаге по времени).
Значения out вне расчетного диапазона индексов обнуляются, поэтому результат совпадает с выделяющей версией.
Операторы без ! выделяют новый массив (используются при записи фреймов и подсчете энергии).

Методы для матриц (AbstractMatrix) - двумерная задача, для векторов (AbstractVector) - одномерная задача вдоль y.
Двумерные методы параллельны по столбцам j (Polyester, @batch), одномерные - последовательные.
=#

####################################################################################################
# --- двумерная задача ---

# --- операторы по x: расчетные столбцы cols, остальные столбцы out обнуляются ---

function sx!(out, array::AbstractMatrix)
    cols = (fictx + 1):(Nx + fictx - 1)
    @batch for j in axes(out, 2)
        if j in cols
            @inbounds for i in axes(array, 1)
                out[i, j] = 0.5 * (array[i, j + 1] + array[i, j])
            end
        else
            @inbounds for i in axes(out, 1)
                out[i, j] = 0.0
            end
        end
    end
    return out
end

function dx!(out, array::AbstractMatrix)
    cols   = (fictx + 1):(Nx + fictx - 1)
    inv_hx = 1 / hx
    @batch for j in axes(out, 2)
        if j in cols
            @inbounds for i in axes(array, 1)
                out[i, j] = inv_hx * (array[i, j + 1] - array[i, j])
            end
        else
            @inbounds for i in axes(out, 1)
                out[i, j] = 0.0
            end
        end
    end
    return out
end

function sStarx!(out, array::AbstractMatrix)
    cols = (fictx + 1):(Nx + fictx)
    @batch for j in axes(out, 2)
        if j in cols
            @inbounds for i in axes(array, 1)
                out[i, j] = 0.5 * (array[i, j] + array[i, j - 1])
            end
        else
            @inbounds for i in axes(out, 1)
                out[i, j] = 0.0
            end
        end
    end
    return out
end

function dStarx!(out, array::AbstractMatrix)
    cols   = (fictx + 1):(Nx + fictx)
    inv_hx = 1 / hx
    @batch for j in axes(out, 2)
        if j in cols
            @inbounds for i in axes(array, 1)
                out[i, j] = inv_hx * (array[i, j] - array[i, j - 1])
            end
        else
            @inbounds for i in axes(out, 1)
                out[i, j] = 0.0
            end
        end
    end
    return out
end

# --- операторы по y: в каждом столбце расчетные строки rows, остальные строки out обнуляются ---

function sy!(out, array::AbstractMatrix)
    rows = (ficty + 1):(Ny + ficty - 1)
    @batch for j in axes(out, 2)
        @inbounds begin
            for i in first(axes(out, 1)):(first(rows) - 1)
                out[i, j] = 0.0
            end
            for i in rows
                out[i, j] = 0.5 * (array[i + 1, j] + array[i, j])
            end
            for i in (last(rows) + 1):last(axes(out, 1))
                out[i, j] = 0.0
            end
        end
    end
    return out
end

function dy!(out, array::AbstractMatrix)
    rows   = (ficty + 1):(Ny + ficty - 1)
    inv_hy = 1 / hy
    @batch for j in axes(out, 2)
        @inbounds begin
            for i in first(axes(out, 1)):(first(rows) - 1)
                out[i, j] = 0.0
            end
            for i in rows
                out[i, j] = inv_hy * (array[i + 1, j] - array[i, j])
            end
            for i in (last(rows) + 1):last(axes(out, 1))
                out[i, j] = 0.0
            end
        end
    end
    return out
end

function sStary!(out, array::AbstractMatrix)
    rows = (ficty + 1):(Ny + ficty)
    @batch for j in axes(out, 2)
        @inbounds begin
            for i in first(axes(out, 1)):(first(rows) - 1)
                out[i, j] = 0.0
            end
            for i in rows
                out[i, j] = 0.5 * (array[i, j] + array[i - 1, j])
            end
            for i in (last(rows) + 1):last(axes(out, 1))
                out[i, j] = 0.0
            end
        end
    end
    return out
end

function dStary!(out, array::AbstractMatrix)
    rows   = (ficty + 1):(Ny + ficty)
    inv_hy = 1 / hy
    @batch for j in axes(out, 2)
        @inbounds begin
            for i in first(axes(out, 1)):(first(rows) - 1)
                out[i, j] = 0.0
            end
            for i in rows
                out[i, j] = inv_hy * (array[i, j] - array[i - 1, j])
            end
            for i in (last(rows) + 1):last(axes(out, 1))
                out[i, j] = 0.0
            end
        end
    end
    return out
end

sx(array::AbstractMatrix)     = sx!(zeros(size(array, 1), Nx + 2 * fictx - 1), array)
dx(array::AbstractMatrix)     = dx!(zeros(size(array, 1), Nx + 2 * fictx - 1), array)
sStarx(array::AbstractMatrix) = sStarx!(zeros(size(array, 1), Nx + 2 * fictx), array)
dStarx(array::AbstractMatrix) = dStarx!(zeros(size(array, 1), Nx + 2 * fictx), array)
sy(array::AbstractMatrix)     = sy!(zeros(Ny + 2 * ficty - 1, size(array, 2)), array)
dy(array::AbstractMatrix)     = dy!(zeros(Ny + 2 * ficty - 1, size(array, 2)), array)
sStary(array::AbstractMatrix) = sStary!(zeros(Ny + 2 * ficty, size(array, 2)), array)
dStary(array::AbstractMatrix) = dStary!(zeros(Ny + 2 * ficty, size(array, 2)), array)

####################################################################################################
# --- одномерная задача (вдоль y) ---

function sy!(out, array::AbstractVector)
    fill!(out, 0.0)
    @inbounds for i in (ficty + 1):(Ny + ficty - 1)
        out[i] = 0.5 * (array[i + 1] + array[i])
    end
    return out
end

function dy!(out, array::AbstractVector)
    fill!(out, 0.0)
    inv_hy = 1 / hy
    @inbounds for i in (ficty + 1):(Ny + ficty - 1)
        out[i] = inv_hy * (array[i + 1] - array[i])
    end
    return out
end

function sStary!(out, array::AbstractVector)
    fill!(out, 0.0)
    @inbounds for i in (ficty + 1):(Ny + ficty)
        out[i] = 0.5 * (array[i] + array[i - 1])
    end
    return out
end

function dStary!(out, array::AbstractVector)
    fill!(out, 0.0)
    inv_hy = 1 / hy
    @inbounds for i in (ficty + 1):(Ny + ficty)
        out[i] = inv_hy * (array[i] - array[i - 1])
    end
    return out
end

sy(array::AbstractVector)     = sy!(zeros(Ny + 2 * ficty - 1), array)
dy(array::AbstractVector)     = dy!(zeros(Ny + 2 * ficty - 1), array)
sStary(array::AbstractVector) = sStary!(zeros(Ny + 2 * ficty), array)
dStary(array::AbstractVector) = dStary!(zeros(Ny + 2 * ficty), array)
