#=
Файл, содержащий функции создания узлов и центров ячеек равномерной сетки по данным из конфиг-файла.
Включен в Main.jl через include.
=#

function CreateGrid2D()
    x = zeros(Nx + 2 * fictx)
    x = [(i - fictx - 1) * hx for i in eachindex(x)]

    xc = zeros(Nx + 2 * fictx - 1)
    xc = [(i - fictx - 1) * hx + 0.5 * hx for i in eachindex(xc)]

    y = zeros(Ny + 2 * ficty)
    y = [(j - ficty - 1) * hy for j in eachindex(y)]

    yc = zeros(Ny + 2 * ficty - 1)
    yc = [(j - ficty - 1) * hy + 0.5 * hy for j in eachindex(yc)]

    return x, xc, y, yc
end

# --- одномерная задача: сетка только вдоль y ---
function CreateGrid1D()
    y = zeros(Ny + 2 * ficty)
    y = [(j - ficty - 1) * hy for j in eachindex(y)]

    yc = zeros(Ny + 2 * ficty - 1)
    yc = [(j - ficty - 1) * hy + 0.5 * hy for j in eachindex(yc)]

    return y, yc
end
