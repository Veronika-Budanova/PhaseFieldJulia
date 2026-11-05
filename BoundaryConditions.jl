#=
Файл, содержащий функции применения граничных условий типа жесткой стенки и периодики для скаляров, векторов и тензоров.
Включен в Main.jl через include.

Методы для матриц - двумерная задача, для векторов - одномерная задача вдоль y.
В двумерной задаче вектор - пара матриц [u, v] (или кортеж (u, v)), в одномерной - вектор v на гранях по y.
Тензорное поле - матрица (2D) или вектор (1D) матриц Tensor3 (неизменяемых, в ячейку записывается целая матрица).
=#

# --- знаки компонент тензора при отражении относительно стенки x = const и y = const ---
const SIGN_X = Tensor3([ 1.0 -1.0 -1.0;
                        -1.0  1.0  1.0;
                        -1.0  1.0  1.0])
const SIGN_Y = Tensor3([ 1.0 -1.0  1.0;
                        -1.0  1.0 -1.0;
                         1.0 -1.0  1.0])

####################################################################################################
# --- двумерная задача ---

function BoundCondScalar(A::AbstractMatrix{<:Real})
    if bound_x == "periodic_x"
        for j in 1:fictx, i in 1:Ny + 2 * ficty - 1
            A[i, j] = A[i, Nx + j - 1]
            A[i, Nx + fictx + j - 1] = A[i, fictx + j]
        end
    elseif bound_x == "wall_x"
        for j in 1:fictx, i in 1:Ny + 2 * ficty - 1
            A[i, fictx - j + 1] = A[i, fictx + j]
            A[i, Nx + fictx + j - 1] = A[i, Nx + fictx - j]
        end
    end
 
    if bound_y == "periodic_y"
        for j in 1:Nx + 2 * fictx - 1, i in 1:ficty
            A[i, j] = A[Ny + i - 1, j]
            A[Ny + ficty + i - 1, j] = A[ficty + i, j]
        end
    elseif bound_y == "wall_y"
        for j in 1:Nx + 2 * fictx - 1, i in 1:ficty
            A[ficty - i + 1, j] = A[ficty + i, j]
            A[Ny + ficty + i - 1, j] = A[Ny + ficty - i, j]
        end
    end
    return A
end

function BoundCondVector(U::Union{AbstractVector{<:AbstractMatrix}, Tuple{AbstractMatrix, AbstractMatrix}})
    u = U[1]
    v = U[2]
 
    if bound_x == "periodic_x"
        for j in 1:fictx, i in 1:Ny + 2 * ficty - 1
            u[i, j] = u[i, Nx + j - 1]
            u[i, Nx + fictx + j] = u[i, fictx + j + 1]
        end
        for i in 1:Ny + 2 * ficty - 1
            u[i, Nx + fictx] = u[i, fictx + 1]
        end
        for j in 1:fictx, i in 1:Ny + 2 * ficty
            v[i, j] = v[i, Nx + j - 1]
            v[i, Nx + fictx + j - 1] = v[i, fictx + j]
        end
    elseif bound_x == "wall_x"
        for j in 1:fictx, i in 1:Ny + 2 * ficty - 1
            u[i, fictx - j + 1] = -u[i, fictx + j + 1]
            u[i, Nx + fictx + j] = -u[i, Nx + fictx - j]
        end
        for i in 1:Ny + 2 * ficty - 1
            u[i, fictx + 1] = 0.0
            u[i, Nx + fictx] = 0.0
        end
        for j in 1:fictx, i in 1:Ny + 2 * ficty
            v[i, fictx - j + 1] = v[i, fictx + j]
            v[i, Nx + fictx + j - 1] = v[i, Nx + fictx - j]
        end
    end
 
    if bound_y == "periodic_y"
        for j in 1:Nx + 2 * fictx, i in 1:ficty
            u[i, j] = u[Ny + i - 1, j]
            u[Ny + ficty + i - 1, j] = u[ficty + i, j]
        end
        for j in 1:Nx + 2 * fictx - 1, i in 1:ficty
            v[i, j] = v[Ny + i - 1, j]
            v[Ny + ficty + i, j] = v[ficty + i + 1, j]
        end
        for j in 1:Nx + 2 * fictx - 1
            v[Ny + ficty, j] = v[ficty + 1, j]
        end
    elseif bound_y == "wall_y"
        for j in 1:Nx + 2 * fictx, i in 1:ficty
            u[ficty - i + 1, j] = u[ficty + i, j]
            u[Ny + ficty + i - 1, j] = u[Ny + ficty - i, j]
        end
        for j in 1:Nx + 2 * fictx - 1, i in 1:ficty
            v[ficty - i + 1, j] = -v[ficty + i + 1, j]
            v[Ny + ficty + i, j] = -v[Ny + ficty - i, j]
        end
        for j in 1:Nx + 2 * fictx - 1
            v[ficty + 1, j] = 0.0
            v[Ny + ficty, j] = 0.0
        end
    end
    return U
end

function BoundCondTensor(A::AbstractMatrix{<:AbstractMatrix})
    if bound_x == "periodic_x"
        for j in 1:fictx, i in 1:Ny + 2 * ficty - 1
            A[i, j] = A[i, Nx + j - 1]
            A[i, Nx + fictx + j - 1] = A[i, fictx + j]
        end
    elseif bound_x == "wall_x"
        for j in 1:fictx, i in 1:Ny + 2 * ficty - 1
            A[i, fictx - j + 1] = SIGN_X .* A[i, fictx + j]
            A[i, Nx + fictx + j - 1] = SIGN_X .* A[i, Nx + fictx - j]
        end
    end
 
    if bound_y == "periodic_y"
        for j in 1:Nx + 2 * fictx - 1, i in 1:ficty
            A[i, j] = A[Ny + i - 1, j]
            A[Ny + ficty + i - 1, j] = A[ficty + i, j]
        end
    elseif bound_y == "wall_y"
        for j in 1:Nx + 2 * fictx - 1, i in 1:ficty
            A[ficty - i + 1, j] = SIGN_Y .* A[ficty + i, j]
            A[Ny + ficty + i - 1, j] = SIGN_Y .* A[Ny + ficty - i, j]
        end
    end
    return A
end

####################################################################################################
# --- одномерная задача (вдоль y) ---

function BoundCondScalar(A::AbstractVector{<:Real})
    if bound_y == "periodic_y"
        for i in 1:ficty
            A[i] = A[Ny + i - 1]
            A[Ny + ficty + i - 1] = A[ficty + i]
        end
    elseif bound_y == "wall_y"
        for i in 1:ficty
            A[ficty - i + 1] = A[ficty + i]
            A[Ny + ficty + i - 1] = A[Ny + ficty - i]
        end
    end
    return A
end

function BoundCondVector(v::AbstractVector{<:Real})
    if bound_y == "periodic_y"
        for i in 1:ficty
            v[i] = v[Ny + i - 1]
            v[Ny + ficty + i] = v[ficty + i + 1]
        end
        v[Ny + ficty] = v[ficty + 1]
    elseif bound_y == "wall_y"
        for i in 1:ficty
            v[ficty - i + 1] = -v[ficty + i + 1]
            v[Ny + ficty + i] = -v[Ny + ficty - i]
        end
        v[ficty + 1] = 0.0
        v[Ny + ficty] = 0.0
    end
    return v
end

function BoundCondTensor(A::AbstractVector{<:AbstractMatrix})
    if bound_y == "periodic_y"
        for i in 1:ficty
            A[i] = A[Ny + i - 1]
            A[Ny + ficty + i - 1] = A[ficty + i]
        end
    elseif bound_y == "wall_y"
        for i in 1:ficty
            A[ficty - i + 1] = SIGN_Y .* A[ficty + i]
            A[Ny + ficty + i - 1] = SIGN_Y .* A[Ny + ficty - i]
        end
    end
    return A
end
