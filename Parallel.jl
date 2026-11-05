#=
Файл, содержащий средства распараллеливания двумерной задачи (пакет Polyester, макрос @batch).
Одномерная задача считается последовательно.
Включен в Main.jl через include (после загрузки конфига, так как используется константа dim).

Число потоков задается при запуске: julia -t 8 Main.jl <config>.jl

Каждая ячейка в параллельных циклах считается независимо и с тем же порядком операций,
поэтому результат совпадает с последовательной версией побитово.
=#

# --- цикл параллельный при dim == 2 и последовательный при dim == 1 (для функций, общих для обеих размерностей) ---
# --- dim - константа, поэтому лишняя ветка убирается компилятором ---
macro batch2d(loop)
    # --- в @batch передается копия цикла, чтобы последовательная ветка не зависела от того, как @batch его разбирает ---
    loop_par = copy(loop)
    return esc(quote
        if dim == 2
            @batch $loop_par
        else
            $loop
        end
    end)
end

# --- поэлементная операция out[k] = f(a[k], ...) для двумерных массивов одного размера, параллельно по k ---
# --- out может совпадать с одним из аргументов (например, a .*= b записывается как Map2D!(*, a, a, b)) ---
function Map2D!(f::F, out::AbstractMatrix, a) where {F}
    @batch for k in eachindex(out)
        @inbounds out[k] = f(a[k])
    end
    return out
end

function Map2D!(f::F, out::AbstractMatrix, a, b) where {F}
    @batch for k in eachindex(out)
        @inbounds out[k] = f(a[k], b[k])
    end
    return out
end

function Map2D!(f::F, out::AbstractMatrix, a, b, c) where {F}
    @batch for k in eachindex(out)
        @inbounds out[k] = f(a[k], b[k], c[k])
    end
    return out
end
