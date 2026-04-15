# BOUNDARY CONDITIONS

function BoundCondScalar(A)
    for i in 1:fict
        if bound == "periodic"
            if isa(A[i], Number)
                A[i] = A[N + i - 1]
                A[N + fict + i - 1] = A[fict + i]
            else
                A[i] .= A[N + i - 1]
                A[N + fict + i - 1] .= A[fict + i]
            end
        elseif bound == "wall"
            if isa(A[i], Number)
                A[fict - i + 1] = A[fict + i]
                A[N + fict + i - 1] = A[N + fict - i]
            else
                A[fict - i + 1] .= A[fict + i]
                A[N + fict + i - 1] .= A[N + fict - i]
            end
        end
    end

    return A

end

function BoundCondVector(A) 
    if bound == "periodic"
        for i in 1:fict
            if isa(A[i], Number)
                A[i] = A[N + i]
                A[N + fict + i] = A[fict + i + 1]
            else 
                A[i] .= A[N + i]
                A[N + fict + i] .= A[fict + i + 1]
            end
        end
        A[N + fict] = A[fict + 1]
    elseif bound == "wall"
        for i in 1:fict
            if isa(A[i], Number)
                A[fict - i + 1] = -A[fict + i + 1]
                A[N + fict + i] = -A[N + fict - i]
            else
                A[fict - i + 1] .= -A[fict + i + 1]
                A[N + fict + i] .= -A[N + fict - i]
            end
        end

        if isa(A[fict + 1], Number)
            A[fict + 1] = 0
        else
            A[fict + 1] .= 0
        end

        if isa(A[N + fict], Number)
            A[N + fict] = 0
        else
            A[N + fict] .= 0
        end
    end      
    return A
end

function BoundCondFlux(A)
    if bound == "periodic"
        for i in 1:fict
            A[i] = A[N + i - 1]
            A[N + fict + i] = A[fict + i + 1]
        end
        A[N + fict] = A[fict + 1]
    elseif bound == "wall"
        for i in 1:fict
            A[fict - i + 1] = A[fict + i + 1]
            A[N + fict + i] = A[N + fict - i]
        end
        A[fict + 1] = 0
        A[N + fict] = 0
    end
    return A
end