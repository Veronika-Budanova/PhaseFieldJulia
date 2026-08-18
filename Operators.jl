# s, d, s*, d*
function s!(out, array)
    @inbounds for i in (fict + 1):(N + fict - 1)
        out[i] = 0.5 * (array[i + 1] + array[i])
    end
    return out
end

function d!(out, array)
    inv_h = 1 / h
    @inbounds for i in (fict + 1):(N + fict - 1)
        out[i] = inv_h * (array[i + 1] - array[i])
    end
    return out
end

function sStar!(out, array)
    @inbounds for i in (fict + 1):(N + fict)
        out[i] = 0.5 * (array[i] + array[i - 1])
    end
    return out
end

function dStar!(out, array)
    inv_h = 1 / h
    @inbounds for i in (fict + 1):(N + fict)
        out[i] = inv_h * (array[i] - array[i - 1])
    end
    return out
end


