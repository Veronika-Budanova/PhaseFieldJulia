##########
# GRID
##########

function CreateGrid()

    x = zeros(N + 2 * fict)
    x = [(i - fict - 1) * h for i in eachindex(x)]

    xc = zeros(N + 2 * fict - 1)
    xc = [(i - fict - 1) * h + 0.5 * h for i in eachindex(xc)]

    return x, xc

end