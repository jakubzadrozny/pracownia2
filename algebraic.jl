# ZAMIANA JEDNOSTEK -- przyjemniejsze obliczenia
function toalgebraic(sat; onedim::Bool=false)
    newsat = deepcopy(sat)
    if(onedim)
        newsat[4] *= 10^3
        newsat[1:3] /= 10^2
    else
        for i = 1:length(sat)
            newsat[i][4] *= 10^3
            newsat[i][1:3] /= 10^2
        end
    end
    newsat
end

# PODSTAWIENIE POD INNA WSPOLRZEDNA
function algebraic2(sats)
    x = [sat[1] for sat in sats]
    leftx = [sat[2:4] for sat in sats]
    cvect = [1, 1, -c^2]
    A = hcat([-cvect .* (leftx[i]-leftx[4]) for i = 1:3]...)'
    B = x[1:3] .- x[4]
    sum4 = sum(cvect .* (leftx[4].^2))
    C = 0.5 * [x[4]^2-x[i]^2 - sum(cvect .* (leftx[i].^2)) + sum4 for i = 1:3]

    a = A \ B
    b = A \ C
    bdiff = b - leftx[4]

    coeffs = Array{Float64}(3)
    coeffs[1] = sum(cvect .* (bdiff .^ 2)) + x[4]^2
    coeffs[2] = 2*(-x[4] + sum(cvect .* a .* bdiff))
    coeffs[3] = 1 + sum(cvect .* (a.^2))
    guesses = filter(x -> isreal(x), roots(Poly(coeffs)))

    location = vcat(guesses[1], [a[i]*guesses[1] + b[i] for i = 1:3])
    if(length(guesses) > 1)
        location2 = vcat(guesses[2], [a[i]*guesses[2] + b[i] for i = 1:3])
        if(sum(abs2, location2) < sum(abs2, location))
            location = location2
        end
    end
    location
#     vcat(location[1:3]*10^2, location[4]/10^3)
end
