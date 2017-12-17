using Polynomials
c = 299792.458;
earth_rad = 6370.0
sat_rad = 20000.0

# --------------- UTILS -------------------------

function randSgn()
    x = rand(0:1)
    if x==0
        return -1
    else
        return 1
    end
end

function rand_position(radius)
    x = randSgn() * rand() * radius
    y = randSgn() * rand() * (radius - abs(x))
    z = randSgn() * sqrt(radius^2 - x^2 - y^2)
    t = randSgn() * rand() / 10
    [x, y, z, t]
end

function prepsat!(x, sat, inaccurate=false)
    for i = 1:length(sat)
        sat[i][4] = (norm(x[1:3]-sat[i][1:3]) / c + x[4]) * inaccuracy(inaccurate)
    end
end

function table(arr::Matrix; pre=true, column_names=nothing)
    if column_names !== nothing
        arr = vcat(reshape(column_names, 1, length(column_names)), arr)
    end
    s =
    """
    <table border="1">
    """
    for r in 1:size(arr, 1)
        s *= "<tr>"
        for c in 1:size(arr, 2)
            s *= "<td>"
            if pre s *= "<pre>" end
            s *= string(arr[r, c])
            if pre s *= "</pre>" end
        end
        s *= "</tr>"
    end
    s *= "</table>"
    return HTML(s)
end;

function path(maxiter, lon0, lat0)
    a = pi/4
    d = 0.00000005
    lon = lon0 + d*cos(a)
    lat = lat0 + d*sin(a)
    t=0
    res = [[], []]
    for i in 1:maxiter
        push!(res[1], lon)
        push!(res[2], lat)
        if t==0
            t = rand(20:50)
            da = rand()*randSgn()
            a +=da
        end
        lon += d*cos(a)
        lat += d*sin(a)
        t-=1
    end
    res
end

function LLAtoXYZ(LL)
    rad = 6378.1
    lon = LL[1]
    lat = LL[2]
    cosLat = cos(lat * pi / 180.0)
    sinLat = sin(lat * pi / 180.0)
    cosLon = cos(lon * pi / 180.0)
    sinLon = sin(lon * pi / 180.0)
    x = rad * cosLat * cosLon
    y = rad * cosLat * sinLon
    z = rad * sinLat
    [x, y, z]
end

function XYZtoLLA(X)
    rad = 6378.1
    x = X[1]
    y = X[2]
    z = X[3]
    lat = asin(z / rad) * 180 / pi
    lon = atan(y / x) * 180 / pi
    return [lon, lat]
end

inacc_coeff = 250000000

function inaccuracy(inaccurate)
    if inaccurate
        return (1 + rand()/inacc_coeff)
    end
    return 1
end

function countTime(mistake, Pos, sat, inaccurate=false)
    C = 299792.458
    dist = norm(Pos[1:3]-sat[1:3])
    t = dist / C
    return [t + mistake]*inaccuracy(inaccurate)
end

function createSats(lon, lat, mistake, inaccuracy=false)
    [ vcat(sat_pos, countTime(mistake, LLAtoXYZ([lon, lat]), sat_pos, inaccuracy) ) for sat_pos in sat_coords ]
end

function GPS_newton(coords, mistake, maxiter=20, inaccuracy=false)
    res = [[],[]]
    for i in 1:length(coords[1])
        sat = createSats(coords[1][i], coords[2][i], mistake, inaccuracy)
        X = newton(sat[1:4], maxiter)
        LL = XYZtoLLA(X[1:3])
        push!(res[1], LL[1])
        push!(res[2], LL[2])
    end
    return res
end

function GPS_leastSquares(coords, mistake, maxiter=20, satcnt=sat_count, inaccuracy=false)
    res = [[],[]]
    for i in 1:length(coords[1])
        sat = createSats(coords[1][i], coords[2][i], mistake, inaccuracy)
        X = leastSquares(sat, maxiter, satcnt)
        LL = XYZtoLLA(X[1:3])
        push!(res[1], LL[1])
        push!(res[2], LL[2])
    end
    return res
end

function GPS_alg(coords, mistake, inaccuracy=false)
    res = [[],[]]
    for i in 1:length(coords[1])
        sat = createSats(coords[1][i], coords[2][i], mistake, inaccuracy)
        X = algebraic(sat[1:4])
        LL = XYZtoLLA(X[1:3])
        push!(res[1], LL[1])
        push!(res[2], LL[2])
    end
    return res
end

function GPS_bancroft(coords, mistake, inaccuracy=false)
    res = [[],[]]
    for i in 1:length(coords[1])
        sat = createSats(coords[1][i], coords[2][i], mistake, inaccuracy)
        X = bancroft(sat[1:4])
        LL = XYZtoLLA(X[1:3])
        push!(res[1], LL[1])
        push!(res[2], LL[2])
    end
    return res
end

function GPS_heura(coords, mistake, inaccuracy=false)
    res = [[],[]]
    for i in 1:length(coords[1])
        sat = createSats(coords[1][i], coords[2][i], mistake, inaccuracy)
        X = heura(sat[1:4])
        LL = XYZtoLLA(X[1:3])
        push!(res[1], LL[1])
        push!(res[2], LL[2])
    end
    return res
end

function MAXdist(coords, gps_coords)
    res = 0.0
    for i in 1:length(coords[1])
        X1 = LLAtoXYZ([coords[1][i], coords[2][i]])
        X2 = LLAtoXYZ([gps_coords[1][i], gps_coords[2][i]])
        d =sqrt(sum([(X1[i]-X2[i])^2   for i in 1:3] ))
        res = max(res, d)
    end
    return res*1000
end

function MINdist(coords, gps_coords)
    res = 100000.0
    for i in 1:length(coords[1])
        X1 = LLAtoXYZ([coords[1][i], coords[2][i]])
        X2 = LLAtoXYZ([gps_coords[1][i], gps_coords[2][i]])
        d =sqrt(sum([(X1[i]-X2[i])^2   for i in 1:3] ))
        res = min(res, d)
    end
    return res*1000
end

function MinMaxMistakeTable(methods, results, coords)
    maxmistake = [ MAXdist(coords, results[i]) for i in 1:length(methods) ]
    minmistake = [ MINdist(coords, results[i]) for i in 1:length(methods) ]
    table(hcat(methods, maxmistake, minmistake), column_names=[:Metoda, :"Największy błąd", :"Najmniejszy błąd"])
end;

# --------------- NEWTON METHOD -----------------------

function product(a, b)
    a[1]*b[1] + a[2]*b[2] + a[3]*b[3] - a[4]*b[4]
end

function f(x, A)
    [ product(x-A[i, :], x-A[i, :]) for i=1:size(A, 1) ]
end

function jacobian(x, A)
    i = vcat(ones(size(A, 2) - 1), -1)
    hcat([ i .* (x - A[j, :]) for j = 1:size(A, 1) ]...)'
end

function newton(sats, maxiter=20)
    A = hcat(sats...)'
    A[:, 4] *= c
    x = zeros(4)
    for i = 1:maxiter
        J = jacobian(x, A)
        b = -0.5 * f(x, A)
        x += pinv(J) * b
    end
    x[4] /= c
    x
end

# -------------- LEAST SQUARES METHOD ---------------------

function leastSquares(sats, maxiter=20, satcnt=7)
    newton(sats[1:satcnt], maxiter)
end

# -------------- ALGEBRAIC METHOD -------------------------

function algebraic(sats, bestfit=false)
    M = hcat(sats...)'
    M[:, 4] *= c
    x = [sat[1:3] for sat in sats]
    t = c * [sat[4] for sat in sats]

    last = length(sats)
    A = hcat([x[i] - x[last] for i = 1:(last-1)]...)'
    B = [ t[i] - t[last] for i = 1:(last -1) ]
    R = [ product(M[i, :], M[i, :]) for i = 1:last ]
    C = 0.5 * [ R[i] - R[last] for i = 1:(last-1)]

    invA = pinv(A)
    a = invA * B
    b = invA * C
    bdiff = b - x[last]

    coeffs = Array{Float64}(3)
    coeffs[1] = sum(abs2, bdiff) - t[last]^2
    coeffs[2] = 2*(sum(a .* bdiff) + t[last])
    coeffs[3] = sum(abs2, a) - 1
    rts = filter(x -> isreal(x), roots(Poly(coeffs)))

    if(bestfit)
        if(length(rts) < 2)
            push!(rts, 0)
        end
        loc1 = vcat(a*rts[1] + b, rts[1])
        loc2 = vcat(a*rts[2] + b, rts[2])
        if(norm(vcat(f(loc1, M), loc1)) < norm(vcat(f(loc2, M), loc2)))
            location = loc1
        else
            location = loc2
        end
        location[4] /= c
        return location
    else
        guess = rts[indmin(abs.(rts))]
        return vcat(a * guess + b, guess / c)
    end
end

# ----------------- BANCROFT's METHOD -------------------

function bancroft(sats)
    A = hcat(sats...)'
    A[:, 4] *= c
    i = ones(length(sats))
    r = 0.5 * [product(A[i, :], A[i, :]) for i=1:size(A, 1)]
    B = pinv(A)

    u = B * i
    v = B * r
    E = product(u, u)
    F = product(u, v) - 1
    G = product(v, v)
    rts = filter(x -> isreal(x), roots(Poly([G, 2*F, E])))
    if(length(rts) < 2)
        push!(rts, 0)
    end

    j = [1, 1, 1, -1]
    loc1 = j .* (rts[1] * u + v)
    loc2 = j .* (rts[2] * u + v)
    if(norm(vcat(f(loc1, A), loc1)) < norm(vcat(f(loc2, A), loc2)))
        location = loc1
    else
        location = loc2
    end
    location[4] /= c
    location
end

# ------------------ HEURISTIC ---------------------

function heura(sats)
    j = [1, 1, 1, -1]
    last = length(sats)
    for i=1:last
        sats[i][4] *= c
    end

    A = hcat([j .* (sats[i] - sats[last]) for i = 1:(last-1)]...)'
    c1 = product(sats[last], sats[last])
    b = -0.5 * [ c1 - product(sats[i], sats[i]) for i = 1:(last-1) ]

    x = pinv(A) * b
    x[4] /= c

    for i=1:last
        sats[i][4] /= c
    end
    x
end
