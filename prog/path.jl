using PyPlot

c = 299792458.0

function fsingle(x, sat)
    a = (x - sat) .^ 2
    a[4] *= -c^2
    sum(a)
end

function fvect(x, sats)
    [fsingle(x, sats[i]) for i = 1:length(sats)]
end

function getcoeffs(x, sats)
    A = [x - sats[i] for i = 1:length(sats)]
    A = hcat(A...)'
    A[:, 4] *= -c^2
    A
end

function newton(sats; maxiter::Int64=3)
    x = [0, 0, 6370, 0]
    for i = 1:maxiter
        A = getcoeffs(x, sats)
        B = fvect(x, sats)
        h = -0.5 * (A \ B)
        x = x + h
    end
    x
end;

function randSgn()
    x = rand(0:1)
    if x==0
        return -1
    else
        return 1
    end
end

function path(maxiter, lon0, lat0)
    a = pi/4
    d = 0.00001
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

function countTime(mistake, Pos, sat)
    dist = sqrt(sum( [ (Pos[i]-sat[i])^2 for i in 1:3]))
    t = dist / c
    return [t + mistake]
end


sat1 = [15600,7540,20140]
sat2 = [18760,2750,18610]
sat3 = [17610,14630,13480]
sat4 = [19170,610,18390]

function GPS_newton(coords, mistake)
    res = [[],[]]
    for i in 1:length(coords[1])
        sat = Array{Array{Float64}}(4)
        sat[1] = vcat(sat1, countTime(mistake, LLAtoXYZ([coords[1][i], coords[2][i]]), sat1) )
        sat[2] = vcat(sat2, countTime(mistake, LLAtoXYZ([coords[1][i], coords[2][i]]), sat2) )
        sat[3] = vcat(sat3, countTime(mistake, LLAtoXYZ([coords[1][i], coords[2][i]]), sat3) )
        sat[4] = vcat(sat4, countTime(mistake, LLAtoXYZ([coords[1][i], coords[2][i]]), sat4) )
        X = newton(sat)
        LL = XYZtoLLA(X[1:3])
        push!(res[1], LL[1])
        push!(res[2], LL[2])
    end
    return res
end

function MAXdist(coords, gps)
    res = 0.0
    for i in 1:length(coords[1])
        X1 = LLAtoXYZ([coords[1][i], coords[2][i]])
        X2 = LLAtoXYZ([gps[1][i], gps[2][i]])
        d =sqrt(sum([(X1[i]-X2[i])^2   for i in 1:3] ))
        res = max(res, d)
    end
    return res*1000
end

function MINdist(coords, gps)
    res = 100000.0
    for i in 1:length(coords[1])
        X1 = LLAtoXYZ([coords[1][i], coords[2][i]])
        X2 = LLAtoXYZ([gps[1][i], gps[2][i]])
        d =sqrt(sum([(X1[i]-X2[i])^2   for i in 1:3] ))
        res = min(res, d)
    end
    return res*1000
end

coords = path(3000, -70, -32.5)
x = coords[1]
y = coords[2]

GPScoords = GPS_newton(coords, 0.00001)
xgps = GPScoords[1]
ygps = GPScoords[2]

println("BIGGEST MISTAKE : ", MAXdist(coords, GPScoords), " meters")
println("SMALLEST MISTAKE : ", MINdist(coords, GPScoords), " meters")

fig, ax = subplots()

ax[:plot](x, y, "-", color="red", linewidth=5, alpha=0.6)
ax[:plot](xgps, ygps, "-", color="black", linewidth=1, alpha=0.6)
show()
