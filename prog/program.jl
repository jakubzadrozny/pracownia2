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

## METODA NEWTON -- opis
# Z uogólnionego twierdzenia Taylora mamy
# \begin{equation}
# f_i(x + h_x, y + h_y, z + h_z, t + h_t) \approx f_i(x, y, z, t) + h_x \frac{\partial f_i}{\partial x}(x, y, z, t)
# + h_y \frac{\partial f_i}{\partial y}(x, y, z, t) + h_z \frac{\partial f_i}{\partial z}(x, y, z, t)
# + h_t \frac{\partial f_i}{\partial t}(x, y, z, t)
# \end{equation}
#
# Ponadto
# \begin{equation}
# \frac{\partial f_i}{\partial x} = 2(x-x_i) \\
# \frac{\partial f_i}{\partial y} = 2(y-y_i) \\
# \frac{\partial f_i}{\partial z} = 2(z-z_i) \\
# \frac{\partial f_i}{\partial t} = -2c^2(t-t_i)
# \end{equation}
#
# Stąd
# \begin{equation}
# f_i(x + h_x, y + h_y, z + h_z, t + h_t) \approx f_i(x, y, z, t) + 2h_x(x-x_i) + 2h_y(y-y_i) + 2h_z(z-z_i) -2c^2h_t(t-t_i)
# \end{equation}
# Zastosujemy metodę Newtona. Chcemy znaleźć $x, y, z, t$, dla których $f_i(x, y, z, t) = 0$ dla $1 \leq i \leq 4$.
# Niech $x_n \in \mathbb{R}^4$ będzie n-tym przybliżeniem metody Newtona. Szukamy $h \in \mathbb{R}^4$, takiego że
# dla $x_{n+1} = x_n + h$ zajdzie $f_i(x_{n+1}) = 0$ dla $1 \leq i \leq 4$.
#
# Niech
# $$
# x_n =
#     \begin{pmatrix}
#     x \\
#     y \\
#     z \\
#     t
#     \end{pmatrix}
# $$
#
# Rozwiążemy układ równań
# \begin{equation}
# 0 = f_1(x, y, z, t) + 2h_x(x-x_1) + 2h_y(y-y_1) + 2h_z(z-z_1) -2c^2h_t(t-t_1) \\
# 0 = f_2(x, y, z, t) + 2h_x(x-x_2) + 2h_y(y-y_2) + 2h_z(z-z_2) -2c^2h_t(t-t_2) \\
# 0 = f_3(x, y, z, t) + 2h_x(x-x_3) + 2h_y(y-y_3) + 2h_z(z-z_3) -2c^2h_t(t-t_3) \\
# 0 = f_4(x, y, z, t) + 2h_x(x-x_4) + 2h_y(y-y_4) + 2h_z(z-z_4) -2c^2h_t(t-t_4)
# \end{equation}
# Inaczej
# \begin{equation}
# \begin{pmatrix}
# x-x_1  & y-y_1 & z-z_1 & -c^2(t-t_1) \\
# x-x_2  & y-y_2 & z-z_2 & -c^2(t-t_2) \\
# x-x_3  & y-y_3 & z-z_3 & -c^2(t-t_3) \\
# x-x_4  & y-y_4 & z-z_4 & -c^2(t-t_4) \\
# \end{pmatrix}
# \begin{pmatrix}
# h_x \\ h_y \\ h_z \\ h_t
# \end{pmatrix}
# =
# -\frac{1}{2}
# \begin{pmatrix}
# f_1(x, y, z, t) \\
# f_2(x, y, z, t) \\
# f_3(x, y, z, t) \\
# f_4(x, y, z, t)
# \end{pmatrix}
# \end{equation}
#
# Wtedy
# $$
# x_{n+1} = x_n + h = x_n + \begin{pmatrix} h_x \\ h_y \\ h_z \\ h_t \end{pmatrix}
# $$
#
# ZAMIANA JEDNOSTEK -- przyjemniejsze obliczenia
# function toalgebraic(sat; onedim::Bool=false)
#     newsat = deepcopy(sat)
#     if(onedim)
#         newsat[4] *= 10^3
#         newsat[1:3] /= 10^2
#     else
#         for i = 1:length(sat)
#             newsat[i][4] *= 10^3
#             newsat[i][1:3] /= 10^2
#         end
#     end
#     newsat
# end
#
# PODSTAWIENIE POD INNA WSPOLRZEDNA
# function algebraic2(sats)
#     x = [sat[1] for sat in sats]
#     leftx = [sat[2:4] for sat in sats]
#     cvect = [1, 1, -c^2]
#     A = hcat([-cvect .* (leftx[i]-leftx[4]) for i = 1:3]...)'
#     B = x[1:3] .- x[4]
#     sum4 = sum(cvect .* (leftx[4].^2))
#     C = 0.5 * [x[4]^2-x[i]^2 - sum(cvect .* (leftx[i].^2)) + sum4 for i = 1:3]
#
#     a = A \ B
#     b = A \ C
#     bdiff = b - leftx[4]
#
#     coeffs = Array{Float64}(3)
#     coeffs[1] = sum(cvect .* (bdiff .^ 2)) + x[4]^2
#     coeffs[2] = 2*(-x[4] + sum(cvect .* a .* bdiff))
#     coeffs[3] = 1 + sum(cvect .* (a.^2))
#     guesses = filter(x -> isreal(x), roots(Poly(coeffs)))
#
#     location = vcat(guesses[1], [a[i]*guesses[1] + b[i] for i = 1:3])
#     if(length(guesses) > 1)
#         location2 = vcat(guesses[2], [a[i]*guesses[2] + b[i] for i = 1:3])
#         if(sum(abs2, location2) < sum(abs2, location))
#             location = location2
#         end
#     end
#     location
# #     vcat(location[1:3]*10^2, location[4]/10^3)
# end
#
# Zastosujemy inne podejście (nieiteracyjne). Rozważany układ równań wygląda następująco
#
# $$
# x^2 + y^2 + z^2 -2xx_1 -2yy_1 -2zz_1 + x_1^2 + y_1^2 + z_1^2 = c^2(t^2 -2tt_1 + t_1^2) \\
# x^2 + y^2 + z^2 -2xx_2 -2yy_2 -2zz_2 + x_2^2 + y_2^2 + z_2^2 = c^2(t^2 -2tt_2 + t_2^2) \\
# x^2 + y^2 + z^2 -2xx_3 -2yy_3 -2zz_3 + x_3^2 + y_3^2 + z_3^2 = c^2(t^2 -2tt_3 + t_3^2) \\
# x^2 + y^2 + z^2 -2xx_4 -2yy_4 -2zz_4 + x_4^2 + y_4^2 + z_4^2 = c^2(t^2 -2tt_4 + t_4^2)
# $$
#
# Po odjęciu czwartego równania stronami od pierwszych trzech otrzymamy
#
# $$
# -2x(x_1-x_4) -2y(y_1-y_4) -2z(z_1-z_4) + x_1^2 + y_1^2 + z_1^2 - (x_4^2 + y_4^2 + z_4^2) = c^2(-2t(t_1-t_4)+t_1^2-t_4^2) \\
# -2x(x_2-x_4) -2y(y_2-y_4) -2z(z_2-z_4) + x_2^2 + y_2^2 + z_2^2 - (x_4^2 + y_4^2 + z_4^2) = c^2(-2t(t_2-t_4)+t_2^2-t_4^2) \\
# -2x(x_3-x_4) -2y(y_3-y_4) -2z(z_3-z_4) + x_3^2 + y_3^2 + z_3^2 - (x_4^2 + y_4^2 + z_4^2) = c^2(-2t(t_3-t_4)+t_3^2-t_4^2)
# $$
#
# Otrzymaliśmy więc układ 3 równań liniowych na 4 zmiennych. Układ taki nie posiada jednozanczego rozwiązania, ale - o ile nie jest sprzeczny - posiada rozwiązania parametryczne. W naszym przypadku (dla realnych danych) powinno to być rozwiązanie zależne od jednego parametru. Oznacza to, że pozostałe zmienne można wyrazić jako kombinacje liniowe tego parametru tak, aby dla dowolnej wartości parametru układ równań był spełniony. Ponieważ dla faktycznych danych żadna zmienna nie powinna być z góry ustalona (na trzech równaniach), to możemy założyć, że parametrem jest $t$. Wtedy dla pewnych rzeczywistych $a_x, a_y, a_z, b_x, b_y, b_z$ mamy
#
# $$
# x = a_xt + b_x \\
# y = a_yt + b_y \\
# z = a_zt + b_z
# $$
# Wówczas dla dowolnego $t$ zachodzi
#
# $$
# -2(a_xt+b_x)(x_1-x_4)-2(a_yt+b_y)(y_1-y_4)-2(a_xt+b_z)(z_1-z_4)+x_1^2+y_1^2+z_1^2-c_0=c^2(-2t(t_1-t_4)+t_1^2-t_4^2) \\
# -2(a_xt+b_x)(x_2-x_4)-2(a_yt+b_y)(y_2-y_4)-2(a_xt+b_z)(z_2-z_4)+x_2^2+y_2^2+z_2^2-c_0=c^2(-2t(t_2-t_4)+t_2^2-t_4^2) \\
# -2(a_xt+b_x)(x_3-x_4)-2(a_yt+b_y)(y_3-y_4)-2(a_xt+b_z)(z_3-z_4)+x_3^2+y_3^2+z_3^2-c_0=c^2(-2t(t_3-t_4)+t_3^2-t_4^2)
# $$
#
# gdzie $c_0 = x_4^2 + y_4^2 + z_4^2$.
#
# Otrzymaliśmy zatem równości trzech funckji liniowych dla każdego argumentu.
#
# Stąd
# $$
# a_x(x_1-x_4) + a_y(y_1-y_4) + a_z(z_1-z_4) = c^2(t_1-t_4) \\
# a_x(x_2-x_4) + a_y(y_2-y_4) + a_z(z_2-z_4) = c^2(t_2-t_4) \\
# a_x(x_3-x_4) + a_y(y_3-y_4) + a_z(z_3-z_4) = c^2(t_3-t_4) \\
# -2(b_x(x_1-x_4)+b_y(y_1-y_4)+b_z(z_1-z_4)) + x_1^2+y_1^2+z_1^2 - c_0 = c^2(t_1^2 - t_4^2) \\
# -2(b_x(x_2-x_4)+b_y(y_2-y_4)+b_z(z_2-z_4)) + x_2^2+y_2^2+z_2^2 - c_0 = c^2(t_2^2 - t_4^2) \\
# -2(b_x(x_3-x_4)+b_y(y_3-y_4)+b_z(z_3-z_4)) + x_3^2+y_3^2+z_3^2 - c_0 = c^2(t_3^2 - t_4^2)
# $$
#
# Inaczej
# $$
# \begin{pmatrix}
# x_1-x_4 & y_1-y_4 & z_1-z_4 \\
# x_2-x_4 & y_2-y_4 & z_2-z_4 \\
# x_3-x_4 & y_3-y_4 & z_3-z_4
# \end{pmatrix}
# \begin{pmatrix}
# a_x \\ a_y \\ a_z
# \end{pmatrix}
# =
# c^2
# \begin{pmatrix}
# t_1-t_4 \\ t_2-t_4 \\ t_3-t_4
# \end{pmatrix}
# \\
# \begin{pmatrix}
# x_1-x_4 & y_1-y_4 & z_1-z_4 \\
# x_2-x_4 & y_2-y_4 & z_2-z_4 \\
# x_3-x_4 & y_3-y_4 & z_3-z_4
# \end{pmatrix}
# \begin{pmatrix}
# b_x \\ b_y \\ b_z
# \end{pmatrix}
# =
# -\frac{1}{2}
# \begin{pmatrix}
# c^2(t_1^2-t_4^2) - x_1^2+y_1^2+z_1^2 + c_0 \\
# c^2(t_2^2-t_4^2) - x_2^2+y_2^2+z_1^2 + c_0 \\
# c^2(t_3^2-t_4^2) - x_3^2+y_3^2+z_1^2 + c_0
# \end{pmatrix}
# $$
#
# Z ostatnich dwóch równości możemy łatwo wyznaczyć współczynniki $a_x, a_y, a_z, b_x, b_y, b_z$.
#
# Aby otrzymać konkretne rozwiązanie, możemy podstawić otrzymane zależności do czwartego równania. Wtedy otrzymamy równanie kwadratowe jednej zmiennej następującej postaci
#
# $$
#   t^2(a_x^2+a_y^2+a_x^2-c^2) + 2t(a_x(b_x-x_4)+a_y(b_y-y_4)+a_z(b_z-z_4)-c^2t_4^2) + (b_x-x_4)^2+(b_y-y_4)^2+(b_z-z_4)-c^2t_4^2 = 0
# $$
#
# Po rozwiązaniu równania otrzymamy dwa kandydaty (?) na $t$, które wyznaczą dwa możliwe rozwiązania $x, y, z, t$. W poniższym kodzie przyjęto, że szukanym rozwiązaniem jest to o mniejszym bezwzględnym błędzie zegara.
