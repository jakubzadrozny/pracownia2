C= BigFloat(1.0)

function f(coefs)
    a = coefs[1]
    b = coefs[2]
    c = coefs[3]
    t = coefs[4]
    return (x, y, z, T) -> BigFloat((x-a)^2 + (y-b)^2 + (z-c)^2 - C*C*(T-t)^2)
end

function df(coefs, der)
    a = coefs[1]
    b = coefs[2]
    c = coefs[3]
    t = coefs[4]
    if der == "x"
        return (x, y, z, T) -> BigFloat(2*(x-a))
    elseif der == "y"
        return (x, y, z, T) -> BigFloat(2*(y-b))
    elseif der == "z"
        return (x, y, z, T) -> BigFloat(2*(z-c))
    elseif der == "T"
        return (x, y, z, T) -> BigFloat(-2*C*C*(T-t))
    end
end

functions = [] # i - fi
derivatives = [] # i j - dfi/dxj
vars = ["x", "y", "z", "T"]

function norm(V)
    s = BigFloat(0.0)
    for i in 1:4
        s += (functions[i](V[1], V[2], V[3], V[4]))^2
    end
    return sqrt(s)
end

function jacobian(X)
    x= X[1]
    y = X[2]
    z = X[3]
    T = X[4]
    J =[]
    for i in 1:4
        push!(J, [])
        for j in 1:4
            push!(J[i], derivatives[i][j](x, y, z, T))
        end
    end
    return J
end

### invert a matrix ############################################################
function A(a, b, c, d, M)
    BigFloat(M[1][a]*M[2][b]*M[3][c]*M[4][d])
end

function B(M, r, c)
    m = []
    for i in 1:4
        if i != r
            push!(m, [])
            for j in 1:4
                if j!=c
                    push!(m[length(m)], M[i][j])
                end
            end
        end
    end
    return BigFloat(m[1][1]*m[2][2]*m[3][3] + m[2][1]*m[3][2]*m[1][3]+
           m[3][1]*m[1][2]*m[2][3] - m[1][3]*m[2][2]*m[3][1] -
           m[2][3]*m[3][2]*m[1][1] - m[3][3]*m[1][2]*m[2][1])
end

function invert(M)
    MR = []
    det =BigFloat( A(1, 2, 3, 4, M) + A(1, 3, 4, 2, M) + A(1, 4, 2, 3, M) +
          A(2, 1, 4, 3, M) + A(2, 3, 1, 4, M) + A(2, 4, 3, 1, M) +
          A(3, 1, 2, 4, M) + A(3, 2, 4, 1, M) + A(3, 4, 1, 2, M) +
          A(4, 1, 3, 2, M) + A(4, 2, 1, 3, M) + A(4, 3, 2, 1, M) -
          A(1, 2, 4, 3, M) - A(1, 3, 2, 4, M) - A(1, 4, 3, 2, M) -
          A(2, 1, 3, 4, M) - A(2, 3, 4, 1, M) - A(2, 4, 1, 3, M) -
          A(3, 1, 4, 2, M) - A(3, 2, 1, 4, M) - A(3, 4, 2, 1, M) -
          A(4, 1, 2, 3, M) - A(4, 2, 3, 1, M) - A(4, 3, 1, 2, M))
    for i in 1:4
        push!(MR, [])
        for j in 1:4
            push!(MR[i], B(M, j, i)/det)
        end
    end
    return MR
end

################################################################################

function fn(X) ### value of F(X), where F(X) = (f1(X), f2(X), ...)
    x = X[1]
    y = X[2]
    z = X[3]
    T = X[4]
    F = []
    for i in 1:4
        push!(F, functions[i](x, y, z, T))
    end
    return F
end

function mult(J, F) ### matrix * vector
    V = []
    for i in 1:length(J)
        s = BigFloat(0.0)
        for j in 1:length(J)
            s += J[i][j]*F[j]
        end
        push!(V, s)
    end
    return V
end

function newton(e)
    X = [BigFloat(0.0), BigFloat(0.02), BigFloat(0.0), BigFloat(0.0)]
    while norm(X) > e
        J = invert(jacobian(X))
        F = fn(X)
        H = mult(J, F)
        X = X .- H
        println(F)
        println("------------------------------")
    end
    return X
end

function solveNewton(coefs)
    for i in 1:4
        # initializing fi and dfi/dxj
        push!(functions, f(coefs[i]))
        push!(derivatives, [])
        for j in 1:4
            push!(derivatives[i], df(coefs[i], vars[j]))
        end
    end
    return newton(BigFloat(0.01))
end

 println(solveNewton(([
        [BigFloat(1.0), BigFloat(0.0), BigFloat(0.0), BigFloat(1.0)],
        [BigFloat(0.0), BigFloat(1.0), BigFloat(0.0), BigFloat(1.0)],
        [BigFloat(0.0), BigFloat(0.0), BigFloat(1.0), BigFloat(1.0)],
        [BigFloat(1.0), BigFloat(1.0), BigFloat(1.0), BigFloat(sqrt(3))]
       ])))

#for i in 1:4
#    println(M[i][1], " ", M[i][2], " ", M[i][3], " ", M[i][4])
#end
