# defines the global precision
CC = AcbField(62)

function get_arc(g::Geodesic)
    c = get_circle_center(g)
    p1_c = g.p1 - c
    p2_c = g.p2 - c
    r = abs(p1_c)
    theta1 = imag(log(CC(p1_c // r)))
    theta2 = imag(log(CC(p2_c // r)))
    return r, (theta1, theta2)
end


# Define the @recipe function for HyperbolicPlane
@recipe function f(H::HyperbolicPlane)
    return Plots.partialcircle(0, 2π, 100, 1)
end

# Define the @recipe function for acb
@recipe function f(::Type{acb}, c::acb)
    return (convert(Float64, real(c)), convert(Float64, imag(c)))
end

@recipe function f(q::qqbar)
    return [CC(q)]
end

@recipe function f(v::Vector{qqbar})
    for q in v
        @series begin
            CC(q)
        end
    end
end

# Define the @recipe function for Vector{acb}
@recipe function f(::Type{Vector{acb}}, v::Vector{acb})
    return [(convert(Float64, real(c)), convert(Float64, imag(c))) for c in v]
end

@recipe function f(g::Geodesic{qqbar})
    c = get_circle_center(g)
    r, (θ1, θ2) = get_arc(g)

    RR = ArbField(64)
    θ1_f = convert(Float64, θ1)
    θ2_f = convert(Float64, θ2)
    
    if θ1_f < 0
        θ1_f += 2 * π
    end
    
    radius = convert(Float64, RR(r))
    points = Plots.partialcircle(θ1_f, θ2_f, 100, radius)
    
    c_real = convert(Float64, real(CC(c)))
    c_imag = convert(Float64, imag(CC(c)))
    [p .+ (c_real, c_imag) for p in points]
end

## Define the @recipe function for Vector{Geodesic}
@recipe function f(v::Vector{Geodesic{qqbar}})
    for g in v
        @series begin
            g
        end
    end
end

@recipe function f(H::HyperbolicPlane)
    @series begin
        Plots.partialcircle(0, 2π, 100, 1)
    end
    D = fundamental_domain(H)
    @series begin
        reduce(vcat, D.geodesics)
    end

    for (A, A_inv) in D.deck_transformations
        @series begin
            [on_geodesic(g, A) for g in D.geodesics]
        end
        @series begin
            [on_geodesic(g, A_inv) for g in D.geodesics]
        end
    end
end
