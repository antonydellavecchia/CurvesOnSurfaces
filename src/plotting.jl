# defines the global precision
CC = AcbField(62)

# Define the @recipe function for HyperbolicPlane
@recipe function f(H::HyperbolicPlane)
    return Plots.partialcircle(0, 2π, 100, 1)
end

# Define the @recipe function for acb
@recipe function f(::Type{acb}, c::acb)
    return (convert(Float64, real(c)), convert(Float64, imag(c)))
end

@recipe function f(q::qqbar)
    @series begin
        [CC(q)]
    end
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
    θ1_f = convert(Float64, RR(θ1)) * π
    θ2_f = convert(Float64, RR(θ2)) * π
    
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

    #for A in D.deck_transformations
    #    @series begin
    #        gs = [on_geodesic(g, A) for g in D.geodesics]
    #        points = [[g.p1, g.p2] for g in gs]
    #    end
    #end
end
