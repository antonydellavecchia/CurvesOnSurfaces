# defines the global precision
const PLOT_CC = AcbField(62)

function get_arc(g::Geodesic)
    c = get_circle_center(g)
    p1_c = g.p1 - c
    p2_c = g.p2 - c
    r = abs(p1_c)

    theta1 = angle(PLOT_CC(p1_c // r))
    theta2 = angle(PLOT_CC(p2_c // r))

    return r, (theta1, theta2)
end

# Define the @recipe function for acb
@recipe function f(::Type{acb}, c::acb)
    return (convert(Float64, real(c)), convert(Float64, imag(c)))
end

@recipe function f(q::qqbar)
    return [PLOT_CC(q)]
end

@recipe function f(v::Vector{qqbar})
    for q in v
        @series begin
            PLOT_CC(q)
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
    
    radius = convert(Float64, RR(r))

    if abs(θ1_f - θ2_f) > π
        if θ1_f < 0
            θ1_f += 2 * π
        else
            θ1_f -= 2 * π
        end
    end
    points = Plots.partialcircle(θ1_f, θ2_f, 100, radius)
    
    c_real = convert(Float64, real(PLOT_CC(c)))
    c_imag = convert(Float64, imag(PLOT_CC(c)))
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
    legend = false
    @series begin
        color := :deepskyblue
        Plots.partialcircle(0, 2π, 100, 1)
    end
    D = fundamental_domain(H)
    @series begin
        color := :deepskyblue
        reduce(vcat, D.geodesics)
    end

    geodesics = D.geodesics
    n_iter = 1
    for i in 1:n_iter
        new_geodesics = Geodesic[]
        for (A, A_inv) in D.deck_transformations
            transformed_geodesics = reduce(vcat,
                               [[on_geodesic(g, A) for g in geodesics],
                                [on_geodesic(g, A_inv) for g in geodesics]])
            @series begin
                color := :deepskyblue
                transformed_geodesics
            end
            new_geodesics = reduce(vcat, [transformed_geodesics, new_geodesics])
        end
        geodesics = new_geodesics
    end
end
