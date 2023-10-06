# defines the global plot precision
const PLOT_CC = AcbField(15)
const PLOT_RR = ArbField(15)

function arc(g::T) where T <: Union{FiniteGeodesic, Geodesic}
    if T <: FiniteGeodesic
        c = circle_center(geodesic(g))
        start_intersection = intersection(geodesic(g), g.g1)
        end_intersection = intersection(geodesic(g), g.g2)
        p1_c = start_intersection - c
        p2_c = end_intersection - c
    else
        c = circle_center(g)
        p1_c = g.p1 - c
        p2_c = g.p2 - c
    end
    r = abs(PLOT_CC(p1_c))
    theta1 = angle(PLOT_CC(p1_c) / r)
    theta2 = angle(PLOT_CC(p2_c) / r)

    theta1_f = convert(Float64, theta1)
    theta2_f = convert(Float64, theta2)
    
    radius = convert(Float64, PLOT_RR(r))

    if abs(theta1_f - theta2_f) > pi
        if theta1_f < 0
            theta1_f += 2 * pi
        else
            theta1_f -= 2 * pi
        end
    end
    points = Plots.partialcircle(theta1_f, theta2_f, 100, radius)

    c_real = convert(Float64, real(PLOT_CC(c)))
    c_imag = convert(Float64, imag(PLOT_CC(c)))
    return [p .+ (c_real, c_imag) for p in points]
end

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

@recipe function f(::Type{Vector{acb}}, v::Vector{acb})
    return [(convert(Float64, real(c)), convert(Float64, imag(c))) for c in v]
end

@recipe function f(g::Geodesic)
    arc(g)
end

@recipe function f(v::Vector{Geodesic})
    for g in v
        @series begin
            g
        end
    end
end

@recipe function f(fg::FiniteGeodesic)
    #start geodesic
    @series begin
        color := :green
        fg.g1
    end

    #end geodesic
    @series begin
        color := :red
        fg.g2
    end

    # main geodesic
    @series begin
        color := :blue
        arc(fg)
    end
end

@recipe function f(dg::DomainGeodesic)
    @series begin
        dg.fg
    end
end

@recipe function f(sg::SurfaceGeodesic)
    @series begin
        universal_cover(surface(sg))
    end

    for dg in domain_geodesics(sg)
        @series begin
            dg
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
    sides = identified_sides(D)
    geodesics = reduce(vcat, [[side.first, side.second] for side in sides])

    # plot main polygon
    @series begin
        color := :deepskyblue
        geodesics
    end

    # plot n_iter of translates1
    transformations = deck_transformations(D)
    invs = inv_transformations(D)
    n_iter = 1
    for i in 1:n_iter
        new_geodesics = Geodesic[]
        for (i, A) in enumerate(transformations)
            if !isdefined(invs, i)
                invs[i] = inv(A)
            end

            A_inv = invs[i]
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

