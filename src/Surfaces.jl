struct Surface
    fundamental_group::FPGroup
    universal_cover::HyperbolicPlane
end

function hyperbolic_surface(genus::Int; with_puncture::Bool = true)
    @req genus > 0 "genus needs to be positive"
    F = free_group(2 * genus)

    if !with_puncture
        @req genus > 1 "surface is not hyperbolic"
        first_part_word = reduce(*, gens(F)[1:end-1])
        surface_word = comm(first_part_word, gens(F)[end])
        fundamental_group = quo(F, surface_word)
        return Surface(fundamental_group, hyperbolic_plane(4 * UInt(genus); padding=QQ(-1//2)))
    end

    return Surface(F, hyperbolic_plane(4 * UInt(genus)))
end

fundamental_group(s::Surface) = s.fundamental_group
universal_cover(s::Surface) = s.universal_cover
fundamental_domain(s::Surface) = fundamental_domain(universal_cover(s))

################################################################################
# FreeHomotopyClass

struct FreeHomotopyClass
    word::FPGroupElem
    S::Surface
    g1::Geodesic
    g2::Geodesic
end

function free_homotopy_class(word::FPGroupElem, S::Surface)
    D = fundamental_domain(universal_cover(S))
    sides = identified_sides(D)
    #acts by last letter in the word first since we use a left action
    init_letter = letters(word)[end]

    if is_negative(init_letter)
        g1 = sides[abs(init_letter)][2]
    else
        g1 = sides[abs(init_letter)][1]
    end

    m_rep = map_word(word, deck_transformations(D); genimgs_inv=inv_transformations(D))

    g2 = on_geodesic(g1, m_rep)

    return FreeHomotopyClass(word, S, g1, g2)
end

surface(c::FreeHomotopyClass) = c.S
init_geodesic(c::FreeHomotopyClass) = c.g1
end_geodesic(c::FreeHomotopyClass) = c.g2

################################################################################
# SurfaceGeodesic

struct SurfaceGeodesic
    class::FreeHomotopyClass
    lift::FiniteGeodesic
    parts::Vector{FiniteGeodesic}
end

function surface_geodesic(c::FreeHomotopyClass, frac::QQFieldElem)
    S = surface(c)
    g1 = init_geodesic(c)
    g2 = end_geodesic(c)

    start_angle = (log_pi_i(g1.p2 // g1.p1)) * QQBar(frac) + log_pi_i(g1.p1)
    start_point = exp_pi_i(start_angle)

    m = moebius([
        start_point => start_point,
        g1.p1 => g2.p2,
        g1.p2 => g2.p1
    ])

    # find other fixed point of moebius transformation knowinga
    # one (start point) by factoring a monic polynomial
    QQBarx, x = QQBar["x"]
    p = x^2 + ((m[2, 2] - m[1, 1]) * x - m[1, 2]) * (QQBar(1) // m[2, 1])
    q = divexact(p, (x - start_point))
    other_fixed_point = -evaluate(q, QQBar(0))

    @req is_zero(evaluate(p, other_fixed_point)) "Couldn't find fixed point"

    main = geodesic(start_point, other_fixed_point)
    lift = FiniteGeodesic(main, g1, g2)
    
    return lift
end
