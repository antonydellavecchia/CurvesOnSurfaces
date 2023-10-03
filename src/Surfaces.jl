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
    word_letters = letters(word)

    # we don't allow conjugations words
    @req word_letters[end] != word_letters[1] "word should be cyclically reduced"
    
    D = fundamental_domain(universal_cover(S))
    sides = identified_sides(D)
    #acts by last letter in the word first since we use a left action
    init_letter = letters(word)[end]

    # start the translate from opposite side so lift goes across fundamental domain
    # but we start the lift from g1, this is so the loop is closed
    if is_positive(init_letter)
        g1, start_side = sides[abs(init_letter)]
    else
        start_side, g1 = sides[abs(init_letter)]
    end

    m_rep = map_word(word, deck_transformations(D); genimgs_inv=inv_transformations(D))
    
    g2 = on_geodesic(start_side, m_rep)
    
    return FreeHomotopyClass(word, S, g1, g2)
end

surface(c::FreeHomotopyClass) = c.S
word(c::FreeHomotopyClass) = c.word
init_geodesic(c::FreeHomotopyClass) = c.g1
end_geodesic(c::FreeHomotopyClass) = c.g2

################################################################################
# SurfaceGeodesic

struct SurfaceGeodesic
    class::FreeHomotopyClass
    lift::FiniteGeodesic
    domain_geodesics::Vector{DomainGeodesic}
end

free_homotopy_class(sg::SurfaceGeodesic) = sg.class
domain_geodesics(sg::SurfaceGeodesic) = sg.domain_geodesics
surface(sg::SurfaceGeodesic) = surface(free_homotopy_class(sg))

function surface_geodesic(c::FreeHomotopyClass, frac::QQFieldElem; check::Bool=false)
    S = surface(c)
    g1 = init_geodesic(c)
    g2 = end_geodesic(c)

    start_angle = (log_pi_i(g1.p2 // g1.p1)) * QQBar(frac) + log_pi_i(g1.p1)
    start_point = exp_pi_i(start_angle)

    @time "finding moebius for lift" m = moebius([
        start_point => start_point,
        g1.p1 => g2.p2,
        g1.p2 => g2.p1
    ])

    # find other fixed point of moebius transformation knowing
    # one (start point) by factoring a monic polynomial
    QQBarx, x = QQBar["x"]
    p = x^2 + ((m[2, 2] - m[1, 1]) * x - m[1, 2]) * (QQBar(1) // m[2, 1])
    @time "dividing" q = divexact(p, (x - start_point))
    other_fixed_point = -evaluate(q, QQBar(0))

    if check
        @req is_zero(evaluate(p, other_fixed_point)) "Couldn't find fixed point"
    end

    main = geodesic(start_point, other_fixed_point)
    lift = FiniteGeodesic(main, g1, g2)

    w = word(c)
    D = fundamental_domain(universal_cover(S))
    transformations = deck_transformations(D)
    invs = inv_transformations(D)
    
    sides = identified_sides(D)
    start_g = g1
    translated_main = main
    domain_fgs = DomainGeodesic[]
    current_translate = Int[]
    
    for l in letters(w)
        # end_side is a pair of geodesics that are identified
        end_side = sides[abs(l)]
        end_g = ispositive(l) ? end_side.second : end_side.first

        fg = FiniteGeodesic(translated_main, start_g, end_g)
        dg = domain_geodesic(fg, copy(current_translate))
        push!(domain_fgs, dg)
        start_g = opposite_side(D, end_g)

        m = ispositive(l) ? invs[l] : transformations[-l]
        translated_main = on_geodesic(translated_main, m)
        current_translate = pushfirst!(current_translate, -l)
    end
    
    end_g = opposite_side(D, g1)
    fg = FiniteGeodesic(translated_main, start_g, end_g)
    dg = DomainGeodesic(fg, copy(current_translate))
    push!(domain_fgs, dg)

    return SurfaceGeodesic(c, lift, domain_fgs)
end

function self_intersections(sg::SurfaceGeodesic)
    intersections = DomainIntersection[]
    dgs = domain_geodesics(sg)
    for (i, dg1) in enumerate(dgs)
        dg1_intersections = Tuple{DomainIntersection, qqbar}[]
        dg1_start = init_geodesic(dg1)
        dg1_end = end_geodesic(dg1)
        dg1_geodesic = geodesic(dg1)

        start_point = intersection(dg1_geodesic, dg1_start)
        end_point = intersection(dg1_geodesic, dg1_end)
        center = circle_center(dg1_geodesic)
        s = start_point - center
        e = end_point - center
        
        for (j, dg2) in enumerate(dgs)
            if i == j
                continue
            end
            
            inter = intersection(dg1_geodesic, geodesic(dg2))

            if isnothing(inter)
                continue
            end

            inter -= center
            
            if compare_angles(s, e) && (compare_angles(inter, s) || compare_angles(e, inter))
                continue
            end
            if compare_angles(e, s) && (compare_angles(inter, e) || compare_angles(s, inter))
                continue
            end

            # intersection lies in fundamental domain
            push!(dg1_intersections, (domain_intersection(dg1, dg2), inter))
        end
        sort!(dg1_intersections; lt=compare_angles, by = x -> x[2], rev = compare_angles(e, s))
        intersections = vcat(intersections, [dg_inter[1] for dg_inter in dg1_intersections])
    end
    return intersections
end
