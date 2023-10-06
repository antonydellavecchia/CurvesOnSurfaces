import Base.length

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
    @req word_letters[end] != - word_letters[1] "word should be cyclically reduced"
    
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
lift(sg::SurfaceGeodesic) = sg.lift

function surface_geodesic(class::FreeHomotopyClass;
                          frac::QQFieldElem = QQ(1//2), check::Bool=false)
    S = surface(class)
    g1 = init_geodesic(class)
    g2 = end_geodesic(class)

    mid_point = (g1.p2 + g1.p1) // abs(g1.p1 + g1.p2)
    to_disc_m = moebius([
        QQBar(1) => g1.p1,
        QQBar(-1) => g1.p2,
        QQBar(0) => mid_point
    ])
    start_point = on_point(exp_pi_i(QQBar(frac)), to_disc_m)

    mid_point = (g2.p2 + g2.p1) // abs(g2.p1 + g2.p2)
    to_disc_m = moebius([
        QQBar(1) => g2.p1,
        QQBar(-1) => g2.p2,
        QQBar(0) => mid_point
    ])
    end_point = on_point(exp_pi_i(QQBar(frac)), to_disc_m)
    
    @time "finding moebius for lift" m = moebius([
        start_point => end_point,
        g1.p1 => g2.p2,
        g1.p2 => g2.p1
    ])

    # find other fixed point of moebius transformation knowing
    # one (start point) by factoring a monic polynomial
    b = (m[2, 2] - m[1, 1]) // m[2, 1]
    c = - m[1, 2] // m[2, 1]
    root1 = (-b + sqrt(b^2 - 4 * c)) // 2
    root2 = (-b - sqrt(b^2 - 4 * c)) // 2

    if check
        @req on_point(root1, m) == root1 "can't find root"
        @req on_point(root2, m) == root2 "can't find root"
    end

    main = geodesic(root1, root2)
    lift = FiniteGeodesic(main, g1, g2)

    w = word(class)
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

    return SurfaceGeodesic(class, lift, domain_fgs)
end

function self_intersections(sg::SurfaceGeodesic)
    intersections = Dict{Int, Vector{Tuple{DomainIntersection, Int}}}()
    dgs = domain_geodesics(sg)
    for (i, dg1) in enumerate(dgs)
        dg1_intersections = Tuple{DomainIntersection, Int, qqbar}[]
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
            push!(dg1_intersections, (domain_intersection(dg1, dg2), j, inter))
        end
        sort!(dg1_intersections; lt=compare_angles, by = x -> x[3], rev = compare_angles(e, s))
        intersections[i] = [dg_inter[1:2] for dg_inter in dg1_intersections]
    end
    return intersections
end

function length(sg::SurfaceGeodesic)
    l = lift(sg)
    p = intersection(geodesic(l), init_geodesic(l))
    q = intersection(geodesic(l), end_geodesic(l))
c
    # since log is monotone increasing, when comparing length we need only compare
    # the following expression
    return (abs(1 - conj(p) * q) + abs(q - p)) // (abs(1 - conj(p) * q) -  abs(q - p))
end

function min_conjugacy_class(c::FreeHomotopyClass)
    # calculates conjugacy class minimizing the unique geodesic
    S = surface(c)
    w = word(c)
    G =  parent(w)
    min_classes = SurfaceGeodesic[]
    conj_class = w
    
    for l in letters(w)
        c = free_homotopy_class(conj_class, S)
        sg = surface_geodesic(c)

        if isempty(min_classes) || length(sg) < length(min_classes[1])
            min_classes = [sg]
        elseif length(sg) == length(min_classes[1])
            push!(min_classes, sg)
        end
        
        conj_class = l > 0 ? conj(conj_class, G[l]) : conj(conj_class, inv(G[-l]))
    end
    return min_classes
end

################################################################################
# Polyhedral Complex

mutable struct HalfEdge
    next::Union{HalfEdge, Nothing}
    flip::Union{HalfEdge, Nothing}
    label:: Pair{Vector{Int}, Vector{Int}}
end

function set_next(h::HalfEdge, next::HalfEdge)
    h.next = next
end

function set_flip(h::HalfEdge, flip::HalfEdge)
    h.flip = flip
end

function half_edge(label::Pair{Vector{Int}, Vector{Int}})
    h = HalfEdge(nothing, nothing, label)
    f = HalfEdge(nothing, h, reverse(label))
    set_flip(h, f)
    return h
end

function polyhedral_subdivision(sg::SurfaceGeodesic)
    # self intersections are ordered by lines following along the curve
    # keys start at 1 and increase, each value is an ordered list of intersections
    # they are tuples where the second entry is the key to the line they also appear on
    self_inters = self_intersections

    for i in 1:length(self_inters)
        for self_inter in self_inters[i]
            # self inter is the intersection of the line i with some other line
            inter_label = label(self_inter[1])
            # j is the index for the self intersection on the other line tha
            j = self_inter[2]
            k = findfirst(x -> x[2] == i, self_inters[j])
        end
    end
end
