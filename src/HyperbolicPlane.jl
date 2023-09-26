import Oscar.QQBar
    
function coordinates(z::qqbar)
    return [real(z), imag(z)]
end

function half_angle_coords(angle_coord::qqbar)
    r = real(angle_coord)
    cos = sqrt((1 + r) // 2)
    sin = sqrt((1 - r) // 2)
    
    return cos + sqrt(QQBar(-1)) * sin
end

function on_point(q::qqbar, M::MatElem{qqbar})
    z = (M[1, 1] * q + M[1, 2]) // (M[2, 1] * q + M[2, 2])
    
    if isone(abs(q))
        @req isone(abs(z)) "point was moved off boundary"
    end
    return z
end

const PointMapType{T} = Pair{T, T} where T <: FieldElem
function moebius(p_maps::Vector{<:PointMapType}; check::Bool=false)
    mat = Vector{qqbar}[]
    for p_map in p_maps
        z1 = p_map.first
        z2 = p_map.second
        eq = qqbar[z1, QQBar(1), - z1 * z2, -z2]
        push!(mat, eq)
    end
    M = matrix(QQBar, mat)

    dim, gen = kernel(M)
    @req dim == 1 "cannot find moebius transformation"
    m = matrix(QQBar, [gen[1] gen[2]; gen[3] gen[4]])
    
    if check
        for p_map in p_maps
            @req on_point(p_map.first, m) == p_map.second "Error mapping $(p_map.first) to $p_map.second"
        end
    end

    return m
end

################################################################################
# Geodesic

struct Geodesic
    # points should be ordered counter clockwise
    p1::qqbar
    p2::qqbar
end

struct FiniteGeodesic
    main::Geodesic 
    g1::Geodesic # start
    g2::Geodesic # end
end

function compare_angles(p1::qqbar, p2::qqbar; precision::Int = 32)
    CC = AcbField(precision)
    RR = ArbField(precision)
    theta1 = angle(CC(p1))
    theta2 = angle(CC(p2))
    pi = const_pi(RR)

    if abs(theta1 - theta2) > pi
        if isnegative(theta1)
            theta1 += 2 * pi
        else
            theta2 += 2 * pi
        end
    end

    if overlaps(theta1, theta2)
        return compare_angles(p1, p2; precision=precision + 1)
    end

    return theta1 < theta2
end

function geodesic(p1::qqbar, p2::qqbar)
    @req isone(abs(p1)) && isone(abs(p2)) && p1 != p2 "Geodesics require 2 unique points on the boundary"
    if compare_angles(p1, p2)
        return Geodesic(p1, p2)
    end
    return Geodesic(p2, p1)
end
    
function circle_center(g::Geodesic; check::Bool = true)
    v1 = g.p1
    v2 = g.p2
    v1_perp = perp(g.p1)
    v2_perp = perp(g.p2)
    M = matrix(QQBar, hcat(v1_perp, v2_perp))
    b = matrix(QQBar, 2, 1, [v2[1] - v1[1], v2[2] - v1[2]])
    t = solve(M, b)
    result = - v1_perp + v1

    if check
        @req result == t[2] .* v1_perp + v1 "Error computing center"
    end
    i = sqrt(QQBar(-1))

    return result[1] + i * result[2]
end

function on_geodesic(g::Geodesic, M::MatElem{qqbar}; copy::Bool = true)
    p1 = on_point(g.p1, M)
    p2 = on_point(g.p2, M)

    if copy
        return geodesic(p1, p2)
    end
    g.p1 = p1
    g.p2 = p2
    return g
end

function perp(q::qqbar)
    return q * sqrt(QQBar(-1))
end

function intersection(g1::Geodesic, g2::Geodesic)
    @time "get center" c1 = circle_center(g1)
    r1 = abs(g1.p1 - c1)
    c2 = circle_center(g2)
    r2 = abs(g2.p1 - c2)
    d =  c2 - c1
    @time "dist" abs_d = abs(d)
    @time "parameter" t = abs_d // 2 + (- r2^2 + r1^2) // (2 * abs_d)

    @time "scaling" d_d_abs = d // abs_d
    # projection of intersection points onto line between centers
    @time "translate from center" x = c1 + d_d_abs * t
    @time "perp" d_perp = perp(d_d_abs)
    @time "dist to inter" dist_to_intersection = sqrt((r1 + t) * (r1 - t))
    @time "direction and dist" direction_and_dist = dist_to_intersection * d_perp
    y = direction_and_dist  + x

    if abs(y) < 1
        return y
    else
        return - direction_and_dist + x
    end
end

################################################################################
# Fundamental Domain
struct FundamentalDomain
    identified_sides::Vector{Pair{Geodesic, Geodesic}}
    deck_transformations::Vector{MatElem{qqbar}}
    inv_transformations::Vector{MatElem{qqbar}}
end

deck_transformations(D::FundamentalDomain) = D.deck_transformations
inv_transformations(D::FundamentalDomain) = D.inv_transformations
identified_sides(D::FundamentalDomain) = D.identified_sides

function opposite_side(D::FundamentalDomain, g::Geodesic)
    for identified_side in identified_sides(D)
        if g == identified_side.first
            return identified_side.second
        elseif g == identified_side.second
            return identified_side.first
        end
    end
    @req true "geodesic $g doesn't correspond to a side of the fundamental domain"
end

function fundamental_domain(n::UInt; padding::QQFieldElem = QQ(0))
    @req n % 4 == 0 "Can only divide boundary into multiples of 4"
    @req abs(padding) <= 1 "padding needs to be between -1 and 1"

    generating_angle = exp_pi_i(QQBar(QQ(2 // n)))
    padding_angle = exp_pi_i(QQBar(padding // QQ(n)))
    boundary_points = qqbar[]
    current_angle = QQBar(1)
    geodesics = Geodesic[]

    for k in 1 : n
        if (k % n) == 0
            new_angle = conj(padding_angle)
        else
            new_angle = generating_angle * current_angle * conj(padding_angle)
        end

        push!(geodesics, geodesic(current_angle, new_angle))
        current_angle = new_angle * padding_angle
    end

    deck_transformations = MatElem{qqbar}[]
    identified_sides = Pair{Geodesic, Geodesic}[]
    
    n_2 = divexact(n, 2)
    for k in 1:n_2
        antipodal_k = k + n_2
        t = deck_transformation(geodesics[k], geodesics[antipodal_k])
        push!(deck_transformations, t)

        identified_side = geodesics[k] => geodesics[antipodal_k]
        push!(identified_sides, identified_side)
    end
    return FundamentalDomain(identified_sides, deck_transformations,
                             inv.(deck_transformations))
end

function deck_transformation(g1::Geodesic, g2::Geodesic)
    fixed_point = (g2.p1 + g2.p2) // abs(g2.p1 + g2.p2)
    p_maps = [
        g1.p1 => g2.p2,
        g1.p2 => g2.p1,
        fixed_point => fixed_point
    ]
    m = moebius(p_maps)
    return m
end

################################################################################
# Hyperbolic Plane
struct HyperbolicPlane
    D::FundamentalDomain
end

function hyperbolic_plane(n::UInt; padding::QQFieldElem = QQ(0))
    return HyperbolicPlane(fundamental_domain(n; padding=padding))
end

fundamental_domain(H::HyperbolicPlane) = H.D
deck_transformations(H::HyperbolicPlane) = deck_transformations(H.D)
