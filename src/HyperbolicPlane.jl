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

const PointMapType{T} = Pair{T, T} where T <: FieldElem
function get_moebius(p_maps::Vector{<:PointMapType})
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

    
    A = matrix(QQBar, [gen[1] gen[2]; gen[3] gen[4]])
    sqrt_det = sqrt(det(A))
    return matrix(QQBar, (1 // sqrt_det) .* A)
end

################################################################################
# Geodesic

mutable struct Geodesic{S <: FieldElem}
    # points should be ordered counter clockwise
    p1::S
    p2::S
end

function compare_angles(p1::qqbar, p2::qqbar; precision::Int = 32)
    CC = AcbField(precision)
    theta1 = angle(CC(p1))
    theta2 = angle(CC(p2))

    if overlaps(theta1, theta2)
        return compare_angles(p1, p2; precision=precision + 1)
    end

    return theta1 < theta2
end

function geodesic(p1::qqbar, p2::qqbar)
    @req isone(abs(p1)) && isone(abs(p2)) && p1 != p2 "Geodesics require 2 unique points on the boundary"
    if compare_angles(p1, p2)
        Geodesic(p1, p2)
    end
    return Geodesic(p2, p1)
end

function get_circle_center(g::Geodesic)
    v1 = coordinates(g.p1)
    v2 = coordinates(g.p2)
    
    v1_perp = [-v1[2], v1[1]]
    v2_perp = [-v2[2], v2[1]]
    M = matrix(QQBar, hcat(v1_perp, v2_perp))
    b = [v2[1] - v1[1], v2[2] - v1[2]]
    t = inv(M) * b
    result = t[1] .* v1_perp + v1

    @req result == t[2] .* v1_perp + v1 "Error computing center"
    i = sqrt(QQBar(-1))

    return result[1] + i * result[2]
end

function on_point(q::qqbar, M::MatElem{qqbar})
    z = (M[1, 1] * q + M[1, 2]) // (M[2, 1] * q + M[2, 2])
    
    if isone(abs(q))
        @req isone(abs(z)) "point was moved off boundary"
    end
    return z
end

function on_geodesic(g::Geodesic{qqbar}, M::MatElem{qqbar}; copy::Bool = true)
    p1 = on_point(g.p1, M)
    p2 = on_point(g.p2, M)

    if copy
        return geodesic(p1, p2)
    end
    g.p1 = p1
    g.p2 = p2
    return g
end
################################################################################
# Fundamental Domain
struct FundamentalDomain
    geodesics::Vector{Geodesic}
    deck_transformations::Vector{Tuple{MatElem{qqbar}, MatElem{qqbar}}}
end

function fundamental_domain(n::UInt)
    @req n % 4 == 0 "Can only divide boundary into multiples of 4"
    
    angle_coord = exp_pi_i(QQBar(QQ(2 // n)))

    boundary_points = [angle_coord^k for k in 0:n - 1]

    geodesics = Geodesic{qqbar}[]

    for k in 1:n
        push!(geodesics, geodesic(boundary_points[k], boundary_points[(k % n) + 1]))
    end

    deck_transformations = Tuple{MatElem{qqbar}, MatElem{qqbar}}[]
    n_2 = divexact(n, 2)
    for k in 1:n_2
        antipodal_k = k + n_2
        t = deck_transformation(geodesics[k], geodesics[antipodal_k])
        push!(deck_transformations, (t, inv(t)))
    end
    return FundamentalDomain(geodesics, deck_transformations)
end

function deck_transformation(g1::Geodesic, g2::Geodesic)
    fixed_point = (g2.p1 + g2.p2) // abs(g2.p1 + g2.p2)
    p_maps = [
        g1.p1 => g2.p2,
        g1.p2 => g2.p1,
        fixed_point => fixed_point
    ]
    m = get_moebius(p_maps)
    @req on_point(fixed_point, m) == fixed_point "Point is not being fixed"
    @req on_point(g1.p1, m) == g2.p2 "Error with mapping first point"
    @req on_point(g1.p2, m) == g2.p1 "Error with mapping second point"

    return m
end

################################################################################
# Hyperbolic Plane
struct HyperbolicPlane
    D::FundamentalDomain
    geodesics::Vector{Geodesic}
end

function hyperbolic_plane(n::UInt)
    return HyperbolicPlane(fundamental_domain(n), Geodesic[])
end

fundamental_domain(H::HyperbolicPlane) = H.D
