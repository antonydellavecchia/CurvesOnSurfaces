
struct HyperbolicPlane{T<:Field, S <:Field}
    base_field::T
    field::S
end

struct Geodesic{S <: FieldElem}
    # points should be ordered counter clockwise
    p1::S
    p2::S
end

struct FundamentalDomain
    H::HyperbolicPlane
    geodesics::Vector{Geodesic}
    deck_transformations::Vector{Matrix{<:FieldElem}}
end

function hyperbolic_plane()
    Qx, x = QQ["x"]
    F, _ = number_field(x^2 + 1)
    return HyperbolicPlane(F, F)
end

base_field(H::HyperbolicPlane) = H.base_field
field(H::HyperbolicPlane) = H.field

function set_field(H::HyperbolicPlane, F::Field)
    return HyperbolicPlane(base_field(H), F)
end

function fundamental_domain(n::UInt, H::Union{HyperbolicPlane, Nothing} = nothing)
    if isnothing(H)
        H = hyperbolic_plane()
    end
        
    @req n % 4 == 0 "Can only divide boundary into multiples of 4"
    n_halves = n // 4 - 1
    i = gen(base_field(H))
    angle_coord = i

    for _ in 1:n_halves
        H, angle_coord = half_angle_coords(H, angle_coord)
    end

    boundary_points = [angle_coord^k for k in 0:n - 1]

    geodesics = Geodesic{elem_type(field(H))}[]

    for k in 1:n
        push!(geodesics, Geodesic(boundary_points[k], boundary_points[(k % n) + 1]))
    end

    deck_transformations = Matrix{<:FieldElem}[]
    n_2 = divexact(n, 2)
    for k in 1:n_2
        antipodal_k = k + n_2
        t = deck_transformation(H, geodesics[k], geodesics[antipodal_k])
        push!(deck_transformations, t)
    end
    return FundamentalDomain(H, geodesics, deck_transformations)
end

function deck_transformation(H::HyperbolicPlane, g1::Geodesic, g2::Geodesic)
    fixed_point = (g2.p1 + g2.p2) // norm(g2.p1 + g2.p2)
    p_maps = [
        g1.p1 => g2.p2,
        g1.p2 => g2.p1,
        fixed_point => fixed_point
    ]
    return get_moebius(H, p_maps)
end

function half_angle_coords(H::HyperbolicPlane, angle_coord::FieldElem)
    base = base_field(H)
    i = gen(base_field(H))
    parent_field = parent(angle_coord)
    Fx, x = parent_field["x"]
    norm = absolute_norm(angle_coord)
    @req isone(norm) "angle coordinate not on unit circle"
    conjugate_coord = parent_field(norm) // angle_coord 
    real = (angle_coord + conjugate_coord) // 2
    q = x^2 - (1 + real) // 2 # cos
    p = x^2 - (1 - real) // 2 # sin

    if p == q
        E, root_half = number_field(q)
        H = set_field(H, E)
        return H, root_half + root_half * i
    elseif iszero(real^2 - 1//2)
        E, a = number_field(q)
        b = a // (1 + 2 * real)
        H = set_field(H, E)
        return H, a + b * i
    end

    E, (c, s) = number_field([q, p], check=false)
    H = set_field(H, E)
    return H, c + i * s
end

const PointMapType{T} = Pair{T, T} where T <: FieldElem
function get_moebius(H::HyperbolicPlane, p_maps::Vector{<:PointMapType})
    F = field(H)
    mat = Vector{elem_type(F)}[]
    for p_map in p_maps
        z1 = p_map.first
        z2 = p_map.second
        eq = elem_type(F)[z1, F(1), - z1 * z2, -z2]
        push!(mat, eq)
    end
    M = matrix(F, mat)
    dim, gen = kernel(M)
    @req dim == 1 "cannot find moebius transformation"
    
    A = [gen[1] gen[2]; gen[3] gen[4]];
    d = (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])
    return [gen[1] // d  gen[2] // d; gen[3] gen[4]];
end
