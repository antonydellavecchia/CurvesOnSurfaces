using Oscar, Test

include("../src/main.jl")

@testset "Hyperbolic Geometry" begin
    S = hyperbolic_surface(2)
    X = universal_cover(S)
    D = fundamental_domain(X)
    id_sides = identified_sides(D)
    F = fundamental_group(S)
    
    @test length(id_sides) == 4

    w1 = F[2]^2 * F[1]^2 * F[3] * F[4]^-1
    c1 = free_homotopy_class(w, S)
    sg1 = surface_geodesic(c, QQ(1//2))

    w2 = conj(w, F[2])
    c2 = free_homotopy_class(w, S)
    sg2 = surface_geodesic(c, QQ(1//2))

    
    @test length(sg1) == length(sg2)
end
