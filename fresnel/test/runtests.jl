using fresnel
using Test
using FresnelIntegrals

@testset "my implementation of Fresnel C(x) function is accurate" begin

    @testset "interval (0,1)" begin
        x_vec = [0.0, 1e-6, 1e-3, 0.01, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999]
        for x in x_vec
            approx = fresnel_C(x) # my implementation
            reference = fresnelc(x) # reference from FresnelIntegrals
            @test isapprox(approx, reference; rtol=5e-11)
        end
    end

    @testset "1 and beyond" begin
        x_vec = [1, 1+1e-6, 1+1e-3, 1.01, 1.3, 1.5, 2, 3, 5, 10, 13, 20, 40, 50, 100, 200, 300, 1000, 5e3, 1e4, 5e4, 1e5]
        for x in x_vec
            approx = fresnel_C(x) # my implementation
            reference = fresnelc(x) # reference from FresnelIntegrals
            @test isapprox(approx, reference; rtol=5e-11)
        end
    end

    @testset "negative inputs" begin
        x_vec = [-1e-6, -1e-3, -0.01, -0.3, -0.99, -0.999, -1, -(1+1e-6), -1.01,  -1.5, -2, -5, -10,-20, -100, -1000, -5e3, -1e4, -5e4, -1e5]
        for x in x_vec
            approx = fresnel_C(x) # my implementation
            reference = fresnelc(x) # reference from FresnelIntegrals
            @test isapprox(approx, reference; rtol=5e-11)
        end
    end

    @testset "exact anti symmetry" begin
        for x in [0.3, 0.99, 1, 1.5, 10, 1e4]
            @test -fresnel_C(-x) == fresnel_C(x)   # exact, not isapprox
        end
    end
    
    @testset "continuity across X_cutoff" begin
        eps = 1e-9
        @test isapprox(fresnel_C(1.0 - eps), fresnel_C(1.0 + eps); rtol=5e-11)
    end
    

end
