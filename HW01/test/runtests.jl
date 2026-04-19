using Test
using HW01

@testset "Cubic spline interpolation is accurate" begin
    
    @testset "sin" begin
        x = -5:1:5
        y = [sin(i) for i in x]

        Z = interpoliraj(x, y)

        for (xi, yi) in zip(x, y)
            @test isapprox(vrednost(Z, float(xi)), yi; atol=1e-8)
        end
    end

    @testset "linear function" begin
        x = 0:1:10
        y = [2*i + 3 for i in x]

        Z = interpoliraj(x, y)

        for xi in 0:0.5:10
            @test isapprox(vrednost(Z, xi), 2*xi + 3; atol=1e-8)
        end
    end

    @testset "constant function" begin
        x = -10:2:10
        y = fill(5.0, length(x))

        Z = interpoliraj(x, y)

        for xi in -10:0.5:10
            @test isapprox(vrednost(Z, xi), 5.0; atol=1e-8)
        end
    end

    @testset "uneven spaces" begin
        x = [0.0, 1.0, 2.0, 10.0, 50.0]
        y = [sin(i) for i in x]

        Z = interpoliraj(x, y)

        for (xi, yi) in zip(x, y)
            @test isapprox(vrednost(Z, float(xi)), yi; atol=1e-8)
        end
    end

end ;

