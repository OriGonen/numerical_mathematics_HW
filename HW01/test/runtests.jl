using Test
using HW01

@testset "Greet should run fine" begin
    HW01.greet()
    @test true
end