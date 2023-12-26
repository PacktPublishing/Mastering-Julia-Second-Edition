using Funky,Test

@testset "Maths Functions" begin
    @test fac(5)  == 120
    @test fib(10) == 55
    @test abs(basel(10^7) - (pi*pi/6)) < 10e-6
    @test size(hailstone(117))[1] == 21
end
