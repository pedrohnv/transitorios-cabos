using Test

@testset "Example tests" begin

    @testset "Cabos" begin
        include("test-cabos.jl")
    end

    @testset "Fórmulas" begin
        include("test-formulas.jl")
    end

end