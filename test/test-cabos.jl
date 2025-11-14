# Testa os construtores dos Cabos.

using Test
include("../src/cabos.jl")

@testset "struct_to_dict and struct_from_dict" begin
    comp = CableComponent(0.0, 1.0, 1.1, 1.72, 1.0, 1.0, 2.0, "core")
    d = struct_to_dict(comp)
    @test typeof(d) == Dict{String,Any}
    comp2 = struct_from_dict(d)
    @test comp2.radius_in ≈ comp.radius_in
    @test comp2.radius_ext ≈ comp.radius_ext
    @test comp2.radius_ext_insulator ≈ comp.radius_ext_insulator
    @test comp2.rho_c ≈ comp.rho_c
    @test comp2.name == comp.name
end

@testset "CableComponent validations" begin
    @test_throws ArgumentError CableComponent(-1.0, 1.0, 1.1, 1.72, 1.0, 1.0, 2.0)
    @test_throws ArgumentError CableComponent(2.0, 1.0, 0.8, 1.72, 1.0, 1.0, 2.0)
end

@testset "CoaxialCable and PipeCable constructors and structure" begin
    inner_comp = CableComponent(0.0, 1.0, 1.1, 1.6, 1.0, 1.0, 2.0, "inner")
    outer_comp = CableComponent(1.1, 2.0, 2.05, 1.7, 1.0, 1.0, 2.0, "outer")
    coax = CoaxialCable([inner_comp, outer_comp], 0.0, 0.0, "CoaxA", 10.0)
    @test outer_radius(coax) ≈ outer_comp.radius_ext_insulator
    @test coax.name == "CoaxA"
    # Test show
    io = IOBuffer()
    show(io, coax)
    str = String(take!(io))
    @test occursin("CoaxialCable", str)
    @test occursin("componentes = 2", str)

    # Invalid connection between components
    bad_comp = CableComponent(1.2, 2.2, 2.3, 1.7, 1.0, 1.0, 2.0)
    @test_throws ArgumentError CoaxialCable([outer_comp, bad_comp])
end

@testset "PipeCable validations and methods" begin
    inner_comp = CableComponent(0.0, 1.0, 1.1, 1.6, 1.0, 1.0, 2.0)
    coax = CoaxialCable([inner_comp], 0.0, 0.0)
    pipe = PipeCable(2.5, 3.0, 3.3, 1.8, 1.2, 1.2, 1.25, 2.2, 2.3; cables=[coax], x=0.0, y=0.0)
    @test outer_radius(pipe) ≈ pipe.radius_ext_insulator
    @test pipe.name == "Pipe"
    io = IOBuffer()
    show(io, pipe)
    str = String(take!(io))
    @test occursin("PipeCable", str)

    # Test invalid nested cable positioning
    coax2 = CoaxialCable([inner_comp], 2.6, 2.6)
    @test_throws ArgumentError PipeCable(2.5, 3.0, 3.2, 1.8, 1.2, 1.2, 1.25, 2.2, 2.3; cables=[coax2], x=0.0, y=0.0)
end

@testset "count_conductors_cable and name_index_dict" begin
    inner_comp = CableComponent(0.0, 1.0, 1.1, 1.6, 1.0, 1.0, 2.0, "center")
    outer_comp = CableComponent(1.1, 2.0, 2.05, 1.7, 1.0, 1.0, 2.0, "shell")
    coax = CoaxialCable([inner_comp, outer_comp], 0.0, 0.0, "Coax", 10.0)
    pipe = PipeCable(2.5, 3.0, 3.05, 2.3, 1.8, 1.15, 1.20, 2.2, 2.3; cables=[coax], x=0.0, y=0.0)
    @test count_conductors_cable(coax) == 2
    @test count_conductors_cable(pipe) == 3

    # Test name_index_dict methods
    coax_comp = deepcopy(coax)
    coax_comp.components[1]._index = 100
    coax_comp.components[2]._index = 101
    d = name_index_dict(coax_comp)
    @test d[100] == coax_comp.components[1].name
    @test d[101] == coax_comp.components[2].name

    pipe._index = 105
    pipe.cables[1].components[1]._index = 110
    d2 = name_index_dict(pipe)
    @test d2[105] == "Pipe"
    @test d2[110] == "center"
end

@testset "move_group" begin
    inner_comp = CableComponent(0.0, 1.0, 1.1, 1.6, 1.0, 1.0, 2.0, "core")
    coax = CoaxialCable([inner_comp], 0.0, 0.0)
    pipe = PipeCable(2.1, 3.0, 3.05, 2.3, 1.8, 1.15, 1.20, 2.2, 2.3; cables=[coax], x=10.0, y=20.0)
    move_group(pipe, 5.0, -2.0)
    @test isapprox(pipe.x, 15.0)
    @test isapprox(pipe.y, 18.0)
    @test isapprox(pipe.cables[1].x, 5.0)
    @test isapprox(pipe.cables[1].y, -2.0)
end
