module TransitoriosCabos

include("cabos/cabos.jl")
include("formulas/formulas.jl")
include("utils.jl")

using .Cabos
using .Formulas

export struct_to_dict, struct_from_dict, count_conductors_cable, outer_radius,
       CableComponent, CoaxialCable, PipeCable,
       zy_cabo, ynodal, ynodal_array, modos_propagacao

end # module
