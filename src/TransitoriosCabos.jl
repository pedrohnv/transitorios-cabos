module TransitoriosCabos

export struct_to_dict, struct_from_dict, count_conductors_cable, CableComponent, CoaxialCable, PipeCable
export zy_cabo, ynodal, ynodal_array, modos_propagacao

include("cabos.jl")
include("formulas.jl")

end  # module
