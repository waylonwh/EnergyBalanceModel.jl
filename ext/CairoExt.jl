module CairoExt

import EnergyBalanceModel as EBM

import CairoMakie

EBM.Plot.isloaded(::Val{:CairoMakie})::Bool = true

function EBM.Plot.init_backend(::Val{:CairoMakie})::Module
    if CairoMakie.Makie.current_backend() !== CairoMakie
        CairoMakie.activate!()
    end # if !==
    return CairoMakie
end # function EBM.Plot.init_backend

for t in (Float64, Int)
    precompile(EBM.Plot.contourf_tiles, (Vector{Float64}, EBM.Vec, EBM.Layout{Matrix{Float64}}))
end # for t

end # module CairoExt
