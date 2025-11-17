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

end # module CairoExt
