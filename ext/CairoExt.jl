module CairoExt

import EnergyBalanceModel as EBM

import CairoMakie

EBM.Plot.isloaded(::Val{:CairoMakie})::Bool = true

function EBM.Plot.init_backend(::Val{:CairoMakie})::Module
    CairoMakie.Makie.current_backend() === CairoMakie || CairoMakie.activate!()
    return CairoMakie
end # function EBM.Plot.init_backend

EBM.Plot.precompile(CairoMakie)

end # module CairoExt
