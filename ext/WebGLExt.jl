module WebGLExt

import EnergyBalanceModel as EBM

import WGLMakie

EBM.Plot.isloaded(::Val{:WGLMakie})::Bool = true

function EBM.Plot.init_backend(::Val{:WGLMakie})::Module
    WGLMakie.Makie.current_backend() === WGLMakie || WGLMakie.activate!()
    return WGLMakie
end # function EBM.Plot.init_backend

EBM.Plot.precompile(WGLMakie)

end # module WebGLExt
