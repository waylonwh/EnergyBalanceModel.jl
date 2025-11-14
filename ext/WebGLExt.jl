module WebGLExt

import EnergyBalanceModel as EBM

import WGLMakie

EBM.Plot.isloaded(::Val{:WGLMakie})::Bool = true

function EBM.Plot.init_backend(::Val{:WGLMakie})::Module
    if WGLMakie.Makie.current_backend() !== WGLMakie
        WGLMakie.activate!()
    end # if !==
    return WGLMakie
end # function EBM.Plot.init_backend

end # module WebGLExt
