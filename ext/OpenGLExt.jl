module OpenGLExt

import EnergyBalanceModel as EBM

import GLMakie

function EBM.Plot.init_backend(::Val{:GLMakie})::Module
    if GLMakie.Makie.current_backend() !== GLMakie
        GLMakie.activate!(; float=true, focus_on_show=true)
    end # if !==
    return GLMakie
end # function EBM.Plot.init_backend

end # module OpenGLExt
