module OpenGLExt

import EnergyBalanceModel as EBM

import GLMakie

EBM.Plot.isloaded(::Val{:GLMakie})::Bool = true

function EBM.Plot.init_backend(::Val{:GLMakie})::Module
    GLMakie.Makie.current_backend() === GLMakie || GLMakie.activate!()
    return GLMakie
end # function EBM.Plot.init_backend

EBM.Plot.precompile(GLMakie)

end # module OpenGLExt
