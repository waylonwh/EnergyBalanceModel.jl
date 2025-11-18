module OpenGLExt

import EnergyBalanceModel as EBM

import GLMakie

function EBM.Plot.init_backend(::Val{:GLMakie})::Module
    if GLMakie.Makie.current_backend() !== GLMakie
        GLMakie.activate!(; float=true, focus_on_show=true)
    end # if !==
    return GLMakie
end # function EBM.Plot.init_backend

for t in (Float64, Int)
    precompile(EBM.Plot.contourf_tiles, (Vector{Float64}, EBM.Vec, EBM.Layout{Matrix{Float64}}))
end # for t

end # module OpenGLExt
