module IO # EnergyBalanceModel.

using ..Utilities, ..Infrastructure
import ..Plot: unsafesave

import Dates, JLD2, TimeZones as TZ

export save, load!

function unsafesave(sols, path::String; spwarn::Bool=false)::String
    if !spwarn
        @warn "`unsafesave` may overwrite existing files. Use `save` instead."
    end # if !
    JLD2.save_object(path, sols)
    return path
end

"""
    save(obj, path::String=joinpath(pwd(), string(reprhex(unique_id()), ".dat"))) -> String

Save `obj` to the specified `path`. If a file already exists at `path`, it is renamed to
include a unique identifier before saving `obj`.

If `obj` is a `Makie.Figure`, additional keyword arguments are passed to `Makie.save`.

# Examples
```julia
julia> save("Hello World", "./greating.jld2")
"./greating.jld2"

julia> save("Hello again", "./greating.jld2")
┌ Warning: File ./greating.jld2 already exists. Last modified on 21 Oct 2025 at 11:03:13. The EXISTING file has been renamed to ./greating_768e991e.jld2.
└ @ EnergyBalanceModel.IO EnergyBalanceModel.jl/src/io.jl:48
"./greating.jld2"
```
"""
function save(obj, path::String=joinpath(pwd(), string(reprhex(unique_id()), ".dat")); kwargs...)::String
    if isfile(path)
        modified = Dates.format(
            TZ.astimezone(
                TZ.ZonedDateTime(Dates.unix2datetime(mtime(path)), TZ.tz"UTC"),
                TZ.localzone()
            ),
            Dates.dateformat"on d u Y at HH:MM:SS"
        ) # Dates.format
        nameext = splitext(path)
        newpath = string(nameext[1], '_', reprhex(unique_id()), nameext[2])
        @warn "File $path already exists. Last modified $modified. The EXISTING file has been renamed to $newpath."
        mv(path, newpath)
    end # if isfile
    return unsafesave(obj, path; spwarn=true, kwargs...)
end # function save

function unsafeload(path::String; spwarn::Bool=false)
    if !spwarn
        @warn "`unsafeload` could overwrite existing variables. Use `load!` instead."
    end # if !
    return JLD2.load_object(path)
end # function unsafeload

"""
    load!(to::Symbol, path::String, modu::Module=Main; house::Symbol=:SAFEHOUSE)

Load the object stored at `path` into the variable `to` in module `modu`. If a variable
named `to` already exists in `modu`, its value is moved to the safehouse specified by
`house` before loading the new value.

# Examples
```julia-repl
julia> save("Hello World", "./greating.jld2");

julia> load!(:greating, "./greating.jld2")
"Hello World"

julia> load!(:greating, "./greating.jld2")
┌ Warning: Variable `greating` already defined in Main. The existing value has been stored in safehouse `Main.safehouse` with ID 2f55a55a.
└ @ EnergyBalanceModel.IO EnergyBalanceModel.jl/src/io.jl:84
"Hello World"

julia> greating
"Hello World"
```
"""
function load!(to::Symbol, path::String, modu::Module=Main; house::Symbol=:SAFEHOUSE)
    if isdefined(modu, to)
        refugee = house!(to, safehouse(modu, house))
        @warn "Variable `$to` already defined in $modu. The existing value has been stored in safehouse `$modu.$safehouse` with ID $(reprhex(refugee.id))."
    end # if isdefined
    loaded = unsafeload(path; spwarn=true)
    @eval modu $to = $loaded
    return loaded
end # function load!

end # module IO
