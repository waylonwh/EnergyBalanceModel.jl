module IO # EnergyBalanceModel.

using ..Utilities, ..Infrastructure

import Dates, JLD2, TimeZones as TZ

export save, load!

function unsafesave(sols, path::String; spwarn::Bool=false)::String
    if !spwarn
        @warn "`unsafesave` may overwrite existing files. Use `save` instead."
    end # if !
    JLD2.save_object(path, sols)
    return path
end

function save(obj, path::String=joinpath(pwd(), string(reprhex(unique_id()), ".dat")))::String
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
    return unsafesave(obj, path; spwarn=true)
end # function save

function unsafeload(path::String; spwarn::Bool=false)
    if !spwarn
        @warn "`unsafeload` could overwrite existing variables. Use `load!` instead."
    end # if !
    return JLD2.load_object(path)
end # function unsafeload

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
