let x = nothing
    global func(a) = a+1
    global function func(a, b)
        return a+b
    end
end

@show func(3)
@show func(3, 4)
