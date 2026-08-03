# Tests for the `EnergyBalanceModel.Utilities` submodule.

# `@persistent` closes over `counter`, so the state survives between calls. Defined here at
# top level because the macro emits a `global const`.
UT.@persistent(
    counter::Int = 0,
    function persistent_counter(step::Int)::Int
        counter += step
        return counter
    end # function persistent_counter
) # UT.@persistent

@testset "Utilities" begin

    @testset "crossmean" begin
        @test UT.crossmean([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]) == [3.0, 4.0]
        @test UT.crossmean([[1.0, -1.0]]) == [1.0, -1.0]
        @test UT.crossmean([[1, 2], [3, 4]]) == [2.0, 3.0]
        # the result has one entry per gridbox, not per snapshot
        @test length(UT.crossmean([collect(1.0:5.0) for _ in 1:3])) == 5
        # ragged input is rejected by the bounds check
        @test_throws BoundsError UT.crossmean([[1.0, 2.0], [3.0]])
    end # @testset "crossmean"

    @testset "condset! and zeroref!" begin
        v = [1.0, 2.0, 3.0, 4.0]
        @test UT.condset!(copy(v), 0.0, >(2.0), v) == [1.0, 2.0, 0.0, 0.0]
        # the condition is evaluated on `ref`, the assignment happens on `to`
        @test UT.condset!(zeros(4), 7.0, >(2.0), v) == [0.0, 0.0, 7.0, 7.0]
        nanned = UT.condset!(copy(v), NaN, iszero, [0.0, 1.0, 0.0, 1.0])
        @test all(isnan, nanned[[1, 3]])
        @test nanned[[2, 4]] == v[[2, 4]]
        # condset! mutates and returns its first argument
        target = copy(v)
        @test UT.condset!(target, 0.0, isone, v) === target
        @test target == [0.0, 2.0, 3.0, 4.0]

        @test UT.zeroref!(copy(v), [1.0, 0.0, 3.0, 0.0]) == [1.0, 0.0, 3.0, 0.0]
        @test UT.zeroref!(copy(v), zeros(4)) == zeros(4)
        @test UT.zeroref!(copy(v), v) == v # nothing is zero in `v`
    end # @testset "condset! and zeroref!"

    @testset "iobuffer" begin
        buffer = UT.iobuffer(stdout)
        @test buffer isa IOContext
        @test get(buffer, :limit, false)
        @test get(buffer, :compact, false)
        @test get(buffer, :color, false)
        @test get(buffer, :displaysize, (0, 0)) == displaysize(stdout)
        @test get(UT.iobuffer(stdout; sizemodifier=(1, -2)), :displaysize, (0, 0)) ==
            displaysize(stdout) .+ (1, -2)
        # the underlying buffer is writable and starts empty
        print(buffer, "abc")
        @test String(take!(buffer.io)) == "abc"
    end # @testset "iobuffer"

    @testset "display_time" begin
        @test UT.display_time(0.0) == "0:00"
        @test UT.display_time(9.0) == "0:09"
        @test UT.display_time(65.0) == "1:05"
        @test UT.display_time(600.0) == "10:00"
        @test UT.display_time(Inf) == "-:--"
        @test UT.display_time(NaN) == "-:--"
    end # @testset "display_time"

    @testset "Progress" begin
        prog = UT.Progress(10, "Counting", 0.0)
        @test prog.total == 10
        @test prog.current == -1 # the first update! moves it to 0
        @test prog.updates == 0
        @test isnan(prog.started)
        @test prog.barwidth == 50 - (ndigits(10) * 2 + 1) - 2 - 5 - 3

        out = capture_stdout() do
            foreach(i -> UT.update!(prog, i), 0:10)
        end # do
        @test prog.current == 10
        @test prog.updates == 11
        @test !isnan(prog.started)
        @test occursin("Counting", out)
        @test occursin("Done", out) # printed once current == total
        @test occursin("10/10", out)

        # a step of one is the default
        prog2 = UT.Progress(3, "Default", 0.0)
        capture_stdout() do
            UT.update!(prog2)
            UT.update!(prog2)
        end # do
        @test prog2.current == 1

        # the custom info feed is appended below the bar
        prog3 = UT.Progress(4, "Feed", 0.0; width=60, infofeed=(t -> "t = $t"))
        out3 = capture_stdout() do
            UT.update!(prog3, 4; feedargs=(1.5,))
        end # do
        @test occursin("t = 1.5", out3)
        @test prog3.barwidth == 60 - (ndigits(4) * 2 + 1) - 2 - 5 - 3

        # with an infinite update frequency only the first and the final update print
        prog4 = UT.Progress(3, "Rare", Inf)
        first_out = capture_stdout() do
            UT.update!(prog4, 1)
        end # do
        quiet_out = capture_stdout() do
            UT.update!(prog4, 2)
        end # do
        @test !isempty(first_out)
        @test isempty(quiet_out)
        @test prog4.current == 2
    end # @testset "Progress"

    @testset "@persistent" begin
        # the state is carried between calls
        @test persistent_counter(2) == 2
        @test persistent_counter(3) == 5
        @test persistent_counter(0) == 5
        @test persistent_counter isa Function
    end # @testset "@persistent"

    @testset "@isdebugging" begin
        @test UT.@isdebugging() isa Bool
    end # @testset "@isdebugging"

end # @testset "Utilities"
