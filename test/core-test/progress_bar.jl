@testitem "Progress Bar" begin
    bar_width = 30
    strLength = 67 + bar_width # including "\r" in the beginning of the string
    prog = ProgressBar(bar_width, enable = true, bar_width = bar_width, interval = 0.2)
    for p in 1:bar_width
        sleep(0.3)

        output = sprint((t, s) -> next!(s, t), prog)

        if p < bar_width
            @test length(output) == strLength
        else # the last output has an extra "\n" in the end
            @test length(output) == strLength + 1
        end
    end

    prog = ProgressBar(bar_width, enable = true, bar_width = bar_width, interval = 0.2)
    # No delay, so no output printed
    for p in 1:bar_width
        output = sprint((t, s) -> next!(s, t), prog)

        if p < bar_width
            @test length(output) == 0
        else
            @test length(output) == strLength + 1
        end
    end
end
