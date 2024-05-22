export ProgressBar, next!

mutable struct ProgressBar{CT,T1<:Integer,T2<:Real}
    counter::CT
    max_counts::T1
    enable::Bool
    bar_width::T1
    start_time::T2
end

function ProgressBar(max_counts::Int; enable::Bool = true, bar_width::Int = 30)
    return ProgressBar(Threads.Atomic{Int}(0), max_counts, enable, bar_width, time())
end

function next!(p::ProgressBar, io::IO = stdout)
    p.counter[] >= p.max_counts && return

    Threads.atomic_add!(p.counter, 1)

    !p.enable && return

    counter = p.counter[]
    max_counts = p.max_counts
    bar_width = p.bar_width
    start_time = p.start_time

    percentage = counter / max_counts
    percentage_100 = lpad(round(100 * percentage, digits = 1), 5, " ")
    progress = floor(Int, bar_width * percentage)

    # Calculate the elapsed time in seconds
    elapsed_time = floor(Int, time() - start_time)
    # Convert the elapsed time into a string in hours, minutes and seconds
    elapsed_time_str = string(
        elapsed_time ÷ 3600,
        "h ",
        lpad((elapsed_time % 3600) ÷ 60, 2, "0"),
        "m ",
        lpad(elapsed_time % 60, 2, "0"),
        "s",
    )

    # Calculate the estimated time of arrival
    eta = floor(Int, elapsed_time / counter * (max_counts - counter))
    # convert eta into a string in hours, minutes and seconds
    eta_str = string(eta ÷ 3600, "h ", lpad((eta % 3600) ÷ 60, 2, "0"), "m ", lpad(eta % 60, 2, "0"), "s")

    # Construct the progress bar string
    bar = "[" * repeat("=", progress) * repeat(" ", bar_width - progress) * "]"

    print(io, "\rProgress: $bar $percentage_100% --- Elapsed Time: $elapsed_time_str (ETA: $eta_str)")

    counter == p.max_counts && print(io, "\n")

    return flush(io)
end
