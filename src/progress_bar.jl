export ProgressBar

struct ProgressBar{CT,T1<:Integer, T2<:Real}
    counter::CT
    max_counts::T1
    enable::Bool
    bar_width::T1
    start_time::T2
    lock::ReentrantLock
end

function ProgressBar(max_counts; enable=true, bar_width=30)
    return ProgressBar(Ref{Int64}(0), max_counts, enable, bar_width, time(), ReentrantLock())
end

function next!(p::ProgressBar)
    p.counter[] >= p.max_counts && return

    lock(p.lock)

    p.counter[] += 1

    !p.enable && return

    counter = p.counter[]
    max_counts = p.max_counts
    bar_width = p.bar_width
    start_time = p.start_time
    
    percentage = counter / max_counts
    precentage_100 = round(100 * percentage, digits=2)
    progress = floor(Int, bar_width * percentage)

    # Calculate the elapsed time is seconds
    elapsed_time = time() - start_time
    # Convert the elapsed time into a string in hours, minutes and seconds
    elapsed_time_str = string(floor(Int, elapsed_time / 3600), "h ", floor(Int, (elapsed_time % 3600) / 60), "m ", floor(Int, elapsed_time % 60), "s")

    # Calculate the estimated time of arrival
    eta = elapsed_time / counter * (max_counts - counter)
    # convert eta into a string in hours, minutes and seconds
    eta_str = string(floor(Int, eta / 3600), "h ", floor(Int, (eta % 3600) / 60), "m ", floor(Int, eta % 60), "s")

    # Construct the progress bar string
    bar = "[" * repeat("=", progress) * repeat(" ", bar_width - progress) * "]"

    print("\rProgress: $bar $precentage_100%  --- Elapsed Time: $elapsed_time_str  (ETA: $eta_str)     ")
    flush(stdout)

    unlock(p.lock)

    return
end