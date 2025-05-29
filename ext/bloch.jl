using QuantumToolbox
using Makie
using LinearAlgebra

function QuantumToolbox.render(b::QuantumToolbox.Bloch; location = nothing)
    fig, location = _getFigAndLocation(location)

    if isnothing(b.fig) || b.fig !== fig
        b.fig = fig
        b.ax = Makie.Axis3(
            location;
            aspect = :data,
            limits = (-1.02, 1.02, -1.02, 1.02, -1.02, 1.02),
            # Disable all axis elements
            xgridvisible = false, ygridvisible = false, zgridvisible = false,
            xticklabelsvisible = false, yticklabelsvisible = false, zticklabelsvisible = false,
            xticksvisible = false, yticksvisible = false, zticksvisible = false,
            xlabel = "", ylabel = "", zlabel = "",
            backgroundcolor = RGBAf(1, 1, 1, 0.0),  # Transparent background
            xypanelvisible = false, xzpanelvisible = false, yzpanelvisible = false,
            xspinesvisible = false, yspinesvisible = false, zspinesvisible = false,
            protrusions = (0, 0, 0, 0),  # Remove margins
            viewmode = :fit
        )
        # Set view angles
        b.ax.azimuth[] = deg2rad(b.view_angles[1])
        b.ax.elevation[] = deg2rad(b.view_angles[2])
    else
        if !(b.ax in contents(location))
            location[] = b.ax
        end
    end

    Makie.empty!(b.ax)

    # ===== BASE VISUALS =====
    # Sphere
    sphere_color = RGBAf(0.9, 0.9, 0.95, 0.15)
    Makie.mesh!(b.ax, Makie.Sphere(Point3f(0, 0, 0), 1.0f0); 
        color = sphere_color, transparency = false)

    # Wireframe (equator + 2 meridians)
    wire_color = RGBAf(0.5, 0.5, 0.5, 0.4)
    φ = range(0, 2π, length=100)
    θ = range(-π, π, length=100)
    Makie.lines!(b.ax, [Point3f(cos(φ), sin(φ), 0) for φ in φ];  # Equator
        color=wire_color, linewidth=1.0)
    Makie.lines!(b.ax, [Point3f(0, sin(θ), cos(θ)) for θ in θ];  # Meridian 1
        color=wire_color, linewidth=1.0)
    Makie.lines!(b.ax, [Point3f(sin(θ), 0, cos(θ)) for θ in θ];  # Meridian 2
        color=wire_color, linewidth=1.0)

    # Axes
    axis_color = RGBAf(0.3, 0.3, 0.3, 0.7)
    Makie.lines!(b.ax, [Point3f(-1.1, 0, 0), Point3f(1.1, 0, 0)]; 
        color=axis_color, linewidth=1.5)
    Makie.lines!(b.ax, [Point3f(0, -1.1, 0), Point3f(0, 1.1, 0)];
        color=axis_color, linewidth=1.5)
    Makie.lines!(b.ax, [Point3f(0, 0, -1.1), Point3f(0, 0, 1.1)];
        color=axis_color, linewidth=1.5)

    # Labels
    label_color = RGBf(0.2, 0.2, 0.2)
    label_size = 20
    Makie.text!(b.ax, "x"; position=Point3f(1.15, 0, 0), color=label_color, fontsize=label_size)
    Makie.text!(b.ax, "y"; position=Point3f(0, 1.15, 0), color=label_color, fontsize=label_size)
    Makie.text!(b.ax, "z"; position=Point3f(0, 0, 1.15), color=label_color, fontsize=label_size)

    # ===== QUANTUM STATE VISUALS =====
    for line_pts in b.lines
        Makie.lines!(b.ax, Makie.Point3f.(line_pts); color = :blue, linewidth = 1.0)
    end

    # Arcs (great circle connections)
   for arc_pts in b.arcs
        if length(arc_pts) >= 2
            v1 = normalize(arc_pts[1])
            v2 = normalize(arc_pts[end])
            n = normalize(cross(v1, v2))
            θ = acos(clamp(dot(v1, v2), -1.0, 1.0))
            if length(arc_pts) == 3
                vm = normalize(arc_pts[2])
                if dot(cross(v1, vm), n) < 0
                    θ = 2π - θ
                end
            end
            t_range = range(0, θ, length=50)
            arc_points = [Makie.Point3f((v1*cos(t) + cross(n, v1)*sin(t))...) for t in t_range]
            Makie.lines!(b.ax, arc_points; color=:orange, linewidth=1.3)
        end
    end


    # Points
    if !isempty(b.points)
        for (i, p) in enumerate(b.points)
            color = get(b.point_color, i, RGBf(0.08, 0.02, 0.02)) # Red default
            marker = get(b.point_marker, i, :circle)
            size = get(b.point_size, i, 5)
            
            Makie.scatter!(b.ax, [Point3f(p)]; 
                color=color,
                marker=marker,
                markersize=size,
                strokecolor=:black,
                strokewidth=1)
        end
    end

    # Vectors
    if !isempty(b.vectors)
        for (i, v) in enumerate(b.vectors)
            color = get(b.vector_color, i, RGBAf(0.2, 0.4, 0.8, 0.9)) # Blue default
            Makie.arrows!(b.ax, 
                [Point3f(0, 0, 0)], 
                [Point3f(normalize(v))];
                color=color,
                linewidth=0.025,
                arrowsize=Vec3f(0.08, 0.08, 0.08))
        end
    end

    return fig, b.ax
end

# Helper function to safely get array elements
get(arr::AbstractVector, i::Int, default) = i ≤ length(arr) ? arr[i] : default
