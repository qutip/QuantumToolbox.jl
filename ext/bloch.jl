using QuantumToolbox
using Makie
using Colors
using LinearAlgebra: normalize
using LaTeXStrings

export add_points!, add_vectors!, add_line!, add_arc!, clear!, render

function QuantumToolbox.render(b::QuantumToolbox.Bloch; location = nothing)
    fig, location = _getFigAndLocation(location)
    if isnothing(b.fig) || b.fig !== fig
        b.fig = fig
        b.ax = Makie.Axis3(
            location;
            aspect = :data,
            limits = (-1.1, 1.1, -1.1, 1.1, -1.1, 1.1),
            xgridvisible = false,
            ygridvisible = false,
            zgridvisible = false,
            xticklabelsvisible = false,
            yticklabelsvisible = false,
            zticklabelsvisible = false,
            xticksvisible = false,
            yticksvisible = false,
            zticksvisible = false,
            xlabel = "",
            ylabel = "",
            zlabel = "",
            backgroundcolor = RGBAf(0, 0, 0, 0),
            xypanelvisible = false,
            xzpanelvisible = false,
            yzpanelvisible = false,
            xspinesvisible = false,
            yspinesvisible = false,
            zspinesvisible = false,
        )
        b.ax.azimuth[] = deg2rad(b.view_angles[1])
        b.ax.elevation[] = deg2rad(b.view_angles[2])
    else
        if !(b.ax in contents(location))
            location[] = b.ax
        end
    end
    Makie.empty!(b.ax)
    c_frame = Makie.to_color(b.frame_color)
    frame_color = Makie.RGBAf(c_frame.r, c_frame.g, c_frame.b, Float32(b.frame_alpha))
    c_sphere = Makie.to_color(b.sphere_color)
    sphere_color = Makie.RGBAf(c_sphere.r, c_sphere.g, c_sphere.b, Float32(b.sphere_alpha))
    sphere_primitive = Makie.Sphere(Makie.Point3f(0, 0, 0), 1.0f0)
    Makie.mesh!(b.ax, sphere_primitive; color = sphere_color, transparency = true)
    Makie.wireframe!(b.ax, sphere_primitive; color = frame_color, linewidth = b.frame_width, transparency = true)

    Makie.lines!(b.ax, [Makie.Point3f(-1, 0, 0), Makie.Point3f(1, 0, 0)]; color = :black, linewidth = 0.5) # X-axis
    Makie.lines!(b.ax, [Makie.Point3f(0, -1, 0), Makie.Point3f(0, 1, 0)]; color = :black, linewidth = 0.5) # Y-axis
    Makie.lines!(b.ax, [Makie.Point3f(0, 0, -1), Makie.Point3f(0, 0, 1)]; color = :black, linewidth = 0.5) # Z-axis

    Makie.text!(
        b.ax,
        Makie.Point3f(b.xlpos[1], 0, 0);
        text = b.xlabel[1],
        align = (:center, :center),
        color = Makie.to_color(b.font_color),
        fontsize = b.font_size,
    )
    Makie.text!(
        b.ax,
        Makie.Point3f(b.xlpos[2], 0, 0);
        text = b.xlabel[2],
        align = (:center, :center),
        color = Makie.to_color(b.font_color),
        fontsize = b.font_size,
    )
    Makie.text!(
        b.ax,
        Makie.Point3f(0, b.ylpos[1], 0);
        text = b.ylabel[1],
        align = (:center, :center),
        color = Makie.to_color(b.font_color),
        fontsize = b.font_size,
    )
    Makie.text!(
        b.ax,
        Makie.Point3f(0, b.ylpos[2], 0);
        text = b.ylabel[2],
        align = (:center, :center),
        color = Makie.to_color(b.font_color),
        fontsize = b.font_size,
    )
    Makie.text!(
        b.ax,
        Makie.Point3f(0, 0, b.zlpos[1]),
        text = b.zlabel[1],
        align = (:center, :center),
        color = Makie.to_color(b.font_color),
        fontsize = b.font_size,
    )
    Makie.text!(
        b.ax,
        Makie.Point3f(0, 0, b.zlpos[2]),
        text = b.zlabel[2],
        align = (:center, :center),
        color = Makie.to_color(b.font_color),
        fontsize = b.font_size,
    )

    for (i, pnt) in enumerate(b.points)
        color_idx = mod1(i, length(b.point_color))
        marker_idx = mod1(i, length(b.point_marker))
        size_idx = mod1(i, length(b.point_size))
        Makie.scatter!(
            b.ax,
            [Makie.Point3f(pnt...)];
            color = Makie.to_color(b.point_color[color_idx]),
            marker = b.point_marker[marker_idx],
            markersize = b.point_size[size_idx],
        )
    end

    for (i, vec) in enumerate(b.vectors)
        color_idx = mod1(i, length(b.vector_color))
        base_color = Makie.to_color(b.vector_color[color_idx])
        vec_rgba_color = Makie.RGBAf(base_color.r, base_color.g, base_color.b, Float32(0.6))
        Makie.arrows!(
            b.ax,
            [Makie.Point3f(0, 0, 0)],
            [Makie.Point3f(vec...)],
            color = vec_rgba_color,
            linewidth = 0.025,
            arrowsize = Makie.Vec3f(0.05, 0.05, 0.05),
            transparency = true,
        )
    end

    for line_pts in b.lines
        Makie.lines!(b.ax, Makie.Point3f.(line_pts); color = :blue, linewidth = 1.0)
    end

    for arc_pts in b.arcs
        Makie.lines!(b.ax, Makie.Point3f.([arc_pts[1], arc_pts[end]]); color = :orange, linewidth = 1.0)
    end

    return fig, b.ax
end
