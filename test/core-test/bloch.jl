
@testset "Bloch" begin
    # 1. Construction
    b = Bloch()
    @test b isa Bloch

    # 2. set_label_convention!
    for conv in
        ["original", "xyz", "sx sy sz", "01", "polarization jones", "polarization jones letters", "polarization stokes"]
        set_label_convention!(b, conv)
        @test length(b.xlabel) == 2
        @test length(b.ylabel) == 2
        @test length(b.zlabel) == 2
    end
    @test_throws ErrorException set_label_convention!(b, "notaconvention")

    # 3. clear!
    add_points!(b, [1.0, 0.0, 0.0])
    add_vectors!(b, [1.0, 0.0, 0.0])
    add_annotation!(b, [1.0, 0.0, 0.0], "test")
    clear!(b)
    @test isempty(b.points)
    @test isempty(b.vectors)
    @test isempty(b.annotations)

    # 4. add_points!
    add_points!(b, [1.0, 0.0, 0.0])
    @test !isempty(b.points)
    clear!(b)
    add_points!(b, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    @test !isempty(b.points)
    clear!(b)
    @test_throws ErrorException add_points!(b, [1.0, 2.0])

    # 5. add_states!
    q = Qobj([1.0, 0.0])
    add_states!(b, q)
    @test !isempty(b.vectors)
    clear!(b)
    add_states!(b, q; kind = "point")
    @test !isempty(b.points)
    clear!(b)
    qs = [Qobj([1.0, 0.0]), Qobj([0.0, 1.0])]
    add_states!(b, qs)
    @test length(b.vectors) == 2
    clear!(b)
    @test_throws ErrorException add_states!(b, q; kind = "badkind")

    # 6. add_vectors!
    add_vectors!(b, [1.0, 0.0, 0.0])
    @test !isempty(b.vectors)
    clear!(b)
    add_vectors!(b, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    @test length(b.vectors) == 2
    clear!(b)
    @test_throws ErrorException add_vectors!(b, [1.0, 2.0])  # invalid

    # 7. add_annotation!
    add_annotation!(b, [1.0, 0.0, 0.0], "label")
    @test !isempty(b.annotations)
    clear!(b)
    add_annotation!(b, Qobj([1.0, 0.0]), "label")
    @test !isempty(b.annotations)
    clear!(b)
    @test_throws ErrorException add_annotation!(b, [1.0, 0.0], "bad")  # not 3D

    # 8. add_line! and add_arc!
    add_line!(b, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    @test !isempty(b._lines)
    clear!(b)
    add_arc!(b, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    @test !isempty(b._arcs)
    clear!(b)
    @test_throws ErrorException add_arc!(b, [0.0, 0.0, 0.0], [0.0, 1.0, 0.0])  # origin

    # 9. render! and show! (smoke test)
    add_points!(b, [1.0, 0.0, 0.0])
    render!(b)
    show!(b)
    clear!(b)

    # 10. save (smoke test)
    add_points!(b, [1.0, 0.0, 0.0])
    save(b; name = "test_bloch.png")
end
