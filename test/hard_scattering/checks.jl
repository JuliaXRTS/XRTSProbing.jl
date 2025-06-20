function _onshell_check(part, mom, atol, rtol)

    return @test isapprox(
        getMass2(mom),
        mass(part)^2,
        atol = atol,
        rtol = rtol,
    )
end
