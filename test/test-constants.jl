using QEDprobing
using Unitful

@testset "Fine Structure Constant" begin
    groundtruth = ELEMENTARY_CHARGE^2 / (4 * pi)
    @test isapprox(groundtruth, ALPHA, atol = 0.0, rtol = 1.0e-6)
end


@testset "Hartree" begin

    # source: https://en.wikipedia.org/wiki/Hartree
    # in eV
    groundtruth = HBARC * ALPHA / (QEDprobing.BOHR_RADIUS_ANG * 1.0e5) * 1.0e6
    @test isapprox(HARTREE, groundtruth, atol = 0.0, rtol = 1.0e-6)
end

@testset "Bohr Radius" begin

    # source: https://en.wikipedia.org/wiki/Bohr_radius
    # in Angstrom
    bohr_inv_me = 1 / ALPHA
    bohr_inv_eV = bohr_inv_me / ELECTRONMASS
    bohr_inv_MeV = bohr_inv_eV * 1.0e6
    bohr_fm = bohr_inv_MeV * HBARC
    bohr_ang = bohr_fm * 1.0e-5

    @test isapprox(BOHR_RADIUS_ANG, bohr_ang, atol = 0.0, rtol = 1.0e-6)

end

# numbers taken from: arXiv:2503.20433v1, Table 1
GROUNDTRUTH_FERMI_ENERGY = [
    # rs(bohr), beta=1/EF (inv Hartree), kF (inv Bohr)
    (2.0, 2.172042885466092, 0.9595791463387565),
    (5.0, 13.575268034163074, 0.3838316585355026),
    (10.0, 54.301072136652294, 0.1919158292677513),
    (20.0, 217.20428854660918, 0.09595791463387565),
    (50.0, 1357.5268034163075, 0.03838316585355026),
    (100.0, 5430.10721366523, 0.01919158292677513),
    (200.0, 21720.42885466092, 0.009595791463387566),
]

ne_from_rs(rs) = inv(4 * pi * rs^3 / 3)

@testset "Fermi Energy/Wavevector" begin
    @testset "$rs" for (rs, gt_beta_inv_hartree, gt_kF_inv_bohr) in GROUNDTRUTH_FERMI_ENERGY

        ne_inv_bohr = ne_from_rs(rs)
        kF_inv_bohr = QEDprobing._fermi_wave_vector(ne_inv_bohr)

        @test isapprox(kF_inv_bohr, gt_kF_inv_bohr, atol = 0.0, rtol = 1.0e-6)

        ne_inv_cm = ne_inv_bohr / ((BOHR_RADIUS_ANG * 1.0e-8)^3)

        ne_inv_me = QEDprobing._internalize_density(ne_inv_cm * 1u"cm^(-3)")
        kF_me = QEDprobing._fermi_wave_vector(ne_inv_me)

        EF_me = QEDprobing._fermi_energy_from_kF(kF_me)

        EF_hartree = EF_me * ELECTRONMASS / HARTREE

        @test isapprox(inv(EF_hartree), gt_beta_inv_hartree, atol = 0.0, rtol = 1.0e-6)
    end
end
