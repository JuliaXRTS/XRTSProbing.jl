function _is_onshell_photon(mom)
    th = getTheta(mom)
    sth, cth = sincos(th)
    phi = getPhi(mom)
    sph, cph = sincos(phi)
    e = getE(mom)
    return isapprox(getX(mom) / (sth * cph), e) &&
        isapprox(getY(mom) / (sth * sph), e) &&
        isapprox(getZ(mom) / cth, e)
end

_axis_index(::XAxis) = 2
_axis_index(::YAxis) = 3
_axis_index(::ZAxis) = 4
