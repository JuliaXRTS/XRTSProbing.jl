include("common.jl")

function _f(xi, A, B)
    return inv(1 + exp(A * xi^2 - B))
end

function _f_prime(xi, A, B)
    term = A * xi^2 - B
    return -2 * A * xi * exp(term) / ((1 + exp(term))^2)
end

function general_integrand2(xi, A, B)
    if isone(xi)
        return -xi * _f_prime(xi, A, B)
    end

    term = (1 - xi^2) * log(abs((xi - 1) / (xi + 1))) / 2

    return _f_prime(xi, A, B) * (-xi + term)
end

function general_integral2(nu, bbar)
    A = bbar * nu^2
    B = _chemical_potential_normalized(bbar) * bbar
    res, _ = quadgk(x -> general_integrand2(x, A, B), 0.0, 1.0, Inf)

    return -nu^2 * res
end

beta_arr = [1.0e-3]

P = plot(
    xscale = :log10,
    xlab = L"\nu",
    ylab = L"g(\nu,\bar\beta)",
    legend = :bottomleft,
    legend_title = L"\bar\beta",
    palette = :seaborn_colorblind,
    title = "\$g(\\nu,\\bar\\beta)\$ (general vs. non-deg)",
)

for beta in reverse(beta_arr)
    @show beta
    #nu_arr = range(1e-3,10000,103)
    nu_arr = 10.0 .^ range(-4, 4, 300)

    gen_res = general_integral2.(nu_arr, beta)
    non_deg_res = non_deg_integal.(nu_arr, beta)

    plot!(P, nu_arr, gen_res, lw = 3, lab = @sprintf("%01.0e", beta))
    plot!(P, nu_arr, non_deg_res, ls = :dash, c = :black, lab = "")
end


file_path = joinpath(PLOTDIR, "g_func_stable.pdf")
savefig(P, file_path)
@info "saved in $file_path"
