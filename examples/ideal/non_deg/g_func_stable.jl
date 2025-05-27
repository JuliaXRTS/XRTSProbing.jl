# TODO:
# - plot negative values of nu as well


include("common.jl")

# this is slighly better than LogExpFunctions version
function _f_prime(xi, A, B)
    #u = A * xi^2 - B
    #sig = logistic(u)
    #return -2 * A * xi * sig*(one(sig)-sig)
    u = A * xi^2 - B
    # Compute sigmoid derivative stably:
    if u ≥ 0
        exp_neg_u = exp(-u)
        s_deriv = exp_neg_u / (1 + exp_neg_u)^2
    else
        exp_u = exp(u)
        s_deriv = exp_u / (1 + exp_u)^2
    end
    return -2 * A * xi * s_deriv
end

function transformed_integrand(xi, A, B)
    f_prime_val = _f_prime(xi, A, B)
    part1 = xi * f_prime_val
    prefac = (xi^2 - 1) / 2
    part2 = -f_prime_val * prefac * log(
        abs(
            (xi - 1) / (xi + 1)
        )
    )

    return part1 + part2
end

function transformed_integral_positive(nu, bbar)
    A = bbar * nu^2
    B = _chemical_potential_normalized(bbar) * bbar
    res, _ = quadgk(x -> transformed_integrand(x, A, B), 0.0, 1.0, Inf)

    return nu^2 * res
end

function transformed_integral(nu, bbar)
    if nu < zero(nu)
        return -transformed_integral_positive(-nu, bbar)
    end
    return transformed_integral_positive(nu, bbar)
end


beta_arr = 10.0 .^ (-3:2)
P = plot(
    xscale = :log10,
    xlab = L"\nu",
    ylab = L"g(\nu,\bar\beta)",
    legend = :topleft,
    legend_title = L"\bar\beta",
    palette = :seaborn_colorblind,
    #title = "\$-g(\\nu,\\bar\\beta)\$ (general vs. non-deg)"
)

for beta in reverse(beta_arr)
    #nu_arr = range(1e-3,10000,103)

    nu_arr = 10.0 .^ range(-4, 4, 300)
    gen_res = -general_integal.(nu_arr, beta)
    non_deg_res = -transformed_integral.(nu_arr, beta)


    plot!(P, nu_arr, gen_res, lw = 3, lab = @sprintf("%01.0e", beta))
    plot!(P, nu_arr, non_deg_res, ls = :dash, c = :black, lab = "")
end

P2 = plot(
    xscale = :log10,
    xlab = L"\nu",
    yscale = :log10,
    ylab = "relative error",
    legend = :outerright,
    legend_title = L"\bar\beta",
    palette = :seaborn_colorblind,
)

for beta in reverse(beta_arr)
    #nu_arr = range(1e-3,10000,103)
    nu_arr = 10.0 .^ range(-4, 4, 103)

    gen_res = general_integal.(nu_arr, beta)
    non_deg_res = transformed_integral.(nu_arr, beta)

    rel_err = @. abs(gen_res - non_deg_res) / abs(gen_res)


    plot!(P2, nu_arr, rel_err, lw = 3, lab = @sprintf("%01.0e", beta))
end

P3 = plot(
    xscale = :log10,
    xlab = L"\nu",
    yscale = :log10,
    ylab = "absolute error",
    legend = :outerright,
    legend_title = L"\bar\beta",
    palette = :seaborn_colorblind,
)

for beta in reverse(beta_arr)
    #nu_arr = range(1e-3,10000,103)
    nu_arr = 10.0 .^ range(-4, 4, 103)

    gen_res = general_integal.(nu_arr, beta)
    non_deg_res = transformed_integral.(nu_arr, beta)

    abs_err = @. abs(gen_res - non_deg_res)


    plot!(P3, nu_arr, abs_err, lw = 3, lab = @sprintf("%01.0e", beta))
end


l = @layout [a{0.2w} b{0.39w} c{0.39w}]
Pall = plot(
    P,
    P2,
    P3,
    plot_title = "\$-g(\\nu,\\bar\\beta)\$ (general vs. transformed)",
    layout = l,
    size = (1400, 400),
    bottom_margin = 40px,
    leftmargin = 30px,
)


file_path = joinpath(PLOTDIR, "g_func_stable.pdf")
savefig(Pall, file_path)
@info "saved in $file_path"
