### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> tags = ["models", "astrostatistics", "Bayesian", "MCMC", "query", "scatter plots", "errorbars", "line plots", "corner plots"]
#> title = "Modeling 3: Bayesian model fitting"
#> description = "Estimate the full uncertainty in your fitted parameters."
#> date = "2026-10-12"
#> layout = "layout.jlhtml"
#> image = "/assets/3_models-bayesian.png"

using Markdown
using InteractiveUtils

# ╔═╡ c4c69a9c-0269-4022-b7cb-2c04e35c3483
begin
    import Pkg
    Pkg.activate(; temp = true)
    Pkg.add(
        [
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(; name = "LinearAlgebra"),
            Pkg.PackageSpec(; name = "Random"),
            Pkg.PackageSpec(; name = "Turing"),
            Pkg.PackageSpec(; name = "PairPlots"),
            Pkg.PackageSpec(; name = "DataFramesMeta"),
            Pkg.PackageSpec(; name = "VirtualObservatory"),
            Pkg.PackageSpec(; name = "CairoMakie"),
            Pkg.PackageSpec(; name = "MathTeXEngine"),
        ]
    )

    # Analysis
    using LinearAlgebra: Diagonal
    using Random: Xoshiro
    using Turing: @model, NUTS, Normal, sample, summarystats, truncated

    # Data handling and visualization
    using DataFramesMeta: @select, DataFrame, dropmissing
    using VirtualObservatory: VizierCatalog, table
    using CairoMakie
    using PairPlots: pairplot
    using MathTeXEngine: set_texfont_family!, FontFamily
    set_texfont_family!(FontFamily("TeXGyreHeros"))

    using PlutoUI: TableOfContents
end;

# ╔═╡ d03f494f-4a91-4ff8-b5c6-cffafe540fcc
md"""
# Modeling 3: Bayesian model fitting

This notebook is adapted from the original JuliaAstro [curve fitting tutorial](https://github.com/JuliaAstro/JuliaAstro.github.io/blob/a8dfa6773ff27ff342653459933785d34e7ced55/docs/src/tutorials/curve-fit.md) by William Thompson.

!!! tip "Learning goals"
    - Build Bayesian models with Turing.jl.
    - Account for intrinsic scatter in our data.
    - Draw samples from a posterior with `NUTS`.
    - Plot posterior samples with the data.
    - Visualize a posterior with a corner plot from PairPlots.jl.

!!! warning "Companion content"
    [learn.JuliaAstro > Modeling 1: Linear model fitting](https://learn.juliaastro.org/tutorials/1_models-linear/)
"""

# ╔═╡ 7f2fc42d-c160-46f9-988d-78b102c7427c
md"""
## Summary

In this tutorial, we will build Bayesian models for the Type II Cepheid and synthetic data sets from Modeling 1, and use them to estimate the full uncertainty in our fitted parameters.
"""

# ╔═╡ 6c2b4069-7015-400e-8df1-201c90bace47
md"""
### Packages 📦
"""

# ╔═╡ e4d65ee5-d7d6-4970-a4bb-12b4df8823d6
md"""
## Introduction

In [Modeling 1](https://learn.juliaastro.org/tutorials/1_models-linear/), we found a single set of best-fit parameters for each of our models. GLM.jl was able to give us confidence intervals for our linear model, but the approximations that it relies on do not hold in general. A Bayesian approach instead gives us the full probability distribution of our model parameters given our data, i.e., the _posterior_.

We will build our models with [Turing.jl](https://turinglang.org/docs/). Instead of writing an objective function, we use the `@model` macro to describe how our data are generated. Each `~` statement either assigns a parameter its prior distribution, or states how our observations are distributed given these parameters (i.e., the likelihood).
"""

# ╔═╡ 789d59f3-eafe-47f4-ae29-a6db10256132
md"""
## Linear model

Let's start with the Type II Cepheid catalog from [Bhardwaj et al. 2017](https://ui.adsabs.harvard.edu/abs/2017A%26A...605A.100B). As in Modeling 1, we use [VirtualObservatory.jl](https://github.com/JuliaAPlavin/VirtualObservatory.jl/tree/master) to download it from Vizier, and then select the periods and Ks-band magnitudes, along with their uncertainties:
"""

# ╔═╡ 2c2faec3-6ddf-44f9-b83b-d3d7bebf3cfc
catalog = VizierCatalog("J/A+A/605/A100"; unitful = false) |> table |> DataFrame |> dropmissing

# ╔═╡ b4149e12-1658-41e6-9d8a-bd67d0250931
df = @select catalog begin
    :log_P = log10.(:Period)
    :Ks = :"<Ksmag>"
    :Ks_err = :"e_<Ksmag>"
end

# ╔═╡ d6215915-fe79-47a1-9582-5a27b89079a9
md"""
Before building our model, it is worth checking how well the weighted least squares fit from Modeling 1 describes these data. We compute it in the same way:
"""

# ╔═╡ 3fe5b684-c32b-4cce-b290-26621cf2c072
X = [ones(length(df.log_P)) df.log_P]  # Design matrix

# ╔═╡ 2456daf2-e526-4b4e-bc1e-cacaaa35f239
W = Diagonal(1 ./ df.Ks_err .^ 2) # Weight matrix

# ╔═╡ d3118358-570a-460f-bd80-2efd5cdbb8b4
β̂_wls = (X' * W * X) \ (X' * W * df.Ks)

# ╔═╡ 163e5ee9-ee91-4b4c-8ab9-74f8048e6fbd
md"""
If our measurement uncertainties accounted for all of the scatter around this line, the reduced ``\chi^2`` would be close to 1:
"""

# ╔═╡ f012fa3f-6314-46cf-91c6-7d5b6e13446d
χ²_red = sum(((df.Ks .- X * β̂_wls) ./ df.Ks_err) .^ 2) / (length(df.Ks) - 2) # Typed \chi<TAB>\^2<TAB>

# ╔═╡ 66756a16-8d51-4dcd-916d-9187fe9d4790
md"""
It is much larger than 1, so our measurement uncertainties alone underestimate the scatter around our line. This is not too surprising, since the period-luminosity relation has some intrinsic spread of its own. We can account for this by adding an intrinsic scatter parameter to our model, ``\sigma_\text{int}``, which adds in quadrature to our measurement uncertainties, ``\sigma_i``:

```math
K_{s, i} \sim \mathcal{N}\left(b + m \log_{10} P_i,\ \sqrt{\sigma_i^2 + \sigma_\text{int}^2}\right)\ .
```
"""

# ╔═╡ e7c5e695-c522-4796-90f8-6c1637e3f496
@model function linear_regression(x, y, σ_y)
    # Priors
    b ~ Normal(0, 20) # y-intercept
    m ~ Normal(0, 10) # Slope
    σ_int ~ truncated(Normal(0, 1); lower = 0) # Intrinsic scatter, must be positive

    # Likelihood
    for i in eachindex(x, y, σ_y)
        y[i] ~ Normal(m * x[i] + b, sqrt(σ_y[i]^2 + σ_int^2))
    end
end

# ╔═╡ 55b0b5ac-c75d-496f-b63c-4d4d41d9a689
md"""
We can now draw samples from our posterior using one of the many samplers available in Turing.jl. Here, we use the No U-Turn Sampler (`NUTS`), and pass a seeded random number generator (`Random.Xoshiro`) so that we draw the same samples each time:
"""

# ╔═╡ dc585840-1c4f-4ca1-95d1-0b1f95d7bbd1
chain_lin = sample(
    Xoshiro(1234),
    linear_regression(df.log_P, df.Ks, df.Ks_err),
    NUTS(),
    2_000;
    progress = false,
)

# ╔═╡ 0a6caf7b-5ec7-4d36-a79f-09cb935c5c7f
md"""
Let's summarize our samples:
"""

# ╔═╡ 4998ab29-b837-4316-adc9-fbe432e2bd04
summarystats(chain_lin)

# ╔═╡ bda1c16c-b788-4c21-9800-343e2c5ca8c4
md"""
Our y-intercept and slope have shifted a bit from our weighted fit, `β̂_wls`. With the intrinsic scatter now dominating over most of our measurement uncertainties, our most precise measurements no longer carry as much weight as before. Let's visualize our results by plotting the line given by each of our posterior samples:
"""

# ╔═╡ 379d0268-642a-4251-942d-cb29790461f4
let
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        title = "Type II Cepheid observations",
        titlesize = 16,
        titlealign = :left,
        subtitle = "Bhardwaj et al. 2017",
        xlabel = L"\mathbf{\log_{10}(\text{Period [days]})}",
        ylabel = "Ks [mag]",
        ylabelfont = :bold,
    )

    # Posterior samples
    ablines!(ax, vec(chain_lin[:b]), vec(chain_lin[:m]); color = (:orange, 0.02))

    # Data
    errorbars!(ax, df.log_P, df.Ks, df.Ks_err; color = :cornflowerblue)
    scatter!(ax, df.log_P, df.Ks; color = :cornflowerblue)

    fig
end

# ╔═╡ 7b681716-9367-49a2-a45a-ecedd7539263
md"""
To examine our model parameters and how they correlate with each other in more detail, we can make a corner plot with [PairPlots.jl](https://sefffal.github.io/PairPlots.jl/dev/):
"""

# ╔═╡ 6d0d38a6-2aca-455a-b36e-e3d7d86cb023
pairplot(chain_lin)

# ╔═╡ 210ec688-2d9e-430c-b226-deea36cee5f7
md"""
Note the strong anticorrelation between our y-intercept and slope. Since our `log_P` values are not centered on zero, a steeper line must also have a higher y-intercept to keep passing through the data.
"""

# ╔═╡ 9608e9c6-49e3-425f-b2b9-b8d190782efb
md"""
## Polynomial model

Next, let's revisit the synthetic data from Modeling 1, which follow a weak parabola with some Gaussian noise added on top. We generate them in the same way, with the same seed:
"""

# ╔═╡ d050bd86-1fd2-42b8-87d6-ff214f2ef125
df_sim = let rng = Xoshiro(1234)
    x = 0.0:5:100
    y = @. (x / 20 - 0.2)^2 + 2
    DataFrame(; x, y = y .+ randn(rng, length(x)))
end

# ╔═╡ 746cec83-3057-47b4-87bb-eb25a5db9135
md"""
We also define the same quadratic model, ``y = \beta_1 + \beta_2 x + \beta_3 x^2``:
"""

# ╔═╡ 68c3f123-ba2b-41c2-9101-63e4611b37a2
@. quadfunc(x, β) = β[1] + β[2] * x + β[3] * x^2

# ╔═╡ ebbcbf01-caa8-49de-8c09-a972b00a0c35
md"""
We do not have any measurement uncertainties this time, so we will fit for the standard deviation of the noise, ``\sigma``, as well:
"""

# ╔═╡ 45c2b4e6-4e3f-4514-b0f4-a3b45223adc0
@model function quad_regression(x, y)
    # Priors, typed \beta<TAB>\_1<TAB>, etc.
    β₁ ~ Normal(0, 5)
    β₂ ~ Normal(0, 0.1)
    β₃ ~ Normal(0, 0.01)
    σ ~ truncated(Normal(0, 10); lower = 0)

    # Likelihood
    for i in eachindex(x, y)
        y[i] ~ Normal(quadfunc(x[i], (β₁, β₂, β₃)), σ)
    end
end

# ╔═╡ 7e0cf119-36f4-44f6-a462-6e3d3cb2a852
md"""
Our polynomial coefficients are strongly correlated with each other, which makes the posterior a bit harder for `NUTS` to explore. We help it out by raising its target acceptance rate from the default of 0.65 to 0.8, which makes it take smaller, more careful steps:
"""

# ╔═╡ bc50ce8f-4f12-4fa5-b4e4-3ee011c21dda
chain_quad = sample(
    Xoshiro(1234),
    quad_regression(df_sim.x, df_sim.y),
    NUTS(0.8),
    2_000;
    progress = false,
)

# ╔═╡ a016b7f2-7927-4e7b-a52a-612dfd7768c3
let
    β₁, β₂, β₃ = vec(chain_quad[:β₁]), vec(chain_quad[:β₂]), vec(chain_quad[:β₃])

    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "y")

    # Every 10th posterior sample
    for i in firstindex(β₁):10:lastindex(β₁)
        lines!(ax, df_sim.x, quadfunc(df_sim.x, (β₁[i], β₂[i], β₃[i])); color = (:orange, 0.1))
    end

    # Data
    scatter!(ax, df_sim.x, df_sim.y; color = :cornflowerblue)

    fig
end

# ╔═╡ 0d86486f-6634-4e11-9a37-3e86727e8d40
pairplot(chain_quad)

# ╔═╡ 3286a9ed-9a71-4319-b1f1-20a0e20aaad2
md"""
For comparison, the expression we used to generate our data expands to ``(x/20 - 0.2)^2 + 2 = 2.04 - 0.02x + 0.0025x^2``, and our noise had a standard deviation of 1. These true values all fall within the bulk of our posterior.

!!! todo
    Mark these true values directly on the corner plot with `PairPlots.Truth` once it supports FlexiChains.jl chains: <https://github.com/sefffal/PairPlots.jl/issues/79>

!!! tip
    See the [Turing.jl documentation](https://turinglang.org/docs/) for many more examples of Bayesian models, and the [PairPlots.jl documentation](https://sefffal.github.io/PairPlots.jl/dev/) for ways to customize these corner plots.
"""

# ╔═╡ 8ba6970b-1e1b-49a2-93bc-34d82a48fa37
md"""
# Notebook setup 🔧
"""

# ╔═╡ a3c1e7d2-5b84-4f0e-9c6a-2d7f18e4b590
TableOfContents(; depth = 4)

# ╔═╡ Cell order:
# ╟─d03f494f-4a91-4ff8-b5c6-cffafe540fcc
# ╟─7f2fc42d-c160-46f9-988d-78b102c7427c
# ╟─6c2b4069-7015-400e-8df1-201c90bace47
# ╠═c4c69a9c-0269-4022-b7cb-2c04e35c3483
# ╟─e4d65ee5-d7d6-4970-a4bb-12b4df8823d6
# ╟─789d59f3-eafe-47f4-ae29-a6db10256132
# ╠═2c2faec3-6ddf-44f9-b83b-d3d7bebf3cfc
# ╠═b4149e12-1658-41e6-9d8a-bd67d0250931
# ╟─d6215915-fe79-47a1-9582-5a27b89079a9
# ╠═3fe5b684-c32b-4cce-b290-26621cf2c072
# ╠═2456daf2-e526-4b4e-bc1e-cacaaa35f239
# ╠═d3118358-570a-460f-bd80-2efd5cdbb8b4
# ╟─163e5ee9-ee91-4b4c-8ab9-74f8048e6fbd
# ╠═f012fa3f-6314-46cf-91c6-7d5b6e13446d
# ╟─66756a16-8d51-4dcd-916d-9187fe9d4790
# ╠═e7c5e695-c522-4796-90f8-6c1637e3f496
# ╟─55b0b5ac-c75d-496f-b63c-4d4d41d9a689
# ╠═dc585840-1c4f-4ca1-95d1-0b1f95d7bbd1
# ╟─0a6caf7b-5ec7-4d36-a79f-09cb935c5c7f
# ╠═4998ab29-b837-4316-adc9-fbe432e2bd04
# ╟─bda1c16c-b788-4c21-9800-343e2c5ca8c4
# ╠═379d0268-642a-4251-942d-cb29790461f4
# ╟─7b681716-9367-49a2-a45a-ecedd7539263
# ╠═6d0d38a6-2aca-455a-b36e-e3d7d86cb023
# ╟─210ec688-2d9e-430c-b226-deea36cee5f7
# ╟─9608e9c6-49e3-425f-b2b9-b8d190782efb
# ╠═d050bd86-1fd2-42b8-87d6-ff214f2ef125
# ╟─746cec83-3057-47b4-87bb-eb25a5db9135
# ╠═68c3f123-ba2b-41c2-9101-63e4611b37a2
# ╟─ebbcbf01-caa8-49de-8c09-a972b00a0c35
# ╠═45c2b4e6-4e3f-4514-b0f4-a3b45223adc0
# ╟─7e0cf119-36f4-44f6-a462-6e3d3cb2a852
# ╠═bc50ce8f-4f12-4fa5-b4e4-3ee011c21dda
# ╠═a016b7f2-7927-4e7b-a52a-612dfd7768c3
# ╠═0d86486f-6634-4e11-9a37-3e86727e8d40
# ╟─3286a9ed-9a71-4319-b1f1-20a0e20aaad2
# ╟─8ba6970b-1e1b-49a2-93bc-34d82a48fa37
# ╠═a3c1e7d2-5b84-4f0e-9c6a-2d7f18e4b590
