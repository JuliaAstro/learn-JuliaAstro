### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> image = "/assets/linear-model-fitting.png"
#> title = "Modeling 1: Linear model fitting"
#> layout = "layout.jlhtml"
#> date = "2025-12-09"
#> description = "Make simple fits to your data."
#> tags = ["models", "model fitting", "astrostatistics", "catalog", "query", "Makie", "plots", "errorbars", "scatter plots"]

using Markdown
using InteractiveUtils

# ╔═╡ b944163c-02c4-4ace-a866-ae6e6f7115ef
begin
    # Can remove this block after AoG #710 is merged:
    # https://github.com/MakieOrg/AlgebraOfGraphics.jl/pull/710
    import Pkg
    Pkg.activate(; temp = true)
    Pkg.add(
        [
            Pkg.PackageSpec(; name = "TOML"),
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(; name = "StatsBase"),
            Pkg.PackageSpec(; name = "LinearAlgebra"),
            Pkg.PackageSpec(; name = "GLM"),
            Pkg.PackageSpec(; name = "Optimization"),
            Pkg.PackageSpec(; name = "OptimizationOptimJL"),
            Pkg.PackageSpec(; name = "DataFramesMeta"),
            Pkg.PackageSpec(; name = "VirtualObservatory"),
            Pkg.PackageSpec(; name = "CairoMakie"),
            Pkg.PackageSpec(; name = "MathTeXEngine"),
            Pkg.PackageSpec(;
                url = "https://github.com/icweaver/AlgebraOfGraphics.jl",
                rev = "glm",
            ),
        ]
    )

    # Analysis
    using StatsBase: coef, predict
    using LinearAlgebra: Diagonal
    using GLM: Normal, @formula, glm, fweights
    using Optimization: OptimizationProblem, solve
    using OptimizationOptimJL: NelderMead

    # Data handling and visualization
    using DataFramesMeta: @select, @transform, DataFrame, Not, disallowmissing, dropmissing
    using VirtualObservatory: VizierCatalog, table
    using CairoMakie
    using AlgebraOfGraphics: aog_theme, data, linear, draw, mapping, visual
    using MathTeXEngine: set_texfont_family!, FontFamily
    set_texfont_family!(FontFamily("TeXGyreHeros"))
end;

# ╔═╡ b880641e-5101-4857-a5bf-558482ca1b21
begin
    using TOML: TOML
    using PlutoUI: TableOfContents
end

# ╔═╡ 101f85ed-9442-4b95-a771-f7516e6d84cb
md"""
## Summary

In this tutorial, we will become familiar with some of the major modeling frameworks available in Julia and learn how to make a quick fit to our data.
"""

# ╔═╡ a68daf19-dae6-4442-8234-d6636ef931c7
md"""
### Packages 📦
"""

# ╔═╡ 6096dc80-a434-4f1a-9af5-40fbbf897d05
md"""
## Data

We are going to start with a linear fit to real data. The data comes from the paper [Bhardwaj et al. 2017](https://ui.adsabs.harvard.edu/abs/2017A%26A...605A.100B). This is a catalog of Type II Cepheids, which is a type of variable stars that pulsate with a period between 1 and 50 days. In this part of the tutorial, we are going to measure the Cepheids Period-Luminosity relation using astropy.modeling. This relation states that if a star has a longer period, the luminosity we measure is higher. We use [VirtualObservatory.jl](https://github.com/JuliaAPlavin/VirtualObservatory.jl/tree/master) to download this data from Vizier:
"""

# ╔═╡ 002fa353-906a-4e3c-bd8e-12d681417240
catalog = VizierCatalog("J/A+A/605/A100"; unitful = false) |> table |> DataFrame |> dropmissing

# ╔═╡ b406eba1-cac7-4d19-854d-f2b2c0eda33f
names(catalog) # List all columns names

# ╔═╡ 05047cc3-b7c7-40eb-86da-fcc6698e41ee
md"""
This catalog has a lot of information, but for this tutorial we are going to work only with periods and magnitudes. Let's grab them using the keywords `Period` and `<Ksmag>`. Note that `e_<Ksmag>` refers to the error bars in the magnitude measurements. We'd also like to rename these to some more convenient labels. We can do this selection and renaming simultaneously with the `DataFramesMeta.@select` macro and `=>` syntax from DataFrames.jl, respectively:
"""

# ╔═╡ 3bfaf85e-7cd8-4fa8-be0c-34afeae5ed06
df = @select catalog begin
    :log_P = log10.(:Period)
    :Ks = :"<Ksmag>"
    :Ks_err = :"e_<Ksmag>"
end

# ╔═╡ ddce8426-9da5-419e-b5f4-65aaf990c7ca
md"""
Let's take a look at the magnitude measurements as a function of period. We'll show a convenient way to do this with AlgebraOfGraphics.jl (AoG.jl), as well as with manual methods if more control is desired.
"""

# ╔═╡ 2132f3ca-cbc2-4b5c-9c45-850cdb14b3b8
md"""
## Automatic plotting with AoG.jl

Starting with AoG.jl, here is a simple one-liner we can start with:
"""

# ╔═╡ a8ae7056-68a4-4bfc-9adb-e52821a1ae49
data(df) * mapping(:log_P, :Ks) * (visual(Scatter) + linear()) |> draw

# ╔═╡ f6ea7c5b-ad96-4f79-a4d4-5c6a1e79fbcb
md"""
A fair bit happened here. In the above line, AoG.jl:

- Made a scatter plot of `log_P` vs. `Ks`.
- Fit a line to the data.
- Plot the fitted line along with its estimated 95% confidence interval.
- Labeled the axes with the appropriate column names used.

While this is convenient for quick visualization to see that there indeed appears to be a linear relationship between the log period of the pulsation period and luminosity (inverse relation to observed magnitude), we really would like to take special care with our statistical analysis.

For example, we would like to weight our fit by the uncertainty in our magnitude measurements, `Ks_err`. It is also important to note that there is a difference between frequency weights, which this package uses by default, and analytic weights, which go like the inverse variance of our measurements ``\left(1 / \sigma_i^2\right)``, where ``\sigma_i \equiv `` `Ks_err` for our purposes. This will also impact how our confidence interval is calculated.

Let's apply these requirements, and also update the styling of our plot a bit:
"""

# ╔═╡ 1d61727f-b373-4394-b486-79eb63a71f37
with_theme(Theme(aog_theme())) do
    # Analytic (inverse-variance) weights, normalized to sum to n so the
    # confidence band uses the true sample size rather than Σw
    w = inv.(df.Ks_err .^ 2)
    df_w = @transform df :weights = w .* (length(w) / sum(w))

    # Common data
    layer_scatter = data(df_w) * mapping(
        :log_P => L"\mathbf{\log_{10}(\text{Period [days]})}",
        :Ks => "Ks [mag]",
    )

    # Errorbars
    layer_errorbars = layer_scatter * mapping(:Ks_err) * visual(Errorbars)

    # Linear model
    layer_model = layer_scatter *
        mapping(weights = :weights) *
        linear()

    # Combined layers
    layer_data = layer_scatter + layer_errorbars
    layers = layer_data * visual(label = "data", color = :cornflowerblue) +
        layer_model * visual(label = "model", color = "orange")

    # Display
    fig = draw(
        layers;
        figure = (;
            title = "Type II Cepheid observations",
            subtitle = "Bhardwaj et al. 2017",
        ),
    )
end

# ╔═╡ 9773d632-f5cd-47d5-b97e-57a7b6ca3bf9
md"""
!!! warning
    AlgebraOfGraphics.jl currently interprets the `weights` column as frequency weights; `linear(; weighttype = :aweights)` will become available once GLM.jl v2 is released. Until then, normalizing inverse-variance weights to sum to the number of data points (as done above) gives results numerically identical to a proper analytic-weights fit.

!!! tip
    Themes can also be set globally with:

    ```julia
    using CairoMakie

    set_theme!(<theme>)
    ```

    See [this section](https://docs.makie.org/stable/explanations/theming/themes) of the Makie.jl documentation for more on theming, and [this tutorial](https://aog.makie.org/stable/tutorials/intro-i) for more on getting started with AoG.
"""

# ╔═╡ 6ef287bb-7b3f-4fbc-aa24-2b317a22b6f9
md"""
Much better! Besides some minor styling instructions given from our end, the entirety of the statistical analysis was handled by this small, but powerful, bit of code:

```julia
# Analytic (inverse-variance) weights, normalized to sum to n
w = inv.(df.Ks_err .^ 2)
df_w = @transform df :weights = w .* (length(w) / sum(w))

# Linear model
layer_model = layer_scatter *
    mapping(weights = :weights) *
    linear()
```

We'll next take a look under the hood to see how these calculations were performed.
"""

# ╔═╡ 4e69c1de-6cf0-448b-8b95-5b6c54795820
md"""
!!! note
    For more on different statistical weightings used in the Julia ecosystem, see [this section](https://juliastats.org/StatsBase.jl/stable/weights/#AnalyticWeights) of the StatsBase.jl documentation.
"""

# ╔═╡ 35f8675d-d8f0-47b0-83b6-156bb75d373d
md"""
## Manual methods

At its core, we are essentially solving the following [linear algebra equation](https://en.wikipedia.org/wiki/Weighted_least_squares):

```math
\mathbf{
    \left( X^\textsf{T} W X \right) \boldsymbol{\hat\beta} =
    X^\textsf{T} W y\
}\ ,
```

where ``\mathbf X`` is our design matrix, ``\mathbf W`` our weight matrix, and ``\boldsymbol{\hat\beta}`` our parameter vector (i.e., the y-intercept and slope of the line we would like to fit). Julia has a very powerful [matrix division operator (`\`)](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#Base.:\\-Tuple{AbstractMatrix,%20AbstractVecOrMat}) built right into the language, which we can use to solve this equation for ``\boldsymbol{\hat\beta}``:
"""

# ╔═╡ bd068a95-2945-4eeb-b2e9-869c03f79e99
md"""
### `Base` Julia
"""

# ╔═╡ e1ad0f97-2d61-43dc-825c-bf7503a708dd
X = [ones(length(df.log_P)) df.log_P]  # Design matrix

# ╔═╡ af209cbf-1144-42cc-8c52-029988b0d5e4
W = Diagonal(1 ./ df.Ks_err .^ 2) # Weight matrix

# ╔═╡ 616d9778-106e-46ae-bf87-92828e98339f
β̂_base = (X' * W * X) \ (X' * W * df.Ks)

# ╔═╡ e6d108e2-6ea1-4875-a3db-239248100bc7
md"""
We now have the y-intercept and slope for our weighted linear model, all in base Julia! We can next conveniently compute the confidence intervals with GLM.jl package, which we show next.
"""

# ╔═╡ d8e6db3b-f1c3-46f5-b092-1ef25d762cdb
md"""
### GLM.jl

[GLM.jl](https://juliastats.org/GLM.jl/dev/) is the linear and generalized linear models package used by AoG.jl to perform model fitting and uncertainty estimation. It is invoked via the `GLM.glm` function, which can be passed our linear formula via the `GLM.@formula` macro, a distribution that our uncertainties are sampled from, and the kind of weights that we are using, e.g., probability, analytic, frequency, etc.
"""

# ╔═╡ 2cf30811-3204-4498-9b73-8a3c3bba28e2
df_glm = @transform df begin
    :log_P = Float64.(:log_P)
    :Ks = Float64.(:Ks)
    :Ks_err = Float64.(:Ks_err)
end;

# ╔═╡ 19df3403-ef9d-4b07-8b14-3b10056475e8
md"""
!!! todo
    Current workaround for <https://github.com/JuliaStats/GLM.jl/issues/260>
"""

# ╔═╡ 4a17877a-6a38-4364-9bd9-91dde011a42a
fit_glm = let w = inv.(df_glm.Ks_err .^ 2)
    # Normalize the inverse-variance weights to sum to n (same as the AoG cell)
    # Can replace with weights = aweights(w) in GLM v2.
    glm(@formula(Ks ~ log_P), df_glm, Normal(); weights = fweights(w .* (length(w) / sum(w))))
end

# ╔═╡ 816426f9-0511-4dde-9ec4-7c9ccb212a45
md"""
!!! note
    See the [GLM.jl documentation](https://juliastats.org/GLM.jl/stable/) for more.
"""

# ╔═╡ 3b7eee06-f730-45d5-aae1-53b7731694d9
md"""
Note that the computed y-intercept and slope (which we can extract with `GLM.coef` for convenience) are quite close to what we computed manually:
"""

# ╔═╡ 8c44bde9-c64a-41c2-8b47-e98d40b3cf42
β̂_glm = coef(fit_glm)

# ╔═╡ 0c73d53c-8b6e-4b13-a253-27046c1bd741
β̂_glm - β̂_base

# ╔═╡ 8641e496-b5f8-4c21-b456-77f14771d0b1
md"""
and we now have associated uncertainty information from `fit_glm` that we can use to estimate our confidence intervals using `GLM.predict`:
"""

# ╔═╡ 851acc22-bf1e-409e-8037-c5c37022fa63
eachcol(df_glm[!, [:log_P]]) |> collect

# ╔═╡ 07193622-e18b-474d-ab54-95cc09cf3bfc
Ks_pred, Ks_lower, Ks_upper = eachcol(predict(fit_glm, df_glm[!, [:log_P]]; interval = :confidence))

# ╔═╡ d03883f3-b488-4176-bb5f-d69dcb9fad77
md"""
Note that `Ks_pred` is the same as computing ``y = \hat β_1 + \hatβ_2x`` ourselves:
"""

# ╔═╡ 7d4b3251-1ce1-43a6-933f-928dd3bf281a
Ks_pred == β̂_glm[1] .+ β̂_glm[2] * df.log_P

# ╔═╡ 38bc6713-6e89-47c0-9e03-b3e90e3182be
md"""
### Optimization.jl

For completeness, we also show how we might accomplish this with Optimization.jl:

!!! note
    This is best suited for nonlinear problems, where the usual linear approximations for estimating confidence intervals we used before [do not hold](https://discourse.julialang.org/t/best-fit-parameter-error-bar-using-optimization-jl/103186/6). At this point, our standard confidence interval estimatation techniques above do not hold, and Bayesian approaches should be used instead. For example, see: <https://juliaastro.org/home/tutorials/curve-fit/#Bayesian-models>
"""

# ╔═╡ e0e8909d-893f-4f35-baf2-c27dafeb23fa
function objective(u, data)
    b, m = u
    x, y, y_err = eachcol(data)

    # Compute weighted residuals
    residuals = @. (y - (m * x + b)) / y_err

    # Sum of squared weighted residuals (χ²)
    return sum(residuals .^ 2)
end

# ╔═╡ a3582cb7-c2a6-414c-bf7b-d764c2311c80
u0 = zeros(2)

# ╔═╡ 3bffbcde-b2b1-445d-b2f7-59a53f33f90a
prob = OptimizationProblem(objective, u0, df)

# ╔═╡ 925eb5ac-a5fa-40d9-8ab1-f10db7dab202
β̂_optim = solve(prob, NelderMead())

# ╔═╡ 4a74192c-c78e-4cea-ac77-7b8bdd486267
md"""
Again, our estimated y-intercept and slope are quite close to our other manual estimates:
"""

# ╔═╡ edf0fc17-b5a9-4475-bdd7-aa6f2673a39d
β̂_glm .- β̂_optim

# ╔═╡ 983cf394-29a5-4363-b07b-ccc91c5475e2
β̂_optim .- β̂_base

# ╔═╡ baeba190-dc5f-42ae-a6ec-516aa49bf612
md"""
We can now use either estimate to produce our linear fit plot in plain Makie.
"""

# ╔═╡ 0f57e294-1c1f-4f73-8e08-ab9c0bba9f0f
md"""
### Plot
"""

# ╔═╡ c3d88d47-93f5-4f39-98ea-eb76f8a4974d
let
    # with_theme(Theme(aog_theme())) do
    log_P, Ks, Ks_err = df.log_P, df.Ks, df.Ks_err

    # Data points
    fig, ax, p = scatter(
        log_P, Ks;
        color = :cornflowerblue,
        label = "data",
        axis = (;
            title = "Type II Cepheid observations",
            titlesize = 16,
            titlealign = :left,
            subtitle = "Bhardwaj et al. 2017",
            xlabel = L"\mathbf{\log_{10}(\text{Period [days]})}",
            ylabel = "Ks [mag]",
            ylabelfont = :bold,
        ),
    )

    # Data uncertainty
    errorbars!(ax, log_P, Ks, Ks_err; color = :cornflowerblue, label = "data")

    # Confidence interval
    band!(
        ax, log_P, disallowmissing(Ks_lower), disallowmissing(Ks_upper);
        color = :orange,
        alpha = 0.15,
        label = "model",
    )

    # Model prediction
    lines!(ax, log_P, Ks_pred; color = :orange, label = "model")

    # Legend
    Legend(fig[1, 2], ax; merge = true)

    fig
end

# ╔═╡ 3ba10da3-1e3c-4b75-9c0c-5d1a2dd4af75
md"""
!!! note
    Note the repeated labeling, styling, and external statistical analysis that AoG saves us from needing to do by hand.
"""

# ╔═╡ b2805e96-5cce-4200-842b-931187007a31
md"""
# Notebook setup 🔧
"""

# ╔═╡ 7cb7ae25-79c2-4138-aa50-fdc27615245b
TableOfContents(; depth = 4)

# ╔═╡ fac4fe53-010d-44e5-956b-76bb2a530011
function frontmatter(path)
    prefix = "#> "
    is_fm = startswith(prefix)
    block = Iterators.takewhile(is_fm, Iterators.dropwhile(!is_fm, eachline(path)))
    toml = TOML.parse(join(chopprefix.(block, prefix), "\n"))
    return toml["frontmatter"]
end

# ╔═╡ 2b57e162-b642-4dfe-88e7-45a6bb3f8447
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ ec1a7344-e375-4847-b4f7-765a53c066d0
md"""
# Modeling 1: Make a quick linear model fit

This notebook is modified from <https://learn.astropy.org/tutorials/1_models-quick-fit.html>

_Original authors: Rocio Kiman, Lia Corrales, Zé Vinícius, Kelle Cruz, Stephanie T. Douglas_

!!! tip "Learning goals"
    - Use VirtualObservatory.jl to download data from Vizier.
    - Use basic models in `Base` Julia, GLM.jl, and Optimization.jl.
    - Learn common functions to fit.
    - Generate a quick fit to data.
    - Plot the model with the data.
    - Compare different models and fitters.

$(keywords())


!!! warning "Companion content"
    Content here.
"""

# ╔═╡ Cell order:
# ╟─ec1a7344-e375-4847-b4f7-765a53c066d0
# ╟─101f85ed-9442-4b95-a771-f7516e6d84cb
# ╟─a68daf19-dae6-4442-8234-d6636ef931c7
# ╠═b944163c-02c4-4ace-a866-ae6e6f7115ef
# ╟─6096dc80-a434-4f1a-9af5-40fbbf897d05
# ╠═002fa353-906a-4e3c-bd8e-12d681417240
# ╠═b406eba1-cac7-4d19-854d-f2b2c0eda33f
# ╟─05047cc3-b7c7-40eb-86da-fcc6698e41ee
# ╠═3bfaf85e-7cd8-4fa8-be0c-34afeae5ed06
# ╟─ddce8426-9da5-419e-b5f4-65aaf990c7ca
# ╟─2132f3ca-cbc2-4b5c-9c45-850cdb14b3b8
# ╠═a8ae7056-68a4-4bfc-9adb-e52821a1ae49
# ╟─f6ea7c5b-ad96-4f79-a4d4-5c6a1e79fbcb
# ╠═1d61727f-b373-4394-b486-79eb63a71f37
# ╟─9773d632-f5cd-47d5-b97e-57a7b6ca3bf9
# ╟─6ef287bb-7b3f-4fbc-aa24-2b317a22b6f9
# ╟─4e69c1de-6cf0-448b-8b95-5b6c54795820
# ╟─35f8675d-d8f0-47b0-83b6-156bb75d373d
# ╟─bd068a95-2945-4eeb-b2e9-869c03f79e99
# ╠═e1ad0f97-2d61-43dc-825c-bf7503a708dd
# ╠═af209cbf-1144-42cc-8c52-029988b0d5e4
# ╠═616d9778-106e-46ae-bf87-92828e98339f
# ╟─e6d108e2-6ea1-4875-a3db-239248100bc7
# ╟─d8e6db3b-f1c3-46f5-b092-1ef25d762cdb
# ╠═2cf30811-3204-4498-9b73-8a3c3bba28e2
# ╟─19df3403-ef9d-4b07-8b14-3b10056475e8
# ╠═4a17877a-6a38-4364-9bd9-91dde011a42a
# ╟─816426f9-0511-4dde-9ec4-7c9ccb212a45
# ╟─3b7eee06-f730-45d5-aae1-53b7731694d9
# ╠═8c44bde9-c64a-41c2-8b47-e98d40b3cf42
# ╠═0c73d53c-8b6e-4b13-a253-27046c1bd741
# ╟─8641e496-b5f8-4c21-b456-77f14771d0b1
# ╠═851acc22-bf1e-409e-8037-c5c37022fa63
# ╠═07193622-e18b-474d-ab54-95cc09cf3bfc
# ╟─d03883f3-b488-4176-bb5f-d69dcb9fad77
# ╠═7d4b3251-1ce1-43a6-933f-928dd3bf281a
# ╟─38bc6713-6e89-47c0-9e03-b3e90e3182be
# ╠═e0e8909d-893f-4f35-baf2-c27dafeb23fa
# ╠═a3582cb7-c2a6-414c-bf7b-d764c2311c80
# ╠═3bffbcde-b2b1-445d-b2f7-59a53f33f90a
# ╠═925eb5ac-a5fa-40d9-8ab1-f10db7dab202
# ╟─4a74192c-c78e-4cea-ac77-7b8bdd486267
# ╠═edf0fc17-b5a9-4475-bdd7-aa6f2673a39d
# ╠═983cf394-29a5-4363-b07b-ccc91c5475e2
# ╟─baeba190-dc5f-42ae-a6ec-516aa49bf612
# ╟─0f57e294-1c1f-4f73-8e08-ab9c0bba9f0f
# ╠═c3d88d47-93f5-4f39-98ea-eb76f8a4974d
# ╟─3ba10da3-1e3c-4b75-9c0c-5d1a2dd4af75
# ╟─b2805e96-5cce-4200-842b-931187007a31
# ╠═7cb7ae25-79c2-4138-aa50-fdc27615245b
# ╟─fac4fe53-010d-44e5-956b-76bb2a530011
# ╟─2b57e162-b642-4dfe-88e7-45a6bb3f8447
# ╠═b880641e-5101-4857-a5bf-558482ca1b21
