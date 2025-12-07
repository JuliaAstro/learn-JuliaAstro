### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ b944163c-02c4-4ace-a866-ae6e6f7115ef
begin
	import Pkg
	Pkg.activate(Base.current_project())
	# Pkg.add(; url = "https://github.com/gragusa/GLM.jl", rev = "JuliaStats-master")
	# Pkg.add(; name = "VirtualObservatory", version = "0.1.14")
	# Pkg.add(["CairoMakie", "InvertedIndices", "DataFramesMeta", "StatsBase", "AlgebraOfGraphics", "LsqFit"])
	using Revise

	using VirtualObservatory: VizierCatalog, table
	using CairoMakie: Errorbars, errorbars, lines!, scatter!
	using InvertedIndices: Not
	using DataFramesMeta
	using StatsBase: coef, predict
	using CairoMakie: Scatter, set_theme!, with_theme, Theme
	using AlgebraOfGraphics: aog_theme, data, linear, draw, mapping, visual
	using LinearAlgebra: Diagonal
	using GLM: GLM, Normal, @formula, glm, aweights
	using Optimization: OptimizationProblem, solve
	using OptimizationOptimJL: NelderMead
end

# ╔═╡ 434b1c97-4bae-49da-8eaf-207a23a6af7d
md"""
# Modeling 1: Make a quick model fit
"""

# ╔═╡ 6096dc80-a434-4f1a-9af5-40fbbf897d05
md"""
## Data
"""

# ╔═╡ 002fa353-906a-4e3c-bd8e-12d681417240
catalog = VizierCatalog("J/A+A/605/A100") |> table |> DataFrame |> dropmissing

# ╔═╡ 3bfaf85e-7cd8-4fa8-be0c-34afeae5ed06
df = @select catalog begin
	:log_P = log10.(:Period)
	:Ks = :"<Ksmag>"
	:Ks_err = :"e_<Ksmag>"
end

# ╔═╡ 3f0abe3e-08cb-4c3c-9997-e654341bce12
md"""
## AoG.jl

Similar: See [this JuliaAstro tutorial](https://learn.juliaastro.org/tutorials/fits-images/#Plotting-with-Makie.jl-+-AoG.jl)
"""

# ╔═╡ 64da9892-c33b-4d00-ac8b-6564dea7b4e1
md"""
One-liner, unweighted (good for quick viz)
"""

# ╔═╡ a8ae7056-68a4-4bfc-9adb-e52821a1ae49
data(df) * mapping(:log_P, :Ks) * (visual(Scatter) + linear()) |> draw

# ╔═╡ 4e69c1de-6cf0-448b-8b95-5b6c54795820
md"""
!!! note
	Assumes frequency weights. Pass ``1 / \sigma_i^2`` to treat as analytic weights, match other libraries like astropy. For a discussion on some of the other weights used in the Julia ecosystem, see [this section](https://juliastats.org/StatsBase.jl/stable/weights/#AnalyticWeights) of the StatsBase.jl documentation.
"""

# ╔═╡ da58364b-78f3-4113-9fc1-1a020a9e15dc
md"""
Easy to customize and weight fit
"""

# ╔═╡ 1d61727f-b373-4394-b486-79eb63a71f37
with_theme(Theme(aog_theme())) do
	# Common data
	layer_scatter = data(df) * mapping(
		:log_P => "log₁₀(Period [days])",
		:Ks => "Ks [mag]",
	)
	
	# Errorbars
	layer_errorbars = layer_scatter * mapping(:Ks_err) * visual(Errorbars)

	# Linear model
	layer_model = layer_scatter *
		mapping(weights = :Ks_err => (x -> inv(x^2))) *
		linear()
	
	# Combined layers
	layer_data = layer_scatter + layer_errorbars
	layers = layer_data * visual(label="data", color=:cornflowerblue) +
		layer_model * visual(label="model", color="orange")

	# Display
	draw(layers;
		figure = (;
			title = "Type II Cepheid observations",
			subtitle = "Bhardwaj et al. 2017",
		),
	)
end

# ╔═╡ 15516d38-b5b0-4dc4-8ce5-8d853dfc2c3d
md"""
!!! tip "To-do"
	See if this can be upstreamed to AoG: <https://github.com/icweaver/AlgebraOfGraphics.jl/tree/glm>
"""

# ╔═╡ 9773d632-f5cd-47d5-b97e-57a7b6ca3bf9
md"""
!!! tip
	Themese can also be set globally with:

	```julia
	using CairoMakie
	
	set_theme!(<theme>)
	```

	See [this section](https://docs.makie.org/stable/explanations/theming/themes) of the Makie.jl documentation for more on themeing, and [this tutorial](https://aog.makie.org/stable/tutorials/intro-i) for more on getting started with AoG.
"""

# ╔═╡ 35f8675d-d8f0-47b0-83b6-156bb75d373d
md"""
## Manual method
"""

# ╔═╡ af209cbf-1144-42cc-8c52-029988b0d5e4
W = Diagonal(1 ./ df.Ks_err .^ 2)

# ╔═╡ e1ad0f97-2d61-43dc-825c-bf7503a708dd
X = [ones(length(df.log_P)) df.log_P]

# ╔═╡ 616d9778-106e-46ae-bf87-92828e98339f
sol_manual = m, b = (X' * W * X) \ (X' * W * df.Ks)

# ╔═╡ e6d108e2-6ea1-4875-a3db-239248100bc7
md"""
We now have the slope and y-intercept for our weighted linear model, all in base Julia!
"""

# ╔═╡ d8e6db3b-f1c3-46f5-b092-1ef25d762cdb
md"""
## GLM.jl

Used under-the-hood by AoG.jl
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
fit_glm = glm(@formula(Ks ~ log_P), df_glm, Normal(); wts = GLM.aweights(inv.(df_glm.Ks_err .^ 2)))

# ╔═╡ 1efa1d00-34da-48d1-90aa-c6b1d45a28ee
fit_glm2 = glm(@formula(Ks ~ log_P), df_glm, Normal(); wts = inv.(df_glm.Ks_err .^ 2))

# ╔═╡ 480d8ef0-42d7-4089-9da2-1543baa2d02b
md"""
!!! tip
	`@formula` is a convenience macro
"""

# ╔═╡ 07193622-e18b-474d-ab54-95cc09cf3bfc
predict(fit_glm, df_glm[!, [:log_P]]; interval = :confidence)

# ╔═╡ 32609df9-b89a-491f-b420-7847e8c177d6
md"""
!!! todo
	Rename wts to weights once this is in: <https://github.com/JuliaStats/GLM.jl/pull/570>
"""

# ╔═╡ 18c7a687-7c3a-4c0d-8928-50b388683ca1
md"""
!!! todo
	Confidence intervals still using frequency waits. See this PR for correctly handling analytic weights <https://github.com/JuliaStats/GLM.jl/pull/487>
"""

# ╔═╡ 8c44bde9-c64a-41c2-8b47-e98d40b3cf42
sol_glm = coef(fit_glm)

# ╔═╡ 38bc6713-6e89-47c0-9e03-b3e90e3182be
md"""
## Optimization.jl
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
sol_optim = solve(prob, NelderMead())

# ╔═╡ 2df6dced-3320-4d1a-a6c8-a46cb4c2bcf1
md"""
!!! note
	Best suited for nonlinear problems, where the usual linear approximations for estimating confidence intervals do not hold. Should use MCMC approaches instead. See <https://juliaastro.org/home/tutorials/curve-fit/#Bayesian-models>

	More discussion here: <https://discourse.julialang.org/t/best-fit-parameter-error-bar-using-optimization-jl/103186/6>
"""

# ╔═╡ 0f57e294-1c1f-4f73-8e08-ab9c0bba9f0f
md"""
## Plot
"""

# ╔═╡ edf0fc17-b5a9-4475-bdd7-aa6f2673a39d
sol_glm .- sol_optim

# ╔═╡ 0c73d53c-8b6e-4b13-a253-27046c1bd741
sol_glm .- sol_manual

# ╔═╡ 983cf394-29a5-4363-b07b-ccc91c5475e2
sol_optim .- sol_manual

# ╔═╡ c3d88d47-93f5-4f39-98ea-eb76f8a4974d
# let
# 	fig, ax, p = errorbars(log_period, k_mag, k_mag_err)
	
# 	scatter!(ax, log_period, k_mag)

# 	lines!(ax, log_period, m .* log_period .+ b; color = :orange)

# 	ax.xlabel = L"\log_{10}(\text{Period [days]})"
# 	ax.ylabel = "Ks"

# 	fig
# end

# ╔═╡ Cell order:
# ╠═b944163c-02c4-4ace-a866-ae6e6f7115ef
# ╟─434b1c97-4bae-49da-8eaf-207a23a6af7d
# ╟─6096dc80-a434-4f1a-9af5-40fbbf897d05
# ╠═002fa353-906a-4e3c-bd8e-12d681417240
# ╠═3bfaf85e-7cd8-4fa8-be0c-34afeae5ed06
# ╟─3f0abe3e-08cb-4c3c-9997-e654341bce12
# ╟─64da9892-c33b-4d00-ac8b-6564dea7b4e1
# ╠═a8ae7056-68a4-4bfc-9adb-e52821a1ae49
# ╟─4e69c1de-6cf0-448b-8b95-5b6c54795820
# ╟─da58364b-78f3-4113-9fc1-1a020a9e15dc
# ╠═1d61727f-b373-4394-b486-79eb63a71f37
# ╟─15516d38-b5b0-4dc4-8ce5-8d853dfc2c3d
# ╟─9773d632-f5cd-47d5-b97e-57a7b6ca3bf9
# ╟─35f8675d-d8f0-47b0-83b6-156bb75d373d
# ╠═af209cbf-1144-42cc-8c52-029988b0d5e4
# ╠═e1ad0f97-2d61-43dc-825c-bf7503a708dd
# ╠═616d9778-106e-46ae-bf87-92828e98339f
# ╟─e6d108e2-6ea1-4875-a3db-239248100bc7
# ╟─d8e6db3b-f1c3-46f5-b092-1ef25d762cdb
# ╠═2cf30811-3204-4498-9b73-8a3c3bba28e2
# ╟─19df3403-ef9d-4b07-8b14-3b10056475e8
# ╠═4a17877a-6a38-4364-9bd9-91dde011a42a
# ╠═1efa1d00-34da-48d1-90aa-c6b1d45a28ee
# ╟─480d8ef0-42d7-4089-9da2-1543baa2d02b
# ╠═07193622-e18b-474d-ab54-95cc09cf3bfc
# ╟─32609df9-b89a-491f-b420-7847e8c177d6
# ╟─18c7a687-7c3a-4c0d-8928-50b388683ca1
# ╠═8c44bde9-c64a-41c2-8b47-e98d40b3cf42
# ╟─38bc6713-6e89-47c0-9e03-b3e90e3182be
# ╠═e0e8909d-893f-4f35-baf2-c27dafeb23fa
# ╠═a3582cb7-c2a6-414c-bf7b-d764c2311c80
# ╠═3bffbcde-b2b1-445d-b2f7-59a53f33f90a
# ╠═925eb5ac-a5fa-40d9-8ab1-f10db7dab202
# ╟─2df6dced-3320-4d1a-a6c8-a46cb4c2bcf1
# ╟─0f57e294-1c1f-4f73-8e08-ab9c0bba9f0f
# ╠═edf0fc17-b5a9-4475-bdd7-aa6f2673a39d
# ╠═0c73d53c-8b6e-4b13-a253-27046c1bd741
# ╠═983cf394-29a5-4363-b07b-ccc91c5475e2
# ╠═c3d88d47-93f5-4f39-98ea-eb76f8a4974d
