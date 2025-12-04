### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> title = "Plotting cosmological redshift and age"
#> layout = "layout.jlhtml"
#> date = "2025-11-21"
#> description = "Work with units in astrophysical calculations."
#> tags = ["cosmology", "physics", "units", "plots"]
#> 
#>     [[frontmatter.author]]
#>     name = "Ian Weaver"
#>     url = "https://github.com/icweaver"

using Markdown
using InteractiveUtils

# ╔═╡ 1d4c2ee6-c6eb-11f0-8669-cd67adc8e577
begin
	using Pkg
	Pkg.add(["PlutoUI", "CairoMakie", "DynamicQuantities", "LaTeXStrings", "MathTeXEngine"])
	Pkg.add(; url = "https://github.com/JuliaAstro/Cosmology.jl", rev = "units")
	
	# Cosmological analysis
	using Cosmology: cosmology, angular_diameter_dist, age
	
	# Plotting
	using CairoMakie:
		Makie,
		Axis,
		Cycled,
		lines,
		lines!,
		hideydecorations!,
		linkxaxes!,
		axislegend
	using LaTeXStrings: @L_str
	using MathTeXEngine: set_texfont_family!, FontFamily
	set_texfont_family!(FontFamily("TeXGyreHeros"))
	
	# Units
	using DynamicQuantities: @u_str, @us_str, ustrip
end

# ╔═╡ f74f122b-2320-45c2-a3e7-ae049f6a897d
using PlutoUI: TableOfContents

# ╔═╡ 33515b2a-5ee4-4eab-9c7e-4aa6780ee369
md"""
# Plotting cosmological redshift and age

*Authors: Ian Weaver*

This notebook is modified from <https://learn.astropy.org/tutorials/redshift-plot.html>

!!! tip ""
	## Learning Goals
	
	* Plot relationships using Makie.jl

    * Add a second axis to the plot

    * Relate distance, redshift, and age for two different types of cosmologies


!!! note ""
	## Keywords

	cosmology, physics, units, plots

!!! warning ""
	## Summary

	Each redshift corresponds to an age of the universe, so if we're plotting some quantity against redshift, it's often useful show the universe age too. Using [Cosmology.jl](https://juliaastro.org/Cosmology), we'll visualize how this relationship between the two changes depending on the type of cosmology we assume.
"""

# ╔═╡ 3f0b8192-6c8d-4fc4-9503-4f5f4876b6f8
md"""
### Imports
"""

# ╔═╡ 3124c523-bfa7-45d5-9cb9-14b26f838ec4
md"""
## Define model

In this tutorial we'll show how to create the following plot:
"""

# ╔═╡ 644f473e-3ca2-4727-b629-66ee83a4d9bf
md"""
We start with a cosmology object representing a flat cosmology with the following parameters:

```math
\begin{align}
h &= 0.7 \\
\Omega_\Lambda &= 0.7 \\
\Omega_M &= 0.3 \\
\Omega_R &= 0
\end{align}
```
"""

# ╔═╡ 8d53ad6d-9a30-4593-b284-3d86e6137bed
cosmo = cosmology(h = 0.7, OmegaM = 0.3, OmegaR = 0)

# ╔═╡ 9be9ec4c-8a13-46a5-aeec-e68c80e9eba1
md"""
## Plot 

Now we need an example quantity to plot versus redshift. Let's use the angular diameter distance, which is the ratio of the proper transverse size of an object at redshift ``z_2`` to its angular size in radians, as seen by an observer at ``z_1``. By default, we assume ``z_1 = 0``.

For a collection of redshifts, `zvals`, we compute the following angular diameter distances using `Cosmology.angular_diameter_dist`:
"""

# ╔═╡ 620d15bc-f613-431b-8e6e-d1dbe100d933
zvals = 0 : 0.1 : 6

# ╔═╡ 5f70dcc4-8aab-44b6-aa82-ffd4dae92b74
d = [angular_diameter_dist(us"Constants.Mpc", cosmo, z) for z in zvals]

# ╔═╡ 20d549be-006c-48df-8e80-1bec80d62b55
md"""
!!! todo
	Trying out experimental DQ unit support here: <https://github.com/JuliaAstro/Cosmology.jl/tree/units>
"""

# ╔═╡ c6c7ee1f-ea7c-47d3-bda4-77e87c83ad75
md"""
Plotting against `zvals`, we produce the following curve:
"""

# ╔═╡ 0ba60166-cbd1-4f10-9b2b-259b7837a796
let
	fig, ax, p = lines(zvals, d;
		color = Cycled(2),
		axis = (
			# Uncomment this line for manual unit control
			# dim2_conversion = Makie.DQConversion(us"Constants.Mpc"),
			xminorticksvisible = true,
			xlabel = "Redshift",
			ylabel = "Angular diameter distance",
		),
	)

	fig
end

# ╔═╡ b8a0440f-3da8-488e-a52f-b9cf98f977b2
md"""
!!! note
	We use the [`Cycled`](https://docs.makie.org/v0.21/explanations/theming/themes#Manual-cycling-using-Cycled) object from Makie.jl to plot the second default color in our colormap series instead of the first. We do this because we are going to plot a second cosmology on this plot soon, and would like it to appear first in the series.

	For more on plotting with DynamicQuantities, see [this section](https://docs.makie.org/dev/explanations/dim-converts#Experimental-DynamicQuantities.jl-support) of the Makie.jl documentation.
"""

# ╔═╡ 0b7a66e1-75d1-4197-90ac-723420f276e0
md"""
### Twin axis

It would be useful to see the corresponding universe ages at each redshift. Let's compute this with `Cosmology.age` and plot this along the top axis:
"""

# ╔═╡ ba65782d-fa00-4e3d-9e49-fdafd3a7e47c
# Just show a few of the ages for clarity
age_ticks, age_vals = let
	n = 15 # Sample every n point
	z_sampled = zvals[begin:n:end]
	age_sampled = [age(us"Gyr", cosmo, z) for z in z_sampled]
	age_sampled_vals = ustrip.(u"Gyr", age_sampled)
	z_sampled, age_sampled_vals
end

# ╔═╡ 5ff49cc2-35b1-463a-ad95-5934069f8420
let
	# Default model
	f, ax1, p = lines(zvals, d;
		color = Cycled(2),
		axis = (
			dim2_conversion = Makie.DQConversion(us"Constants.Mpc"),
		),
	)
	
	# Age axis
	ax2 = Axis(f[1, 1];
		xaxisposition = :top,
		xticks = age_ticks,
		xtickformat = x -> string.(round.((age_vals); digits = 2)),
	)

	# Labels
	ax1.xlabel = "Redshift"
	ax1.ylabel = "Angular diameter distance"
	ax2.xlabel = "Time since Big Bang (Gyr)"

	# Grid lines + ticks
	ax1.xminorticksvisible = true
	ax1.xgridvisible = false
	ax2.xgridvisible = false
	hideydecorations!(ax2)
	linkxaxes!(ax1, ax2)

	f
end

# ╔═╡ dd63c9c9-50e9-4d7b-86e3-8874134d502b
md"""
!!! note
	For clarity, we suppress the vertical grid lines coming from each axis.
"""

# ╔═╡ 4d29f960-6f15-4d18-aa17-e80be8a2770a
md"""
## Adding another cosmology

Finally, let's add a second cosmology for comparison. For this example, we will use the [Planck 2013](https://ui.adsabs.harvard.edu/abs/2014A%26A...571A..16P/abstract) model:
"""

# ╔═╡ 8a54c1e4-abca-4da6-a738-0b2d0d2ac3da
cosmo_planck = cosmology(h = 0.6777, OmegaM = 0.30712)

# ╔═╡ 9f01bb56-455c-45c2-98b2-20a38f141ebd
d_planck = [angular_diameter_dist(us"Constants.Mpc", cosmo_planck, z) for z in zvals]

# ╔═╡ 5879e1c3-b6b0-4216-b00c-af0bf819ab47
md"""
!!! todo
	See if we can upstream Chris's work <https://github.com/JuliaAstro/Cosmology.jl/issues/43#issuecomment-2790915234>
"""

# ╔═╡ dbcc87a3-fe96-42aa-87b7-e8147ac1a48c
fig = let
	# Plank 2013 model
	f, ax1, p = lines(zvals, d_planck;
		label = "Planck 2013",
		axis = (;
			dim2_conversion = Makie.DQConversion(us"Constants.Mpc"),
			xlabel = "Redshift",
			ylabel = "Angular diameter distance",
		)
	)

	# Default model
	lines!(ax1, zvals, d; label = L"h = 0.7,\ \Omega_M = 0.3,\ \Omega_\Lambda = 0.7")
		
	# Age axis
	ax2 = Axis(f[1, 1];
		xaxisposition = :top,
		xticks = age_ticks,
		xtickformat = x -> string.(round.((age_vals); digits = 2)),
	)
	ax2.xlabel = "Time since Big Bang (Gyr)"

	# Grid lines + ticks
	ax1.xminorticksvisible = true
	ax1.xgridvisible = false
	ax2.xgridvisible = false
	hideydecorations!(ax2)
	linkxaxes!(ax1, ax2)

	# Legend
	axislegend(ax1; position = :rb)

	# Display
	f
end

# ╔═╡ 6734b49a-5a4b-4654-9fd7-e94df8872e88
fig

# ╔═╡ d7d4e075-0445-46c2-b88e-56a3f1df0bba
md"""
!!! tip
	 This figure can be saved with:

	```julia
	using CairoMakie

	save("my_plot.pdf", fig)
	```
"""

# ╔═╡ 6afad2fa-0555-400a-9de7-e341e5956955
md"""
# Notebook setup 🔧
"""

# ╔═╡ aea0ba58-a37f-4b55-8fe6-55cbb79e058a
TableOfContents()

# ╔═╡ Cell order:
# ╟─33515b2a-5ee4-4eab-9c7e-4aa6780ee369
# ╟─3f0b8192-6c8d-4fc4-9503-4f5f4876b6f8
# ╠═1d4c2ee6-c6eb-11f0-8669-cd67adc8e577
# ╟─3124c523-bfa7-45d5-9cb9-14b26f838ec4
# ╟─6734b49a-5a4b-4654-9fd7-e94df8872e88
# ╟─644f473e-3ca2-4727-b629-66ee83a4d9bf
# ╠═8d53ad6d-9a30-4593-b284-3d86e6137bed
# ╟─9be9ec4c-8a13-46a5-aeec-e68c80e9eba1
# ╠═620d15bc-f613-431b-8e6e-d1dbe100d933
# ╠═5f70dcc4-8aab-44b6-aa82-ffd4dae92b74
# ╟─20d549be-006c-48df-8e80-1bec80d62b55
# ╟─c6c7ee1f-ea7c-47d3-bda4-77e87c83ad75
# ╠═0ba60166-cbd1-4f10-9b2b-259b7837a796
# ╟─b8a0440f-3da8-488e-a52f-b9cf98f977b2
# ╟─0b7a66e1-75d1-4197-90ac-723420f276e0
# ╠═ba65782d-fa00-4e3d-9e49-fdafd3a7e47c
# ╠═5ff49cc2-35b1-463a-ad95-5934069f8420
# ╟─dd63c9c9-50e9-4d7b-86e3-8874134d502b
# ╟─4d29f960-6f15-4d18-aa17-e80be8a2770a
# ╠═8a54c1e4-abca-4da6-a738-0b2d0d2ac3da
# ╠═9f01bb56-455c-45c2-98b2-20a38f141ebd
# ╟─5879e1c3-b6b0-4216-b00c-af0bf819ab47
# ╠═dbcc87a3-fe96-42aa-87b7-e8147ac1a48c
# ╟─d7d4e075-0445-46c2-b88e-56a3f1df0bba
# ╟─6afad2fa-0555-400a-9de7-e341e5956955
# ╠═f74f122b-2320-45c2-a3e7-ae049f6a897d
# ╠═aea0ba58-a37f-4b55-8fe6-55cbb79e058a
