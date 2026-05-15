### A Pluto.jl notebook ###
# v0.20.25

#> [frontmatter]
#> title = "Modeling 2: Create a User Defined Model using astropy.modeling"
#> layout = "layout.jlhtml"
#> date = "2025-12-31"
#> description = "TODO"
#> tags = ["models", "model fitting", "astrostatistics", "catalog", "query", "Makie", "plots", "errorbars", "scatter plots"]

using Markdown
using InteractiveUtils

# ╔═╡ 33e2fa03-a93d-400e-b3aa-87683e5bbd5a
begin
    import Pkg
    Pkg.activate(; temp = true)

    Pkg.add(
        [
            Pkg.PackageSpec(; name = "Downloads"),
            Pkg.PackageSpec(; name = "FITSFiles"),
            Pkg.PackageSpec(; name = "DynamicQuantities"),
            Pkg.PackageSpec(; name = "LsqFit"),
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(;
                url = "https://github.com/MakieOrg/Makie.jl",
                subdir = "Makie",
                rev = "ff/breaking-0.25",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/MakieOrg/Makie.jl",
                subdir = "ComputePipeline",
                rev = "ff/breaking-0.25",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/MakieOrg/Makie.jl",
                subdir = "CairoMakie",
                rev = "ff/breaking-0.25",
            ),
        ]
    )

    Pkg.add("TOML")

    using Downloads: download
    using FITSFiles: FITSFiles
    using DynamicQuantities
    using CairoMakie
    using LsqFit

    deps_ready = true
end

# ╔═╡ 9d88288c-6415-4851-a52f-0008fecacf0e
begin
    deps_ready

    using TOML: TOML
    using PlutoUI: TableOfContents
    using Test: @test
end

# ╔═╡ 4027b761-1f9c-4fbc-9cf7-8d0d479e8a33
md"""
## Summary

In this tutorial, we will learn how to define a new model in two ways: with a compound model and with a custom model.
"""

# ╔═╡ f5b35c02-47ce-4f61-8371-c3033f1c8017
md"""
### Packages 📦
"""

# ╔═╡ 6f86435d-3dd2-465d-a707-02fee87f4179
md"""
## Introduction

M dwarfs are low mass stars (less than half of the mass of the sun). Currently we do not understand completely the physics inside low mass stars because they do not behave the same way higher mass stars do. For example, they stay magnetically active longer than higher mass stars. One way to measure magnetic activity is the height of the [Hα emission line](https://en.wikipedia.org/wiki/Hydrogen-alpha). It is located at 6563 Angstroms at the spectrum.

Let's search for a spectrum of an M dwarf in the Sloan Digital Sky Survey (SDSS). First, we are going to look for the spectrum in the [SDSS database](https://dr12.sdss.org/basicSpectra). SDSS has a particular way to identify the stars it observes: it uses three numbers: Plate, Fiber and MJD (Modified Julian Date). The star we are going to use has:

- Plate: 1349
- MJD: 52797
- Fiber: 216

So go ahead, put these numbers in the website and click on Plot to visualize the spectrum. Try to localize the Hα line.

We will download the spectrum by hand from their website.

!!! todo
	Add SDSS data querying functionality, maybe in VirtualObservatory.jl. Could possibly pull some functionality from the old [SDSS.jl package](https://github.com/kbarbary/SDSS.jl) too.
"""

# ╔═╡ ea94f687-01d6-4b6d-b436-eafa14f933c0
spectrum_filename = download("https://data.sdss.org/sas/dr17/sdss/spectro/redux/26/spectra/1349/spec-1349-52797-0216.fits")

# ╔═╡ 21120560-b6e7-47ff-9834-592463a942b9
spectrum = FITSFiles.fits(spectrum_filename)

# ╔═╡ fa4d872a-a906-4a01-bf68-1185fce80d4e
md"""
Let's take a look at the `data` and `cards` fields for our `COADD` table:
"""

# ╔═╡ e52f5fa7-8f40-4292-8c38-91d326092e30
spectrum[:COADD].data

# ╔═╡ a30cf18a-8dc2-4ee4-a238-792ff0e71e93
spectrum[:COADD].cards # Can also access by index: spectrum[2].cards

# ╔═╡ 860f0297-5a3d-46ea-b7d5-940da5635dbb
md"""
To plot the spectrum from this table, we need the flux as a function of wavelength (usually called lambda or λ). Note that the wavelength is in log scale: loglam, so we calculate ``10^λ`` to remove this scale:
"""

# ╔═╡ 75ff1bd7-7390-41ef-805e-da7e9fcc5ff7
lam = exp10.(spectrum[:COADD].data["loglam"]) * us"Å" # Typed \Angstrom<TAB>

# ╔═╡ c12338da-d598-4dc3-93ba-56f09cf2aaf6
flux = spectrum[:COADD].data["flux"] * us"erg/cm^2/s/Å"

# ╔═╡ c3e3d957-92f5-4e51-a266-d557689a5d0a
md"""
!!! note "A note on units"
	To find the units for flux and wavelength, we look in `spectrum[:COADD].cards`.

	**flux:** FITS standard requires that the header keyword `bunit` or `BUNIT` contains the physical units of the array values. That's where we'll find the flux units.

	**wavelength:** Different sources will define wavelength information differently, so we need to check the documentation. For example, this [SDSS tutorial](https://www.sdss.org/dr12/tutorials/quicklook/#python) tells us what header keyword to look at:

"""

# ╔═╡ e5724917-10e1-4895-b4b5-957bb71b20f8
spectrum[1].cards["BUNIT"]

# ╔═╡ b589e768-59a2-48e2-81c2-88183173406a
spectrum[1].cards["WAT1_001"]

# ╔═╡ 43ef1e98-2a7a-46f7-8b78-a9cf830d4313
md"""
Now we are ready to plot the spectrum with all the information.
"""

# ╔═╡ f3ea8e3e-eb59-423e-8fb4-eaee84aae1dd
md"""
## Observed data
"""

# ╔═╡ a6cce4f0-d0f8-4955-a75c-8e5ef2fccb5e
let
    fig, ax, p = lines(
        lam, flux;
        color = :black,
        axis = (; xlabel = "Wavelength", ylabel = "Flux 1e-17"),
    )

    vlines!(ax, 6563u"Å", linestyle = :dash)

    xlims!(ax, 6300u"Å", 6700u"Å")

    fig
end

# ╔═╡ 904d7aec-03c5-4e64-9e28-ffd5326d8163
md"""
We just plotted our spectrum! Check different ranges of wavelength to see how the full spectrum looks like in comparison to the one we saw before.
"""

# ╔═╡ bd754db1-3e42-4986-8c0d-bb481383106b
md"""
## Fit a Gaussian Model

The blue dashed line marks the Hα emission line. We can tell this is an active star because it has a strong emission line.

Now, we would like to measure the height of this line. Let's fit a gaussian to the Hα line, initialized at its center. The idea is that the gaussian amplitude will tell us the height of the line.
"""

# ╔═╡ 62cdedc8-f522-4a10-b166-94889d96346f
@. gaussian(x, p) = p[1] * exp(-(x - p[2])^2 / (2 * p[3]^2))

# ╔═╡ 96322bfd-0d14-4bda-bd3e-eff954eb63e6
p_fit = let
    p0 = [1.0, 6563.0, 10.0] # [amplitude, μ, σ]

    model = curve_fit(gaussian, lam.value, flux.value, p0)

    @info model

    coef(model)
end

# ╔═╡ 4b4beab1-172e-4b16-bbd4-4b4a62c381b8
let
    # Data
    fig, ax, p = lines(
        lam, flux;
        color = :black,
        axis = (; xlabel = "Wavelength", ylabel = "Flux 1e-17"),
    )

    # Model
    model = gaussian(lam.value, p_fit)u"erg/cm^2/s/Å"
    lines!(ax, lam, model)

    vlines!(ax, 6563u"Å", linestyle = :dash)
    xlims!(ax, 6300u"Å", 6700u"Å")

    fig
end

# ╔═╡ e4c5ec2a-b3de-4dc9-9536-717cf6cca7f5
md"""
We can see the fit is not doing a good job. Inspect the model parameters of this fit printed above the plot.
"""

# ╔═╡ 4e9a8448-86a2-46ad-807e-09ab4f56e40d
md"""
### Exercise

Go back to the previous plot and try to make the fit work. Note: Do not spend more than 10 minutes in this exercise. A couple of ideas to try:

- Is it not working because of the model we chose to fit?
- Is it not working because of the fitter we chose?
- Is it not working because of the range of data we are fitting?
- Is it not working because of how we are plotting the data?
"""

# ╔═╡ 2589fab7-3440-4d03-bc61-90b91f3f8fde
md"""
## Compound models

One model is not enough to make this fit work. We need to combine a couple of models to make a compound model. The idea is that we can add, divide, or multiply models that already exist and fit the compound model to our data.

For our problem we are going to combine the gaussian with a polynomial of degree 1 to account for the background spectrum close to the line. Take a look at the plot we made before to convince yourself that this is the case.

Now let’s make our compound model!
"""

# ╔═╡ fa9061a3-25fe-4467-9cb1-c50921e5e480
@. gaussian_compound(x, p) = gaussian(x, [p[1:3]]) + p[4] + p[5] * x

# ╔═╡ 7f03cf99-1cd5-41bd-a6eb-a072b71c0f87
md"""
After this point, we fit the data in exactly the same way as before, except we use a compound model instead of the gaussian model:
"""

# ╔═╡ 8d578108-7bbb-40c7-a64e-0677a30e41b9
p_fit_compound = let
    p0 = [1.0, 6563.0, 0.1, 0.0, 0.0]

    model = curve_fit(
        gaussian_compound,
        lam.value,
        flux.value,
        p0;
        lower = [0.0, 6563 - 0.5, 0, -Inf, -Inf],
        upper = [Inf, 6563 + 0.5, 10, Inf, Inf],
    )

    @info model

    coef(model)
end

# ╔═╡ 0cb42da2-3b6b-4ef2-868d-fc25dd5832da
let
    # Data
    fig, ax, p = lines(
        lam, flux;
        color = :black,
        axis = (; xlabel = "Wavelength", ylabel = "Flux 1e-17"),
    )

    # Model
    model = gaussian_compound(lam.value, p_fit_compound)u"erg/cm^2/s/Å"
    lines!(ax, lam, model)

    vlines!(ax, 6563u"Å", linestyle = :dash)
    xlims!(ax, 6300u"Å", 6700u"Å")

    fig
end

# ╔═╡ 8879c1a4-accd-4afc-b636-a975c6cf929e
md"""
### Exercise
"""

# ╔═╡ 3c0d99d6-ae1a-414b-af86-ead5c8211543
md"""
# Notebook setup 🔧
"""

# ╔═╡ 5d9e2890-00fb-4370-9337-f6e93d49e5ed
TableOfContents(; depth = 4)

# ╔═╡ b04e9eeb-e476-49ff-8bbf-5a7ae199df31
function frontmatter(path)
	prefix = "#> "
	is_fm = startswith(prefix)
	block = Iterators.takewhile(is_fm, Iterators.dropwhile(!is_fm, eachline(path)))
	toml = TOML.parse(join(chopprefix.(block, prefix), "\n"))
	return toml["frontmatter"]
end

# ╔═╡ a42a6e27-f5a9-4958-9b26-905fbb3dad9d
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ 2d19d3b7-096c-44e5-bb86-7551095e0df9
md"""
# Modeling 2: Create a User Defined Model

This notebook is modified from <https://learn.astropy.org/tutorials/2_user-defined-model.html>

_Original authors: Rocio Kiman, Lia Corrales, Zé Vinícius, Stephanie T. Douglas_

!!! tip "Learning goals"
	- Define a new model
	- Identify cases were a user-defined model could be useful
	- Define models in two different ways:
	  - Compound models
	  - Custom models

$(keywords())

!!! warning "Companion content"
	[learn.JuliaAstro > Modeling 1: Linear model fitting](https://learn.juliaastro.org/tutorials/models-1_linear_fitting/)
"""

# ╔═╡ Cell order:
# ╟─2d19d3b7-096c-44e5-bb86-7551095e0df9
# ╟─4027b761-1f9c-4fbc-9cf7-8d0d479e8a33
# ╟─f5b35c02-47ce-4f61-8371-c3033f1c8017
# ╠═33e2fa03-a93d-400e-b3aa-87683e5bbd5a
# ╟─6f86435d-3dd2-465d-a707-02fee87f4179
# ╠═ea94f687-01d6-4b6d-b436-eafa14f933c0
# ╠═21120560-b6e7-47ff-9834-592463a942b9
# ╟─fa4d872a-a906-4a01-bf68-1185fce80d4e
# ╠═e52f5fa7-8f40-4292-8c38-91d326092e30
# ╠═a30cf18a-8dc2-4ee4-a238-792ff0e71e93
# ╟─860f0297-5a3d-46ea-b7d5-940da5635dbb
# ╠═75ff1bd7-7390-41ef-805e-da7e9fcc5ff7
# ╠═c12338da-d598-4dc3-93ba-56f09cf2aaf6
# ╟─c3e3d957-92f5-4e51-a266-d557689a5d0a
# ╠═e5724917-10e1-4895-b4b5-957bb71b20f8
# ╠═b589e768-59a2-48e2-81c2-88183173406a
# ╟─43ef1e98-2a7a-46f7-8b78-a9cf830d4313
# ╟─f3ea8e3e-eb59-423e-8fb4-eaee84aae1dd
# ╠═a6cce4f0-d0f8-4955-a75c-8e5ef2fccb5e
# ╟─904d7aec-03c5-4e64-9e28-ffd5326d8163
# ╟─bd754db1-3e42-4986-8c0d-bb481383106b
# ╠═62cdedc8-f522-4a10-b166-94889d96346f
# ╠═96322bfd-0d14-4bda-bd3e-eff954eb63e6
# ╠═4b4beab1-172e-4b16-bbd4-4b4a62c381b8
# ╟─e4c5ec2a-b3de-4dc9-9536-717cf6cca7f5
# ╟─4e9a8448-86a2-46ad-807e-09ab4f56e40d
# ╟─2589fab7-3440-4d03-bc61-90b91f3f8fde
# ╠═fa9061a3-25fe-4467-9cb1-c50921e5e480
# ╟─7f03cf99-1cd5-41bd-a6eb-a072b71c0f87
# ╠═8d578108-7bbb-40c7-a64e-0677a30e41b9
# ╠═0cb42da2-3b6b-4ef2-868d-fc25dd5832da
# ╟─8879c1a4-accd-4afc-b636-a975c6cf929e
# ╟─3c0d99d6-ae1a-414b-af86-ead5c8211543
# ╠═5d9e2890-00fb-4370-9337-f6e93d49e5ed
# ╟─b04e9eeb-e476-49ff-8bbf-5a7ae199df31
# ╟─a42a6e27-f5a9-4958-9b26-905fbb3dad9d
# ╠═9d88288c-6415-4851-a52f-0008fecacf0e
