### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> title = "Modeling 2: Create a User Defined Model using astropy.modeling"
#> layout = "layout.jlhtml"
#> date = "2025-12-31"
#> description = "TODO"
#> tags = ["TODO"]

using Markdown
using InteractiveUtils

# ╔═╡ cedb4512-0c81-47e7-8555-967e0067cd9c
using PlutoUI

# ╔═╡ 3a8e4e8d-9a90-40b8-beee-3e8ed55041c8
using PythonCall

# ╔═╡ d53e27ed-173f-43f2-997c-1633e363f005
using Unitful

# ╔═╡ de204217-aef5-436d-91f4-6254a46eda2b
using CairoMakie

# ╔═╡ 8cf89505-6bc3-4834-8caa-84d0ae04a2fd
using Distributions

# ╔═╡ 2d19d3b7-096c-44e5-bb86-7551095e0df9
md"""
# Modeling 2: Create a User Defined Model using astropy.modeling
"""

# ╔═╡ 36180cec-6175-11ef-337f-47bbbc831921
md"""
Original notebook at <https://learn.astropy.org/tutorials/User-Defined-Model.html>
"""

# ╔═╡ f85102a6-5dd2-46fc-b742-b4c327201576
TableOfContents()

# ╔═╡ 6d59a7dd-498a-4d3d-b5e3-7bebe3bbaf3b
md"""
## Authors

Rocio Kiman, Lia Corrales, Zé Vinícius, Stephanie T. Douglas
"""

# ╔═╡ 280d1a3e-59fa-4ad9-932e-40c468530048
md"""
## Learning Goals

- Define a new model with astropy
- Identify cases were a user-defined model could be useful
- Define models in two different ways:
  - Compound models
  - Custom models

This tutorial assumes the student knows how to fit data using astropy.modeling. This topic is covered in the [Models-Quick-Fit tutorial](https://learn.astropy.org/tutorials/Models-Quick-Fit.html).
"""

# ╔═╡ 15f69500-8188-4f8b-ba5a-956d559e442d
md"""
## Keywords

modeling, FITS, astrostatistics, matplotlib, model fitting, error bars, scatter plots
"""

# ╔═╡ 5f53f67f-e951-4561-8582-5308b24a24b7
md"""
## Summary
"""

# ╔═╡ 4027b761-1f9c-4fbc-9cf7-8d0d479e8a33
md"""
In this tutorial, we will learn how to define a new model in two ways: with a compound model and with a custom model.
"""

# ╔═╡ 4950cdb4-cdb6-4450-9083-ff6aebd0801e
md"""
Required packages for this notebook:
"""

# ╔═╡ 552f6911-9483-4e86-af21-ac0c016645c3
Pkg.status |> with_terminal

# ╔═╡ 2ca99127-b24b-4a5f-ad45-a6840bad2b5e
md"""
## Fit an emission line in a stellar spectrum
"""

# ╔═╡ 6f86435d-3dd2-465d-a707-02fee87f4179
md"""
M dwarfs are low mass stars (less than half of the mass of the sun). Currently we do not understand completely the physics inside low mass stars because they do not behave the same way higher mass stars do. For example, they stay magnetically active longer than higher mass stars. One way to measure magnetic activity is the height of the Hα emission line. It is located at 6563 Angstroms at the spectrum.
"""

# ╔═╡ 209a0f49-5405-4853-9771-aecc0512aba7
md"""
Let's search for a spectrum of an M dwarf in the Sloan Digital Sky Survey (SDSS). First, we are going to look for the spectrum in the SDSS database. SDSS has a particular way to identify the stars it observes: it uses three numbers: Plate, Fiber and MJD (Modified Julian Date). The star we are going to use has:

- Plate: 1349
- Fiber: 216
- MJD: 52797

So go ahead, put this numbers in the website and click on Plot to visualize the spectrum. Try to localize the Hα line.
"""

# ╔═╡ d47c0afa-d70b-4a75-b527-02a387e07a69
md"""
We could download the spectrum by hand from this website, but we are going to import it using the SDSSClass from astroquery.sdss. We can get the spectrum using the plate, fiber and mjd in the following way:
"""

# ╔═╡ d665e15e-7500-4116-b384-05eb451a6ef8
SDSS = pyimport("astroquery.sdss" => "SDSS")

# ╔═╡ ae9ac725-d5ee-4e10-adec-aabc8b73dd5b
spectrum = SDSS.get_spectra(plate = 1349, fiberID = 216, mjd = 52797) |> only |> collect

# ╔═╡ fa4d872a-a906-4a01-bf68-1185fce80d4e
md"""
One way to check what is inside the fits file spectrum is the following:
"""

# ╔═╡ e52f5fa7-8f40-4292-8c38-91d326092e30
spectrum[2].columns

# ╔═╡ 860f0297-5a3d-46ea-b7d5-940da5635dbb
md"""
To plot the spectrum we need the flux as a function of wavelength (usually called lambda or λ). Note that the wavelength is in log scale: loglam, so we calculate ``10^λ`` to remove this scale.
"""

# ╔═╡ e9f7aac6-5780-43d4-9eb9-720b03c7f9ea
flux = pyconvert(Vector, spectrum[2].data["flux"])

# ╔═╡ 34a8d7cb-64b6-4813-ae83-3acd55d40f68
lam = exp10.(pyconvert(Vector, spectrum[2].data["loglam"]))

# ╔═╡ e63c8f64-d1b6-4306-aefb-6133052c67d6
md"""
To find the units for flux and wavelength, we look in `fitsfile[0].header`.

FITS standard requires that the header keyword 'bunit' or 'BUNIT' contains the physical units of the array values. That's where we'll find the flux units.
"""

# ╔═╡ e5724917-10e1-4895-b4b5-957bb71b20f8
units_flux = spectrum[1].header["bunit"]

# ╔═╡ 398423be-e873-40e2-82ca-255c2c26ceb6
md"""
Different sources will definite wavelength information differently, so we need to check the documentation. For example, this [SDSS tutorial](https://www.sdss.org/dr12/tutorials/quicklook/#python) tells us what header keyword to look at.
"""

# ╔═╡ b589e768-59a2-48e2-81c2-88183173406a
units_wavelength_full = spectrum[1].header["WAT1_001"]

# ╔═╡ 43ef1e98-2a7a-46f7-8b78-a9cf830d4313
md"""
Now we are ready to plot the spectrum with all the information.
"""

# ╔═╡ a6cce4f0-d0f8-4955-a75c-8e5ef2fccb5e
let
    f = Figure()
    ax = Axis(f[1,1], xlabel = "Wavelength (Å)", ylabel = "Flux ($units_flux)")

    lines!(ax, lam, flux, color = :black)
    vlines!(ax, [6563], linestyle = :dash)

    xlims!(ax, 6300, 6700)

    f
end

# ╔═╡ 904d7aec-03c5-4e64-9e28-ffd5326d8163
md"""
We just plotted our spectrum! Check different ranges of wavelength to see how the full spectrum looks like in comparison to the one we saw before.
"""

# ╔═╡ 91f4040e-e38f-4804-8a64-434eb734e8b6
md"""
## Fit an Emission Line with a Gaussian Model
"""

# ╔═╡ bd754db1-3e42-4986-8c0d-bb481383106b
md"""
The blue dashed line marks the Hα emission line. We can tell this is an active star because it has a strong emission line.
"""

# ╔═╡ 30a9a40b-4271-4251-a65b-a9acc723efea
md"""
Now, we would like to measure the height of this line. Let's use astropy.modeling to fit a gaussian to the Hα line. We are going to initialize a gaussian model at the position of the Hα line. The idea is that the gaussian amplitude will tell us the height of the line.
"""

# ╔═╡ 62cdedc8-f522-4a10-b166-94889d96346f
model = fit(Normal, flux)

# ╔═╡ 4b4beab1-172e-4b16-bbd4-4b4a62c381b8
let
    f = Figure()
    ax = Axis(f[1,1], xlabel = "Wavelength (Å)", ylabel = "Flux ($units_flux)")

    max_flux = maximum(flux)

    lines!(ax, lam, flux, color = :black)
    vlines!(ax, [6563], linestyle = :dash)
    #lines!(ax, lam, pdf(model, ).* flux*max_flux)

    xlims!(ax, 6300, 6700)

    f
end

# ╔═╡ 8879c1a4-accd-4afc-b636-a975c6cf929e
md"""
## Exercise
"""

# ╔═╡ Cell order:
# ╟─2d19d3b7-096c-44e5-bb86-7551095e0df9
# ╟─36180cec-6175-11ef-337f-47bbbc831921
# ╠═cedb4512-0c81-47e7-8555-967e0067cd9c
# ╟─f85102a6-5dd2-46fc-b742-b4c327201576
# ╟─6d59a7dd-498a-4d3d-b5e3-7bebe3bbaf3b
# ╟─280d1a3e-59fa-4ad9-932e-40c468530048
# ╟─15f69500-8188-4f8b-ba5a-956d559e442d
# ╟─5f53f67f-e951-4561-8582-5308b24a24b7
# ╟─4027b761-1f9c-4fbc-9cf7-8d0d479e8a33
# ╟─4950cdb4-cdb6-4450-9083-ff6aebd0801e
# ╟─552f6911-9483-4e86-af21-ac0c016645c3
# ╟─2ca99127-b24b-4a5f-ad45-a6840bad2b5e
# ╟─6f86435d-3dd2-465d-a707-02fee87f4179
# ╟─209a0f49-5405-4853-9771-aecc0512aba7
# ╟─d47c0afa-d70b-4a75-b527-02a387e07a69
# ╠═3a8e4e8d-9a90-40b8-beee-3e8ed55041c8
# ╠═d665e15e-7500-4116-b384-05eb451a6ef8
# ╠═ae9ac725-d5ee-4e10-adec-aabc8b73dd5b
# ╟─fa4d872a-a906-4a01-bf68-1185fce80d4e
# ╠═e52f5fa7-8f40-4292-8c38-91d326092e30
# ╟─860f0297-5a3d-46ea-b7d5-940da5635dbb
# ╠═d53e27ed-173f-43f2-997c-1633e363f005
# ╠═e9f7aac6-5780-43d4-9eb9-720b03c7f9ea
# ╠═34a8d7cb-64b6-4813-ae83-3acd55d40f68
# ╟─e63c8f64-d1b6-4306-aefb-6133052c67d6
# ╠═e5724917-10e1-4895-b4b5-957bb71b20f8
# ╟─398423be-e873-40e2-82ca-255c2c26ceb6
# ╠═b589e768-59a2-48e2-81c2-88183173406a
# ╟─43ef1e98-2a7a-46f7-8b78-a9cf830d4313
# ╠═de204217-aef5-436d-91f4-6254a46eda2b
# ╠═a6cce4f0-d0f8-4955-a75c-8e5ef2fccb5e
# ╟─904d7aec-03c5-4e64-9e28-ffd5326d8163
# ╟─91f4040e-e38f-4804-8a64-434eb734e8b6
# ╟─bd754db1-3e42-4986-8c0d-bb481383106b
# ╟─30a9a40b-4271-4251-a65b-a9acc723efea
# ╠═8cf89505-6bc3-4834-8caa-84d0ae04a2fd
# ╠═62cdedc8-f522-4a10-b166-94889d96346f
# ╠═4b4beab1-172e-4b16-bbd4-4b4a62c381b8
# ╟─8879c1a4-accd-4afc-b636-a975c6cf929e
