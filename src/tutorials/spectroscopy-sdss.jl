### A Pluto.jl notebook ###
# v0.20.25

#> [frontmatter]
#> layout = "layout.jlhtml"
#> title = "Spectroscopy with SDSS"
#> tags = ["spectroscopy", "SDSS", "units", "dust exctinction", "cosmology", "plots"]
#> date = "2026-05-16"
#> description = "Explore common tasks used for analyzing SDSS spectroscopic data"
#>
#>     [[frontmatter.author]]
#>     name = "Aditya Kumar Pandey"
#>     [[frontmatter.author]]
#>     name = "Chris Garling"
#>     [[frontmatter.author]]
#>     name = "Ian Weaver"

using Markdown
using InteractiveUtils

# ╔═╡ d7ed3be9-00c1-444d-a49d-bb7ce7c5bf03
begin
    import Pkg
    Pkg.activate(; temp = true)

    # TODO: Can remove when https://github.com/MakieOrg/Makie.jl/pull/5623
    # is upstreamed

    # Data viz
    Pkg.add(
        [
            Pkg.PackageSpec(;
                url = "https://github.com/icweaver/Makie.jl",
                subdir = "Makie",
                rev = "units-matrix",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/MakieOrg/Makie.jl",
                subdir = "CairoMakie",
                rev = "ff/breaking-0.25",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/MakieOrg/Makie.jl",
                subdir = "ComputePipeline",
                rev = "ff/breaking-0.25",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/icweaver/Measurements.jl",
                rev = "makie-v0.25.0",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaAstro/SpectrumBase.jl",
            ),
        ]
    )

    Pkg.add(["TOML", "PlutoUI", "DustExtinction", "Cosmology", "SkyCoords", "FITSFiles", "Downloads", "Unitful", "UnitfulAstro", "UnitfulEquivalences", "MathTeXEngine", "DataInterpolations"])

    # Analysis
    using DustExtinction, Cosmology, SkyCoords, FITSFiles, SpectrumBase

    using Downloads

    # Plotting
    using CairoMakie
    using MathTeXEngine: set_texfont_family!, FontFamily
    set_texfont_family!(FontFamily("TeXGyreHeros"))

    # Units
    using Unitful, UnitfulAstro, UnitfulEquivalences, Measurements

    deps_ready = true # Can drop this when above patches upstreamed
end;

# ╔═╡ e4c5b28c-c4fb-4726-ab7a-8f1cc04280bd
begin
    deps_ready

    using TOML: TOML
    using PlutoUI
    using Test: @test
end

# ╔═╡ 29945a91-d0c5-470b-81e4-99329e874ee0
md"""
## Summary

Spectroscopy is one of the most powerful tools in an astronomer's toolkit. By splitting light into its component wavelengths, we can learn about a source's chemical composition, temperature, velocity, and much more.

In this tutorial we walk through a realistic end-to-end spectroscopy workflow using:

- [FITSFiles.jl](https://juliaastro.org/FITSFiles): FITS file I/O
- [SpectrumBase.jl](https://juliaastro.org/SpectrumBase/stable/): core spectral types
- [Unitful.jl](https://juliaphysics.github.io/Unitful.jl/stable/) + [UnitfulAstro.jl](https://juliaastro.org/UnitfulAstro/stable/): physical units
- [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl): uncertainty propagation
- [DustExtinction.jl](https://juliaastro.org/DustExtinction/stable/): dust dereddening
- [Cosmology.jl](https://juliaastro.org/Cosmology/stable/): cosmological distances
- [Makie.jl](https://docs.makie.org/stable/): plotting
"""

# ╔═╡ 5d3a54c7-97d8-48c8-af80-8b7e2ad13ff8
md"""
### Packages 📦
"""

# ╔═╡ e35fd5df-e01e-43af-99cf-546ca8c8f342
md"""
## Part 1: Loading a real SDSS spectrum

We use a publicly available spectrum from SDSS DR14 -- a galaxy on plate 1323, MJD 52797, fiber 12. The FITS file stores flux in units of 10⁻¹⁷ erg s⁻¹ cm⁻² Å⁻¹ and an `ivar` (inverse-variance) column that encodes real per-pixel measurement uncertainties:
"""

# ╔═╡ 645822c8-d735-4ec5-8db7-7cdbb3fd191b
sdss_file = let
    sdss_url = "https://dr14.sdss.org/optical/spectrum/view/data/format=fits/spec=lite?plateid=1323&mjd=52797&fiberid=12"
    sdss_file = joinpath(@__DIR__, "sdss_example.fits")
    isfile(sdss_file) ? sdss_file : Downloads.download(sdss_url, sdss_file)
end

# ╔═╡ 8b96f755-9cfd-479d-a32b-6a664004793b
md"""
We use FITSFiles.jl to load this file and explore its contents:
"""

# ╔═╡ 7758319f-471a-4d9d-9e31-103ee81fe9ab
hdus = fits(sdss_file)

# ╔═╡ 9a8dac0b-1ca0-43ff-b822-27690eae157c
md"""
We next load the wavelength and flux information from the `COADD` table into a spectrum object with SpectrumBase.jl, complete with unit support from Unitful.jl and uncertainty support with Measurements.jl:
"""

# ╔═╡ ae4d28da-8c69-481c-8589-149e51ddd9fa
spec = let
    # Unpack fields
    hdu = hdus[2].data
    loglam = hdu["loglam"] # log wavelength
    ivar = hdu["ivar"] # Inverse variance (σ_flux = 1 / √ivar)
    flux_err = inv.(sqrt.(ivar))

    # Store result
    spectrum(
        exp10.(loglam) * u"Å", # wavelength
        (hdu["flux"] .± flux_err) * u"erg/s/cm^2/Å" # flux
    )
end

# ╔═╡ df8b47f8-7692-4032-a9e5-e99397f7808a
md"""
We then can use Makie.jl to plot this spectrum directly:

!!! todo
	Upstream a SpectrumBase.jl recipe
"""

# ╔═╡ 143e0860-0704-437e-b2d3-22f5298a9198
let
    wav, flux = spec.spectral_axis, spec.flux_axis

    fig, ax, p = band(
        wav,
        flux;
        axis = (;
            xlabel = "Wavelength",
            ylabel = "Flux density * 1e-17",
            title = "SDSS Galaxy — plate 1323, fiber 12",
        ),
        color = :orange,
    )

    lines!(ax, wav, flux)

    xlims!(ax, 6450u"Å", 6800u"Å")

    fig
end

# ╔═╡ 366d72e2-8f5f-4b8a-986b-467f662aaf5b
md"""
## Part 2: Physical unit conversions

[`Unitful.jl`](https://painterqubits.github.io/Unitful.jl/stable/) also makes it easy to convert between wavelength, frequency, and energy representations:
"""

# ╔═╡ f0833976-0cd9-45e3-ac8e-0a2ed1073d1f
H_alpha = 6563.0u"Å"

# ╔═╡ 24edd839-07bb-477b-b059-730bfe57db5f
H_alpha |> u"nm"

# ╔═╡ c0b77309-56c6-46df-b4ab-360c179c1934
H_alpha |> u"μm"

# ╔═╡ 0d5c2f5c-f57d-43ef-9633-94f3411aa977
uconvert(u"erg", H_alpha, Spectral())

# ╔═╡ eae70e1e-7182-4257-aaa0-4b8f88c901ce
md"""
## Part 3: Dust extinction and dereddening

We can use DustExctinction.jl to deredden our spectrum.
"""

# ╔═╡ 43fad909-2733-48db-8c0b-cefd532ed267
eq = ICRSCoords(deg2rad(178.90417), deg2rad(0.66278))

# ╔═╡ c8529662-9746-4d0f-bef5-796c33d09ac5
md"""
!!! warning
	Where did these coords come from? They don't seem to match from the header I think?
"""

# ╔═╡ 6e6972a6-1648-40e4-a1e6-beed98bbc93e
hdus[1].cards["RA"], hdus[1].cards["DEC"]

# ╔═╡ 45849da1-b473-4d19-b140-d9d3fca39799
gal = convert(GalCoords, eq)

# ╔═╡ e58714fd-8432-43c2-be11-8c9a22e2de85
ENV["DATADEPS_ALWAYS_ACCEPT"] = true

# ╔═╡ d42d14f8-0e3a-46d0-b04f-c7e735591aee
dustmap = SFD98Map()

# ╔═╡ eff1ff04-b875-4b01-b7c6-67b07c682396
ebv = dustmap(gal.l, gal.b)

# ╔═╡ 1ba509c8-a4fa-4aeb-b811-6e08259a6d18
Av = 3.1 * ebv

# ╔═╡ 9866a403-9374-42db-8beb-57e9673e71b8
spec_dered = SpectrumBase.deredden(spec, Av)

# ╔═╡ 9e824cb8-ff28-40e0-be4f-36e579fae709
let
    fig, ax, p = lines(
        spec.spectral_axis, spec.flux_axis;
        axis = (;
            xlabel = "Wavelength",
            ylabel = "Flux density",
            title = "CCM89 Dust Dereddening",
        ),
        label = "original",
    )

    lines!(ax, spec_dered.spectral_axis, spec_dered.flux_axis; label = "dereddened")

    axislegend()

    fig
end

# ╔═╡ 8d6e0d8b-bbea-4144-bf22-9a94c38b33e6
md"""
## Part 4: Blackbody spectra and stellar luminosity density

`blackbody` returns spectral radiance ``B(\lambda, T)``. The stellar luminosity density is:

```math
L_\lambda(\lambda, T) = 4\pi^2 R^2 B(\lambda, T)
```
"""

# ╔═╡ 07baa347-0ac7-43e6-ab67-4e99c3774f65
wav_bb = range(3_000, 20_000; length = 500)u"Å"

# ╔═╡ 57f114fb-3098-47d0-bb95-3e0f9fe41790
spec_solar = blackbody(wav_bb, 5778.0u"K")

# ╔═╡ 97bf0849-d92e-47cb-914a-6f9a1d02ce94
B_solar = spec_solar.flux_axis

# ╔═╡ 84052581-f53d-4bed-9776-1ee31d7c38fe
L(λ, T; R = u"Rsun") = 4 * π^2 * R^2 * blackbody(λ, T).flux_axis

# ╔═╡ 3d996ca3-26c6-4ebe-896e-c6b0378ae85c
L_solar = L(wav_bb, 5778.0u"K")

# ╔═╡ 48a59ba2-73bc-42dd-bec4-7be6b2aaff90
maximum(B_solar) |> u"erg/s/cm^2/sr/Å"

# ╔═╡ b2ac804d-2ebd-42e7-8625-99a91e0458f6
maximum(L_solar) |> u"erg/s/Å"

# ╔═╡ 4147f9f6-df33-4ac1-b31f-17dab96c951f
md"""
Or in a loop for different temperatures:
"""

# ╔═╡ b742966f-a179-4a6a-93ed-fbabb3e3803e
let
    # Set up figure
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        yscale = log10,
        yticks = LogTicks(LinearTicks(5)),
        dim2_conversion = Makie.UnitfulConversion(u"erg/s/cm^2/Å"),
        xlabel = "Wavelength",
        ylabel = "Flux density",
        title = "Blackbody Spectra",
    )

    # Add plots
    temps_bb = [30_000, 10_000, 5778.0, 3_000]u"K"
    colors_bb = [:purple, :steelblue, :orange, :red]

    for (T, color) in zip(temps_bb, colors_bb)
        B = blackbody(wav_bb, T).flux_axis
        lines!(ax, wav_bb, B; color, label = string(T))
    end

    axislegend()

    fig
end

# ╔═╡ 5b3a311f-ad28-44cf-bd9a-53bc8f6842b9
md"""
## Part 5: Cosmological redshift and surface-brightness dimming

For a source at redshift ``z``, the observed flux density is:

```math
F_\lambda^\text{obs}(\lambda_\text{obs}, z) =
	\frac{L_\lambda(\lambda_\text{rest})}{(1 + z) 4\pi d_L^2(z)}
```
"""

# ╔═╡ 61bfc289-e129-4bca-b9ec-36e39d91c151
cosmo = cosmology()

# ╔═╡ 15601673-c10d-48aa-98c7-5fd389f35aa3
function F_obs(λ_obs, z; T = 5778.0u"K", cosmo = cosmo)
    λ_rest = λ_obs * inv(1 + z)
    d_L = luminosity_dist(cosmo, z)
    return L(λ_rest, T) / ((1 + z) * 4 * π * d_L^2)
end

# ╔═╡ cbf45b3f-8adf-4a40-bfc1-f73a4596e113
let
    fig = Figure()

    ax = Axis(
        fig[1, 1];
        yticks = LinearTicks(5),
        dim2_conversion = Makie.UnitfulConversion(u"erg/s/cm^2/Å"),
        xlabel = "Wavelength (observed)",
        ylabel = "Flux density",
        title = "Cosmological Surface-Brightness Dimming",
    )

    redshifts = [0.4 => :black, 0.5 => :blue, 1.0 => :green, 2.0 => :red]
    for (z, color) in redshifts
        wav_obs = wav_bb * (1 + z)
        lines!(ax, wav_obs, F_obs(wav_obs, z); label = "z = $(z)", color)
    end

    axislegend(; position = :rt)

    fig
end

# ╔═╡ aeba92bf-15de-425a-8f7e-b2ea0518756c
md"""
# Notebook setup 🔧
"""

# ╔═╡ d97007f0-5c6c-4644-bc2d-e035c9d851ac
function frontmatter()
    path = split(@__FILE__, "#==#") |> first |> string
    prefix = "#>"
    is_fm = startswith(prefix)
    block = Iterators.takewhile(is_fm, Iterators.dropwhile(!is_fm, eachline(path)))
    toml = TOML.parse(join((chopprefix(chopprefix(l, prefix), " ") for l in block), "\n"))
    return toml["frontmatter"]
end

# ╔═╡ 4668b8af-8ca2-4a51-a7f1-567aa73ee67d
title(_fm) = _fm["title"]

# ╔═╡ 2aa83581-6911-4fdc-b19e-95196c4651e0
function authors(_fm)
    names = [d["name"] for d in _fm["author"]]
    body = join(names, ", ")
    return md"_Authors: $(body)_"
end

# ╔═╡ 48a43023-5a00-4f93-bf97-fdf534f04e65
function keywords(_fm; kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = _fm["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ 4ca2f579-4240-40b4-a07c-29896c3684b4
begin
    _fm = frontmatter()

    md"""
    # $(title(_fm))

    $(authors(_fm))

    !!! tip "Learning goals"
    	Goals here.

    $(keywords(_fm))

    !!! warning "Companion content"
    	Content here.
    """
end

# ╔═╡ 89579410-1058-4b49-932c-c1715ce662a8
TableOfContents(; depth = 4)

# ╔═╡ Cell order:
# ╟─4ca2f579-4240-40b4-a07c-29896c3684b4
# ╟─29945a91-d0c5-470b-81e4-99329e874ee0
# ╟─5d3a54c7-97d8-48c8-af80-8b7e2ad13ff8
# ╠═d7ed3be9-00c1-444d-a49d-bb7ce7c5bf03
# ╟─e35fd5df-e01e-43af-99cf-546ca8c8f342
# ╠═645822c8-d735-4ec5-8db7-7cdbb3fd191b
# ╟─8b96f755-9cfd-479d-a32b-6a664004793b
# ╠═7758319f-471a-4d9d-9e31-103ee81fe9ab
# ╟─9a8dac0b-1ca0-43ff-b822-27690eae157c
# ╠═ae4d28da-8c69-481c-8589-149e51ddd9fa
# ╟─df8b47f8-7692-4032-a9e5-e99397f7808a
# ╠═143e0860-0704-437e-b2d3-22f5298a9198
# ╟─366d72e2-8f5f-4b8a-986b-467f662aaf5b
# ╠═f0833976-0cd9-45e3-ac8e-0a2ed1073d1f
# ╠═24edd839-07bb-477b-b059-730bfe57db5f
# ╠═c0b77309-56c6-46df-b4ab-360c179c1934
# ╠═0d5c2f5c-f57d-43ef-9633-94f3411aa977
# ╟─eae70e1e-7182-4257-aaa0-4b8f88c901ce
# ╠═43fad909-2733-48db-8c0b-cefd532ed267
# ╟─c8529662-9746-4d0f-bef5-796c33d09ac5
# ╠═6e6972a6-1648-40e4-a1e6-beed98bbc93e
# ╠═45849da1-b473-4d19-b140-d9d3fca39799
# ╠═e58714fd-8432-43c2-be11-8c9a22e2de85
# ╠═d42d14f8-0e3a-46d0-b04f-c7e735591aee
# ╠═eff1ff04-b875-4b01-b7c6-67b07c682396
# ╠═1ba509c8-a4fa-4aeb-b811-6e08259a6d18
# ╠═9866a403-9374-42db-8beb-57e9673e71b8
# ╠═9e824cb8-ff28-40e0-be4f-36e579fae709
# ╟─8d6e0d8b-bbea-4144-bf22-9a94c38b33e6
# ╠═07baa347-0ac7-43e6-ab67-4e99c3774f65
# ╠═57f114fb-3098-47d0-bb95-3e0f9fe41790
# ╠═97bf0849-d92e-47cb-914a-6f9a1d02ce94
# ╠═3d996ca3-26c6-4ebe-896e-c6b0378ae85c
# ╠═84052581-f53d-4bed-9776-1ee31d7c38fe
# ╠═48a59ba2-73bc-42dd-bec4-7be6b2aaff90
# ╠═b2ac804d-2ebd-42e7-8625-99a91e0458f6
# ╟─4147f9f6-df33-4ac1-b31f-17dab96c951f
# ╠═b742966f-a179-4a6a-93ed-fbabb3e3803e
# ╟─5b3a311f-ad28-44cf-bd9a-53bc8f6842b9
# ╠═61bfc289-e129-4bca-b9ec-36e39d91c151
# ╠═15601673-c10d-48aa-98c7-5fd389f35aa3
# ╠═cbf45b3f-8adf-4a40-bfc1-f73a4596e113
# ╟─aeba92bf-15de-425a-8f7e-b2ea0518756c
# ╟─d97007f0-5c6c-4644-bc2d-e035c9d851ac
# ╟─4668b8af-8ca2-4a51-a7f1-567aa73ee67d
# ╟─2aa83581-6911-4fdc-b19e-95196c4651e0
# ╟─48a43023-5a00-4f93-bf97-fdf534f04e65
# ╠═89579410-1058-4b49-932c-c1715ce662a8
# ╠═e4c5b28c-c4fb-4726-ab7a-8f1cc04280bd
