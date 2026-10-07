### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> tags = ["file I/O", "FITS", "images", "image processing", "plots", "colorbars", "data cubes", "radio astronomy", "WCS"]
#> title = "FITS Cubes"
#> description = "View and manipulate FITS cubes"
#> date = "2026-07-12"
#> layout = "layout.jlhtml"
#> order = 2

using Markdown
using InteractiveUtils

# ╔═╡ 61c0bf34-302b-4732-a44d-4c2da611eb74
begin
    import Pkg
    Pkg.activate(; temp = true)
    Pkg.add(
        [
            Pkg.PackageSpec(; name = "Downloads"),
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(; rev = "main", url = "https://github.com/JuliaAstro/FITSFiles.jl"),
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
            Pkg.PackageSpec(;
                rev = "docs-spectral",
                url = "https://github.com/JuliaAstro/AstroImages.jl",
            ),
            Pkg.PackageSpec(;
                rev = "compat/makie-v0.25",
                # url = "https://github.com/JuliaAstro/SpectrumBase.jl",
                path = "../../../SpectrumBase.jl/",
            ),

        ]
    )

    using Downloads: download
    using FITSFiles
    using SpectrumBase
    using CairoMakie
    using AstroImages

    deps_ready = true
end

# ╔═╡ 7d07caf5-e203-4152-8bb9-c1f396c4f80c
begin
    deps_ready

    using PlutoUI: TableOfContents
end

# ╔═╡ 3c48207e-ae5d-4597-8010-587d6ed8736b
md"""
# Working with FITS cubes

This notebook is modified from <https://learn.astropy.org/tutorials/FITS-cubes.html>

_Original authors: Dhanesh Krishnarao (DK), Shravan Shetty, Diego Gonzalez-Casanova, Audra Hernandez, Kris Stern, Kelle Cruz, Stephanie Douglas_

!!! tip "Learning goals"
    - Find and download data
    - Read and plot slices across different dimensions of a data cube
    - Compare different data sets (2D and 3D) by overploting contours
    - Transform coordinate projections and match data resolutions
    - Create intensity moment maps / velocity maps

!!! warning "Companion content"
    - <https://juliaastro.org/AstroImages.jl/previews/PR125/manual/spec/>
"""

# ╔═╡ 91f00e98-e69c-4435-b9d0-10d30006efef
md"""
## Summary

In this tutorial we will visualize 2D and 3D data sets in Galactic and equatorial coordinates.

The tutorial will walk you though a visual analysis of the Small Magellanic Cloud (SMC) using HI 21cm emission and a Herschel 250 micron map. We will learn how to read in data from a file, query and download matching data from Herschel, and plot the resulting images.
"""

# ╔═╡ a6e33cf8-1fe3-4810-a66b-adc07166871e
md"""
### Packages 📦
"""

# ╔═╡ 8cc6c458-5737-495b-bc04-895aa148661d
hi_datafile = download("http://data.astropy.org/tutorials/FITS-cubes/reduced_TAN_C14.fits")

# ╔═╡ 433fbef7-91de-4a59-aa3c-503141f858ed
md"""
!!! note "Data management"
    For caching, auto-updates, and more, see this (non-exhaustive) list of packages:

    - [RemoteFiles.jl](https://github.com/helgee/RemoteFiles.jl)
    - [DataDeps.jl](https://github.com/oxinabox/DataDeps.jl)
    - [DataToolkit.jl](https://github.com/tecosaur/DataToolkit.jl)
"""

# ╔═╡ 93125bb2-efe7-4869-a459-10a7760bd353
hi_data = fits(hi_datafile)

# ╔═╡ 012eeb56-1442-4f13-876e-5d0b193c88fc
cube = load(hi_datafile)

# ╔═╡ 8986633d-e8b3-4ebb-a04a-c547267ab994
cube[Z = 300] # Equivalent to cube[:, :, 300]

# ╔═╡ 016c4e70-121c-45e0-b984-bde2413288ce
cube[Z = 300] |> implotview

# ╔═╡ 0614a265-5d93-4016-9676-5c107bf48891
vels = [pixel_to_world(cube, [75, 75, z])[3] for z in 1:size(cube, 3)] ./ 1.0e3  # km/s

# ╔═╡ 18b5864a-b945-4b8f-8518-36db9dc18f4c
lines(vels, collect(cube[X = 75, Y = 75])) # km/s vs. K

# ╔═╡ 56640bd0-9302-4bfa-a398-33cceda0fb54
md"""
!!! todo
    Find a nicer way to do this, maybe with SpectrumBase.jl
"""

# ╔═╡ c65018aa-e30a-4727-ad4e-b853a1479a40
md"""
# Notebook setup 🔧
"""

# ╔═╡ 89e8f2a6-9d3b-44b8-8805-91daa24124c3
TableOfContents(; depth = 4)

# ╔═╡ Cell order:
# ╟─3c48207e-ae5d-4597-8010-587d6ed8736b
# ╟─91f00e98-e69c-4435-b9d0-10d30006efef
# ╟─a6e33cf8-1fe3-4810-a66b-adc07166871e
# ╠═61c0bf34-302b-4732-a44d-4c2da611eb74
# ╠═8cc6c458-5737-495b-bc04-895aa148661d
# ╟─433fbef7-91de-4a59-aa3c-503141f858ed
# ╠═93125bb2-efe7-4869-a459-10a7760bd353
# ╠═012eeb56-1442-4f13-876e-5d0b193c88fc
# ╠═8986633d-e8b3-4ebb-a04a-c547267ab994
# ╠═016c4e70-121c-45e0-b984-bde2413288ce
# ╠═0614a265-5d93-4016-9676-5c107bf48891
# ╠═18b5864a-b945-4b8f-8518-36db9dc18f4c
# ╟─56640bd0-9302-4bfa-a398-33cceda0fb54
# ╟─c65018aa-e30a-4727-ad4e-b853a1479a40
# ╠═89e8f2a6-9d3b-44b8-8805-91daa24124c3
# ╟─7d07caf5-e203-4152-8bb9-c1f396c4f80c
