### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> image = "/assets/fits-images.png"
#> order = 2
#> title = "FITS images"
#> layout = "layout.jlhtml"
#> date = "2025-11-19"
#> description = "View and manipulate data from FITS images."
#> tags = ["file I/O", "FITS", "images", "image processing", "plots", "histograms", "colorbars"]

using Markdown
using InteractiveUtils

# ╔═╡ 61c0bf34-302b-4732-a44d-4c2da611eb74
begin
    # Can remove this block after Makie v0.25 is merged:
    # https://github.com/MakieOrg/Makie.jl/pull/5484
    import Pkg
    Pkg.activate(; temp = true)
    Pkg.add(
        [
            Pkg.PackageSpec(; name = "Downloads"),
            Pkg.PackageSpec(; name = "TOML"),
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(; name = "DataFramesMeta"),
            Pkg.PackageSpec(; name = "StatsBase"),
            Pkg.PackageSpec(; name = "Printf"),
            Pkg.PackageSpec(; name = "MathTeXEngine"),
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
                rev = "makie",
                url = "https://github.com/JuliaAstro/AstroImages.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaAstro/AstroAngles.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaAstro/FITSFiles.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaAstro/FITSWCS.jl",
            ),
            Pkg.PackageSpec(;
                rev = "makie-0.25",
                url = "https://github.com/icweaver/DimensionalData.jl",
            ),
        ]
    )

    # Analysis
    using StatsBase: Histogram, fit, mean, median, std

    # Data handling and visualization
    using DataFramesMeta: DataFrame, @rsubset!
    using Downloads: download
    using Printf: @sprintf
    using AstroImages
    using CairoMakie

    deps_ready = true # Can remove this too after upstream merge
end

# ╔═╡ 7d07caf5-e203-4152-8bb9-c1f396c4f80c
begin
    deps_ready

    using TOML: TOML
    using PlutoUI: TableOfContents
end

# ╔═╡ 91f00e98-e69c-4435-b9d0-10d30006efef
md"""
## Summary

This tutorial demonstrates how to use [AstroImages.jl](https://juliaastro.org/AstroImages) to view image data from FITS files and create publication-ready visualizations with different color scales and stretches. We will also demonstrate how to perform simple image stacking.
"""

# ╔═╡ a6e33cf8-1fe3-4810-a66b-adc07166871e
md"""
### Packages 📦
"""

# ╔═╡ 74f155a2-fe8a-404f-b0c6-7a4f2724c408
md"""
!!! tip
    Try adjusting the plot options above, or try adding your own! See the [Getting started](https://docs.makie.org/stable/tutorials/getting-started) section of the Makie.jl documentation for a comprehensive tutorial.

You may notice that there are still a few things missing from our plot that would be nice to have by default, e.g., labeled axes and a formatted colorbar. We will show an ergonomic way to do this next.
"""

# ╔═╡ 0598bfc5-0900-4b04-ad84-9638b89869b9
md"""
## Opening FITS files and viewing the image data

For this tutorial, we will work with an astronomical image of the Horsehead Nebula taken with a photographic plate. The image has been digitized, i.e., scanned by a computer and converted to a 2D array. Each position in the array corresponds with the projected position on the sky, and bright areas of the image have high values while dark areas have low values in the array.

Images taken with astronomical instruments called CCDs or "[charge-coupled devices](https://en.wikipedia.org/wiki/Charge-coupled_device)" are organized similarly. When illuminated by light, CCDs accumulate electrons, converting brightness values to electron counts. A CCD image is essentially a 2D array, where each position on the array represents a single CCD pixel, and the values in that array represent the number of counts registered in that pixel.
"""

# ╔═╡ 622ae936-36d2-4154-a82b-e8fe33773f92
img_file = download("http://data.astropy.org/tutorials/FITS-images/HorseHead.fits");

# ╔═╡ 3b813a3a-b24d-4801-bd1a-d7ac3106518b
img = AstroImage(img_file; scale = false)

# ╔═╡ d9881179-49b8-40af-8aa0-233d7bb1496c
md"""
!!! tip
     If you are in an interactive environment like VSCode, Jupyter, or Pluto (this tutorial), instead of a REPL, AstroImages are automatically rendered to images and displayed.
"""

# ╔═╡ a9ec676e-8d14-48af-b3ea-6482dce86832
md"""
Here is some information about our image:
"""

# ╔═╡ d413263f-8896-4060-bc6a-e849ef683a19
size(img) # NAXIS1 x NAXIS2 dimension

# ╔═╡ 798590cd-b1c9-44d7-90b7-71d470abce5a
eltype(img) # Data type of each pixel in our image

# ╔═╡ 3a5ec935-41f8-4845-bfeb-d6e5836353bf
header(img) # The stored header file

# ╔═╡ fc3001ac-dd95-467a-9be8-a581d9fef628
img.data # Access the underlying matrix data

# ╔═╡ a9548ae8-ab0e-42be-825d-8146b95196a2
md"""
Under the hood, AstroImages.jl calls [FITSFiles.jl](https://juliaastro.org/FITSFiles) to load our data for us. By default, this takes the first HDU to be our image data, but we could also load all HDUs by passing `:` to the above `load` function:
"""

# ╔═╡ 4cb7062f-4d86-4ac6-884e-cd7546ef5b90
load(img_file, :)

# ╔═╡ ee2c3044-9c44-4a49-a0c8-2ff77ee5a83d
md"""
We see that in addition to the first block holding our image data, we also have a second block containing table data. AstroImages.jl is specialized for image data, so for more information on accessing the table information, it is recommended to use the lower-level FITSFiles.jl package as shown in the first tutorial ([learn.JuliaAstro > Working with FITS tables](/tutorials/fits-tables/)).
"""

# ╔═╡ 722db014-6b52-40d9-8183-3fc6bc067937
md"""
### `imview`

Before producing a final image for publication, we can use the `imview` function provided by AstroImages.jl to apply different transformations to our image data, displayed directly in the notebook. This is what is called automatically when displaying AstroImage objects in this notebook. Try adjusting the different available options below to modify the image:
"""

# ╔═╡ 647115c6-856a-47a8-996d-29c9e2f2cde5
imview(
    img;
    clims = Percent(99.5),
    stretch = identity,
    cmap = :magma,
    contrast = 1.0,
    bias = 0.5,
)

# ╔═╡ 8cd4bbda-c0d0-4694-ae9c-f3282d1ccc9c
md"""
See the [`imview`](https://juliaastro.org/AstroImages.jl/dev/api/#AstroImages.imview) documentation for more.
"""

# ╔═╡ b3c6961a-b8a6-4597-aaf7-a97cae793670
md"""
### Getting basic statistics

`AstroImage` objects behave much like regular arrays, which allows them to support the usual array options and statistical computations:
"""

# ╔═╡ 010ddf8b-7355-41cf-8d92-0d217e382978
img[1:10, 1:10]

# ╔═╡ 22824a6b-9a68-4519-86af-2c5b60da5288
extrema(img) # Or (minimum(img), maximum(img))

# ╔═╡ 2f151bab-e310-46a2-b277-3082e880360b
mean(img)

# ╔═╡ 45b56d0c-71b7-4151-95b5-bcd67239f40d
median(img)

# ╔═╡ c3011e6e-6590-4b88-b3b3-a03c768a82c6
std(img)

# ╔═╡ d0d4bea0-24d8-490d-850f-59aa6b760e3d
md"""
### Plotting a histogram

Let's flatten our image with `vec` next to take a look at its histogram:
"""

# ╔═╡ 0fa6d28a-4608-444b-ae2c-5845ec8a21b1
stephist(
    vec(img);
    bins = 50,
    axis = (; xlabel = "Pixel value", ylabel = "Counts")
)

# ╔═╡ 5dd43020-71e9-47e7-9b57-295e57a98bce
md"""
### Plotting with Makie.jl

AstroImages.jl uses Makie.jl as its plotting backend, which we used directly to produce the histogram plot above. Makie.jl is a modern plotting ecosystem written in Julia. Its [set of backends](https://docs.makie.org/stable/explanations/backends/backends#What-is-a-backend) allows us to produce plots for a wide range of contexts. For this tutorial, we use the CairoMakie.jl backend to produce publication-quality vector graphic plots.

Let's start with the highlevel plotting function `implotview` provided by `AstroImages`.jl:
"""

# ╔═╡ cd852bfc-5996-4ac3-83f4-19457a9f0825
implotview(img)

# ╔═╡ 2dddbd57-e15d-47ad-ab5d-282cf970f75a
md"""
!!! note
    See <https://juliaastro.org/AstroImages/stable/manual/dimensions-and-world-coordinates/> for more
"""

# ╔═╡ 9b1ec23d-e773-4893-83ba-c465fd510778
md"""
`AstroImage` objects use the [DimensionalData.jl](https://rafaqz.github.io/DimensionalData.jl/) interface, which allows them to participate in extended array handling and plotting functionality:
"""

# ╔═╡ 0e690360-cdb1-47c9-b637-b0184508ac03
img[X = 500, Y = Near(500.1)] == img[500, 500]

# ╔═╡ 6500f709-d767-44dd-911f-20fb29740671
fig, ax, p = plot(
    img;
    colorscale = log10, # log scale the colors
    colormap = :greys,
    colorbar = (
        ticks = [4.0e3, 5.0e3, 6.0e4, 1.0e4, 2.0e4], # Based on histogram
        minorticksvisible = true,
        minorticks = IntervalsBetween(9),
    )
)

# ╔═╡ e2e9b225-88c2-4297-bcc7-fe31b1a8ca9f
md"""
For more on working with AstroImage data, see the [Getting started](https://juliaastro.org/AstroImages/stable/manual/getting-started/) section of the AstroImages.jl manual.

We end by looking at a brief image stacking example.
"""

# ╔═╡ 5fb6d8fa-d890-45c5-afa3-944e96bc818e
md"""
## Image stacking

For this example, we'll stack several images of M13 taken with a ~10" telescope.
"""

# ╔═╡ 4e9d0b96-358d-4baa-8ea9-7eb5aeff36f8
md"""
### Load data

Let's start by opening a series of FITS files and storing the data in a vector called `imgs`:
"""

# ╔═╡ e2defd90-7781-4c64-98ad-4aabdfa99d63
# We use the @sprintf macro from the base Printf.jl Julia module
# to format our strings
fpaths = map(1:5) do i
    @sprintf("http://data.astropy.org/tutorials/FITS-images/M13_blue_%04d.fits", i)
end

# ╔═╡ 215e0183-adfa-4da0-80d7-108894f73f23
md"""
!!! tip
    This is just an alternative syntax to array comprehensions. We could have just as easily done:

    ```julia
    [
        @sprintf(
            "http://data.astropy.org/tutorials/FITS-images/M13_blue_%04d.fits",
            i
        )
        for i in 1:5
    ]
    ```

    See [this section of the Julia manual](https://docs.julialang.org/en/v1/manual/functions/#Do-Block-Syntax-for-Function-Arguments) for more on do-block syntax.
"""

# ╔═╡ 4dc4074e-b3d3-4deb-a657-8a7847a8c156
imgs = [load(download(fpath)) for fpath in fpaths];

# ╔═╡ 51f3adcf-4e9a-49b5-a3a1-162a132e7f68
md"""
### Visualize

We can then directly sum up this vector of image data to produce a stacked image:
"""

# ╔═╡ 1cc2e29e-e9b3-4313-86b2-9eb7032c0ba4
img_stacked = sum(imgs)

# ╔═╡ 134e86e5-995e-432c-87a5-bf6a496fe7f6
md"""
As in the examples above, this does the usual image transformations defined in `imview` to produce a nice image by default. We can also directly control this in the same way:
"""

# ╔═╡ 9c5450c2-807f-438d-a907-d78090a631ab
stephist(vec(img_stacked); bins = 50)

# ╔═╡ 21bd0738-97cb-499d-ad8d-b6301c6a0bdc
md"""
The pixel values looks to be mostly around [2000, 3000] counts, so we set our colorbar limits there:
"""

# ╔═╡ 1be1886c-122f-4ebf-bcbf-0bd197aadd98
implotview(img_stacked; clims = (2.0e3, 3.0e3))

# ╔═╡ eccd1973-1082-4fb6-920c-0039ddf5482a
md"""
### Save

Finally, we can save our underlying stacked image + header data using the `save` function exported from AstroImages.jl:

```julia
using AstroImages

save("test.fits", img_stacked)
```
"""

# ╔═╡ c65018aa-e30a-4727-ad4e-b853a1479a40
md"""
# Notebook setup 🔧
"""

# ╔═╡ 89e8f2a6-9d3b-44b8-8805-91daa24124c3
TableOfContents(; depth = 4)

# ╔═╡ c61c61a3-9582-4680-95e5-59a6b818ef9d
function frontmatter(path)
    prefix = "#> "
    is_fm = startswith(prefix)
    block = Iterators.takewhile(is_fm, Iterators.dropwhile(!is_fm, eachline(path)))
    toml = TOML.parse(join(chopprefix.(block, prefix), "\n"))
    return toml["frontmatter"]
end

# ╔═╡ e892e25d-06a0-496a-bf63-40a8a988d089
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ 3c48207e-ae5d-4597-8010-587d6ed8736b
md"""
# Working with FITS images

This notebook is modified from <https://learn.astropy.org/tutorials/FITS-images.html>

_Original authors: Lia Corrales, Kris Stern, Stephanie T. Douglas, Kelle Cruz, Lúthien Liu, Zihao Chen, Saima Siddiqui_

!!! tip "Learning goals"
    1. Open FITS files and load image data
    1. Make a 2D histogram with image data
    1. Stack several images into a single image
    1. Write image data to a FITS fil

$(keywords())

!!! warning "Companion content"
    [learn.JuliaAstro > Working with FITS tables](/tutorials/fits-tables/)
"""

# ╔═╡ Cell order:
# ╟─3c48207e-ae5d-4597-8010-587d6ed8736b
# ╟─91f00e98-e69c-4435-b9d0-10d30006efef
# ╟─a6e33cf8-1fe3-4810-a66b-adc07166871e
# ╠═61c0bf34-302b-4732-a44d-4c2da611eb74
# ╟─74f155a2-fe8a-404f-b0c6-7a4f2724c408
# ╟─0598bfc5-0900-4b04-ad84-9638b89869b9
# ╠═622ae936-36d2-4154-a82b-e8fe33773f92
# ╠═3b813a3a-b24d-4801-bd1a-d7ac3106518b
# ╟─d9881179-49b8-40af-8aa0-233d7bb1496c
# ╟─a9ec676e-8d14-48af-b3ea-6482dce86832
# ╠═d413263f-8896-4060-bc6a-e849ef683a19
# ╠═798590cd-b1c9-44d7-90b7-71d470abce5a
# ╠═3a5ec935-41f8-4845-bfeb-d6e5836353bf
# ╠═fc3001ac-dd95-467a-9be8-a581d9fef628
# ╟─a9548ae8-ab0e-42be-825d-8146b95196a2
# ╠═4cb7062f-4d86-4ac6-884e-cd7546ef5b90
# ╟─ee2c3044-9c44-4a49-a0c8-2ff77ee5a83d
# ╟─722db014-6b52-40d9-8183-3fc6bc067937
# ╠═647115c6-856a-47a8-996d-29c9e2f2cde5
# ╟─8cd4bbda-c0d0-4694-ae9c-f3282d1ccc9c
# ╟─b3c6961a-b8a6-4597-aaf7-a97cae793670
# ╠═010ddf8b-7355-41cf-8d92-0d217e382978
# ╠═22824a6b-9a68-4519-86af-2c5b60da5288
# ╠═2f151bab-e310-46a2-b277-3082e880360b
# ╠═45b56d0c-71b7-4151-95b5-bcd67239f40d
# ╠═c3011e6e-6590-4b88-b3b3-a03c768a82c6
# ╟─d0d4bea0-24d8-490d-850f-59aa6b760e3d
# ╠═0fa6d28a-4608-444b-ae2c-5845ec8a21b1
# ╟─5dd43020-71e9-47e7-9b57-295e57a98bce
# ╠═cd852bfc-5996-4ac3-83f4-19457a9f0825
# ╟─2dddbd57-e15d-47ad-ab5d-282cf970f75a
# ╟─9b1ec23d-e773-4893-83ba-c465fd510778
# ╠═0e690360-cdb1-47c9-b637-b0184508ac03
# ╠═6500f709-d767-44dd-911f-20fb29740671
# ╟─e2e9b225-88c2-4297-bcc7-fe31b1a8ca9f
# ╟─5fb6d8fa-d890-45c5-afa3-944e96bc818e
# ╟─4e9d0b96-358d-4baa-8ea9-7eb5aeff36f8
# ╠═e2defd90-7781-4c64-98ad-4aabdfa99d63
# ╟─215e0183-adfa-4da0-80d7-108894f73f23
# ╠═4dc4074e-b3d3-4deb-a657-8a7847a8c156
# ╟─51f3adcf-4e9a-49b5-a3a1-162a132e7f68
# ╠═1cc2e29e-e9b3-4313-86b2-9eb7032c0ba4
# ╟─134e86e5-995e-432c-87a5-bf6a496fe7f6
# ╠═9c5450c2-807f-438d-a907-d78090a631ab
# ╟─21bd0738-97cb-499d-ad8d-b6301c6a0bdc
# ╠═1be1886c-122f-4ebf-bcbf-0bd197aadd98
# ╟─eccd1973-1082-4fb6-920c-0039ddf5482a
# ╟─c65018aa-e30a-4727-ad4e-b853a1479a40
# ╠═89e8f2a6-9d3b-44b8-8805-91daa24124c3
# ╟─c61c61a3-9582-4680-95e5-59a6b818ef9d
# ╟─e892e25d-06a0-496a-bf63-40a8a988d089
# ╠═7d07caf5-e203-4152-8bb9-c1f396c4f80c
