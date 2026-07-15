### A Pluto.jl notebook ###
# v1.0.3

#> [frontmatter]
#> order = 2
#> title = "Working with Celestial Coordinates in WCS 1"
#> layout = "layout.jlhtml"
#> date = "2026-07-08"
#> description = "Work with images containing WCS information."
#> tags = ["wcs", "coordinates", "units", "file I/O", "plotting"]

using Markdown
using InteractiveUtils

# ╔═╡ b7045a18-7aa0-11f1-a432-3f9f2d689dd8
begin
    import Pkg

    Pkg.activate(; temp = true)

    Pkg.add(
        [
            Pkg.PackageSpec(; name = "Downloads"),
            Pkg.PackageSpec(; name = "TOML"),
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(; url = "https://github.com/JuliaAstro/FITSFiles.jl"),
            Pkg.PackageSpec(; url = "https://github.com/JuliaAstro/FITSWCS.jl"),
            Pkg.PackageSpec(; url = "https://github.com/JuliaAstro/AstroAngles.jl"),
            Pkg.PackageSpec(;
                rev = "makie",
                url = "https://github.com/JuliaAstro/AstroImages.jl",
            ),
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
                url = "https://github.com/icweaver/DimensionalData.jl",
                rev = "makie-0.25",
            ),
        ]
    )

    # Analysis
    using FITSWCS

    # # Data handling and visualization
    using Downloads: download
    using AstroImages, CairoMakie

    AstroImages.set_cmap!(:cividis)

    deps_ready = true
end;

# ╔═╡ 4d5ec068-f75d-459e-97c6-12f2ebdc337c
begin
    deps_ready

    using TOML: TOML
    using PlutoUI: TableOfContents, details
end

# ╔═╡ fa245c6e-9a32-4600-b579-6a42c0f0fc3a
md"""
## Summary

This tutorial series aims to show how the content of Chapter 1 of "An Introduction to Modern Astrophysics" by Carroll and Ostlie can be applied to real life astrophysics research situations, using tools in the JuliaAstro ecosystem. We will introduce two different approaches to building a `FITSWCS.WCSTransform` object, which contains meta-data that (in this case) defines a mapping between image coordinates and sky coordinates.

The FITSWCS.jl package conforms to the standards of the FITS World Coordinate System (WCS) used extensively by the astronomy research community. We will created a 2D WCS for an image of the iconic the Helix nebula (a planetary nebula) and display an image of the nebula with sky coordinates (here, equatorial, ICRS RA and Dec.) labeled. Finally, we will over-plot a scale bar on the Helix nebula image using WCS to give the reader a sense of the angular size of the image.
"""

# ╔═╡ c6e9643c-7b8e-4de8-bcf0-baf9a51910f2
md"""
## Packages 📦
"""

# ╔═╡ 163c21e5-ad94-4d25-aac0-c06c213dc95d
md"""
## Section 1: Two ways to create a WCS object

*World coordinates* serve to locate a measurement in some multi-dimensional parameter space. A World Coordinate System (WCS) specifies the physical, or world, coordinates to be attached to each pixel or voxel of an N-dimensional image or array. An [elaborate set of standards and conventions](https://fits.gsfc.nasa.gov/fits_wcs.html) have been developed for the Flexible Image Transport System (FITS) format ([Wells et al. 1981](https://ui.adsabs.harvard.edu/abs/1981A&AS...44..363W/abstract)). A typical WCS example is to specify the Right Ascension (RA) and Declination (Dec) on the sky associated with a given pixel or spaxel location in a 2-dimensional celestial image ([Greisen & Calabretta 2002](https://ui.adsabs.harvard.edu/abs/2002A&A...395.1061G/abstract); [Calabretta and Greisen 2002](https://ui.adsabs.harvard.edu/abs/2002A&A...395.1077C/abstract)).
"""

# ╔═╡ 09b72d31-629c-4be8-a8ef-911ffd23c4d5
md"""
The [FITSWCS.jl package](https://github.com/JuliaAstro/FITSWCS.jl) implements FITS standards and conventions for World Coordinate Systems. Using the `FITSWCS.WCSTransform` object and Makie.jl, we can generate images of the sky that have axes labeled with coordinates such as right ascension (RA) and declination (Dec).
"""

# ╔═╡ f6ca5f4a-1d3c-4f5c-91ba-93d2b95e2cc0
md"""
There are two main ways to initialize a `WCS` object: with a dictionary (or dictionary-like object, like a FITS file header) or with vectors. In this set of examples, we will initialize a `FITSWCS.WCSTransform` object with two dimensions, as would be needed to represent an image.
"""

# ╔═╡ c0f4fcb8-dcbb-40a4-b57f-54c9ccfd8c5c
md"""
The WCS standard defines a set of keywords that are used to represent the world coordinate system for a given set of data (e.g., image). Here is a list of the essential WCS keywords and their uses. In each case, the integer ``n`` denotes the dimensional axis (starting with 1) to which the keyword is being applied. In our examples below, we will have two image dimensions (axes), so ``n`` will either be 1 or 2.
"""

# ╔═╡ 15a26fa5-24fe-4664-8300-ad8f806e3d94
md"""
- **CRVALn**: the coordinate value at a reference point (e.g., RA and DEC value in degrees)
- **CRPIXn**: the pixel location of the reference point (e.g., CRPIX1=1, CRPIX2=1 describes the center of a corner pixel)
- **CDELTn**: the coordinate increment at the reference point (e.g., the difference in 'RA' value from the reference pixel to its neighbor along the RA axis)
- **CTYPEn**: an 8-character string describing the axis type (e.g., 'RA---TAN' and 'DEC---TAN' describe the typical tangent-plane sky projection that astronomers use)
- **CUNITn**: a string describing the unit for each axis (if not specified, the default unit is degrees.)
- **NAXISn**: an integer defining the number of pixels in each axis


Some good references of the WCS standard can be found [here](https://fits.gsfc.nasa.gov/fits_wcs.html).
"""

# ╔═╡ ff7e4376-68df-45c8-899d-ebddbba9e708
md"""
### Method 1: Building a WCS object with a dictionary

One way to define a WCS object is to construct a dictionary containing all of the essential information (i.e., specifying values for the WCS keywords listed above) that map the pixel coordinate space to the world coordinate space.

In this example, we define two coordinate axes with:

- A two dimensional image (**NAXIS = 2**)
- A Gnomonic (tangent-plane, **CTYPEn**) projection, which corresponds to the RA/Dec coordinate system
- A reference location of (RA,DEC) = (337.52, -20.83), as defined by the **CRVALn** keys
- The pixel at coordinate value (1,1) as the reference location (**CRPIXn** keys)
- Units of degrees (**CUNITn = "deg"**)
- Pixel sizes of 1 x 1 arcsec (**CDELTn = 0.002778** in degrees)
"""

# ╔═╡ c0047c6e-28c1-41e7-861e-12e7bd4a7062
wcs_helix_dict = Dict(
    "NAXIS" => 2,
    "CTYPE1" => "RA---TAN",
    "CUNIT1" => "deg",
    "CDELT1" => -0.0002777777778,
    "CRPIX1" => 1,
    "CRVAL1" => 337.5202808,
    "CTYPE2" => "DEC--TAN",
    "CUNIT2" => "deg",
    "CDELT2" => 0.0002777777778,
    "CRPIX2" => 1,
    "CRVAL2" => -20.833333059999998,
) |> WCS

# ╔═╡ 134f9f2a-63c2-4ab2-ad45-2bdaae1eaac3
md"""
!!! note
	Instead of storing the size with `NAXIS1 = ...`, `NAXIS2 = ...`, we store the dimensionality with `NAXIS => 2`
"""

# ╔═╡ ebe11f46-8f9b-498a-83ed-ad73c868f463
md"""
### Method 2: Create an empty WCS object before assigning values

Alternatively, we could construct the `FITSWCS.WCSTransform` object with keyword values corresponding to each respective axis:
"""

# ╔═╡ b4e2cc78-6327-47ee-86c6-8869ec1c2a88
WCS(
    2;
    crpix = [1, 1],
    crval = [337.5202808, -20.833333059999998],
    cunit = ["deg", "deg"],
    ctype = ["RA---TAN", "DEC--TAN"],
    cdelt = [-0.0002777777778, 0.0002777777778],
)

# ╔═╡ ca0aa269-b6ac-4804-86da-408f988e7d99
md"""
## Section 2: Show an image of the Helix nebula with RA and Dec labeled

Most of the time we can obtain the required `FITSWCS.WCSTransform` object from the header of the FITS file from a telescope or astronomical database. This process is described below.
"""

# ╔═╡ e12e893b-7ff6-4c25-89d1-050c20eeb1c0
md"""
### Step 1: Read in the FITS file

We will read the FITS file containing an image of the Helix nebula from the astropy-data GitHub repository using the [AstroImages.jl](https://juliaastro.org/AstroImages/) package:
"""

# ╔═╡ 4eea74ab-53a4-4e7c-8505-b1aefe678ef5
img = let
    url = "https://github.com/astropy/astropy-data/raw/6d92878d18e970ce6497b70a9253f65c925978bf/tutorials/celestial-coords1/tailored_dss.22.29.38.50-20.50.13_60arcmin.fits"
    f = joinpath(@__DIR__, "data", "tailored_dss.22.29.38.50-20.50.13_60arcmin.fits")
    isfile(f) ? load(f) : load(download(url, f))
end;

# ╔═╡ 1fd4e3ce-412a-4bd2-b649-f35f090bfa60
md"""
!!! note
	This image (FITS file) was originally accessed from the Digitized Sky Survey but is provided in the astropy-data repository for convenience.
"""

# ╔═╡ f4b32c4e-f289-4ee0-bde3-44dbb0387e32
md"""
FITS files are a binary file format that is mainly used by astronomers and can contain information arranged in many “extensions,” which contain both header information (e.g., metadata) and data (e.g., image data). We can check how many extensions there are in a FITS file, as well as view a summary of the contents in each extension, by printing the FITS object information:
"""

# ╔═╡ 91d60255-18d6-47f7-8f66-09b9cd397602
h = header(img)

# ╔═╡ 7bae002b-d32f-4c36-a83c-020b530f7431
md"""
Please note that the original header (as downloaded from the DSS) violates the FITS WCS standards (because it includes both CDELTn keywords and a matrix of CD values; including deprecated PC-matrix keywords). The header has been cleaned up to conform to the existing standards.
"""

# ╔═╡ 5448300e-1726-4fa4-a390-df5cac27549f
md"""
### Step 2: Read in the FITS image coordinate system

Because the header contains WCS information, we can construct a `WCSTransform` object directly from the FITS header:
"""

# ╔═╡ afa923fe-7023-43b4-9c10-bc8c2d21acc1
wcs_helix = WCS(h)

# ╔═╡ 8df60e64-fa50-4bf7-8eb5-4563ad15b016
md"""
!!! note
    This differs slightly from astropy and WCS.jl, which both use the underlying wcslib to detect the legacy DSS keywords in the header and rebuild the WCS from the plate polynomial.

    In practice, these both map to virtually the same pixel (within about 1e-7°), only the reference point differs.
"""

# ╔═╡ 17dfadfb-1c0e-4233-90ae-1c1f592ac564
md"""
### Step 3: Plot the Helix nebula with sky coordinate axes (RA and Dec)

The image data, `img`, is an `AstroImages.AstroImage`, containing WCS information within its header that can automatically be plotted along with its bundled 2D image data:
"""

# ╔═╡ c167fb7f-5ab6-4b32-a2c5-beeefa1ae2c8
implotview(img; gridcolor = :coral)

# ╔═╡ 5a3f62b9-0d39-46f5-b763-49244c13d4f0
md"""
!!! todo
	Switch to Makie.jl

	```julia
    fig, ax, p = image(img; colormap = :cividis)

    colsize!(fig.layout, 1, Aspect(1, size(img, 1) / size(img, 2)))

    resize_to_layout!(fig)

    fig
	```

    and add a recipe that replaces `implot`.
"""

# ╔═╡ 656bf1f9-eaee-46a4-9f44-422a384a19ee
md"""
### Exercise

Copy the code block above and instead overlay a coordinate grid in Galactic coordinates.
"""

# ╔═╡ 3b3ceefb-8abf-4bf5-857e-f0ca0780d960
md"""
!!! todo
	Implement Reproject.jl: <https://juliaastro.org/AstroImages/stable/guide/reproject/>
"""

# ╔═╡ bfe22652-e591-43db-9a2d-0de983e18f5e
md"""
## Section 3: Plot a scale marker on an image with WCS

!!! todo
	Implement this once we've switched over to Makie.jl as backend.
"""

# ╔═╡ a7df5632-ec33-41d7-92eb-e68fcb6d1f4c
md"""
### Exercise

Make a horizontal bar with the same length. Keep in mind that 1 hour angle = 15 degrees.
"""

# ╔═╡ 1321b8db-291c-4ebf-a5d0-ac254629e950
md"""
# Notebook setup 🔧
"""

# ╔═╡ eee39591-e860-4c91-bed4-e05e5dfaf0de
TableOfContents()

# ╔═╡ a8cb2e7a-790b-4896-b49a-e9fa6b98bffa
function frontmatter(path)
    prefix = "#> "
    is_fm = startswith(prefix)
    block = Iterators.takewhile(is_fm, Iterators.dropwhile(!is_fm, eachline(path)))
    toml = TOML.parse(join(chopprefix.(block, prefix), "\n"))
    return toml["frontmatter"]
end

# ╔═╡ b18be7a0-a053-4d56-8713-3eb473044ac2
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ bb93ae12-1f64-474d-9762-a4c9813ab801
md"""
# Working with Celestial Coordinates in WCS 1: Specifying, reading, and plotting

This notebook is modified from <https://learn.astropy.org/tutorials/FITS-images.html>

_Original authors: Kris Stern, Kelle Cruz, Lia Corrales, David Shupe, Adrian Price-Whelan_

!!! tip "Learning goals"
    1. Demonstrate two ways to build a `FITSWCS.WCSTransform` object
    1. Show an image of the Helix nebula with RA and Dec labeled
    1. Plot a scale bar on an image with WCS information

$(keywords())

!!! warning "Companion content"
    1. "An Introduction to Modern Astrophysics" ([Carroll & Ostlie](https://ui.adsabs.harvard.edu/abs/2006ima..book.....C/abstract))
    1. [FITS WCS page at GSFC](https://fits.gsfc.nasa.gov/fits_wcs.html)
    1. [learn.JuliaAstro > FITS Images](https://learn.juliaastro.org/tutorials/fileio-fits_images/)
"""

# ╔═╡ Cell order:
# ╟─bb93ae12-1f64-474d-9762-a4c9813ab801
# ╟─fa245c6e-9a32-4600-b579-6a42c0f0fc3a
# ╟─c6e9643c-7b8e-4de8-bcf0-baf9a51910f2
# ╠═b7045a18-7aa0-11f1-a432-3f9f2d689dd8
# ╟─163c21e5-ad94-4d25-aac0-c06c213dc95d
# ╟─09b72d31-629c-4be8-a8ef-911ffd23c4d5
# ╟─f6ca5f4a-1d3c-4f5c-91ba-93d2b95e2cc0
# ╟─c0f4fcb8-dcbb-40a4-b57f-54c9ccfd8c5c
# ╟─15a26fa5-24fe-4664-8300-ad8f806e3d94
# ╟─ff7e4376-68df-45c8-899d-ebddbba9e708
# ╠═c0047c6e-28c1-41e7-861e-12e7bd4a7062
# ╟─134f9f2a-63c2-4ab2-ad45-2bdaae1eaac3
# ╟─ebe11f46-8f9b-498a-83ed-ad73c868f463
# ╠═b4e2cc78-6327-47ee-86c6-8869ec1c2a88
# ╟─ca0aa269-b6ac-4804-86da-408f988e7d99
# ╟─e12e893b-7ff6-4c25-89d1-050c20eeb1c0
# ╠═4eea74ab-53a4-4e7c-8505-b1aefe678ef5
# ╟─1fd4e3ce-412a-4bd2-b649-f35f090bfa60
# ╟─f4b32c4e-f289-4ee0-bde3-44dbb0387e32
# ╠═91d60255-18d6-47f7-8f66-09b9cd397602
# ╟─7bae002b-d32f-4c36-a83c-020b530f7431
# ╟─5448300e-1726-4fa4-a390-df5cac27549f
# ╠═afa923fe-7023-43b4-9c10-bc8c2d21acc1
# ╟─8df60e64-fa50-4bf7-8eb5-4563ad15b016
# ╟─17dfadfb-1c0e-4233-90ae-1c1f592ac564
# ╠═c167fb7f-5ab6-4b32-a2c5-beeefa1ae2c8
# ╟─5a3f62b9-0d39-46f5-b763-49244c13d4f0
# ╟─656bf1f9-eaee-46a4-9f44-422a384a19ee
# ╟─3b3ceefb-8abf-4bf5-857e-f0ca0780d960
# ╟─bfe22652-e591-43db-9a2d-0de983e18f5e
# ╟─a7df5632-ec33-41d7-92eb-e68fcb6d1f4c
# ╟─1321b8db-291c-4ebf-a5d0-ac254629e950
# ╠═eee39591-e860-4c91-bed4-e05e5dfaf0de
# ╟─b18be7a0-a053-4d56-8713-3eb473044ac2
# ╟─a8cb2e7a-790b-4896-b49a-e9fa6b98bffa
# ╠═4d5ec068-f75d-459e-97c6-12f2ebdc337c
