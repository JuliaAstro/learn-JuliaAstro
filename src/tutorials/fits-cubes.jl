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
            Pkg.PackageSpec(; name = "CodecZlib"),
            Pkg.PackageSpec(; name = "Contour"),
            Pkg.PackageSpec(; name = "CSV"),
            Pkg.PackageSpec(; name = "DataFrames"),
            Pkg.PackageSpec(; name = "Downloads"),
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(; name = "Tar"),
            Pkg.PackageSpec(; name = "URIs"),
            Pkg.PackageSpec(;
                rev = "hdu-size",
                url = "https://github.com/JuliaAstro/FITSFiles.jl",
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
                rev = "reproject",
                url = "https://github.com/JuliaAstro/AstroImages.jl",
            ),
            Pkg.PackageSpec(;
                rev = "compat/makie-v0.25",
                url = "https://github.com/JuliaAstro/SpectrumBase.jl",
            ),
            Pkg.PackageSpec(;
                rev = "makie-v0.25",
                url = "https://github.com/JuliaAstro/SkyCoords.jl/",
            ),
            Pkg.PackageSpec(;
                rev = "fits-wcs",
                url = "https://github.com/JuliaAstro/Reproject.jl",
            ),
        ]
    )

    # Analysis
    using Reproject: reproject
    using SkyCoords: FK5Coords, GalCoords
    using SpectrumBase
    import Contour

    # Data handling and visualization
    using CodecZlib: GzipDecompressor
    using CSV
    using DataFrames: DataFrame
    using Downloads: download
    using FITSFiles
    using URIs: escapeuri
    import Tar
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
    - Compare different data sets (2D and 3D) by overplotting contours
    - Transform coordinate projections and match data resolutions
    - Create intensity moment maps / velocity maps

!!! warning "Companion content"
    - <https://juliaastro.org/AstroImages.jl/previews/PR125/manual/spec/>
"""

# ╔═╡ 91f00e98-e69c-4435-b9d0-10d30006efef
md"""
## Summary

In this tutorial we will visualize 2D and 3D data sets in Galactic and equatorial coordinates.

The tutorial will walk you through a visual analysis of the Small Magellanic Cloud (SMC) using HI 21 cm emission and a Herschel 350 micron map. We will learn how to read in data from a file, query and download matching data from Herschel, and plot the resulting images in a multitude of ways.

The primary packages we will be using are [AstroImages.jl](https://juliaastro.org/AstroImages/) for working with the data and its World Coordinate System (WCS) information, and [Makie.jl](https://docs.makie.org/stable/) for plotting.
"""

# ╔═╡ a6e33cf8-1fe3-4810-a66b-adc07166871e
md"""
### Packages 📦
"""

# ╔═╡ 90111f7f-3a00-4167-8237-4f98036f4a25
md"""
## Download the HI data

We'll be using HI 21 cm emission data from the [HI4PI survey](https://ui.adsabs.harvard.edu/abs/2016A&A...594A.116H/abstract). We want to look at neutral gas emission from the Magellanic Clouds and learn about the kinematics of the system and column densities. Using the VizieR catalog, we've found a relevant data cube to use that covers this region of the sky. You can also download an all-sky data cube, but this is a very large file, so picking out sub-sections can be useful!

For us, the [relevant file is available from CDS Strasbourg](http://cdsarc.u-strasbg.fr/vizier/ftp/cats/J/A+A/594/A116/CUBES/GAL/TAN/TAN_C14.fits). We have a reduced version of it, which is a FITS data cube in Galactic coordinates using the tangential sky projection:
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

# ╔═╡ 5f267771-585c-442e-8884-4d9b3b254904
md"""
## Load the HI data cube

Awesome, so now we have a copy of the data file (a FITS file). So how do we do anything with it?

At the lowest level, [FITSFiles.jl](https://juliaastro.org/FITSFiles/) reads the file into a vector of header-data units (HDUs). This file has a single primary HDU, which holds a 3D array:
"""

# ╔═╡ 93125bb2-efe7-4869-a459-10a7760bd353
hi_data = fits(hi_datafile)

# ╔═╡ 4f95ad6c-7532-45e8-8a20-c4cb20b92e7b
md"""
To work with the data, we use [AstroImages.jl](https://juliaastro.org/AstroImages/) instead. Its `load` function reads the first image HDU into an `AstroImage`: an array that also carries along the FITS header, and the WCS information stored in it.
"""

# ╔═╡ 012eeb56-1442-4f13-876e-5d0b193c88fc
cube = load(hi_datafile)

# ╔═╡ 08b0bac6-3de1-48a6-b1ee-106aa3edcc99
md"""
### Some things to pay attention to here

A data cube has three axes. In this case, there is Galactic longitude (`X`), Galactic latitude (`Y`), and a spectral axis in terms of a local standard of rest (LSR) radial velocity (`Z`).

Julia arrays are column-major, so the data are stored with shape `(n_x, n_y, n_s)`: the first axis of the array is the longitude axis, the second is latitude, and the third is the spectral axis. This is the same order as the `NAXISn` keywords in the FITS header (NumPy reverses it to `(n_s, n_y, n_x)`).

The header also tells us the units of everything. The data are temperatures (`BUNIT = 'K'`), the spatial axes are in degrees, and the spectral axis is in meters per second:
"""

# ╔═╡ d382158a-1e90-4012-b955-99b7c4f76142
header(cube)

# ╔═╡ 7a5f834e-ee79-41f6-92b9-b38f7bf67f58
md"""
AstroImages.jl uses the WCS information in the header to convert between pixel coordinates and world coordinates. For example, here are the longitude, latitude, and velocity of the first and the last voxel of the cube:
"""

# ╔═╡ f3f4a5f8-cb78-4efd-b5b2-fca926a85816
pixel_to_world(cube, [1, 1, 1]), pixel_to_world(cube, [size(cube)...])

# ╔═╡ de560272-e35a-4b56-84cb-2da7d3ab1030
md"""
### Taking a quick look

So let's see what our data actually looks like! Indexing a single channel along the spectral axis gives a 2D image, which is automatically displayed in notebook environments like this one:
"""

# ╔═╡ 8986633d-e8b3-4ebb-a04a-c547267ab994
cube[Z = 300] # Equivalent to cube[:, :, 300]

# ╔═╡ 1c04dd86-56af-4ed1-ab61-530c416dadc3
md"""
To label the axes with world coordinates and add a colorbar, we can plot it with `implotview` instead. The velocity of the channel appears in the title:
"""

# ╔═╡ 016c4e70-121c-45e0-b984-bde2413288ce
cube[Z = 300] |> implotview

# ╔═╡ d0e0220b-248b-48b1-8b5c-d471b2f4d080
md"""
Indexing down to a single spatial position instead gives a 1D spectrum along the spectral axis. To plot it against velocity rather than channel number, we convert each channel to its velocity with `pixel_to_world`, which accepts a matrix with one column per position. We will need this again later, so let's wrap it in a small function:
"""

# ╔═╡ 0614a265-5d93-4016-9676-5c107bf48891
"""
    velocity_axis(img)

Return the radial velocity (in km/s) of each channel along the spectral axis (`Z`) of the data cube `img`.
"""
function velocity_axis(img)
    spec = img[X = 1, Y = 1] # The spectral axis is the same at every spatial position
    return vec(pixel_to_world(spec, axes(spec, 1)')) ./ 1.0e3 # m/s → km/s
end

# ╔═╡ 4ed104bc-56cb-4b8b-8650-b5f0a9ccc330
md"""
We can then bundle these velocities together with the brightness temperatures at a given position into a spectrum object with [SpectrumBase.jl](https://juliaastro.org/SpectrumBase/stable/):
"""

# ╔═╡ 7630e988-72da-4d0b-b48a-de6298781b45
hi_spectrum = spectrum(velocity_axis(cube), collect(cube[X = 75, Y = 75]))

# ╔═╡ 18b5864a-b945-4b8f-8518-36db9dc18f4c
lines(
    spectral_axis(hi_spectrum),
    flux_axis(hi_spectrum);
    axis = (; xlabel = "Radial velocity (km/s)", ylabel = "Brightness temperature (K)"),
)

# ╔═╡ 56640bd0-9302-4bfa-a398-33cceda0fb54
md"""
The narrow line near 0 km/s is HI in the disk of the Milky Way, while the faint bump near +170 km/s is gas associated with the Magellanic system.

!!! todo
    Upstream a SpectrumBase.jl Makie recipe so that this can just be `lines(hi_spectrum)`
"""

# ╔═╡ 2a7650ee-eb8f-43ea-8be5-c549518dee51
md"""
### Exercise

Try messing around with slicing the cube along different axes, or picking out different spectra.
"""

# ╔═╡ 2d942063-8f66-446f-987f-4932497c9c84
md"""
## Make a smaller cube, focusing on the Magellanic Clouds

The HI data cube we downloaded is bigger than we actually need it to be. Let's try zooming in on just the part we need and make a new `sub_cube`. We describe the region that we want in Galactic coordinates, and list the longitude range from left to right as it appears in the images above (Galactic longitude increases to the left):
"""

# ╔═╡ ab31dc8f-d4c2-4053-8f43-f1ae2e4ff990
lon_range = (306, 295) # deg

# ╔═╡ 58fbe8e0-65ee-4e09-9382-c162ce4a7c4d
lat_range = (-46, -40) # deg

# ╔═╡ 2f4ef3b6-e3bf-4196-9de9-d7b2a080d400
md"""
`world_to_pixel` converts the corners of this region to (fractional) pixel coordinates. The cube has three axes, so we also need to pass a velocity, but its value does not matter for the spatial pixel coordinates that we keep here:
"""

# ╔═╡ 14d655d6-b9a9-4402-aad9-f7d5a5c2a322
x_lo, y_lo = world_to_pixel(cube, [lon_range[1], lat_range[1], 0])

# ╔═╡ 94338395-4093-4692-84dc-fceb726afdca
x_hi, y_hi = world_to_pixel(cube, [lon_range[2], lat_range[2], 0])

# ╔═╡ 9723e375-8f24-4a99-835a-f71fb15be1b6
md"""
We then select all of the pixels within these bounds with the `..` interval selector:
"""

# ╔═╡ af5d6c14-e7bc-4bd1-ba1e-235196b7f52a
sub_cube = cube[X = x_lo .. x_hi, Y = y_lo .. y_hi]

# ╔═╡ b9c01fdd-a7e2-4d19-8766-ad8921218c43
md"""
!!! note
    Selectors like `..` act on the lookup values of each dimension, rather than on positions in the array. For an `AstroImage`, these lookup values are the pixel indices of the original cube (e.g., `X` now runs over pixels 66 through 120), so the sub-cube remembers where it came from and its world coordinates stay correct. Its first pixel is the corner of our region:
"""

# ╔═╡ bd8b3061-bf65-4e51-8ee9-b6518063c560
pixel_to_world(sub_cube, [1, 1, 1])

# ╔═╡ 0653fc2f-8e99-452f-aa58-18335f6bcd07
md"""
## Cut along the spectral axis

We don't really need data from such a large velocity range, so let's just extract a little slab between -300 and 300 km/s. The velocity increases with channel number, so this is the range of channels from the first one above our lower limit to the last one below our upper limit:
"""

# ╔═╡ 2507f06c-969e-4b0c-9ba7-7b808cefbeaf
vel_range = (-300, 300) # km/s

# ╔═╡ c6b34b59-a242-4f4c-be4e-912dc8b355ea
sub_cube_slab = let
    vels = velocity_axis(sub_cube)
    sub_cube[Z = findfirst(≥(vel_range[1]), vels):findlast(≤(vel_range[2]), vels)]
end

# ╔═╡ 4189eb5d-7904-4634-ad44-068d00612253
extrema(velocity_axis(sub_cube_slab))

# ╔═╡ 750f7aa9-d2d9-4467-8f94-783959eae06e
md"""
## Moment maps

Moment maps are a useful analysis tool to study data cubes. In short, a moment is a weighted integral along an axis (typically the spectral axis) that can give information about the total intensity (or column density), mean velocity, or velocity dispersion along lines of sight. The zeroth and first moments of the brightness temperature ``T_B`` are:

```math
M_0 = \int T_B \, \mathrm{d}v
\qquad
M_1 = \frac{\int T_B \, v \, \mathrm{d}v}{\int T_B \, \mathrm{d}v}
```

On our grid of channels, these integrals become sums over the spectral axis (`Z`) multiplied by the channel width. The cube stores blank pixels as `NaN`, so we first replace them with zeros to keep them from contributing to the sums:
"""

# ╔═╡ 37aec29e-0b13-4f62-b930-40889a4c3e07
slab = replace(sub_cube_slab, NaN => 0)

# ╔═╡ e4af7841-9e97-4a8b-aca2-190887ee0ca6
dv = slab["CDELT3"] / 1.0e3 # Channel width (km/s)

# ╔═╡ 7aff970f-ee9a-4481-a944-4b8d6fa7f6bb
moment_0 = dropdims(sum(slab; dims = Z); dims = Z) .* dv # K km/s

# ╔═╡ 563f208e-1977-43c7-8dd6-6c2ba771eb2c
md"""
For the first moment, we weight each channel by its velocity. Reshaping the velocities into a ``1 \times 1 \times n_s`` array lines them up with the spectral axis, so that broadcasting multiplies each channel by its own velocity:
"""

# ╔═╡ e4b42310-e09a-4bf7-895f-ba900f868f80
moment_1 = let
    vels = reshape(velocity_axis(slab), 1, 1, :)
    dropdims(sum(slab .* vels; dims = Z); dims = Z) .* dv ./ moment_0 # km/s
end

# ╔═╡ 113a38d1-a114-4b0e-9311-a62e7fbd1068
md"""
Assuming that the HI is optically thin, the zeroth moment converts directly to a column density, ``N_\mathrm{HI} = 1.82 \times 10^{18} \, \mathrm{cm^{-2}} \, (\mathrm{K \, km \, s^{-1}})^{-1} \, M_0``:
"""

# ╔═╡ 94bb0060-6e86-4162-badf-ec6986cd2607
hi_column_density = 1.82e18 .* moment_0 # cm⁻²

# ╔═╡ 74fd7100-3dec-4894-b368-e1cb880233d2
md"""
## Display the moment maps

Each of these moment maps is again an `AstroImage`, which keeps the WCS information of the two spatial axes. Just like for the channel maps above, `implotview` uses it to label the axes with Galactic coordinates.

Collapsing the spectral axis leaves behind a reference to the velocity at its midpoint, which `implotview` would otherwise show in the title, so we set our own title and colorbar label. We then overplot contours of the HI column density on the same axis with Makie's `contour!`:
"""

# ╔═╡ 6dfd219b-0050-457d-92cc-ff05a2c5a6fe
let
    fig, iv = implotview(
        moment_1;
        cmap = Reverse(:RdBu),
        clims = (0, 200),
        colorbar_label = "Velocity (km/s)",
        axis = (; title = "HI 21 cm mean velocity"),
    )

    contour!(
        iv.ax,
        hi_column_density;
        levels = [1.0e20, 5.0e20, 1.0e21, 3.0e21, 5.0e21, 7.0e21, 1.0e22],
        colormap = :Greys,
    )

    fig
end

# ╔═╡ 07125257-e76a-450e-82a2-6e16360dfdde
md"""
The SMC stands out with mean velocities of roughly 100 to 200 km/s, and the column density contours trace its main body.
"""

# ╔═╡ 10ea0fd2-57d7-4216-aae1-674278cc708e
md"""
## Display a longitude-velocity slice

It is often useful to slice along a single latitude and display an image showing longitude and velocity information only (a position-velocity, or longitude-velocity, diagram). Indexing a single position of the slab along `Y` does just that, and `implotview` labels the two remaining axes using the WCS information, noting the latitude of the slice in the title:
"""

# ╔═╡ fc135a57-2bff-499f-bf8e-13a694926cf8
lat_slice = 19 # Position along the latitude axis of the slab

# ╔═╡ 5c212fff-0025-4947-a725-a0ce0963a4d6
implotview(sub_cube_slab[Y = lat_slice])

# ╔═╡ 1fdfc95b-8ae3-411b-bd1a-b2ebdf53bf4c
md"""
As we can see, the SMC seems to be only along positive velocities, well separated from the HI in the disk of the Milky Way near 0 km/s.
"""

# ╔═╡ 19805ca2-9e69-4bfe-af21-8660a98b625c
md"""
### Exercise

Create a new spectral slab isolating just the SMC, and slice along a different dimension to create a latitude-velocity diagram.
"""

# ╔═╡ bfcc7e5c-d68a-454e-855c-d8d81f48efa4
md"""
## Find and download a Herschel image

This is great, but we want to compare the HI emission data with Herschel 350 micron emission to trace some dust. The astropy version of this tutorial uses [astroquery](https://astroquery.readthedocs.io/) to search the [ESASky](https://sky.esa.int/) archive. ESASky also provides a [Table Access Protocol (TAP)](https://www.ivoa.net/documents/TAP/) service, so we can query it directly with a small helper function that submits an [ADQL](https://www.ivoa.net/documents/ADQL/) query and reads the result into a `DataFrame`:
"""

# ╔═╡ 14c01dc5-4aeb-4c4a-9009-7d04efd94c3a
"""
    tap_query(service_url, adql)

Synchronously submit the ADQL query `adql` to the TAP service at `service_url` and return the result as a `DataFrame`.
"""
tap_query(service_url, adql) = CSV.read(
    download("$service_url/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=" * escapeuri(adql)),
    DataFrame,
)

# ╔═╡ ca4e5e99-b5fc-40a6-9963-563e6b96953c
md"""
In the astropy version, ESASky resolves the name "SMC" to a position for us. Here, we look it up in [SIMBAD](https://simbad.u-strasbg.fr/simbad/) with the same helper:
"""

# ╔═╡ f968ed82-fdba-4fcf-8c51-a644f912d3f1
smc = only(
    eachrow(
        tap_query(
            "https://simbad.u-strasbg.fr/simbad/sim-tap",
            "SELECT ra, dec FROM basic JOIN ident ON oidref = oid WHERE id = 'SMC'",
        )
    )
)

# ╔═╡ 5203a827-81b1-4ace-b7e2-dc857a7ec8bd
md"""
ESASky keeps a table of observations for each mission. Herschel's imaging observations live in `observations.mv_v_v_hsa_esasky_photo_fdw_fdw`, where the `fov` column holds the footprint of each observation on the sky. We search for all footprints that intersect a circle with a radius of 1 degree around the SMC:
"""

# ╔═╡ c2416082-4101-4970-89c8-231c23eb6254
herschel_results = tap_query(
    "https://sky.esa.int/esasky-tap/tap",
    """
    SELECT * FROM observations.mv_v_v_hsa_esasky_photo_fdw_fdw
    WHERE 1 = INTERSECTS(CIRCLE('ICRS', $(smc.ra), $(smc.dec), 1), fov)
    """,
)

# ╔═╡ b75307ba-4107-4026-a24b-4d6c797735a7
md"""
Here, the result is a `DataFrame` with 28 Herschel data products that can be downloaded. We can see what information is available by examining its column names:
"""

# ╔═╡ aa689a58-f4fd-4988-af87-88469feec31b
names(herschel_results)

# ╔═╡ 5cdea77f-202a-4668-be2c-3ffad9e6e80a
md"""
We want to find a 350 micron image, so we need to look closer at the filters used for these observations:
"""

# ╔═╡ beec8d4b-cc78-4bdc-ba28-b4478890dcf2
herschel_results.filter

# ╔═╡ eed6e672-757d-4047-986a-8bfd4fab99f3
md"""
Luckily for us, several observations were made with the three filters of Herschel's SPIRE instrument: 250, 350, and 500 microns. These are the ones we want:
"""

# ╔═╡ 5d442b9c-3a66-4738-9272-4c6c7fd0e04b
spire_obs = herschel_results[
    herschel_results.filter .== "250, 350, 500",
    [:observation_id, :instrument, :ra_deg, :dec_deg, :product_url],
]

# ╔═╡ 8ef6816f-9daa-4741-a650-3d463a3b6137
md"""
Each `product_url` points to a tar archive with the processed maps of that observation in all three filters. To match the astropy version of this tutorial, we download the one for observation 1342198565, which covers the SMC bar:

!!! note
    This downloads a ~220 MB file, which can take several minutes. We save it to the `data` folder next to this notebook, so that it only needs to be downloaded once.
"""

# ╔═╡ 6c73a7d4-e223-4490-9e7a-73993ed6e2a8
herschel_tarball = let
    f = joinpath(@__DIR__, "data", "herschel_1342198565.tar")
    if !isfile(f)
        obs = only(filter(row -> row.observation_id == 1342198565, spire_obs))
        download(obs.product_url, f)
    end
    f
end

# ╔═╡ 56105fc3-320c-4c2b-a632-6fabcd4dee9b
md"""
SPIRE map files encode their filter in their names: `psw` (250 μm), `pmw` (350 μm), and `plw` (500 μm). We list the contents of the archive with Julia's [Tar.jl](https://github.com/JuliaIO/Tar.jl) standard library:
"""

# ╔═╡ 429bfeae-fd30-448c-a91e-e425783436ea
Tar.list(herschel_tarball)

# ╔═╡ f9a85248-4828-4413-9921-591b2cfb7c18
md"""
and extract just the 350 μm map:
"""

# ╔═╡ 25812e52-e60f-426b-9d57-b616ae384681
herschel_file = let
    path = only(h.path for h in Tar.list(herschel_tarball) if contains(h.path, "pmw"))
    joinpath(Tar.extract(h -> h.path == path, herschel_tarball), path)
end

# ╔═╡ 161e600b-f470-4e7c-8466-7f1c85aa818a
md"""
This is a gzip-compressed FITS file with several HDUs: the `image` itself, followed by its `error` and `coverage` maps and some processing history. Since we are just doing some qualitative analysis, we only need the image. We read it with the `load_spire_map` helper function defined at the bottom of this notebook:

!!! todo
    Replace `load_spire_map` once FITSFiles.jl can read SPIRE map files directly. It currently errors on their nonstandard `HIERARCH  key.…` header cards, and on their processing history tables.
"""

# ╔═╡ e19ab554-7f48-4b21-996b-846af3eb66df
md"""
Its header holds the WCS information, this time in equatorial coordinates (RA and Dec):
"""

# ╔═╡ 0f8f9c6a-59a6-4470-b523-7232981188e4
md"""
### Display the Herschel image

Just like for our HI maps, `implotview` labels the axes using the WCS information in the header. We use a logarithmic stretch (`logstretch`) to bring out the fainter emission:
"""

# ╔═╡ 26370aef-5c5b-4fd7-98db-89f4365c96ac
md"""
## Overlay HI 21 cm contours on the IR 350 micron image

To visually compare the neutral gas and dust as traced by HI 21 cm emission and IR 350 micron emission, we can draw the HI data on top of the Herschel image. In Makie.jl, a plot can carry a `transformation` that is applied to its data before it is drawn. AstroImages.jl provides `world_transform(iv)`, which maps the world coordinates of an image (here, RA and Dec) to the plotted pixel positions of the image view `iv`. This is the equivalent of `ax.get_transform("world")` in astropy's WCSAxes.

Our HI maps are in Galactic coordinates though, so we first need to convert them to equatorial coordinates.
"""

# ╔═╡ 1356b637-666a-4398-bf64-352be309a5c4
md"""
### From Galactic to equatorial coordinates

[SkyCoords.jl](https://juliaastro.org/SkyCoords/stable/) represents positions on the sky in different coordinate frames, and converts between them with `convert`. The Herschel map uses the FK5 frame at the J2000 equinox, so that is the frame we convert to. SkyCoords.jl stores angles in radians, so we wrap the conversion in a small function that works in degrees, like our WCS information:
"""

# ╔═╡ a1007cb3-13e1-4c9b-88fe-5a5db1dae9d9
"""
    gal_to_eq(l, b)

Convert Galactic longitude `l` and latitude `b` to FK5 (J2000) right ascension and declination, all in degrees.
"""
function gal_to_eq(l, b)
    c = convert(FK5Coords{2000}, GalCoords(deg2rad(l), deg2rad(b)))
    return rad2deg(c.ra), rad2deg(c.dec)
end

# ╔═╡ fafe1a24-7bed-4323-a352-d900633602ef
md"""
As a quick check, the Galactic center should land near Sagittarius A* at an RA and Dec of about (266.4°, -28.9°):
"""

# ╔═╡ a6419723-ea7c-4604-876a-868307f1b134
gal_to_eq(0, 0)

# ╔═╡ b72dad58-268e-427a-bf09-3d54a33d8a82
md"""
### Contours and velocity map

Makie draws contours, but does not expose the computed contour lines that we would need to transform. As in the [contours guide](https://juliaastro.org/AstroImages/stable/guide/contours/) of AstroImages.jl, we compute them with [Contour.jl](https://juliageometry.github.io/Contour.jl/stable/) instead, and convert each of their vertices from HI pixel coordinates to Galactic coordinates (with `pixel_to_world`), and then to equatorial coordinates (with `gal_to_eq`):
"""

# ╔═╡ a057e565-6b43-415e-9584-6dfb95c2a705
"""
    equatorial_contours(img, levels)

Compute the contour lines of the Galactic image `img` at the given `levels`, and return them as `level => vertices` pairs, with the vertices in equatorial (RA, Dec) coordinates.
"""
function equatorial_contours(img, levels)
    result = []
    for cl in Contour.levels(Contour.contours(dims(img)..., img, levels))
        for line in Contour.lines(cl)
            xs, ys = Contour.coordinates(line) # Pixel coordinates
            radec = map(xs, ys) do x, y
                gal_to_eq(pixel_to_world(img, [x, y]; parent = true)...)
            end
            push!(result, Contour.level(cl) => radec)
        end
    end
    return result
end

# ╔═╡ 8b4e79c0-b2c3-4bfd-9e4e-ea5371b393ad
md"""
!!! note
    `dims(img)` holds the pixel indices of the original cube (see above), so Contour.jl returns positions in that frame, and we pass `parent = true` to `pixel_to_world` to interpret them that way.
"""

# ╔═╡ b1ee51cc-153f-4ab5-926a-655bf435c561
hi_levels = [2.0e21, 3.0e21, 5.0e21, 7.0e21, 8.0e21, 1.0e22] # cm⁻²

# ╔═╡ f19cc0d9-1a33-4e51-9b98-b78c7991a7a3
hi_contours = equatorial_contours(hi_column_density, hi_levels);

# ╔═╡ 19410be3-1f24-455f-bd07-7378d5d702ac
md"""
We would also like to see the HI mean velocity map. Makie's `surface` accepts a grid of arbitrary x and y positions, so we can place the center of each HI pixel at its RA and Dec:
"""

# ╔═╡ d073a68e-69ef-44f5-ba76-a6bb8d4ba4fe
hi_radec = [
    gal_to_eq(pixel_to_world(moment_1, [i, j])...)
        for i in axes(moment_1, 1), j in axes(moment_1, 2)
];

# ╔═╡ 8a331faa-8d53-438d-b7b2-4f7f8ab5a3a6
md"""
Passing the same `world_transform` of the Herschel image view to both plots then puts them in the right place:
"""

# ╔═╡ eea69306-6861-49ab-81d2-6a9406200e08
md"""
## Reprojecting to match image resolutions

Reprojection transforms image data from one map projection and resolution to another. A simple use of it is to scale down (or up) the resolution of an image artificially. This could be a useful step if you are trying to get emission line ratios, or directly compare the intensity or flux from a tracer to that of another tracer in the same physical point of the sky.

From our previous images, we can see that the IR Herschel image has a higher spatial resolution than that of the HI data cube. We can look into this more by comparing their pixel sizes (in degrees) from their headers:
"""

# ╔═╡ a0f018b6-ba79-4345-ac21-8a7be9fbb90a
hi_column_density["CDELT1"], hi_column_density["CDELT2"] # HI resolution (dx, dy)

# ╔═╡ 90a34557-df49-4883-87f1-1354851c3f63
md"""
As we can see, the IR pixels are more than 50 times smaller. [Reproject.jl](https://github.com/JuliaAstro/Reproject.jl) can transform an image into a new projection, given the WCS information to transform into. Since we want to match our resolutions, we can just "steal" the WCS information of our HI maps. Like `reproject_interp` in the [reproject](https://reproject.readthedocs.io/) Python package used by the astropy version of this tutorial, Reproject.jl's `reproject` function evaluates the input image at the center of each output pixel by interpolation.

An `AstroImage` keeps track of where a slice or sub-region of a larger image comes from, so we can pass our HI column density map directly, even though it is a 2D slice of a sub-cube of the HI cube. Along with the reprojected image, `reproject` returns a footprint marking the pixels that fall within the Herschel image:
"""

# ╔═╡ b5a1b186-3cb8-48f2-a7a2-449af70c723d
md"""
`rescaled_herschel` shares the header of our HI maps, so it behaves just like them, now with a degraded resolution and in Galactic rather than equatorial coordinates. In particular, we can overlay the HI maps with regular Makie plotting functions, without any transformations:
"""

# ╔═╡ 5aa12cf4-9c6b-4c62-a14b-cf2a734dd40b
md"""
### Challenge

The real power of reprojection is in actually changing the map projection used to display the data. This is done by creating a WCS that contains a different projection type, such as `CTYPE = ["RA---CAR", "DEC--CAR"]`, as opposed to `["RA---TAN", "DEC--TAN"]`.

Use FITSWCS.jl to create a new WCS in a different map projection (see [Working with Celestial Coordinates in WCS 1](/tutorials/wcs-coords/) for how to build one), and reproject the Herschel image onto it with `reproject` to see how distortions in the image change.
"""

# ╔═╡ c65018aa-e30a-4727-ad4e-b853a1479a40
md"""
# Notebook setup 🔧
"""

# ╔═╡ 89e8f2a6-9d3b-44b8-8805-91daa24124c3
TableOfContents(; depth = 4)

# ╔═╡ cc4ba913-ded0-4f71-9171-19bdb9ecf256
"""
    load_spire_map(path)

Load the `image` HDU of the gzip-compressed Herschel SPIRE map at `path` as an `AstroImage`.

FITSFiles.jl cannot parse the nonstandard `HIERARCH  key.…` header cards in these files yet, so we first turn them into `COMMENT` cards. We also stop reading after the image HDU, which skips the processing history tables that follow it.
"""
function load_spire_map(path)
    bytes = transcode(GzipDecompressor, read(path))

    # Header cards are 80 bytes long and aligned to the start of the file
    for i in 1:80:(length(bytes) - 79)
        if view(bytes, i:(i + 13)) == codeunits("HIERARCH  key.")
            bytes[i:(i + 7)] .= codeunits("COMMENT ")
        end
    end

    io = IOBuffer(bytes)
    read(io, HDU) # Primary HDU (header only)
    return AstroImage(read(io, HDU)) # Image HDU
end

# ╔═╡ 73859e33-f129-46f0-95c5-51f87e0d1852
herschel = load_spire_map(herschel_file);

# ╔═╡ 8d3300b1-3c33-4bab-aa52-4239c59b911c
header(herschel)

# ╔═╡ 27e787ed-387e-475a-a10e-4fbd98bb2e2a
herschel_label = "Herschel 350 μm ($(strip(herschel["BUNIT"])))"

# ╔═╡ 9b3619fc-9e9a-4934-a87d-1b3ed25af5fd
implotview(
    herschel;
    cmap = :viridis,
    clims = (2, 50),
    stretch = logstretch,
    colorbar_label = herschel_label,
)

# ╔═╡ 5dacd0b9-a6ab-47e6-8864-c7744c0f5915
let
    fig, iv = implotview(
        herschel;
        cmap = :viridis,
        clims = (5, 50),
        stretch = logstretch,
        colorbar_label = herschel_label,
        figure = (; size = (750, 750)),
    )
    transformation = world_transform(iv)

    # HI mean velocity
    hi_vel = surface!(
        iv.ax,
        first.(hi_radec),
        last.(hi_radec),
        zeros(size(hi_radec));
        color = parent(moment_1),
        colormap = Reverse(:RdBu),
        colorrange = (0, 200),
        alpha = 0.4,
        shading = NoShading,
        transformation,
    )
    Colorbar(
        fig[2, 1],
        hi_vel;
        vertical = false,
        flipaxis = false,
        label = "HI 21 cm mean velocity (km/s)",
    )

    # HI column density
    for (level, radec) in hi_contours
        lines!(
            iv.ax,
            first.(radec),
            last.(radec);
            color = level,
            colormap = :Greys,
            colorrange = extrema(hi_levels),
            transformation,
        )
    end

    fig
end

# ╔═╡ edd0a5bc-2b77-45e5-b645-7fe64439b69b
herschel["CDELT1"], herschel["CDELT2"] # IR resolution (dx, dy)

# ╔═╡ 448f8ef2-a435-4a34-9bb6-fdd37d8bdb95
rescaled_herschel, herschel_footprint = reproject(herschel, hi_column_density);

# ╔═╡ 4c79397f-5380-4662-84c9-d0657f0ead97
let
    fig, iv = implotview(
        rescaled_herschel;
        cmap = :viridis,
        clims = (5, 50),
        stretch = logstretch,
        colorbar_label = herschel_label,
        axis = (; title = "Herschel 350 μm at the HI resolution"),
    )

    # HI mean velocity
    hi_vel = image!(
        iv.ax,
        moment_1;
        colormap = Reverse(:RdBu),
        colorrange = (0, 200),
        alpha = 0.4,
    )
    Colorbar(
        fig[2, 1],
        hi_vel;
        vertical = false,
        flipaxis = false,
        label = "HI 21 cm mean velocity (km/s)",
    )

    # HI column density
    contour!(iv.ax, hi_column_density; levels = hi_levels, colormap = :Greys)

    fig
end

# ╔═╡ Cell order:
# ╟─3c48207e-ae5d-4597-8010-587d6ed8736b
# ╟─91f00e98-e69c-4435-b9d0-10d30006efef
# ╟─a6e33cf8-1fe3-4810-a66b-adc07166871e
# ╠═61c0bf34-302b-4732-a44d-4c2da611eb74
# ╟─90111f7f-3a00-4167-8237-4f98036f4a25
# ╠═8cc6c458-5737-495b-bc04-895aa148661d
# ╟─433fbef7-91de-4a59-aa3c-503141f858ed
# ╟─5f267771-585c-442e-8884-4d9b3b254904
# ╠═93125bb2-efe7-4869-a459-10a7760bd353
# ╟─4f95ad6c-7532-45e8-8a20-c4cb20b92e7b
# ╠═012eeb56-1442-4f13-876e-5d0b193c88fc
# ╟─08b0bac6-3de1-48a6-b1ee-106aa3edcc99
# ╠═d382158a-1e90-4012-b955-99b7c4f76142
# ╟─7a5f834e-ee79-41f6-92b9-b38f7bf67f58
# ╠═f3f4a5f8-cb78-4efd-b5b2-fca926a85816
# ╟─de560272-e35a-4b56-84cb-2da7d3ab1030
# ╠═8986633d-e8b3-4ebb-a04a-c547267ab994
# ╟─1c04dd86-56af-4ed1-ab61-530c416dadc3
# ╠═016c4e70-121c-45e0-b984-bde2413288ce
# ╟─d0e0220b-248b-48b1-8b5c-d471b2f4d080
# ╠═0614a265-5d93-4016-9676-5c107bf48891
# ╟─4ed104bc-56cb-4b8b-8650-b5f0a9ccc330
# ╠═7630e988-72da-4d0b-b48a-de6298781b45
# ╠═18b5864a-b945-4b8f-8518-36db9dc18f4c
# ╟─56640bd0-9302-4bfa-a398-33cceda0fb54
# ╟─2a7650ee-eb8f-43ea-8be5-c549518dee51
# ╟─2d942063-8f66-446f-987f-4932497c9c84
# ╠═ab31dc8f-d4c2-4053-8f43-f1ae2e4ff990
# ╠═58fbe8e0-65ee-4e09-9382-c162ce4a7c4d
# ╟─2f4ef3b6-e3bf-4196-9de9-d7b2a080d400
# ╠═14d655d6-b9a9-4402-aad9-f7d5a5c2a322
# ╠═94338395-4093-4692-84dc-fceb726afdca
# ╟─9723e375-8f24-4a99-835a-f71fb15be1b6
# ╠═af5d6c14-e7bc-4bd1-ba1e-235196b7f52a
# ╟─b9c01fdd-a7e2-4d19-8766-ad8921218c43
# ╠═bd8b3061-bf65-4e51-8ee9-b6518063c560
# ╟─0653fc2f-8e99-452f-aa58-18335f6bcd07
# ╠═2507f06c-969e-4b0c-9ba7-7b808cefbeaf
# ╠═c6b34b59-a242-4f4c-be4e-912dc8b355ea
# ╠═4189eb5d-7904-4634-ad44-068d00612253
# ╟─750f7aa9-d2d9-4467-8f94-783959eae06e
# ╠═37aec29e-0b13-4f62-b930-40889a4c3e07
# ╠═e4af7841-9e97-4a8b-aca2-190887ee0ca6
# ╠═7aff970f-ee9a-4481-a944-4b8d6fa7f6bb
# ╟─563f208e-1977-43c7-8dd6-6c2ba771eb2c
# ╠═e4b42310-e09a-4bf7-895f-ba900f868f80
# ╟─113a38d1-a114-4b0e-9311-a62e7fbd1068
# ╠═94bb0060-6e86-4162-badf-ec6986cd2607
# ╟─74fd7100-3dec-4894-b368-e1cb880233d2
# ╠═6dfd219b-0050-457d-92cc-ff05a2c5a6fe
# ╟─07125257-e76a-450e-82a2-6e16360dfdde
# ╟─10ea0fd2-57d7-4216-aae1-674278cc708e
# ╠═fc135a57-2bff-499f-bf8e-13a694926cf8
# ╠═5c212fff-0025-4947-a725-a0ce0963a4d6
# ╟─1fdfc95b-8ae3-411b-bd1a-b2ebdf53bf4c
# ╟─19805ca2-9e69-4bfe-af21-8660a98b625c
# ╟─bfcc7e5c-d68a-454e-855c-d8d81f48efa4
# ╠═14c01dc5-4aeb-4c4a-9009-7d04efd94c3a
# ╟─ca4e5e99-b5fc-40a6-9963-563e6b96953c
# ╠═f968ed82-fdba-4fcf-8c51-a644f912d3f1
# ╟─5203a827-81b1-4ace-b7e2-dc857a7ec8bd
# ╠═c2416082-4101-4970-89c8-231c23eb6254
# ╟─b75307ba-4107-4026-a24b-4d6c797735a7
# ╠═aa689a58-f4fd-4988-af87-88469feec31b
# ╟─5cdea77f-202a-4668-be2c-3ffad9e6e80a
# ╠═beec8d4b-cc78-4bdc-ba28-b4478890dcf2
# ╟─eed6e672-757d-4047-986a-8bfd4fab99f3
# ╠═5d442b9c-3a66-4738-9272-4c6c7fd0e04b
# ╟─8ef6816f-9daa-4741-a650-3d463a3b6137
# ╠═6c73a7d4-e223-4490-9e7a-73993ed6e2a8
# ╟─56105fc3-320c-4c2b-a632-6fabcd4dee9b
# ╠═429bfeae-fd30-448c-a91e-e425783436ea
# ╟─f9a85248-4828-4413-9921-591b2cfb7c18
# ╠═25812e52-e60f-426b-9d57-b616ae384681
# ╟─161e600b-f470-4e7c-8466-7f1c85aa818a
# ╠═73859e33-f129-46f0-95c5-51f87e0d1852
# ╟─e19ab554-7f48-4b21-996b-846af3eb66df
# ╠═8d3300b1-3c33-4bab-aa52-4239c59b911c
# ╟─0f8f9c6a-59a6-4470-b523-7232981188e4
# ╠═27e787ed-387e-475a-a10e-4fbd98bb2e2a
# ╠═9b3619fc-9e9a-4934-a87d-1b3ed25af5fd
# ╟─26370aef-5c5b-4fd7-98db-89f4365c96ac
# ╟─1356b637-666a-4398-bf64-352be309a5c4
# ╠═a1007cb3-13e1-4c9b-88fe-5a5db1dae9d9
# ╟─fafe1a24-7bed-4323-a352-d900633602ef
# ╠═a6419723-ea7c-4604-876a-868307f1b134
# ╟─b72dad58-268e-427a-bf09-3d54a33d8a82
# ╠═a057e565-6b43-415e-9584-6dfb95c2a705
# ╟─8b4e79c0-b2c3-4bfd-9e4e-ea5371b393ad
# ╠═b1ee51cc-153f-4ab5-926a-655bf435c561
# ╠═f19cc0d9-1a33-4e51-9b98-b78c7991a7a3
# ╟─19410be3-1f24-455f-bd07-7378d5d702ac
# ╠═d073a68e-69ef-44f5-ba76-a6bb8d4ba4fe
# ╟─8a331faa-8d53-438d-b7b2-4f7f8ab5a3a6
# ╠═5dacd0b9-a6ab-47e6-8864-c7744c0f5915
# ╟─eea69306-6861-49ab-81d2-6a9406200e08
# ╠═edd0a5bc-2b77-45e5-b645-7fe64439b69b
# ╠═a0f018b6-ba79-4345-ac21-8a7be9fbb90a
# ╟─90a34557-df49-4883-87f1-1354851c3f63
# ╠═448f8ef2-a435-4a34-9bb6-fdd37d8bdb95
# ╟─b5a1b186-3cb8-48f2-a7a2-449af70c723d
# ╠═4c79397f-5380-4662-84c9-d0657f0ead97
# ╟─5aa12cf4-9c6b-4c62-a14b-cf2a734dd40b
# ╟─c65018aa-e30a-4727-ad4e-b853a1479a40
# ╠═89e8f2a6-9d3b-44b8-8805-91daa24124c3
# ╠═cc4ba913-ded0-4f71-9171-19bdb9ecf256
# ╟─7d07caf5-e203-4152-8bb9-c1f396c4f80c
