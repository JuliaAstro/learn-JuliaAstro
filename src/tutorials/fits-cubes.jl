### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> tags = ["file I/O", "FITS", "images", "contours", "line plots", "colorbars", "data cubes", "radio astronomy", "WCS"]
#> title = "FITS cubes"
#> description = "View and manipulate FITS cubes"
#> date = "2026-10-13"
#> layout = "layout.jlhtml"
#> image = "/assets/fits-cubes.png"
#> order = 2

using Markdown
using InteractiveUtils

# ╔═╡ 61c0bf34-302b-4732-a44d-4c2da611eb74
begin
    # Not used directly, but imported so that Pluto keeps them pinned:
    import ComputePipeline, DimensionalData, DustExtinction, FITSWCS, Makie

    # Analysis
    using Reproject: reproject
    using SkyCoords: FK5Coords, GalCoords
    using SpectrumBase
    import Contour

    # Data handling and visualization
    using CodecZlib: GzipDecompressor
    using CSV
    using DataFrames: DataFrame
    using Downloads: download, RequestError
    using FITSFiles
    using URIs: escapeuri
    import Tar
    using CairoMakie
    using AstroImages
end

# ╔═╡ 7d07caf5-e203-4152-8bb9-c1f396c4f80c
using PlutoUI: TableOfContents

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
        # The archive is sometimes briefly unavailable, so try a few times before giving up.
        # Only server errors and dropped connections are retried, not client errors like 404.
        is_transient(_, e) = e isa RequestError && !(400 ≤ e.response.status < 500)
        patient_download = retry(download; delays = ExponentialBackOff(n = 3, first_delay = 30, max_delay = 300), check = is_transient)
        # Download next to `f` first, so an interrupted download never leaves a partial file at `f`
        mv(patient_download(obs.product_url, f * ".part"), f; force = true)
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AstroImages = "fe3fc30c-9b16-11e9-1c73-17dabf39f4ad"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
CodecZlib = "944b1d66-785c-5afd-91f1-9de20f533193"
ComputePipeline = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
Contour = "d38c429a-6771-53c6-b99e-75d170b6e991"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DimensionalData = "0703355e-b756-11e9-17c0-8b28908087d0"
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
DustExtinction = "fb44c06c-c62f-5397-83f5-69249e0a3c8e"
FITSFiles = "358a0a88-3548-4ad6-b652-8bdbf64af8e5"
FITSWCS = "b9a10f0d-c86c-4e01-8e2e-1f1d7c5fb3f2"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Reproject = "d1dcc2e6-806e-11e9-2897-3f99785db2ae"
SkyCoords = "fc659fc5-75a3-5475-a2ea-3da92c065361"
SpectrumBase = "391af1a9-06f1-59d3-8d21-0be089654739"
Tar = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[sources]
AstroImages = {rev = "reproject", url = "https://github.com/JuliaAstro/AstroImages.jl"}
CairoMakie = {rev = "ff/breaking-0.25", subdir = "CairoMakie", url = "https://github.com/MakieOrg/Makie.jl"}
ComputePipeline = {rev = "ff/breaking-0.25", subdir = "ComputePipeline", url = "https://github.com/MakieOrg/Makie.jl"}
DimensionalData = {rev = "compat-makie-0.25", url = "https://github.com/icweaver/DimensionalData.jl"}
DustExtinction = {rev = "makie-v0.25", url = "https://github.com/JuliaAstro/DustExtinction.jl"}
FITSFiles = {rev = "hdu-size", url = "https://github.com/JuliaAstro/FITSFiles.jl"}
FITSWCS = {rev = "main", url = "https://github.com/JuliaAstro/FITSWCS.jl"}
Makie = {rev = "ff/breaking-0.25", subdir = "Makie", url = "https://github.com/MakieOrg/Makie.jl"}
Reproject = {rev = "fits-wcs", url = "https://github.com/JuliaAstro/Reproject.jl"}
SkyCoords = {rev = "makie-v0.25", url = "https://github.com/JuliaAstro/SkyCoords.jl/"}
SpectrumBase = {rev = "compat/makie-v0.25", url = "https://github.com/JuliaAstro/SpectrumBase.jl"}

[compat]
AstroImages = "~0.6.0"
CSV = "~1.1.0"
CairoMakie = "~0.16.0"
CodecZlib = "~0.7.9"
ComputePipeline = "~0.2.0"
Contour = "~0.6.3"
DataFrames = "~1.8.2"
DimensionalData = "~0.30.3"
DustExtinction = "~0.11.7"
FITSFiles = "~0.4.0"
FITSWCS = "~0.1.0"
Makie = "~0.25.0"
PlutoUI = "~0.7.83"
Reproject = "~0.4.0"
SkyCoords = "~2.0.0"
SpectrumBase = "~0.1.0"
URIs = "~1.7.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.13.2"
manifest_format = "2.1"
project_hash = "8458d15226a22e4e08040c3398ff0d566982cce6"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
registries = "General"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
git-tree-sha1 = "e71ee7b4aa06b045259a7d6101e1cb45ad140bce"
registries = "General"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.4.1"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
registries = "General"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "7063ad1083578215c7c4bf410368150abe8d5524"
registries = "General"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.45"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "dee8388a9bec89c8713ad176f96960ae2f98598a"
registries = "General"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.7.4"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdaptiveArrayPools]]
deps = ["Preferences", "Printf", "Random"]
git-tree-sha1 = "5feaa5055572e809b69b62f06826b1c1c8d3ebf3"
registries = "General"
uuid = "4f381ef7-9af0-4cbe-99d4-cf36d7b0f233"
version = "0.3.11"

    [deps.AdaptiveArrayPools.extensions]
    AdaptiveArrayPoolsCUDAExt = "CUDA"
    AdaptiveArrayPoolsMetalExt = "Metal"

    [deps.AdaptiveArrayPools.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"

[[deps.AdaptivePredicates]]
git-tree-sha1 = "7e651ea8d262d2d74ce75fdf47c4d63c07dba7a6"
registries = "General"
uuid = "35492f91-a3bd-45ad-95db-fcad7dcfedb7"
version = "1.2.0"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
registries = "General"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e092fa223bf66a3c41f9c022bd074d916dc303e7"
registries = "General"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "b4a2d2dfcaedf8f0cb93a970c4ccba6f632f14cc"
registries = "General"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.15.1"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AstroAngles]]
git-tree-sha1 = "e8cc2067b12add82352ad1fd557bc72464af2f6d"
registries = "General"
uuid = "5c4adb95-c1fc-4c53-b4ea-2a94080c53d2"
version = "0.3.0"

[[deps.AstroImages]]
deps = ["AbstractFFTs", "AstroAngles", "ColorSchemes", "DimensionalData", "FITSFiles", "FITSWCS", "FileIO", "ImageAxes", "ImageBase", "ImageIO", "ImageShow", "MappedArrays", "PlotUtils", "PrecompileTools", "Printf", "Statistics", "Tables", "UUIDs"]
git-tree-sha1 = "c9b67254dc7f81b448d7eac937901c54d1c08a07"
repo-rev = "reproject"
repo-url = "https://github.com/JuliaAstro/AstroImages.jl"
uuid = "fe3fc30c-9b16-11e9-1c73-17dabf39f4ad"
version = "0.6.0"
weakdeps = ["Makie", "Reproject"]

    [deps.AstroImages.extensions]
    AstroImagesMakieExt = "Makie"
    AstroImagesReprojectExt = "Reproject"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "94eab0b3ccdcac361188cc661daf69d4433c1818"
registries = "General"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.2.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
registries = "General"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
registries = "General"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.BSplineKit]]
deps = ["ArrayLayouts", "BandedMatrices", "FastGaussQuadrature", "ForwardDiff", "LinearAlgebra", "PrecompileTools", "Random", "Reexport", "SparseArrays", "Static", "StaticArrays", "StaticArraysCore", "StatsAPI"]
git-tree-sha1 = "97ff0843861ef3d10a6f07a4d2830aed50f6145e"
registries = "General"
uuid = "093aae92-e908-43d7-9660-e50ee39d5a0a"
version = "0.19.3"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "54e3d0ce61d4c4863911fe2823d18b27d8d32988"
registries = "General"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.14.0"

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"
    CliqueTreesExt = "CliqueTrees"

    [deps.BandedMatrices.weakdeps]
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BaseDirs]]
git-tree-sha1 = "8c290a1b223deaeea9aea44b235d24546da8eb98"
registries = "General"
uuid = "18cc8868-cbac-4acf-b575-c8ff214dc66f"
version = "1.4.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
registries = "General"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
registries = "General"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CFITSIO]]
deps = ["CFITSIO_jll"]
git-tree-sha1 = "8c6b984c3928736d455eb53a6adf881457825269"
registries = "General"
uuid = "3b1b4be9-1499-4b22-8d78-7db3344d1961"
version = "1.7.2"

[[deps.CFITSIO_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "57ae3edf22a0c6549a239ee798f2eca27cae0754"
registries = "General"
uuid = "b3e40c51-02ae-5482-8a39-3ace5868dcf4"
version = "4.7.0+0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"
version = "1.11.0"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "66188d9d103b92b6cd705214242e27f5737a1e5e"
registries = "General"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.2"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
registries = "General"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.CSV]]
deps = ["CodecZlib", "DataStrings", "Dates", "Downloads", "Durations", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "Printf", "Tables", "Unicode"]
git-tree-sha1 = "ed05d41d96b1a815d52e9de1ef472971bc509832"
registries = "General"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "1.1.0"

    [deps.CSV.extensions]
    CSVDataDecimalsExt = "DataDecimals"
    CSVFilePathsBaseExt = "FilePathsBase"
    CSVInlineStringsExt = "InlineStrings"

    [deps.CSV.weakdeps]
    DataDecimals = "3e2245cb-6932-498d-a7cf-39d39de8bde3"
    FilePathsBase = "48062228-2e41-5def-b9a4-89aafe57970f"
    InlineStrings = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "71aa551c5c33f1a4415867fe06b7844faadb0ae9"
registries = "General"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.1.1"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Cairo_jll", "Colors", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "3ba7cd968eff71e88dcfe1831bd2015bb6e8c534"
repo-rev = "ff/breaking-0.25"
repo-subdir = "CairoMakie"
repo-url = "https://github.com/MakieOrg/Makie.jl"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.16.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "0cfb779e24fdced77a05df1dbe99200821a9f563"
registries = "General"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.8+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9cb23bbb1127eefb022b022481466c0f1127d430"
registries = "General"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.2"

[[deps.ChainRulesCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "4378c3544d8fcc3c8a521782db0e93e2c32d1838"
registries = "General"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.4"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "970758a3d591a2a5c2a907c53f2e2f8c1b1d3537"
registries = "General"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.9"

[[deps.CodecZstd]]
deps = ["TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "da54a6cd93c54950c15adf1d336cfd7d71f51a56"
registries = "General"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.8.7"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON"]
git-tree-sha1 = "07da79661b919001e6863b81fc572497daa58349"
registries = "General"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
registries = "General"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "a7cf827a79b378a381e31e5c7feca853ec19cb94"
registries = "General"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.4"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
registries = "General"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Reexport"]
git-tree-sha1 = "a633156e74096f90b3b5e833529292b651b5671a"
registries = "General"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.3"

[[deps.CommonSolve]]
deps = ["PrecompileTools"]
git-tree-sha1 = "6c389fa857f6ca5a95474b52a52023fd77f24cb7"
registries = "General"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.14"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
registries = "General"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
deps = ["PrecompileTools"]
git-tree-sha1 = "7c4ea0adfb7238bf4678aa36325863fab6163f6f"
registries = "General"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.2.2"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
registries = "General"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.5.5+2"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
registries = "General"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ComputePipeline]]
deps = ["Observables", "Preferences"]
git-tree-sha1 = "e4e5c117d02ec2b62af4852d15283b2abaf3b55a"
repo-rev = "ff/breaking-0.25"
repo-subdir = "ComputePipeline"
repo-url = "https://github.com/MakieOrg/Makie.jl"
uuid = "95dc2771-c249-4cd0-9c9f-1f3b4330693c"
version = "0.2.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
registries = "General"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
registries = "General"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.CoreMath]]
deps = ["CoreMath_jll"]
git-tree-sha1 = "8c0480f92b1b1796239156a1b9b1bfb1b39499b4"
registries = "General"
uuid = "b7a15901-be09-4a0e-87d2-2e66b0e09b5a"
version = "0.1.0"

[[deps.CoreMath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a692a4c1dc59a4b8bc0b6403876eb3250fde2bc3"
registries = "General"
uuid = "a38c48d9-6df1-5ac9-9223-b6ada3b5572b"
version = "0.1.0+0"

[[deps.Crayons]]
git-tree-sha1 = "54b76cbb40d9a0f5368c880725b2f141da77c94f"
registries = "General"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.2.0"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
registries = "General"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataDeps]]
deps = ["Downloads", "Libdl", "Reexport", "SHA", "Scratch", "p7zip_jll"]
git-tree-sha1 = "bc4c8efd0ec8c627028d880b7e584090ff1830ab"
registries = "General"
uuid = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
version = "1.0.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "5fab31e2e01e70ad66e3e24c968c264d1cf166d6"
registries = "General"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.8.2"

[[deps.DataStrings]]
git-tree-sha1 = "b8032263d364177a0b192ea15f277df4e6ab07e2"
registries = "General"
uuid = "48204bd6-5611-42ea-b167-fc1d713f9be2"
version = "1.1.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "b0bc6d2cad1fed8b7fd59a1551a991cb3d2809e6"
registries = "General"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.6"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
registries = "General"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.DelaunayTriangulation]]
deps = ["AdaptivePredicates", "EnumX", "ExactPredicates", "Random"]
git-tree-sha1 = "e1e28f37cbdca4b594ad86ecfe2ee8f25796ebff"
registries = "General"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "1.7.0"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
registries = "General"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "79a2aca180a85c690c58a020d47b426954b590f8"
registries = "General"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.16.0"

[[deps.DimensionalData]]
deps = ["ConstructionBase", "DataAPI", "Dates", "Extents", "Interfaces", "IntervalSets", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "PrecompileTools", "Random", "Statistics", "TableTraits", "Tables"]
git-tree-sha1 = "dfedabc33168c112c713c032e0e1fabe518e4650"
repo-rev = "compat-makie-0.25"
repo-url = "https://github.com/icweaver/DimensionalData.jl"
uuid = "0703355e-b756-11e9-17c0-8b28908087d0"
version = "0.30.3"

    [deps.DimensionalData.extensions]
    DimensionalDataAbstractFFTsExt = "AbstractFFTs"
    DimensionalDataAdaptExt = "Adapt"
    DimensionalDataAlgebraOfGraphicsExt = "AlgebraOfGraphics"
    DimensionalDataArrayInterfaceExt = "ArrayInterface"
    DimensionalDataCategoricalArraysExt = "CategoricalArrays"
    DimensionalDataChainRulesCoreExt = "ChainRulesCore"
    DimensionalDataDiskArraysExt = "DiskArrays"
    DimensionalDataMakieExt = "Makie"
    DimensionalDataNearestNeighborsExt = "NearestNeighbors"
    DimensionalDataPythonCallExt = "PythonCall"
    DimensionalDataRecipesBaseExt = "RecipesBase"
    DimensionalDataSparseArraysExt = "SparseArrays"
    DimensionalDataStaticArraysCoreExt = "StaticArraysCore"
    DimensionalDataStatsBaseExt = "StatsBase"

    [deps.DimensionalData.weakdeps]
    AbstractFFTs = "621f4979-c628-5d54-868e-fcf4e3e8185c"
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
    ArrayInterface = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
    CategoricalArrays = "324d7699-5711-5eae-9e2f-1d82baa6b597"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiskArrays = "3c3547ce-8d99-4f5e-a174-61eb10b00ae3"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    NearestNeighbors = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "Roots", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "2ddbc341251ed713d39e503b504593a5d6789551"
registries = "General"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.132"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
registries = "General"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.7.0"

[[deps.Durations]]
deps = ["Dates"]
git-tree-sha1 = "407b209d76168dcc4355f9dcc4d3b3aaa60dfb3e"
registries = "General"
uuid = "41ac09f5-0a95-4810-aa18-c88bf48ee130"
version = "1.4.1"

    [deps.Durations.extensions]
    DurationsArrowExt = "Arrow"
    DurationsTimeZonesExt = "TimeZones"

    [deps.Durations.weakdeps]
    Arrow = "69666777-d1a9-59fb-9406-91d4454c9d45"
    TimeZones = "f269a46b-ccf7-5d73-abea-4c690281aa53"

[[deps.DustExtinction]]
deps = ["BSplineKit", "DataDeps", "FITSIO", "Healpix", "Unitful", "UnitfulAstro"]
git-tree-sha1 = "8f1df635c4577dccdc7c4204790551ebfc07e94c"
repo-rev = "makie-v0.25"
repo-url = "https://github.com/JuliaAstro/DustExtinction.jl"
uuid = "fb44c06c-c62f-5397-83f5-69249e0a3c8e"
version = "0.11.7"
weakdeps = ["Makie", "SkyCoords"]

    [deps.DustExtinction.extensions]
    MakieExt = "Makie"
    SkyCoordsExt = "SkyCoords"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
registries = "General"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
registries = "General"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "83231673ea4d3d6008ac74dc5079e77ab2209d8f"
registries = "General"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0f666deaf8099a7d7ac229a74402ada88898be1"
registries = "General"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.9.0+0"

[[deps.Extents]]
git-tree-sha1 = "b309b36a9e02fe7be71270dd8c0fd873625332b4"
registries = "General"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.6"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
registries = "General"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FITSFiles]]
deps = ["CodecZlib", "Printf", "Unitful", "UnitfulAngles", "UnitfulAstro", "UnitfulAtomic"]
git-tree-sha1 = "68e710b310233c43c9bd220895e495d00dc8eb57"
repo-rev = "hdu-size"
repo-url = "https://github.com/JuliaAstro/FITSFiles.jl"
uuid = "358a0a88-3548-4ad6-b652-8bdbf64af8e5"
version = "0.4.0"

[[deps.FITSIO]]
deps = ["CFITSIO", "Printf", "Reexport", "Tables"]
git-tree-sha1 = "f57de3f533590c785210893030736dc11c4a4afb"
registries = "General"
uuid = "525bcba6-941b-5504-bd06-fd0dc1a4d2eb"
version = "0.17.5"

[[deps.FITSWCS]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "569f09c3db4eb42a32a97aa10f7400a46349d236"
repo-rev = "main"
repo-url = "https://github.com/JuliaAstro/FITSWCS.jl"
uuid = "b9a10f0d-c86c-4e01-8e2e-1f1d7c5fb3f2"
version = "0.1.0"
weakdeps = ["FITSFiles", "FITSIO"]

    [deps.FITSWCS.extensions]
    FITSWCSFITSFilesExt = "FITSFiles"
    FITSWCSFITSIOExt = "FITSIO"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "4916117dd032ec5959b7633aedbbac408ca5ddeb"
registries = "General"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "1.3.0"

[[deps.FastInterpolations]]
deps = ["AdaptiveArrayPools", "HelpPlots", "LinearAlgebra", "Preferences", "Printf"]
git-tree-sha1 = "f7bcd496f0a2e4e1174e7e1f667f27fb671fd7df"
registries = "General"
uuid = "9ea80cae-fc13-4c00-8066-6eaedb12f34b"
version = "0.4.20"

    [deps.FastInterpolations.extensions]
    FastInterpolationsChainRulesCoreExt = "ChainRulesCore"
    FastInterpolationsColorVectorSpaceExt = ["ColorTypes", "ColorVectorSpace"]
    FastInterpolationsEnzymeExt = "Enzyme"
    FastInterpolationsForwardDiffExt = "ForwardDiff"
    FastInterpolationsRecipesBaseExt = "RecipesBase"
    FastInterpolationsSymbolicsExt = ["Symbolics", "SymbolicUtils"]

    [deps.FastInterpolations.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ColorTypes = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
    ColorVectorSpace = "c3611d14-8923-5661-9e6a-0046d554d3a4"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SymbolicUtils = "d1185830-fcd6-423d-90d6-eec64667417b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "6621fef488e496356c9c9625d0562c12a6070819"
registries = "General"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.20.0"

    [deps.FileIO.extensions]
    HTTPExt = "HTTP"

    [deps.FileIO.weakdeps]
    HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport"]
git-tree-sha1 = "a1b2fbfe98503f15b665ed45b3d149e5d8895e4c"
registries = "General"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.9.0"

    [deps.FilePaths.extensions]
    FilePathsGlobExt = "Glob"
    FilePathsURIParserExt = "URIParser"
    FilePathsURIsExt = "URIs"

    [deps.FilePaths.weakdeps]
    Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
    URIParser = "30578b45-9adc-5946-b283-645ec420af67"
    URIs = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
registries = "General"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "417289ac840d201005d9ed16789a46defccbe36d"
registries = "General"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.18.1"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Random", "Statistics"]
git-tree-sha1 = "59af96b98217c6ef4ae0dfe065ac7c20831d1a84"
registries = "General"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.6"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "e1e9ba3a8f5d71e6ae698f138a1362aeabb42e9f"
registries = "General"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.18.3+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
registries = "General"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "3b0f72e2ffef1a139ac450725c9ecd01c0a3c050"
registries = "General"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.4.6"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
registries = "General"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d95a7a3a7e0f2fffc6291bd33035830aab01310b"
registries = "General"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.15.0+0"

[[deps.FreeTypeAbstraction]]
deps = ["BaseDirs", "ColorVectorSpace", "Colors", "FreeType", "GeometryBasics", "Mmap"]
git-tree-sha1 = "4ebb930ef4a43817991ba35db6317a05e59abd11"
registries = "General"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.8"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f3f50c8001bdfa84f312234fa86715b7c9243771"
registries = "General"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.18+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.Gamma]]
deps = ["LogExpFunctions"]
git-tree-sha1 = "becc397f7cfb06e343496ae6ffb04818a851da51"
registries = "General"
uuid = "a0844989-3bd2-4988-8bea-c9407ab0941b"
version = "1.2.0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "LinearAlgebra", "PrecompileTools", "Random", "StaticArrays"]
git-tree-sha1 = "ec46c5825710fa1a15d468acb2d93cc939a7a5fe"
registries = "General"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.5.13"

    [deps.GeometryBasics.extensions]
    ExtentsExt = "Extents"
    GeometryBasicsGeoInterfaceExt = "GeoInterface"
    IntervalSetsExt = "IntervalSets"

    [deps.GeometryBasics.weakdeps]
    Extents = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
    GeoInterface = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
registries = "General"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
registries = "General"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Giflib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3efbbc027441271444dcd0c0a46f2d119dc4329"
registries = "General"
uuid = "59f7168a-df46-5410-90c8-f2779963d0ec"
version = "6.1.3+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "ea9f1fb37ff10597353b56f3ba0ebd525acef6e0"
registries = "General"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.90.1+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "a641238db938fff9b2f60d08ed9030387daf428c"
registries = "General"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.3"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "69ffb934a5c5b7e086a0b4fee3427db2556fba6e"
registries = "General"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.16+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "ef70da5e123a06a29e2d6ddff0f09985bc226491"
registries = "General"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.11.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "9e544f0758a2341e7e782fde6f3f67610366661c"
registries = "General"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "100.14006.0+0"

[[deps.Healpix]]
deps = ["CFITSIO", "LazyArtifacts", "Libsharp", "LinearAlgebra", "Pkg", "Printf", "Random", "RecipesBase", "StaticArrays", "Test"]
git-tree-sha1 = "611bfa6aa0a90e6f8188ddbb3de1a4b9cf18ca6e"
registries = "General"
uuid = "9f4e344d-96bc-545a-84a3-ae6b9e1b672b"
version = "4.2.4"

[[deps.HelpPlots]]
git-tree-sha1 = "07d389d81ca31e42d68bd88cfa1e900a63b7d378"
registries = "General"
uuid = "5880b2e8-ad0a-48c9-bacd-30dc9c1357e0"
version = "1.1.0"

    [deps.HelpPlots.extensions]
    HelpPlotsExt = "Plots"

    [deps.HelpPlots.weakdeps]
    Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[[deps.HypergeometricFunctions]]
deps = ["Gamma", "LinearAlgebra"]
git-tree-sha1 = "31bb6c92405c084617facc1d7ed9eb6c402d061e"
registries = "General"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.30"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
registries = "General"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
registries = "General"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
registries = "General"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
registries = "General"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "e12629406c6c4442539436581041d372d69c55ba"
registries = "General"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.12"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
registries = "General"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "8c193230235bbcee22c8066b0374f63b5683c2d3"
registries = "General"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.5"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs", "WebP"]
git-tree-sha1 = "f0f005f997dfb8c5fe23920d99458a9619873893"
registries = "General"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.10"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "2a81c3897be6fbcde0802a0ebe6796d0562f63ec"
registries = "General"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.10"

[[deps.ImageShow]]
deps = ["Base64", "ColorSchemes", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "3b5344bcdbdc11ad58f3b1956709b5b9345355de"
registries = "General"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.8"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dcc8d0cd653e55213df9b75ebc6fe4a8d3254c65"
registries = "General"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.2.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
registries = "General"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
registries = "General"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InlineStrings]]
git-tree-sha1 = "06b65886c7577a3784d616e29f1302c2e36e389d"
registries = "General"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.6"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "c72458f1962faeb003bf23cbdb75164fe6280906"
registries = "General"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interfaces]]
git-tree-sha1 = "331ff37738aea1a3cf841ddf085442f31b84324f"
registries = "General"
uuid = "85a1e053-f937-4924-92a5-1367d23b7b87"
version = "0.3.2"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "48922d06068130f87e43edef52382e6a94305ae6"
registries = "General"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.3"
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "CoreMath", "MacroTools", "OpenBLASConsistentFPCSR_jll", "Printf", "Random", "RoundingEmulator"]
git-tree-sha1 = "ff294afb9a15d31d8d7422da138844641a73135f"
registries = "General"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "1.0.11"

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticArblibExt = "Arblib"
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticForwardDiffExt = "ForwardDiff"
    IntervalArithmeticIntervalSetsExt = "IntervalSets"
    IntervalArithmeticIrrationalConstantsExt = "IrrationalConstants"
    IntervalArithmeticLinearAlgebraExt = "LinearAlgebra"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"
    IntervalArithmeticSparseArraysExt = "SparseArrays"

    [deps.IntervalArithmetic.weakdeps]
    Arblib = "fb37089c-8514-4489-9461-98f9c8763369"
    DiffRules = "b552c78f-8df3-52c6-915a-8e097449b14b"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    IrrationalConstants = "92d709cd-6900-40b7-9082-c6be49f344b6"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.IntervalSets]]
git-tree-sha1 = "79d6bd28c8d9bccc2229784f1bd637689b256377"
registries = "General"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.14"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "5cc6ad3e2bfc75a11eed92632ee341d37f4a4536"
registries = "General"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.18"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
registries = "General"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "d22e8bdb3f03c300a103ff23c0dd5a654c6701de"
registries = "General"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.7"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
registries = "General"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
registries = "General"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
registries = "General"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7204148362dafe5fe6a273f855b8ccbe4df8173e"
registries = "General"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.8.0"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "9d8bead1c9ed9f2dc4eb5b51947262163b38ffce"
registries = "General"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.12.1"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "9496de8fb52c224a2e3f9ff403947674517317d9"
registries = "General"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.6"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "037babc10853eeb8e585418922246cb97b8e5b74"
registries = "General"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.2.0+1"

[[deps.JuliaSyntaxHighlighting]]
deps = ["StyledStrings"]
uuid = "ac6e5ff7-fb65-4e79-a425-ec3bc9c03011"
version = "1.12.0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "9eda8292dd3268b3b7ec9df21bbfac24e177ec52"
registries = "General"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.12"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "39bca05343661c347aae0bca57a5994a0bf4f08d"
registries = "General"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.2.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e5b100780d4d30d63b4618d7930d48af409c1772"
registries = "General"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "23.1.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f88f3ccef05a6a72a0cf0ed417c8fd68530f4ab2"
registries = "General"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.1"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "df7566479bd64f20bd16b09960145e70160ffb3b"
registries = "General"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.12"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
registries = "General"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "1.0.0"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "Zlib_jll", "Zstd_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.22.0+0"

[[deps.LibGit2]]
deps = ["LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "LibSSH2_jll", "Libdl", "OpenSSL_jll", "PCRE2_jll", "Zlib_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.9.7+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "OpenSSL_jll", "Zlib_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.104+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
registries = "General"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
registries = "General"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
registries = "General"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b8a91037379e1fd661e9b5ff3d7a770424b1f9f2"
registries = "General"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.4+0"

[[deps.Libsharp]]
deps = ["Libdl", "libsharp2_jll"]
git-tree-sha1 = "e09051a3f95b83091fc9b7a26e6c585c6691b5bc"
registries = "General"
uuid = "ac8d63fe-4615-43ae-9860-9cd4a3820532"
version = "0.2.0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "aebd334d06cee9f24cea70bd19a39749daf73881"
registries = "General"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.3+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9e887184f306b28dbfcff88cefbb0ec80f1aa9e0"
registries = "General"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.4+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.13.0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "28d620bcdb37a3fb02de9f521928ee453690f414"
registries = "General"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "1.0.3"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
registries = "General"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
registries = "General"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "ComputePipeline", "Contour", "Dates", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FastInterpolations", "FileIO", "FilePaths", "FixedPointNumbers", "Format", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageBase", "ImageIO", "InteractiveUtils", "Interpolations", "IntervalSets", "InverseFunctions", "Isoband", "KernelDensity", "LaTeXStrings", "Latexify", "LinearAlgebra", "MacroTools", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "PNGFiles", "Packing", "Pkg", "PlotUtils", "PolygonOps", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun", "Unitful"]
git-tree-sha1 = "4130711b84d0b93bc2b70f744c3e369952eac891"
repo-rev = "ff/breaking-0.25"
repo-subdir = "Makie"
repo-url = "https://github.com/MakieOrg/Makie.jl"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.25.0"

    [deps.Makie.extensions]
    MakieDynamicQuantitiesExt = "DynamicQuantities"
    MakieFFMPEGExt = "FFMPEG_jll"

    [deps.Makie.weakdeps]
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
    FFMPEG_jll = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
registries = "General"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64", "JuliaSyntaxHighlighting", "StyledStrings"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "aa1078778be5a8e5259ff04fbc3d258b3e78d464"
registries = "General"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.6.9"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "Requires"]
git-tree-sha1 = "bdcde8ec04ca84aef5b124a17684bf3b302de00e"
registries = "General"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.11.0"

    [deps.Measurements.extensions]
    MeasurementsBaseTypeExt = "BaseType"
    MeasurementsJunoExt = "Juno"
    MeasurementsRecipesBaseExt = "RecipesBase"
    MeasurementsSpecialFunctionsExt = "SpecialFunctions"
    MeasurementsUnitfulExt = "Unitful"

    [deps.Measurements.weakdeps]
    BaseType = "7fbed51b-1ef5-4d67-9085-a4a9b26f478c"
    Juno = "e5e0dc1b-0480-54bc-9374-aad01c23163d"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
registries = "General"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
registries = "General"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2026.8.13"

[[deps.MuladdMacro]]
deps = ["PrecompileTools"]
git-tree-sha1 = "283bf85d4a767481dd924dff0eee1735e95f449e"
registries = "General"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.7"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "676c0f1c1c34224d6869675047a54b42fdc7b8c5"
registries = "General"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.5"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
registries = "General"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.3.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
registries = "General"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
registries = "General"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.OpenBLASConsistentFPCSR_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "38a93f17e431141c6470bb67a88952a7c4f0e928"
registries = "General"
uuid = "6cdc7f73-28fd-5e50-80fb-958a8875b1af"
version = "0.3.34+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.30+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "97db9e07fe2091882c765380ef58ec553074e9c7"
registries = "General"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.3"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "OpenJPH_jll", "Zlib_jll", "Zstd_jll", "libdeflate_jll"]
git-tree-sha1 = "eb47260eaf61dca5c0e17bd126f4426edac35975"
registries = "General"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.5.2+0"

[[deps.OpenJPH_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "1b1f3cabb1130670812b848f4623d69d578ea805"
registries = "General"
uuid = "95196869-6d61-5868-8eee-c38c33804e0c"
version = "0.32.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.7+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.8+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
registries = "General"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "b079a6affebde167dad8ffe0d670b53917fbcb68"
registries = "General"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "2.0.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.46.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "123266c25174ef6c8d4718920abc206452cf8de6"
registries = "General"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.41"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "32b657a0d57c310a1a172bfc8c8cf68c5e674323"
registries = "General"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.5"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "bc5bf2ea3d5351edf285a06b0016788a121ce92c"
registries = "General"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.1"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
registries = "General"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1912a9f1b9ca55005b03ba075f8e19993583e237"
registries = "General"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.58.2+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools"]
git-tree-sha1 = "663e8b48b789916221e0765393b289ca6c88f24e"
registries = "General"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "3.0.0"

[[deps.PhysicalConstants]]
deps = ["Measurements", "Roots", "Unitful"]
git-tree-sha1 = "3260db49a2599ed83d649b2c3043f1733c3d2895"
registries = "General"
uuid = "5ad8b20f-a522-5ce9-bfc9-ddf1d5bda6ab"
version = "0.2.5"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "e4a6721aa89e62e5d4217c0b21bd714263779dda"
registries = "General"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.46.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "Zstd_jll", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.13.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
registries = "General"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Reexport", "Statistics"]
git-tree-sha1 = "ddcc36fd83fb4cd459a53c891981697bde0d3de6"
registries = "General"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.5.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "e189d0623e7ce9c37389bac17e80aac3b0302e75"
registries = "General"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.83"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
registries = "General"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
registries = "General"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "edbeefc7a4889f528644251bdb5fc9ab5348bc2c"
registries = "General"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.3.4"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "aef1987d90dcc2e2cbd919ad10b2d25bcce4f418"
registries = "General"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.6.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "StyledStrings", "Tables"]
git-tree-sha1 = "99fd8f32ab6728ff205cddd24e2d9ea6cbb498d3"
registries = "General"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.5.0"

    [deps.PrettyTables.extensions]
    PrettyTablesExcelExt = "XLSX"
    PrettyTablesTypstryExt = "Typstry"
    PrettyTablesWriteDocxExt = "WriteDocx"

    [deps.PrettyTables.weakdeps]
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"
    WriteDocx = "d049ceea-54ee-41d7-a26f-ba29db3b6599"
    XLSX = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
registries = "General"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Printf"]
git-tree-sha1 = "ce5ba3edfdeba107a721748d53b5f5887f415853"
registries = "General"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.13.0"

    [deps.ProgressMeter.extensions]
    ProgressMeterDistributedExt = "Distributed"
    ProgressMeterIJuliaExt = "IJulia"

    [deps.ProgressMeter.weakdeps]
    Distributed = "8ba89e20-285c-5b6f-9357-94700520ee1b"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
registries = "General"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "472daaa816895cb7aee81658d4e7aec901fa1106"
registries = "General"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
registries = "General"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "4d8c1b7c3329c1885b857abb50d08fa3f4d9e3c8"
registries = "General"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.7"

[[deps.REPL]]
deps = ["Base64", "Dates", "FileWatching", "InteractiveUtils", "JuliaSyntaxHighlighting", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
registries = "General"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
registries = "General"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
registries = "General"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "0c0dcffc464cda2ff33fc8ffc1dd1bdd53cd98be"
registries = "General"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.4.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
registries = "General"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
registries = "General"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Reproject]]
deps = ["FITSFiles", "FITSWCS", "Interpolations", "SkyCoords"]
git-tree-sha1 = "67bc0be2c91eab3fd719fe76c89e9ac85b093ec2"
repo-rev = "fits-wcs"
repo-url = "https://github.com/JuliaAstro/Reproject.jl"
uuid = "d1dcc2e6-806e-11e9-2897-3f99785db2ae"
version = "0.4.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
registries = "General"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
registries = "General"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d40b2fe70437b01397d2a4d5b020008da4e7019"
registries = "General"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.2+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "e6ffa7bd4f8d0b8744d18f219304738d952f9a37"
registries = "General"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "3.0.11"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "5680a9276685d392c87407df00d57c9924d9f11e"
registries = "General"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.7.1"
weakdeps = ["RecipesBase"]

    [deps.Rotations.extensions]
    RotationsRecipesBaseExt = "RecipesBase"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
registries = "General"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "1.0.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "37ac518b3c924f110d82d50173eaf4ebeaf27ad7"
registries = "General"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.8.1"

[[deps.SciMLPublic]]
deps = ["PrecompileTools"]
git-tree-sha1 = "74685afb51732a464fbce79a72708f7c4203ceb7"
registries = "General"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.3.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
registries = "General"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "084c47c7c5ce5cfecefa0a98dff69eb3646b5a80"
registries = "General"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.10"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays"]
git-tree-sha1 = "57aa595158717ef165e6f5ab639fe2e3178c0a2b"
registries = "General"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.5.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.SignedDistanceFields]]
deps = ["Statistics"]
git-tree-sha1 = "3949ad92e1c9d2ff0cd4a1317d5ecbba682f4b92"
registries = "General"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "7ddb0b49c109481b046972c0e4ab02b2127d6a75"
registries = "General"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.6"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2c79185e5261b159474903598c571ddd9f12ff7b"
registries = "General"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.6"

[[deps.SkyCoords]]
deps = ["ConstructionBase", "LinearAlgebra", "Rotations", "StaticArrays"]
git-tree-sha1 = "47b9e82c678edc5f2267be4fd85856f06a29314c"
repo-rev = "makie-v0.25"
repo-url = "https://github.com/JuliaAstro/SkyCoords.jl/"
uuid = "fc659fc5-75a3-5475-a2ea-3da92c065361"
version = "2.0.0"

    [deps.SkyCoords.extensions]
    AccessorsExt = "Accessors"
    DynamicQuantitiesExt = "DynamicQuantities"
    MakieExt = "Makie"
    NearestNeighborsExt = "NearestNeighbors"
    UnitfulExt = "Unitful"

    [deps.SkyCoords.weakdeps]
    Accessors = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    NearestNeighbors = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "13cd91cc9be159e3f4d95b857fa2aa383b53772a"
registries = "General"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.3"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.13.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "88ef0d98204515bc53dd4aaea5f7788882d1b66c"
registries = "General"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.9.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SpectrumBase]]
deps = ["DustExtinction", "LinearAlgebra", "Measurements", "PhysicalConstants", "RecipesBase", "Unitful", "UnitfulAstro", "UnitfulEquivalences"]
git-tree-sha1 = "5b65a9369188cdd1ccfb04b2d657a4e5801e36e7"
repo-rev = "compat/makie-v0.25"
repo-url = "https://github.com/JuliaAstro/SpectrumBase.jl"
uuid = "391af1a9-06f1-59d3-8d21-0be089654739"
version = "0.1.0"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be1cf4eb0ac528d96f5115b4ed80c26a8d8ae621"
registries = "General"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.2"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "474a5283ad435618090122872eea6a8165ea6bcf"
registries = "General"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.4.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "8a8a52f887f38d44014e24ee46b8a7311c625f6d"
registries = "General"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.24"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
registries = "General"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "e2b53ce13a53367e96601081e33d34746b571bad"
registries = "General"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.5"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
registries = "General"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "adb9da019510162e67a4493fc235c23203d8b09e"
registries = "General"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.13"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91a5737baed20ee31f3faea0e51f57461f6a689e"
registries = "General"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "2.2.1"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StringManipulation]]
deps = ["PrecompileTools", "StyledStrings"]
git-tree-sha1 = "9e3a02d73a8f0a9be04e677dbef758aab40f8bea"
registries = "General"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.6.1"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "4daaddbe37aa80b4551bf7c56adf53e7cfcdf35a"
registries = "General"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.4"

    [deps.StructArrays.extensions]
    AdaptExt = "Adapt"
    GPUArraysCoreExt = "GPUArraysCore"
    KernelAbstractionsExt = "KernelAbstractions"
    LinearAlgebraExt = "LinearAlgebra"
    SparseArraysExt = "SparseArrays"
    StaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "044ed182312a8ad0f28437f048f3e1eae9093a80"
registries = "General"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.11.1"

    [deps.StructUtils.extensions]
    StructUtilsDataFramesExt = ["DataFrames"]
    StructUtilsLazilyInitializedFieldsExt = ["LazilyInitializedFields"]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazilyInitializedFields = "0e77f7df-68c5-4e49-93ce-4cd80f5598bf"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.10.1+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
registries = "General"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "a94d9bdda1b7bed0046cea645639ab3f62196fac"
registries = "General"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.14.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
registries = "General"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["CodecZstd", "ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "PrecompileTools", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "9ca5f1f2d42f80df4b8c9f6ab5a64f438bbd9976"
registries = "General"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.11.9"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
registries = "General"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
registries = "General"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
registries = "General"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.URIs]]
deps = ["Serialization"]
git-tree-sha1 = "50eda8011f098869a76da69a524ff59e0298711a"
registries = "General"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.7.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
registries = "General"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "1f0f9f401753701a7e4113b5056ca38d33875b55"
registries = "General"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.29.0"
weakdeps = ["ConstructionBase", "ForwardDiff", "InverseFunctions", "LaTeXStrings", "Latexify", "NaNMath", "Printf"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "79875b1f2e4bf918f0702a5980816955066d9ae2"
registries = "General"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.7.2"

[[deps.UnitfulAstro]]
deps = ["Unitful", "UnitfulAngles"]
git-tree-sha1 = "fbe44a0ade62ae5ed0240ad314dfdd5482b90b40"
registries = "General"
uuid = "6112ee07-acf9-5e0f-b108-d242c714bf9f"
version = "1.2.2"

[[deps.UnitfulAtomic]]
deps = ["Unitful"]
git-tree-sha1 = "903be579194534af1c4b4778d1ace676ca042238"
registries = "General"
uuid = "a7773ee8-282e-5fa2-be4e-bd808c38a91a"
version = "1.0.0"

[[deps.UnitfulEquivalences]]
deps = ["Unitful"]
git-tree-sha1 = "249cafa9e88c0c28dbd7b55e258a53a01761f974"
registries = "General"
uuid = "da9c4bc3-91c8-4f02-8a40-6b990d2a7e0c"
version = "0.2.1"

[[deps.WebP]]
deps = ["CEnum", "ColorTypes", "FileIO", "FixedPointNumbers", "ImageCore", "libwebp_jll"]
git-tree-sha1 = "aa1ca3c47f119fbdae8770c29820e5e6119b83f2"
registries = "General"
uuid = "e3aaa7dc-3e4b-44e0-be63-ffb868ccd7c1"
version = "0.1.3"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
registries = "General"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e52eca002a11c30a858185efdfb15311e1c7a6bf"
registries = "General"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "ce91b54a5f90c4f27df6950706234b34e49cf79e"
registries = "General"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.14+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
registries = "General"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
registries = "General"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
registries = "General"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a54cf54677c9f6eac0d97e9d757c233b56bc21b9"
registries = "General"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.13+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
registries = "General"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
registries = "General"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.3.2+0"

[[deps.Zstd_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl"]
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
registries = "General"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.15.0+0"

[[deps.libdeflate_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7ee221abc4f125ce18e73d48367cde0e16bc93b9"
registries = "General"
uuid = "46979653-d7f6-5232-b59e-dd310c4598de"
version = "1.26.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "1fcd6cf05d23adab0b5ffc496d020537215751dc"
registries = "General"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.59+1"

[[deps.libsharp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl", "Pkg"]
git-tree-sha1 = "1029b8ace2e7cb7c6bafde9300380acd6faee75a"
registries = "General"
uuid = "180207a7-b08e-5162-af94-7d62a04fe081"
version = "1.0.2+2"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "libpng_jll"]
git-tree-sha1 = "e067c8bae65bb40866552296a2d5b4dc65b8928f"
registries = "General"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.80.702+0"

[[deps.libwebp_jll]]
deps = ["Artifacts", "Giflib_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libglvnd_jll", "Libtiff_jll", "libpng_jll"]
git-tree-sha1 = "52d3b9475133c3bc8c0a7f90f18bfc8cdb5a443a"
registries = "General"
uuid = "c5f90fcd-3b7e-5836-afba-fc50a0988cb2"
version = "1.6.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.70.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.8.3+0"

[registries.General]
url = "https://github.com/JuliaRegistries/General.git"
uuid = "23338594-aafe-5451-b93e-139f81909106"
"""

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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
