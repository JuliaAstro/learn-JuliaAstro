### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Astronomical Coordinates 3: Working with Velocity Data"
#> layout = "layout.jlhtml"
#> date = "2025-12-31"
#> description = "Represent and transform velocity data alongside sky coordinates in Julia, and use proper motions to predict the past position of a fast-moving star on a DSS image."
#> tags = ["coordinates", "gaia"]

using Markdown
using InteractiveUtils

# ╔═╡ 955ad972-038b-4519-9eeb-b5a297b4d185
begin
    # Can remove this block after Makie v0.25 is merged:
    # https://github.com/MakieOrg/Makie.jl/pull/5484
    import Pkg
    Pkg.activate(; temp = true)
    Pkg.add(
        [
            Pkg.PackageSpec(; name = "Downloads"),
            Pkg.PackageSpec(; name = "TOML"),
            Pkg.PackageSpec(; name = "CSV"),
            Pkg.PackageSpec(; name = "DataFramesMeta"),
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(; name = "URIs"),
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
                url = "https://github.com/JuliaAstro/FITSFiles.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaAstro/FITSWCS.jl",
            ),
            Pkg.PackageSpec(;
                rev = "makie",
                url = "https://github.com/JuliaAstro/AstroImages.jl",
            ),
            Pkg.PackageSpec(;
                rev = "makie-0.25",
                url = "https://github.com/icweaver/DimensionalData.jl",
            ),
            Pkg.PackageSpec(;
                rev = "altaz",
                url = "https://github.com/activexray/SkyCoords.jl",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaPhysics/DynamicQuantities.jl",
            ),
        ]
    )

    using Downloads: download
    using CSV
    using DataFramesMeta
    using CairoMakie
    using DynamicQuantities
    using SkyCoords
    using AstroImages
    using FITSWCS
    using URIs: escapeuri
    using LinearAlgebra: ⋅

    using DynamicQuantities.Units: °, yr
    using DynamicQuantities.Constants: pc
end

# ╔═╡ bfd71b9c-ca7e-4480-92d3-c16fc03358d2
begin
    using TOML: TOML
    using PlutoUI: TableOfContents
end

# ╔═╡ 3b782148-7df3-42e9-a30d-add6ecac4e76
md"""
## Summary
In the previous tutorials in this series, we showed how astronomical *positional* coordinates can be represented and transformed in Julia with SkyCoords.jl ([docs](https://juliaastro.org/SkyCoords/stable/)). Many sources, especially stars (thanks to the *Gaia* mission), also have measured velocities or measured components of their velocity (e.g., just proper motion, or just radial velocity).

In this tutorial, we will explore how astronomical coordinates that have associated velocity data can be represented and transformed in Julia. You may find it helpful to keep the documentation linked above open alongside this tutorial for reference or additional reading.

*Note: This is the 3rd tutorial in a series of tutorials about coordinates in Julia. If you are new to this series, you may want to start from the beginning or an earlier tutorial.*
"""

# ╔═╡ d3db4bd5-c2f5-4be9-97ea-a9f4d1fbfd56
md"""
### Packages 📦

We start by importing some general packages that we will need below:
"""

# ╔═╡ cc4fab2a-8a85-4dd8-833f-523d675597e3
md"""
## More Than Just Sky Positions: Including Velocity Data

As we have seen in the previous tutorials, `AbstractSkyCoords` objects can be used to store scalars or vectors of positional coordinate information, and they support transforming between different coordinate frames. Astropy additionally supports representing and transforming *velocity* information along with positions inside its `SkyCoord` objects ([docs](https://docs.astropy.org/en/latest/coordinates/velocities.html)). SkyCoords.jl does not support this (yet!), but as we will see below, composing the same functionality ourselves is only a few lines of Julia:

!!! todo
    Add velocity component support (proper motion and radial velocity) to `AbstractSkyCoords` objects in SkyCoords.jl

## Passing Velocity Data

In Julia, a natural lightweight way to bundle velocity components together with a sky position is a `NamedTuple`. For example, to represent a sky position and a proper motion in the ICRS coordinate frame, in addition to a coordinate object for the position components `ra`, `dec`, we can store values for the proper motion components `pm_ra_cosdec` and `pm_dec` ("pm" for "proper motion") as unitful DynamicQuantities.jl quantities:
"""

# ╔═╡ e5e0904a-eae9-47b0-8758-4b0a7d5b142f
(;
    coord = ICRSCoords(10°, 20°),
    pm_ra_cosdec = 1us"mas/yr",
    pm_dec = 2us"mas/yr",
)

# ╔═╡ e2e15684-5dcf-4889-b966-ef5fc59a982f
md"""
Here, you may notice that the proper motion in right ascension has "cosdec" in the name: This is to explicitly note that the value is expected to be the proper motion scaled by the cosine of the declination, which accounts for the fact that a change in longitude (right ascension) has different physical length at different latitudes (declinations).

Like the examples in previous tutorials demonstrated for positional coordinates, we can also bundle arrays of data for all of the components. In this case, each value in the array represents a quantity for one object in a larger data set. This is beneficial when dealing with a large collection of stars, like a star cluster:
"""

# ╔═╡ f2d636e1-547f-4f96-aa86-ae3c16c11166
(;
    coord = ICRSCoords.(range(0, 10, 5) * °, range(5, 20, 5) * °),
    pm_ra_cosdec = range(-5, 5, 5) * us"mas/yr",
    pm_dec = range(-5, 5, 5) * us"mas/yr",
)

# ╔═╡ cdc27d8a-3150-4269-befc-9c7002a2a39c
md"""
However, for some of the examples below we will continue to use scalar values for brevity.

We can also include radial velocity data with a `radial_velocity` field:
"""

# ╔═╡ ae4aa69a-73b4-4d95-b8bb-ee3765132c80
velocity_coord = (;
    coord = ICRSCoords(10°, 20°),
    pm_ra_cosdec = 1us"mas/yr",
    pm_dec = 2us"mas/yr",
    radial_velocity = 100us"km/s",
)

# ╔═╡ a8f2133d-e970-47d2-b244-5b665f21210d
md"""
The component data can then be accessed using the same names used to store the velocity components:
"""

# ╔═╡ 8291ed3e-4f11-481c-be48-3c50ccaa3839
velocity_coord.pm_ra_cosdec

# ╔═╡ 7b3581dd-01c3-4621-a988-750a8ec4feb8
velocity_coord.radial_velocity

# ╔═╡ 08193019-6c81-4f83-a89f-c2d288192004
md"""
## Transforming Velocity Data

The *position* in our bundle can be transformed to other frames just like the position-only coordinate objects we used in the previous tutorials — for example, with the `GalCoords` constructor. But what about the proper motion components? As discussed in the previous tutorial, the Galactic frame is a fixed 3D *rotation* away from the ICRS frame, with no change of origin. A velocity vector transforms with exactly the same rotation as the position vector, so we can transform the velocity data ourselves using the rotation matrix that SkyCoords.jl already uses internally for its position transformations (a non-exported constant):
"""

# ╔═╡ 911b0b50-ea1f-4d70-bc40-73a21b94bd2c
SkyCoords.ICRS_TO_GAL

# ╔═╡ 84ea56f6-e90a-4744-8136-f0d35cabfa0b
md"""
A proper motion is a rate of change of position *on the sky*, so it lives in the plane tangent to the celestial sphere at the position of the source: the `pm_ra_cosdec` component points along the local "east" direction (increasing longitude), and `pm_dec` along the local "north" direction (increasing latitude). At a spherical position with longitude ``\lambda`` and latitude ``\phi``, these two directions are given by the Cartesian unit vectors:

```math
\hat{e}_\mathrm{east} = (-\sin\lambda, \cos\lambda, 0), \quad
\hat{e}_\mathrm{north} = (-\cos\lambda \sin\phi, -\sin\lambda \sin\phi, \cos\phi)
```
"""

# ╔═╡ f37980ea-708f-42f2-b247-5d31ba033fb1
begin
    ê_east(lon, lat) = [-sin(lon), cos(lon), 0]
    ê_north(lon, lat) = [-cos(lon) * sin(lat), -sin(lon) * sin(lat), cos(lat)]
end

# ╔═╡ 87d2b737-f1c8-4182-8603-b46e97858b43
md"""
To transform a proper motion from ICRS to Galactic, we then: i) assemble the proper motion vector in Cartesian coordinates, ii) rotate it with the frame rotation matrix, and iii) read off its components along the east/north directions at the new (Galactic) sky position:
"""

# ╔═╡ e0331849-286b-4284-9c74-b9d3a3afd122
"""
    galactic_pm(c::ICRSCoords, pm_ra_cosdec, pm_dec)

Transform the given ICRS proper motion components of the source at `c` to the Galactic
frame, returned as a `NamedTuple` with the component names `pm_l_cosb` and `pm_b`.
"""
function galactic_pm(c::ICRSCoords, pm_ra_cosdec, pm_dec)
    g = GalCoords(c)
    v_icrs = ustrip(us"mas/yr", pm_ra_cosdec) * ê_east(c.ra, c.dec) +
        ustrip(us"mas/yr", pm_dec) * ê_north(c.ra, c.dec)
    v_gal = SkyCoords.ICRS_TO_GAL * v_icrs
    return (;
        pm_l_cosb = (v_gal ⋅ ê_east(g.l, g.b))us"mas/yr",
        pm_b = (v_gal ⋅ ê_north(g.l, g.b))us"mas/yr",
    )
end

# ╔═╡ a3cd666d-cf0f-4d88-875e-e66cb1d8a99f
velocity_coord_gal = (;
    coord = GalCoords(velocity_coord.coord),
    galactic_pm(
        velocity_coord.coord,
        velocity_coord.pm_ra_cosdec,
        velocity_coord.pm_dec,
    )...,
    velocity_coord.radial_velocity,
)

# ╔═╡ ecd79657-1c56-42fd-b63b-39d9f248dcc4
md"""
Note that, like the position components, which change from `ra`, `dec` to `l`, `b`, the proper motion component names have changed to reflect naming conventions for the component names in the new frame: `pm_ra_cosdec` and `pm_dec` have become `pm_l_cosb` and `pm_b`:
"""

# ╔═╡ e6dcd0f7-9833-4f98-9197-feea74baad3e
velocity_coord_gal.pm_l_cosb

# ╔═╡ cb37c0d1-7c22-4467-a69e-b984beff4ee6
velocity_coord_gal.pm_b

# ╔═╡ 408ba1ad-8856-4b92-989e-9fc885762944
md"""
The radial velocity, on the other hand, is carried over unchanged: it points along the line of sight, and a pure rotation does not change the line-of-sight component of the velocity.

An important caveat to note when transforming velocity data is that some reference frames require knowing the distance and/or the full velocity vector (i.e., the proper motion components *and* the radial velocity) in order to transform correctly. For example, the transformation to the Galactocentric frame involves both a shift of origin (from the solar system barycenter to the Galactic center) and a velocity offset (the motion of the Sun around the Galactic center), so a source with only sky position and proper motion data **cannot** be transformed to it — attempting this in astropy raises a `ConvertError`. The Galactocentric frame is not yet implemented in SkyCoords.jl:

!!! todo
    Implement Galactocentric coords
"""

# ╔═╡ 3892c3dc-c2eb-42a4-b2b4-46f76b1bc475
md"""
## Evolving Coordinate Positions Between Epochs

For nearby or fast-moving stars, a star's position could change appreciably between two well-spaced observations of the source. For such cases, it might be necessary to compute the position of the star at a given time using its proper motion or velocity. Let's demonstrate this idea by comparing the sky position of a source as measured by [*Gaia* Data Release 2](https://www.cosmos.esa.int/web/gaia/dr2) (given at the epoch J2015.5) to an image near this source from the Digitized Sky Survey (DSS; digital scans of photographic plates observed in the 1950s).

From previous astrometric measurements, we know that the star HD 219829 has a very large proper motion: close to 0.5 arcsec/year! Between the DSS and *Gaia*, we therefore expect that the position of the star has changed by about 0.5 arcmin. Let's see if this is the case!

To start, we will query the *Gaia* catalog to retrieve data for this star, re-using the `from_name` SIMBAD lookup helper that we wrote in the first tutorial of this series to get a sky coordinate to search around:
"""

# ╔═╡ 0d4d2809-9580-4fec-8f3a-766e6e7b9a25
md"""
In the first tutorial, we submitted our archive queries with VirtualObservatory.jl. Unfortunately, that package cannot currently be installed side-by-side with the AstroImages.jl development branch that we will use to display images later in this notebook (they require incompatible versions of AstroAngles.jl), so for this notebook we roll the same one-line TAP query ourselves with a plain HTTP download:

!!! todo
    Switch these queries (back) to VirtualObservatory.jl after the AstroAngles.jl compat bump in its VOTables.jl dependency is merged:
    <https://github.com/JuliaAPlavin/VOTables.jl/pull/18>
"""

# ╔═╡ 6a7ef5e6-18d5-4c9c-9562-bbf32c90e411
"""
    tap_query(service_url, adql)

Synchronously submit the ADQL query `adql` to the TAP service at `service_url` and return the result as a `DataFrame`.
"""
tap_query(service_url, adql) = CSV.read(
    download("$service_url/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=" * escapeuri(adql)),
    DataFrame,
)

# ╔═╡ fc06f86b-2286-4aa1-a11a-2b0238416e8f
"""
    from_name(name)

Look up the ICRS coordinates of `name` in [SIMBAD](https://simbad.u-strasbg.fr/simbad/) and return them as an `ICRSCoords` object.
"""
function from_name(name)
    result = tap_query(
        "https://simbad.u-strasbg.fr/simbad/sim-tap",
        "SELECT ra, dec FROM basic JOIN ident ON oidref = oid WHERE id = '$name'",
    )
    isempty(result) && error("Could not resolve \"$name\" with SIMBAD")
    obj = only(eachrow(result))
    return ICRSCoords(deg2rad(obj.ra), deg2rad(obj.dec))
end

# ╔═╡ acc6f854-eaeb-419b-8bec-51382e06d8b6
hd219829_center = from_name("HD 219829")

# ╔═╡ 61b57f57-3529-4bd3-9ae4-716793858266
md"""
As in the first tutorial, we express the *Gaia* archive cone search in ADQL, and submit it to the archive's TAP service — this time with our own `tap_query` helper. We use a large-ish search radius (1 arcmin) so many sources will be returned:
"""

# ╔═╡ 75a973b2-cb8b-45fb-bdc3-b66c32b1bc50
hd219829_query = """
SELECT
    source_id, ra, dec, parallax, pmra, pmdec,
    radial_velocity, phot_g_mean_mag, ref_epoch,
    DISTANCE(
        POINT('ICRS', ra, dec),
        POINT('ICRS', $(rad2deg(hd219829_center.ra)), $(rad2deg(hd219829_center.dec)))
    ) AS dist
FROM gaiadr2.gaia_source
WHERE 1 = CONTAINS(
    POINT('ICRS', ra, dec),
    CIRCLE('ICRS', $(rad2deg(hd219829_center.ra)), $(rad2deg(hd219829_center.dec)), $(1 / 60))
)
ORDER BY dist ASC
"""

# ╔═╡ 0990f7a6-c401-4de9-841e-d240b460265e
md"""
Running this query requires an internet connection, so we have included the results file next to this notebook. The cell below only submits the query (and saves the results to `data/HD_219829_query_results.csv`) if that file is not already present:
"""

# ╔═╡ d67675df-7291-466e-b40d-7e90ff4d99e0
hd219829_table = let
    f = joinpath(@__DIR__, "data", "HD_219829_query_results.csv")
    if !isfile(f)
        CSV.write(f, tap_query("https://gea.esac.esa.int/tap-server/tap", hd219829_query))
    end
    CSV.read(f, DataFrame)
end

# ╔═╡ d5ff195c-9e28-4e90-8250-4d1d9ea2f136
md"""
We know that HD 219829 will be the brightest source in this small region, so we can extract the row with the smallest G-band magnitude. Let's check the proper motion values for this source to make sure that they are large:
"""

# ╔═╡ a8e57417-9a4e-49e9-9d9d-73d8d8938071
hd219829_row = hd219829_table[argmin(hd219829_table.phot_g_mean_mag), :]

# ╔═╡ 1e3179fe-05b3-4da6-ac4e-1991f4866158
hd219829_row[[:source_id, :pmra, :pmdec]]

# ╔═╡ bb23a77f-2996-49f3-8bca-a7aa25b55e8e
md"""
Indeed, it looks like this is our source! Let's construct a velocity-data bundle for it using the data from the *Gaia* archive. We also convert the measured parallax to a distance like we did in the first tutorial, and record the reference epoch of the measurements as `obstime` (more on this below):

*Note about the Gaia catalog proper motion column names: The names in the Gaia archive and other repositories containing Gaia data give right ascension proper motion values simply as `pmra`. These components implicitly contain the ``\cos(\mathrm{dec})`` term, so we do **not** have to modify these values in order to store them as `pm_ra_cosdec`.*
"""

# ╔═╡ ef024267-f69b-4e12-9aa5-710e12c5d520
parallax_to_distance(ϖ_mas) = 1000 / ϖ_mas * pc

# ╔═╡ cdbff55e-ae03-44d0-8e7e-c60d6f224368
hd219829_coord = (;
    coord = ICRSCoords(hd219829_row.ra * °, hd219829_row.dec * °),
    distance = parallax_to_distance(hd219829_row.parallax) |> us"Constants.pc",
    pm_ra_cosdec = hd219829_row.pmra * us"mas/yr",
    pm_dec = hd219829_row.pmdec * us"mas/yr",
    obstime = hd219829_row.ref_epoch,
)

# ╔═╡ 5b2c80cc-e869-4e95-b301-3e83b5b061db
md"""
We now have a representation of the position and proper motion of the star HD 219829 as measured by *Gaia* and reported at the epoch J2015.5. What does this mean exactly? *Gaia* actually measures the (time-dependent) position of a star every time it scans the part of the sky that contains the source, and this is how *Gaia* is able to measure proper motions of stars. However, if every star is moving and changing its sky position, how do we ever talk about "the sky position" of a star as opposed to "the sky trajectory of a star"?! The key is that catalogs often only report the position of a source at some reference epoch. For a survey that only observes the sky once or a few times (e.g., SDSS or 2MASS), this reference epoch might be "the time that the star was observed." But for a mission like *Gaia*, which scans the sky many times, they perform astrometric fits to the individual position measurements, which allow them to measure the parallax, proper motion, and the reference position at a reference time for each source. For *Gaia* data release 2, the reference time is J2015.5, and the sky positions (and other quantities) reported in the catalog for each source are at this epoch — this is the Julian year that we stored in the `obstime` field above.

Now that we have all of this information for HD 219829, let's compare the position of the star as measured by *Gaia* to its apparent position in an image from the DSS. We can query the [STScI DSS image cutout service](https://archive.stsci.edu/cgi-bin/dss_form) to retrieve a FITS image of the field around this star. As with the catalog query above, the cell below only downloads the image (to `data/dss_hd219829.fits`) if it is not already present:
"""

# ╔═╡ 4f8bd7bc-e292-48c7-a4f2-1c8798f7b4e8
dss_cutout_filename = let
    c = hd219829_coord.coord
    # width/height in arcmin
    url = "http://archive.stsci.edu/cgi-bin/dss_search?f=FITS&ra=$(rad2deg(c.ra))&dec=$(rad2deg(c.dec))&width=4&height=4"
    f = joinpath(@__DIR__, "data", "dss_hd219829.fits")
    isfile(f) ? f : download(url, f)
end;

# ╔═╡ f36ca985-f94b-4a5d-a0c4-5e1c13e09c4f
md"""
We can now load the FITS image with AstroImages.jl ([docs](https://juliaastro.org/AstroImages/)), which bundles the image data together with its header metadata and automatically displays the image in the notebook:
"""

# ╔═╡ 5cfa1e48-91a3-433d-8119-10de94c94e0b
dss_img = load(dss_cutout_filename)

# ╔═╡ 8f1d2372-7ece-4b72-8dc8-29ba3ef5f5d2
md"""
The FITS header contains World Coordinate System (WCS) information — metadata that defines the mapping between pixel coordinates in the image and sky coordinates. (If you are unfamiliar with FITS files or WCS, check out the [FITS Images tutorial](https://learn.juliaastro.org/tutorials/fileio-fits_images/), which explains both in more detail.) We can construct a `WCSTransform` object for this mapping directly from the header with FITSWCS.jl:

!!! note
    The DSS returns its astrometry as a legacy plate solution instead of standard WCS keywords; FITSWCS.jl detects this and rebuilds an equivalent standard (tangent-plane) WCS from it.
"""

# ╔═╡ 4c56f658-6735-416d-892a-8aad7ed9a733
dss_wcs = WCS(header(dss_img))

# ╔═╡ f1f08c80-0603-4c92-b66d-c7e0a426d413
md"""
By converting the *Gaia*-measured sky position of HD 219829 to pixel coordinates with `world_to_pixel`, we can over-plot a marker for it on the DSS image:

!!! todo
    Add a Makie.jl recipe to AstroImages.jl that displays world coordinate axis labels and grid lines (an equivalent of astropy's [WCSAxes](https://docs.astropy.org/en/latest/visualization/wcsaxes/index.html)). For now, we plot in pixel coordinates.
"""

# ╔═╡ 51f6ac42-b6dd-4e50-bd43-f1215ddc1ff0
let
    fig = Figure(size = (600, 600))

    ax = Axis(
        fig[1, 1];
        xlabel = "x [pix]",
        ylabel = "y [pix]",
        title = "DSS image around HD 219829",
        aspect = DataAspect(),
    )

    image!(ax, dss_img)

    c = hd219829_coord.coord
    x, y = world_to_pixel(dss_wcs, rad2deg.([c.ra, c.dec]))
    scatter!(
        ax, x, y;
        markersize = 30,
        color = :transparent,
        strokecolor = :red,
        strokewidth = 2,
    )

    fig
end

# ╔═╡ 41bcfaae-0dd6-4dcf-a8d2-98b8df1b96d0
md"""
The brightest star (as observed by DSS) in this image is our target, and the red circle is where *Gaia* observed this star. As we expected, it has moved quite a bit since the 1950s! We can account for this motion and predict the position of the star at around the time the DSS plate was observed. Let's assume that this plate was observed in 1950 exactly (this is not strictly correct, but should get us close enough).

To account for the proper motion of the source and evolve the position to a new time, astropy provides a `SkyCoord.apply_space_motion()` method. SkyCoords.jl does not have an equivalent (yet!), but for the small sky motions involved here, linearly propagating the sky position along the proper motion direction is an excellent approximation — and only takes a few lines to write ourselves:

!!! todo
    Add an `apply_space_motion`-style epoch propagation function to SkyCoords.jl that properly propagates the full 3D space motion (including the distance and radial velocity).
"""

# ╔═╡ 0a87c508-7149-48db-afe1-2517b0beff59
"""
    apply_space_motion(c::ICRSCoords, pm_ra_cosdec, pm_dec, Δt)

Linearly propagate the sky position `c` along its proper motion over the time interval
`Δt`, returning the new position as an `ICRSCoords` object.
"""
function apply_space_motion(c::ICRSCoords, pm_ra_cosdec, pm_dec, Δt)
    Δra = pm_ra_cosdec * Δt / cos(c.dec)
    Δdec = pm_dec * Δt
    return ICRSCoords(c.ra + ustrip(u"rad", Δra), c.dec + ustrip(u"rad", Δdec))
end

# ╔═╡ 6af46369-2d72-4d47-afd9-3a8b4863d018
md"""
Because we recorded the reference epoch when we defined our bundle for HD 219829:
"""

# ╔═╡ 1aeaf43f-50dd-402a-a1a3-c6996aeeea80
hd219829_coord.obstime

# ╔═╡ ea446b7c-da7e-469d-9d03-4b70ab58efe2
md"""
we can now propagate its position from this epoch back to (the Julian year) 1950:
"""

# ╔═╡ c73921c3-60f7-48a8-aad5-24682a798b95
hd219829_coord_1950 = apply_space_motion(
    hd219829_coord.coord,
    hd219829_coord.pm_ra_cosdec,
    hd219829_coord.pm_dec,
    (1950.0 - hd219829_coord.obstime)yr,
)

# ╔═╡ 5cdfccb8-397e-4e5d-b811-cd1b518126b9
md"""
Let's now plot our predicted position for this source as it would appear in 1950 based on the *Gaia* position and proper motion:
"""

# ╔═╡ 64b1b0aa-5e61-480d-a9b7-bdad0570b2f2
let
    fig = Figure(size = (600, 600))

    ax = Axis(
        fig[1, 1];
        xlabel = "x [pix]",
        ylabel = "y [pix]",
        title = "DSS image around HD 219829",
        aspect = DataAspect(),
    )

    image!(ax, dss_img)

    c = hd219829_coord.coord
    x, y = world_to_pixel(dss_wcs, rad2deg.([c.ra, c.dec]))
    scatter!(
        ax, x, y;
        markersize = 30,
        color = :transparent,
        strokecolor = :red,
        strokewidth = 2,
    )

    # Plot the predicted (past) position:
    x_1950, y_1950 = world_to_pixel(
        dss_wcs,
        rad2deg.([hd219829_coord_1950.ra, hd219829_coord_1950.dec]),
    )
    scatter!(
        ax, x_1950, y_1950;
        markersize = 30,
        color = :transparent,
        strokecolor = :dodgerblue,
        strokewidth = 2,
    )

    fig
end

# ╔═╡ a58365fa-b732-4dbf-870f-0faeff92adec
md"""
The red circle is the same as in the previous image and shows the position of the source in the *Gaia* catalog (in 2015.5). The blue circle shows our prediction for the position of the source in 1950 — this looks much closer to where the star is in the DSS image!
"""

# ╔═╡ 36b7063d-b54a-479f-9b0e-77190f3b62b9
md"""
In this tutorial, we have introduced how to store and transform velocity data along with positional data using SkyCoords.jl objects and a bit of our own Julia. We also demonstrated how to use the velocity of a source to predict its position at an earlier or later time.
"""

# ╔═╡ 057c9c59-13d7-4aea-a2ec-1bdc83548d7d
md"""
## Exercises

Lalande 21185 is the brightest red dwarf star in the northern hemisphere and has a pretty high proper motion. Use the [*Gaia* archive](https://gea.esac.esa.int/archive/) to find values and create a velocity-data bundle (like `hd219829_coord` above) for Lalande 21185. (Hint: earlier in the tutorial, we extracted information from a *Gaia* table and mentioned which *Gaia* column names match with our position components.)

Two things to watch out for here: Lalande 21185 is one of a handful of very bright, fast-moving stars that did not make it into *Gaia* DR2, so we query the DR3 catalog (`gaiadr3.gaia_source`, reference epoch J2016.0) instead. Also, its proper motion is so large (almost 5 arcsec/year!) that by 2016 the star sat over an arcminute away from its J2000 catalog position, so we use a generous 3 arcmin search radius:
"""

# ╔═╡ 71d8d300-f8c5-45d5-9f51-bbf87a8153fe
lalande_row = let
    c = from_name("Lalande 21185")
    query = """
    SELECT TOP 1
        source_id, ra, dec, parallax, pmra, pmdec, phot_g_mean_mag, ref_epoch
    FROM gaiadr3.gaia_source
    WHERE 1 = CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', $(rad2deg(c.ra)), $(rad2deg(c.dec)), $(3 / 60))
    )
    ORDER BY phot_g_mean_mag ASC
    """
    only(eachrow(tap_query("https://gea.esac.esa.int/tap-server/tap", query)))
end

# ╔═╡ 92f0656d-a2a9-4b1a-9b58-d30f73c17c95
lalande_coord = (;
    coord = ICRSCoords(lalande_row.ra * °, lalande_row.dec * °),
    distance = parallax_to_distance(lalande_row.parallax) |> us"Constants.pc",
    pm_ra_cosdec = lalande_row.pmra * us"mas/yr",
    pm_dec = lalande_row.pmdec * us"mas/yr",
    obstime = lalande_row.ref_epoch,
)

# ╔═╡ 7f41fafe-793c-4a28-be64-1587cfa62e02
md"""
# Notebook setup 🔧
"""

# ╔═╡ 2858ead6-0a96-4b82-b721-8d5154660bc7
TableOfContents(; title = "On this page", depth = 4)

# ╔═╡ e43f434d-57c7-447b-b92d-55600797ddad
function frontmatter(path)
    prefix = "#> "
    is_fm = startswith(prefix)
    block = Iterators.takewhile(is_fm, Iterators.dropwhile(!is_fm, eachline(path)))
    toml = TOML.parse(join(chopprefix.(block, prefix), "\n"))
    return toml["frontmatter"]
end

# ╔═╡ c7dc7a08-783e-41ae-8ecd-13c0bdcc6bb1
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ 2afe550d-7e72-49ae-825d-4888a497a62f
md"""
# Astronomical Coordinates 3: Working with Velocity Data

This notebook is modified from <https://learn.astropy.org/tutorials/3_Coordinates-Velocities.html>

_Original authors: Adrian Price-Whelan, Saima Siddiqui, Luthien Liu, Zihao Chen_

!!! tip "Learning Goals"
    - Introduce how to represent and transform velocity data along with sky coordinates
    - Demonstrate how to predict the position of a star at a different time using its proper motion

$(keywords())
"""

# ╔═╡ Cell order:
# ╟─2afe550d-7e72-49ae-825d-4888a497a62f
# ╟─3b782148-7df3-42e9-a30d-add6ecac4e76
# ╟─d3db4bd5-c2f5-4be9-97ea-a9f4d1fbfd56
# ╠═955ad972-038b-4519-9eeb-b5a297b4d185
# ╟─cc4fab2a-8a85-4dd8-833f-523d675597e3
# ╠═e5e0904a-eae9-47b0-8758-4b0a7d5b142f
# ╟─e2e15684-5dcf-4889-b966-ef5fc59a982f
# ╠═f2d636e1-547f-4f96-aa86-ae3c16c11166
# ╟─cdc27d8a-3150-4269-befc-9c7002a2a39c
# ╠═ae4aa69a-73b4-4d95-b8bb-ee3765132c80
# ╟─a8f2133d-e970-47d2-b244-5b665f21210d
# ╠═8291ed3e-4f11-481c-be48-3c50ccaa3839
# ╠═7b3581dd-01c3-4621-a988-750a8ec4feb8
# ╟─08193019-6c81-4f83-a89f-c2d288192004
# ╠═911b0b50-ea1f-4d70-bc40-73a21b94bd2c
# ╟─84ea56f6-e90a-4744-8136-f0d35cabfa0b
# ╠═f37980ea-708f-42f2-b247-5d31ba033fb1
# ╟─87d2b737-f1c8-4182-8603-b46e97858b43
# ╠═e0331849-286b-4284-9c74-b9d3a3afd122
# ╠═a3cd666d-cf0f-4d88-875e-e66cb1d8a99f
# ╟─ecd79657-1c56-42fd-b63b-39d9f248dcc4
# ╠═e6dcd0f7-9833-4f98-9197-feea74baad3e
# ╠═cb37c0d1-7c22-4467-a69e-b984beff4ee6
# ╟─408ba1ad-8856-4b92-989e-9fc885762944
# ╟─3892c3dc-c2eb-42a4-b2b4-46f76b1bc475
# ╟─0d4d2809-9580-4fec-8f3a-766e6e7b9a25
# ╠═6a7ef5e6-18d5-4c9c-9562-bbf32c90e411
# ╠═fc06f86b-2286-4aa1-a11a-2b0238416e8f
# ╠═acc6f854-eaeb-419b-8bec-51382e06d8b6
# ╟─61b57f57-3529-4bd3-9ae4-716793858266
# ╠═75a973b2-cb8b-45fb-bdc3-b66c32b1bc50
# ╟─0990f7a6-c401-4de9-841e-d240b460265e
# ╠═d67675df-7291-466e-b40d-7e90ff4d99e0
# ╟─d5ff195c-9e28-4e90-8250-4d1d9ea2f136
# ╠═a8e57417-9a4e-49e9-9d9d-73d8d8938071
# ╠═1e3179fe-05b3-4da6-ac4e-1991f4866158
# ╟─bb23a77f-2996-49f3-8bca-a7aa25b55e8e
# ╠═ef024267-f69b-4e12-9aa5-710e12c5d520
# ╠═cdbff55e-ae03-44d0-8e7e-c60d6f224368
# ╟─5b2c80cc-e869-4e95-b301-3e83b5b061db
# ╠═4f8bd7bc-e292-48c7-a4f2-1c8798f7b4e8
# ╟─f36ca985-f94b-4a5d-a0c4-5e1c13e09c4f
# ╠═5cfa1e48-91a3-433d-8119-10de94c94e0b
# ╟─8f1d2372-7ece-4b72-8dc8-29ba3ef5f5d2
# ╠═4c56f658-6735-416d-892a-8aad7ed9a733
# ╟─f1f08c80-0603-4c92-b66d-c7e0a426d413
# ╠═51f6ac42-b6dd-4e50-bd43-f1215ddc1ff0
# ╟─41bcfaae-0dd6-4dcf-a8d2-98b8df1b96d0
# ╠═0a87c508-7149-48db-afe1-2517b0beff59
# ╟─6af46369-2d72-4d47-afd9-3a8b4863d018
# ╠═1aeaf43f-50dd-402a-a1a3-c6996aeeea80
# ╟─ea446b7c-da7e-469d-9d03-4b70ab58efe2
# ╠═c73921c3-60f7-48a8-aad5-24682a798b95
# ╟─5cdfccb8-397e-4e5d-b811-cd1b518126b9
# ╠═64b1b0aa-5e61-480d-a9b7-bdad0570b2f2
# ╟─a58365fa-b732-4dbf-870f-0faeff92adec
# ╟─36b7063d-b54a-479f-9b0e-77190f3b62b9
# ╟─057c9c59-13d7-4aea-a2ec-1bdc83548d7d
# ╠═71d8d300-f8c5-45d5-9f51-bbf87a8153fe
# ╠═92f0656d-a2a9-4b1a-9b58-d30f73c17c95
# ╟─7f41fafe-793c-4a28-be64-1587cfa62e02
# ╠═2858ead6-0a96-4b82-b721-8d5154660bc7
# ╟─e43f434d-57c7-447b-b92d-55600797ddad
# ╟─c7dc7a08-783e-41ae-8ecd-13c0bdcc6bb1
# ╠═bfd71b9c-ca7e-4480-92d3-c16fc03358d2
