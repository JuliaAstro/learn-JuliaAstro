### A Pluto.jl notebook ###
# v1.0.4

#> [frontmatter]
#> tags = ["spectroscopy", "radio astronomy", "query", "units", "images", "stairs", "colorbars"]
#> title = "Cube Reprojection Tutorial"
#> description = "Put two spectral cubes observed at different frequencies onto a common spectral and spatial grid with matching resolution"
#> date = "2026-10-22"
#> layout = "layout.jlhtml"
#> image = "/assets/spectral_cube_reproject.png"

using Markdown
using InteractiveUtils

# ╔═╡ 20362a1a-c083-4767-8d01-e231eefbd335
begin
    import Pkg
    Pkg.activate(; temp = true)
    Pkg.add(
        [
            Pkg.PackageSpec(; name = "DynamicQuantities"),
            Pkg.PackageSpec(; name = "Interpolations"),
            Pkg.PackageSpec(; name = "ImageFiltering"),
            Pkg.PackageSpec(; name = "LinearAlgebra"),
            Pkg.PackageSpec(; name = "PlutoUI"),
            Pkg.PackageSpec(; name = "ProgressLogging"),
            Pkg.PackageSpec(; name = "Downloads"),
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
                rev = "fits-wcs",
                url = "https://github.com/JuliaAstro/Reproject.jl/",
            ),
            Pkg.PackageSpec(;
                rev = "hdu-size",
                url = "https://github.com/JuliaAstro/FITSFiles.jl/",
            ),
            Pkg.PackageSpec(;
                url = "https://github.com/JuliaAstro/FITSWCS.jl/",
            ),
        ]
    )

    using AstroImages
    using CairoMakie
    using DynamicQuantities
    using Interpolations: interpolate, extrapolate, Gridded, Linear, NoInterp
    using ImageFiltering: imfilter, centered, NA
    using LinearAlgebra: Diagonal, Symmetric, cholesky, diag, eigen
    using ProgressLogging
    using Downloads: download
    using Reproject: reproject

    deps_ready = true
end

# ╔═╡ e3cc0abc-8bb5-4b60-9fc0-bddb46ca9cdb
begin
    deps_ready
    using PlutoUI: TableOfContents
end

# ╔═╡ ea3a159a-6174-11ef-0cb1-b73af72665d4
md"""
# Cube Reprojection Tutorial

This notebook is modified from <https://learn.astropy.org/tutorials/SpectralCubeReprojectExample.html>

_Original authors: Adam Ginsburg, Eric Koch_

!!! tip "Learning Goals"
    - reproject a cube spectrally
    - smooth it spectrally
    - reproject it spatially
"""

# ╔═╡ 7abaa07b-c97e-4ae1-9eaf-650fb529fe14
md"""
## Summary

Spectroscopic cube observations taken at different wavelengths can trace the motion of gas or stars using spectral lines, but often lines at different wavelengths give different information. For example, one might observe a galaxy in the 21cm line of HI and the 115 GHz line of CO, or a protoplanetary disk in a line of N2H+ and a line of CO, or a galactic disk in the Hα and Hβ lines (in absorption or emission). In order to compare these data sets pixel-by-pixel, they must be placed onto a common grid with common resolution.

This tutorial shows how to take two spectral cubes observed toward the same part of the sky, but at different frequencies, and put them onto the same grid. We load the cubes and handle their world coordinates with [AstroImages.jl](https://juliaastro.org/AstroImages/), resample their spectral axes with [Interpolations.jl](https://juliamath.github.io/Interpolations.jl/stable/), and resample their spatial axes with [Reproject.jl](https://github.com/JuliaAstro/Reproject.jl). A few helper functions at the bottom of this notebook stand in for the [spectral-cube](https://spectral-cube.readthedocs.io/) and [radio-beam](https://radio-beam.readthedocs.io/) methods that the original tutorial uses.

In this example, we do spectral smoothing and interpolation (step 4) before spatial smoothing and reprojection (steps 5 and 6), but if you have a varying-resolution cube (with a different beam size for each channel), you have to do spatial smoothing first. For more information see the [spectral-cube documentation](https://spectral-cube.readthedocs.io/en/latest/smoothing.html).
"""

# ╔═╡ 721cc524-00cf-423c-b1eb-5406a60a0912
md"""
### Packages
"""

# ╔═╡ dd726d8d-5823-4553-8216-3c1e5065219e
md"""
## Step 1: Download the data

We download two example spectral cubes of a point in the Galactic center from a permalink on the ALMA archives. These are moderately large files, with sizes 18 MB and 337 MB, so we save them to the `data` folder next to this notebook. That way, they only need to be downloaded once. The `pluto_download` function defined at the bottom of this notebook shows a progress bar while they download.
"""

# ╔═╡ 12bffd04-b8e1-4c02-bf1d-e577b3db92a1
md"""
## Step 2: Load the cubes

Besides their two spatial axes and their spectral axis, these ALMA cubes have a fourth axis of length one for the Stokes parameter (only Stokes I here), which we drop:
"""

# ╔═╡ 2289fee5-e9f5-4e6e-9a3d-3dfc6f856bb3
md"""
The cubes are at different frequencies:
"""

# ╔═╡ fa6738ff-92c7-4121-939e-a82bde55051d
md"""
The first cube covers the H₂CS 4(1,3)–3(1,2) line at 139.483699 GHz.

The second covers the SiO 5–4 line at 217.104984 GHz.
"""

# ╔═╡ fe301445-fa04-41d6-99aa-8e85a2f9d9c2
md"""
!!! todo
    We use the `find_lines` tool to query [splatalogue](https://splatalogue.online/) with TBD over the spectral range covered by the cube. It returns a table of matching lines. Note that some line names will be repeated because Splatalogue includes several different databases and most chemical species are present in all of these.
"""

# ╔═╡ 1536e606-0ae3-4d10-80c9-336d0adc3c3e
md"""
## Step 3: Convert cubes from frequency to velocity

To compare the kinematic structure of the target, we need to convert from the observed frequency (which must be in a common reference frame; in this case, it already is) to the doppler velocity.

We use the radio convention, ``v = c\,(1 - \nu/\nu_0)``, with the rest frequency ``\nu_0`` of each line. The `with_velocity_axis` function defined at the bottom of this notebook writes this velocity axis into a copy of each cube's header, like spectral-cube's `with_spectral_unit`:
"""

# ╔═╡ 25f62fe1-2985-4b52-8c9f-56953440ab1d
md"""
`pixel_to_world` now returns the radio velocity of a channel, so we can collect the velocities of all channels with the `velocity_axis` helper function:
"""

# ╔═╡ d4706c21-fa59-4533-9461-4ef26dd0f000
md"""
From these velocity ranges, we can see the H₂CS cube (`cube1vel`) is narrower in velocity, so we'll use that as the target spectral reprojection. However, the SiO cube (`cube2vel`) is the smaller footprint on the sky.
"""

# ╔═╡ 6fed040b-fb90-433f-be1b-af4d43f742b8
md"""
### Create spatial maps of the peak intensity to quickly explore the cubes

One way to quickly explore the structure in the data cubes is to produce a peak intensity map, or the maximum along the spectral axis (`Z`):
"""

# ╔═╡ 8473c04b-df19-4233-9cab-28ea3dfed108
md"""
We can do the same thing for the other cube:
"""

# ╔═╡ 5738bf09-cd54-493f-a2b3-9ece1c64a3e3
md"""
## Step 4: Spectral Interpolation
"""

# ╔═╡ 619c2904-39a2-40e1-8c15-5211eb0112cf
md"""
We can choose to do either the spatial or spectral step first.
In this case, we choose the spectral step first because the H₂CS cube is narrower in velocity (cube1vel) and this will reduce the number of channels we need to spatially interpolate over in the next step.
"""

# ╔═╡ 00cd85d9-9af0-4a46-8710-1f8314a68bc8
md"""
We need to match resolution to the cube with the largest channel width:
"""

# ╔═╡ 87429d88-58dd-43c6-a9ec-c53bacef0805
md"""
Next, we will reduce `cube2vel` to have the same spectral range as `cube1vel`:
"""

# ╔═╡ 735320ad-dd59-47bb-9208-54e9f410346c
md"""
Note that it is important for the to-be-interpolated cube, in this case `cube2vel`, to have pixels bounding the spectral axis of `cube1vel`, but in this case it does not, at either end. If the pixel range doesn't overlap perfectly, it may blank out one of the edge pixels. So, to fix this, we add a buffer of one channel at each end:
"""

# ╔═╡ 22f11580-de20-420a-ab19-13881633ecc5
md"""
Our H₂CS cube (`cube1vel`) has broader channels. We need to first smooth `cube2vel` to the broader channel width before resampling it onto the channels of `cube1vel`. To do this, we spectrally smooth with a Gaussian whose width is chosen such that smoothing `cube2vel` results in the same width as `cube1vel`:
"""

# ╔═╡ 628cafa2-9959-4025-8986-ba09df1a43e8
md"""
We smooth each spectrum channel by channel, so we convert this FWHM to the standard deviation of the Gaussian, in channels of `cube2vel`:

!!! note
    The original tutorial passes the standard deviation in km/s directly to `Gaussian1DKernel`, which expects it in channels, so it smooths `cube2vel` less (σ ≈ 0.80 channels).
"""

# ╔═╡ b237961f-f3e8-4332-a381-07c5f7dae8a8
md"""
We then smooth with the `spectral_smooth` function defined at the bottom of this notebook. This does 420 × 420 = 176400 smoothing operations on a 222-channel spectrum:
"""

# ╔═╡ 05e7ef56-3164-48eb-a09f-7ba0fef2e3df
md"""
Now that we've done spectral smoothing, we can resample the spectral axis of `cube2vel_smooth` to match `cube1vel` by linearly interpolating `cube2vel_smooth` onto the channels of `cube1vel`:
"""

# ╔═╡ cfa20b5d-b172-4422-9e0d-81fc8028f6d9
md"""
## Step 5: Spatial Smoothing
"""

# ╔═╡ 590877e7-d2f3-41ab-955f-abcbeed6b7f2
md"""
Now that we’ve done spectral smoothing, we also need to follow a similar procedure of smoothing then resampling for the spatial axes.
"""

# ╔═╡ cb169320-4ca3-43a6-af86-a9432ad8dc02
md"""
The `beam` is the resolution element of our cubes. The `beam` function defined at the bottom of this notebook reads it from the `BMAJ`, `BMIN`, and `BPA` header keywords:
"""

# ╔═╡ dba82ba0-60db-45f0-90c6-738ad6692913
md"""
`cube1` again has the larger beam, so we'll smooth `cube2` to its resolution.
"""

# ╔═╡ 7ea93647-2a92-45cf-91e0-e8504457bf0f
md"""
### Aside: mixed beams
"""

# ╔═╡ 8e6564ce-5869-434c-9034-bf0dddece65f
md"""
If cube1 and cube2 had different sized beams, but neither was clearly larger, we would have to convolve _both_ to a [common beam](https://radio-beam.readthedocs.io/en/latest/commonbeam.html#finding-the-smallest-common-beam).

In this case, it's redundant and we could have just used `cube1`'s beam, but this is the more general approach. The `common_beam` function defined at the bottom of this notebook finds the smallest beam that both beams can be convolved to:
"""

# ╔═╡ 356c9582-bd46-44a6-bcaf-f2edf68057e0
md"""
This is just `beam1`, since `beam2` fits inside it. Two elongated beams at right angles to each other, on the other hand, need a round common beam:
"""

# ╔═╡ 439357d1-0cad-45c5-8285-0e9e38b7a501
md"""
We then convolve `cube2` to this beam with the `convolve_to` function defined at the bottom of this notebook. Since the cubes are in units of Jy/beam, it also scales the intensities by the ratio of the beam areas, like spectral-cube does:
"""

# ╔═╡ 780d4bc8-3e27-4fb2-a8c1-bec8df15b6c8
md"""
## Step 6: Reprojection
"""

# ╔═╡ 65cb16b6-9ace-4f9e-9c7c-4aac7877cffd
md"""
Now we can do the spatial resampling as the final step for producing two cubes matched to the same spatial and spectral pixel grid.

Reproject.jl reprojects two-dimensional images, so the `reproject_cube` function defined at the bottom of this notebook reprojects `cube2` one channel at a time onto the matching channel of `cube1vel`:
"""

# ╔═╡ 63c1cc87-1a58-4347-9112-b3edd4ae99cc
md"""
These two cubes are now on an identical grid, and can be directly compared:
"""

# ╔═╡ 0e2daa33-648e-4271-a779-a3143e872d08
md"""
Their spectra can now be overplotted, as they are in the same unit with the same beam. Here we compare them at the center of the map:
"""

# ╔═╡ a126afc1-14d0-4e6c-a109-dff801e992e7
md"""
At this position, both lines peak near 40 km/s, but the SiO line is much broader than the H₂CS line.
"""

# ╔═╡ 9ce43d1c-16f5-4ee0-ab7c-f15ee48dd630
md"""
# Notebook setup 🔧
"""

# ╔═╡ 9aa6d05f-4a00-412e-ac5d-5bb973ca112d
TableOfContents()

# ╔═╡ 003b15c8-7e25-4a6c-8849-13b637410213
"""
    pluto_download(url, path)

Download `url` to `path` with `Downloads.download`, and show a progress bar in Pluto.

The file is downloaded to a temporary location first, so an interrupted download never leaves a partial file at `path`.
"""
function pluto_download(url, path)
    return @withprogress name = basename(path) begin
        # Pluto keeps every log message, so only log whole-percent steps
        last_percent = -1
        function progress(total, now)
            percent = total > 0 ? 100 * now ÷ total : 0
            if percent > last_percent
                last_percent = percent
                @logprogress percent / 100
            end
        end

        mv(download(url; progress), path)
    end
end

# ╔═╡ d8f76e52-3670-4287-85f4-c99e18f71dc4
filename1 = let
    fname = "member.uid___A001_X1465_X3a33.BrickMaser_sci.spw71.cube.I.manual.image.pbcor.fits"
    local_path = joinpath(@__DIR__, "data", fname)
    isfile(local_path) ? local_path : pluto_download("https://almascience.eso.org/dataPortal/$(fname)", local_path)
end

# ╔═╡ 1534acf5-1b39-4a02-819b-0ba304f70d84
cube1 = dropdims(load(filename1); dims = 4)

# ╔═╡ 772d38fd-20a2-4043-83ab-344b64d74d64
header(cube1)

# ╔═╡ 0d811af2-14aa-4f68-8892-50be1938372d
filename2 = let
    fname = "member.uid___A001_X87d_X141.a_sma1_sci.spw27.cube.I.pbcor.fits"
    local_path = joinpath(@__DIR__, "data", fname)
    isfile(local_path) ? local_path : pluto_download("https://almascience.eso.org/dataPortal/$(fname)", local_path)
end

# ╔═╡ 56ff05c6-365b-4cf0-b16f-c0b9288aaf75
cube2 = dropdims(load(filename2); dims = 4)

# ╔═╡ 5c84410d-8c15-4eed-b389-48fb497ff8fe
header(cube2)

# ╔═╡ e6b7b27e-fb00-4b31-b6df-3e1fb778a2d1
(cube1[:CRVAL3], cube2[:CRVAL3]) .* u"Hz" .|> us"GHz"

# ╔═╡ 4d4255bc-f211-4228-80f8-678c6371f821
md"""
!!! todo
    Find a home for these helper functions upstream. They stand in for these methods that the original tutorial uses:
    - `with_velocity_axis`: spectral-cube's `with_spectral_unit` with `velocity_convention = "radio"`
    - `spectral_smooth` and `spectral_interpolate`: spectral-cube's `spectral_smooth` and `spectral_interpolate`
    - `beam`, `common_beam`, and `deconvolve`: radio-beam's `Beam`, `Beam.commonbeam_with`, and `Beam.deconvolve`
    - `convolve_to` and `reproject_cube`: spectral-cube's `convolve_to` and `reproject`

    Unlike those methods, they assume a single beam per cube, and spatial axes that are aligned with RA and Dec.
"""

# ╔═╡ e65709a7-99ff-4698-8169-4cfe2862aafe
begin
    """
        with_velocity_axis(cube, rest_frequency)

    Copy of `cube` whose spectral axis is the radio velocity, v = c (1 - ν / ν₀), of the line at `rest_frequency`. `cube` must have a linear frequency axis in Hz.
    """
    function with_velocity_axis(cube, rest_frequency)
        ν₀ = ustrip(rest_frequency |> us"Hz")
        c = ustrip(Constants.c |> us"m/s")
        out = copyheader(cube, parent(cube))
        out["CTYPE3"] = "VRAD"
        out["CRVAL3"] = c * (1 - cube["CRVAL3"] / ν₀)
        out["CDELT3"] = -c * cube["CDELT3"] / ν₀
        out["CUNIT3"] = "m/s"
        out["RESTFRQ"] = ν₀
        return out
    end

    """
        velocity_axis(cube)

    Radio velocity of each channel of `cube`, which must have a velocity axis (see `with_velocity_axis`).
    """
    velocity_axis(cube) = [pixel_to_world(cube, [1, 1, k])[3] for k in 1:size(cube, Z)] * u"m/s" .|> us"km/s"

    """
        nanconvolve(A, kernel)

    Convolve `A` with the symmetric `kernel` (odd size along each axis), skipping NaNs and pixels beyond the edges of `A`. Pixels
    that are NaN in `A` stay NaN.
    """
    function nanconvolve(A, kernel)
        out = imfilter(eltype(A), A, centered(kernel), NA())
        out[isnan.(A)] .= NaN
        return out
    end

    """
        spectral_smooth(cube, σ)

    Smooth each spectrum of `cube` with a Gaussian of standard deviation `σ` channels, cut off at 4σ.
    """
    function spectral_smooth(cube, σ)
        kernel = [exp(-k^2 / 2σ^2) for k in -ceil(Int, 4σ):ceil(Int, 4σ)]
        return copyheader(cube, nanconvolve(parent(cube), reshape(kernel, 1, 1, :)))
    end

    """
        spectral_interpolate(cube, target)

    Linearly interpolate the spectra of `cube` onto the velocity channels of `target`.
    """
    function spectral_interpolate(cube, target)
        A = parent(cube)
        v_in, v_out = ustrip.(velocity_axis(cube)), ustrip.(velocity_axis(target))
        itp = extrapolate(
            interpolate((axes(A, 1), axes(A, 2), v_in), A, (NoInterp(), NoInterp(), Gridded(Linear()))),
            NaN,
        )
        data = eltype(A).([itp(i, j, v) for i in axes(A, 1), j in axes(A, 2), v in v_out])
        out = AstroImage(data, deepcopy(header(cube)))
        for key in ("CTYPE3", "CRVAL3", "CDELT3", "CRPIX3", "CUNIT3")
            out[key] = target[key]
        end
        return out
    end
end

# ╔═╡ 12c1c71f-355f-41e9-a4e9-41353e48f649
cube1vel = with_velocity_axis(cube1, 139.483699u"GHz");

# ╔═╡ a50b4db4-2f21-4343-a747-fb0bbb5fb481
let
    peak = dropdims(maximum(cube1vel; dims = Z); dims = Z)
    fig, _ = implotview(
        peak;
        cmap = :viridis,
        colorbar_label = "Peak intensity (Jy/beam)",
        axis = (; title = "H₂CS 4(1,3)–3(1,2)"),
    )
    fig
end

# ╔═╡ 9f1656ab-a779-4b9d-9367-067b352cc174
cube2vel = with_velocity_axis(cube2, 217.104984u"GHz");

# ╔═╡ bbbcb89b-88d9-461e-87c4-cfd51089f106
let
    peak = dropdims(maximum(cube2vel; dims = Z); dims = Z)
    fig, _ = implotview(
        peak;
        cmap = :viridis,
        colorbar_label = "Peak intensity (Jy/beam)",
        axis = (; title = "SiO 5–4"),
    )
    fig
end

# ╔═╡ 1b64c447-4e6c-41cb-82d5-5f6441fc3371
vel1 = velocity_axis(cube1vel)

# ╔═╡ 34f7bd96-be78-44de-8493-e71012fe90aa
vel2 = velocity_axis(cube2vel)

# ╔═╡ 9dbf33c3-6b06-40d4-922e-7b95044c911a
extrema(vel1), extrema(vel2)

# ╔═╡ aa4af422-79b5-489f-be59-406aa899f798
velocity_res_1, velocity_res_2 = vel1[2] - vel1[1], vel2[2] - vel2[1]

# ╔═╡ 2dba98a0-f827-4089-a80d-3146349479d9
fwhm_gaussian = sqrt(velocity_res_1^2 - velocity_res_2^2)

# ╔═╡ 3861bfc0-e44b-4327-998d-f887c14d7da4
σ_channels = let
    fwhm_to_sigma = sqrt(8 * log(2))
    ustrip(fwhm_gaussian / abs(velocity_res_2)) / fwhm_to_sigma
end

# ╔═╡ 8e6269ed-306b-4e15-b7ab-4a58e3406325
let
    channels = findall(v -> minimum(vel1) ≤ v ≤ maximum(vel1), vel2)
    extrema(vel1), extrema(vel2[channels])
end

# ╔═╡ 6cd61a39-b384-4bd4-988a-3f17a8cbce5d
cutout_channels = let
    inside = findall(v -> minimum(vel1) - velocity_res_2 ≤ v ≤ maximum(vel1) + velocity_res_2, vel2)
    first(inside):last(inside)
end

# ╔═╡ 4d6f1df3-3300-45d3-8e6d-35928f2e0749
cube2vel_cutout = cube2vel[Z = cutout_channels];

# ╔═╡ 885c0a19-82d4-4575-968b-9edf8d4b73cf
extrema(velocity_axis(cube2vel_cutout))

# ╔═╡ c4b42c3f-b0d7-4d83-8efa-6fce52df205e
cube2vel_smooth = spectral_smooth(cube2vel_cutout, σ_channels);

# ╔═╡ f544ef92-f259-4e18-b7e5-3645de3c6fa4
cube2vel_spectralresample = spectral_interpolate(cube2vel_smooth, cube1vel);

# ╔═╡ 38bbdb00-d2cd-4a04-880e-f5f88c12e181
size(cube2vel_spectralresample)

# ╔═╡ e2922e65-12cb-4cce-9728-44444f10ebb2
begin
    """
        beam(img)

    The restoring beam in the header of `img`: the FWHM of its major and minor axes, and its position angle east of north.
    """
    beam(img) = (;
        major = img["BMAJ"] * u"deg" |> us"arcsec",
        minor = img["BMIN"] * u"deg" |> us"arcsec",
        pa = img["BPA"] * u"deg" |> us"deg",
    )

    # Covariance matrix (in arcsec², along east and north) of the Gaussian `beam`
    function beam_covariance(beam)
        σ_major, σ_minor = ustrip.((beam.major, beam.minor) .|> us"arcsec") ./ sqrt(8 * log(2))
        θ = ustrip(beam.pa |> us"deg")
        major_axis, minor_axis = [sind(θ), cosd(θ)], [cosd(θ), -sind(θ)]
        return σ_major^2 * major_axis * major_axis' + σ_minor^2 * minor_axis * minor_axis'
    end

    # Inverse of `beam_covariance`
    function beam_from_covariance(Σ)
        (λ_minor, λ_major), V = eigen(Symmetric(Σ))
        east, north = V[:, 2]
        return (;
            major = sqrt(λ_major) * sqrt(8 * log(2)) * us"arcsec",
            minor = sqrt(max(λ_minor, 0)) * sqrt(8 * log(2)) * us"arcsec",
            pa = (mod(atand(east, north) + 90, 180) - 90) * us"deg",
        )
    end

    """
        common_beam(beam1, beam2)

    The smallest beam that both `beam1` and `beam2` can be convolved to.
    """
    function common_beam(beam1, beam2)
        # In coordinates where `beam1` is a circle, `beam2` is an ellipse with semi-axes √λ
        # (relative to the circle). The smallest ellipse containing both has semi-axes √max(λ, 1).
        L = cholesky(Symmetric(beam_covariance(beam1))).L
        λ, Q = eigen(Symmetric(L \ beam_covariance(beam2) / L'))
        return beam_from_covariance(L * Q * Diagonal(max.(λ, 1)) * Q' * L')
    end

    """
        deconvolve(beam1, beam2)

    The Gaussian that `beam2` must be convolved with to give `beam1`.
    """
    function deconvolve(beam1, beam2)
        Σ = beam_covariance(beam1) - beam_covariance(beam2)
        minimum(eigen(Symmetric(Σ)).values) > -1.0e-10 ||
            throw(ArgumentError("$beam2 cannot be convolved to $beam1"))
        return beam_from_covariance(Σ)
    end

    """
        convolve_to(cube, new_beam)

    Smooth `cube` (in Jy/beam) spatially, from the beam in its header to the larger `new_beam`.
    """
    function convolve_to(cube, new_beam)
        old_beam = beam(cube)
        # From offsets (east, north) in arcsec to offsets in pixels
        J = Diagonal(1 ./ (3600 .* [cube["CDELT1"], cube["CDELT2"]]))
        Σ = J * beam_covariance(deconvolve(new_beam, old_beam)) * J
        r = ceil.(Int, 4 .* sqrt.(diag(Σ)))
        kernel = [exp(-[x, y]' * (Σ \ [x, y]) / 2) for x in -r[1]:r[1], y in -r[2]:r[2], _ in 1:1]
        data = nanconvolve(parent(cube), kernel)
        # A larger beam collects more flux
        data .*= ustrip((new_beam.major * new_beam.minor) / (old_beam.major * old_beam.minor))
        out = copyheader(cube, data)
        out["BMAJ"], out["BMIN"], out["BPA"] = ustrip.((new_beam.major, new_beam.minor, new_beam.pa) .|> us"deg")
        return out
    end
end

# ╔═╡ 6facfe35-2d72-45a4-88a8-351070440ce3
beam1, beam2 = beam(cube1vel), beam(cube2vel_spectralresample)

# ╔═╡ f07f71a0-c5bb-4fdf-bdb2-b2cfee9fd214
target_beam = common_beam(beam1, beam2)

# ╔═╡ dea7a814-8631-4851-a771-59dcdf6be112
common_beam(
    (; major = 2.0us"arcsec", minor = 1.0us"arcsec", pa = 0.0us"deg"),
    (; major = 2.0us"arcsec", minor = 1.0us"arcsec", pa = 90.0us"deg"),
)

# ╔═╡ f2f259a1-2d1f-4ca5-896d-9308619b5ef8
cube2vel_spatialspectralsmooth = convolve_to(cube2vel_spectralresample, target_beam);

# ╔═╡ 39a53be4-0e6b-4224-a4d3-50de05e686f4
beam(cube2vel_spatialspectralsmooth)

# ╔═╡ e5fb7fc9-4753-45a0-b884-5da29206339f
"""
    reproject_cube(cube, target)

Reproject each channel of `cube` onto the pixel grid of the same channel of `target`.
"""
function reproject_cube(cube, target)
    data = similar(parent(target))
    for k in 1:size(target, Z)
        data[:, :, k] = first(reproject(cube[Z = k], target[Z = k]))
    end
    out = copyheader(target, data)
    out["RESTFRQ"] = cube["RESTFRQ"] # Same velocities, but a different line
    return out
end

# ╔═╡ f2adf9fd-a826-4dc8-80b0-8bf897eb1ac7
cube2vel_reproj = reproject_cube(cube2vel_spatialspectralsmooth, cube1vel);

# ╔═╡ c41c098c-1f2d-4e64-8473-917acea5cd58
size(cube2vel_reproj), size(cube1vel)

# ╔═╡ 32b3f31e-e4de-4ded-b989-582fb39cdbc3
let
    fig = Figure(size = (900, 400))
    for (pos, cube, title) in (
            (fig[1, 1], cube1vel, "H₂CS 4(1,3)–3(1,2)"),
            (fig[1, 2], cube2vel_reproj, "SiO 5–4, reprojected"),
        )
        peak = dropdims(maximum(cube; dims = Z); dims = Z)
        implotview(
            pos, peak;
            cmap = :viridis,
            colorbar_label = "Peak intensity (Jy/beam)",
            axis = (; title),
        )
    end
    fig
end

# ╔═╡ 0e0c0d60-7f99-442d-956b-868d0dec825d
let
    fig = Figure()

    ax = Axis(fig[1, 1]; xlabel = "Radio velocity", ylabel = "Intensity (Jy/beam)")

    stairs!(ax, vel1, collect(cube1vel[X = 126, Y = 126]); label = "H₂CS 4(1,3)–3(1,2)")

    stairs!(ax, vel1, collect(cube2vel_reproj[X = 126, Y = 126]); label = "SiO 5–4")

    axislegend(ax)

    fig
end

# ╔═╡ Cell order:
# ╟─ea3a159a-6174-11ef-0cb1-b73af72665d4
# ╟─7abaa07b-c97e-4ae1-9eaf-650fb529fe14
# ╟─721cc524-00cf-423c-b1eb-5406a60a0912
# ╠═20362a1a-c083-4767-8d01-e231eefbd335
# ╟─dd726d8d-5823-4553-8216-3c1e5065219e
# ╠═d8f76e52-3670-4287-85f4-c99e18f71dc4
# ╠═0d811af2-14aa-4f68-8892-50be1938372d
# ╟─12bffd04-b8e1-4c02-bf1d-e577b3db92a1
# ╠═1534acf5-1b39-4a02-819b-0ba304f70d84
# ╠═772d38fd-20a2-4043-83ab-344b64d74d64
# ╠═56ff05c6-365b-4cf0-b16f-c0b9288aaf75
# ╠═5c84410d-8c15-4eed-b389-48fb497ff8fe
# ╟─2289fee5-e9f5-4e6e-9a3d-3dfc6f856bb3
# ╠═e6b7b27e-fb00-4b31-b6df-3e1fb778a2d1
# ╟─fa6738ff-92c7-4121-939e-a82bde55051d
# ╟─fe301445-fa04-41d6-99aa-8e85a2f9d9c2
# ╟─1536e606-0ae3-4d10-80c9-336d0adc3c3e
# ╠═12c1c71f-355f-41e9-a4e9-41353e48f649
# ╠═9f1656ab-a779-4b9d-9367-067b352cc174
# ╟─25f62fe1-2985-4b52-8c9f-56953440ab1d
# ╠═1b64c447-4e6c-41cb-82d5-5f6441fc3371
# ╠═34f7bd96-be78-44de-8493-e71012fe90aa
# ╠═9dbf33c3-6b06-40d4-922e-7b95044c911a
# ╟─d4706c21-fa59-4533-9461-4ef26dd0f000
# ╟─6fed040b-fb90-433f-be1b-af4d43f742b8
# ╠═a50b4db4-2f21-4343-a747-fb0bbb5fb481
# ╟─8473c04b-df19-4233-9cab-28ea3dfed108
# ╠═bbbcb89b-88d9-461e-87c4-cfd51089f106
# ╟─5738bf09-cd54-493f-a2b3-9ece1c64a3e3
# ╟─619c2904-39a2-40e1-8c15-5211eb0112cf
# ╟─00cd85d9-9af0-4a46-8710-1f8314a68bc8
# ╠═aa4af422-79b5-489f-be59-406aa899f798
# ╟─87429d88-58dd-43c6-a9ec-c53bacef0805
# ╠═8e6269ed-306b-4e15-b7ab-4a58e3406325
# ╟─735320ad-dd59-47bb-9208-54e9f410346c
# ╠═6cd61a39-b384-4bd4-988a-3f17a8cbce5d
# ╠═4d6f1df3-3300-45d3-8e6d-35928f2e0749
# ╠═885c0a19-82d4-4575-968b-9edf8d4b73cf
# ╟─22f11580-de20-420a-ab19-13881633ecc5
# ╠═2dba98a0-f827-4089-a80d-3146349479d9
# ╟─628cafa2-9959-4025-8986-ba09df1a43e8
# ╠═3861bfc0-e44b-4327-998d-f887c14d7da4
# ╟─b237961f-f3e8-4332-a381-07c5f7dae8a8
# ╠═c4b42c3f-b0d7-4d83-8efa-6fce52df205e
# ╟─05e7ef56-3164-48eb-a09f-7ba0fef2e3df
# ╠═f544ef92-f259-4e18-b7e5-3645de3c6fa4
# ╠═38bbdb00-d2cd-4a04-880e-f5f88c12e181
# ╟─cfa20b5d-b172-4422-9e0d-81fc8028f6d9
# ╟─590877e7-d2f3-41ab-955f-abcbeed6b7f2
# ╟─cb169320-4ca3-43a6-af86-a9432ad8dc02
# ╠═6facfe35-2d72-45a4-88a8-351070440ce3
# ╟─dba82ba0-60db-45f0-90c6-738ad6692913
# ╟─7ea93647-2a92-45cf-91e0-e8504457bf0f
# ╟─8e6564ce-5869-434c-9034-bf0dddece65f
# ╠═f07f71a0-c5bb-4fdf-bdb2-b2cfee9fd214
# ╟─356c9582-bd46-44a6-bcaf-f2edf68057e0
# ╠═dea7a814-8631-4851-a771-59dcdf6be112
# ╟─439357d1-0cad-45c5-8285-0e9e38b7a501
# ╠═f2f259a1-2d1f-4ca5-896d-9308619b5ef8
# ╠═39a53be4-0e6b-4224-a4d3-50de05e686f4
# ╟─780d4bc8-3e27-4fb2-a8c1-bec8df15b6c8
# ╟─65cb16b6-9ace-4f9e-9c7c-4aac7877cffd
# ╠═f2adf9fd-a826-4dc8-80b0-8bf897eb1ac7
# ╟─63c1cc87-1a58-4347-9112-b3edd4ae99cc
# ╠═c41c098c-1f2d-4e64-8473-917acea5cd58
# ╠═32b3f31e-e4de-4ded-b989-582fb39cdbc3
# ╟─0e2daa33-648e-4271-a779-a3143e872d08
# ╠═0e0c0d60-7f99-442d-956b-868d0dec825d
# ╟─a126afc1-14d0-4e6c-a109-dff801e992e7
# ╟─9ce43d1c-16f5-4ee0-ab7c-f15ee48dd630
# ╠═9aa6d05f-4a00-412e-ac5d-5bb973ca112d
# ╠═e3cc0abc-8bb5-4b60-9fc0-bddb46ca9cdb
# ╠═003b15c8-7e25-4a6c-8849-13b637410213
# ╟─4d4255bc-f211-4228-80f8-678c6371f821
# ╠═e65709a7-99ff-4698-8169-4cfe2862aafe
# ╠═e2922e65-12cb-4cce-9728-44444f10ebb2
# ╠═e5fb7fc9-4753-45a0-b884-5da29206339f
