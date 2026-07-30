### A Pluto.jl notebook ###
# v0.2.6

#> [frontmatter]
#> title = "Spectroscopic Data Reduction Part 1: Tracing"
#> layout = "layout.jlhtml"
#> date = "2025-12-31"
#> description = "Derive a spectroscopic trace model and extract a 1D spectrum."
#> tags = ["spectroscopy"]

using Markdown
using InteractiveUtils

# ╔═╡ 865165c9-ff3d-4f61-9fd1-4e682be885ed
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
            Pkg.PackageSpec(; name = "Images"),
            Pkg.PackageSpec(; name = "StatsBase"),
            Pkg.PackageSpec(; name = "Optim"),
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

    using CairoMakie
    using CairoMakie: Axis
    using Downloads: download
    using Images
    using StatsBase: mean, median, weights
    using Optim: optimize, minimizer, NelderMead

    using Makie: ClosedInterval as CI
    set_theme!(Image = (; colormap = :cividis))

    deps_ready = true
end

# ╔═╡ 77cf408c-81de-483f-90ee-619842b67d73
begin
    deps_ready

    using TOML: TOML
    using PlutoUI: TableOfContents, Dump
end

# ╔═╡ 794e3ad1-c003-4016-87b8-4a21539373f3
md"""
## Summary
This tutorial will walk through the derivation of a spectroscopic trace model and extraction.

A spectroscopic trace is the path of a point source (star) spectrum through a two-dimensional dispersed spectrum. The trace is needed because, in general, spectra are not perfectly aligned with the axes of a detector.
"""

# ╔═╡ 0b11ec87-3d67-4aed-81d8-23f711fda033
md"""
### Packages 📦
"""

# ╔═╡ 0a4d0e77-864d-4b67-883e-754341a220a6
md"""
## Step 1: Examine the spectrum
"""

# ╔═╡ 1bc961ee-6383-418f-ad59-1f32a4073a24
md"""
We'll work with a 2D spectrum that contains no attached metadata, so we have to infer many of the features ourselves.

All we know is that this is a spectrum of a star, Aldebaran.

Our data are in the form of `.bmp` (bitmap) files, so we need a package ([Images.jl](http://juliaimages.org/)) to open them.

While `.bmp` files are not astronomical standard FITS files, as are commonly delivered from professional observatories, image formats like `.bmp`, `.jpg`, `.raw`, `.png`, etc. produced by consumer cameras may also be used for spectroscopy.

In this case, our images are monochromatic, which is similar to standard FITS images.
"""

# ╔═╡ 07f1c80c-f8ce-4db9-b581-f8cf4dc00382
img_data = load(download("https://raw.githubusercontent.com/astropy/astropy-tutorials/b2fb15153efc44c14bbc7deab36aec52f58d1367/tutorials/SpectroscopicDataReductionBasics/aldebaran_3s_1.bmp"))

# ╔═╡ e0ddd7dd-0db9-4219-9154-43f149eb785d
md"""
Below are some details about our image:
"""

# ╔═╡ 71cb2b8b-48a6-4339-b096-3a9666e702d6
typeof(img_data)

# ╔═╡ 96bbbbb3-4fd5-4a1b-8285-89dd79dc8220
eltype(img_data) |> Dump

# ╔═╡ 3a3290ba-35f2-47be-9597-7ca7b696bfee
# Image is monochrome
red.(img_data) == green.(img_data) == blue.(img_data)

# ╔═╡ 06d967b8-cc56-4edc-91ce-670a2362c4be
red.(img_data)

# ╔═╡ ff8d460e-04c6-474f-b975-b50d04572e79
md"""
Our data are unsigned 8-bit integers (0-255) representing a monochromatic image.

Since all three color channels are identical, we can work with just one of them. `rawview` from Images.jl lets us view the underlying 8-bit integer values of that channel as a plain array, verifying that it is indeed 2-dimensional:
"""

# ╔═╡ 41bd8731-4c7f-43a7-967f-50731b3133b3
img_array = rawview(blue.(img_data))

# ╔═╡ 79fd7099-5b53-4b7e-8b5a-f43aea125aa7
let
    fig, ax, p = image(img_array')
    Colorbar(fig[1, 2], p)
    fig
end

# ╔═╡ f22689b3-3f91-47eb-9eec-79a7f2b19810
md"""
!!! note "Why does the plot seem "flipped"?"
    Our plot above appears reflected about the x-axis relative to `img_data` because it follows the default convention of placing the origin in the lower left-hand corner, while the image displayed in the notebook places it in the top left corner.

    For more on data vs. image orientation in Makie.jl plots, see the [discussion here](https://github.com/MakieOrg/Makie.jl/issues/205).
"""

# ╔═╡ b2e89b4d-358a-4dd4-a71e-007d01043953
md"""
The main goal of the trace is to obtain a model ``f(x)`` defining the vertical position of the light (the signal) along the detector.

We're going to start by assuming that wavelength dispersion is in the X-direction and the Y-direction is entirely spatial.

This is an approximation made by inspecting the image by eye.
"""

# ╔═╡ d69c30fb-b9e1-4dab-9429-3e5daaba74fa
md"""
## Step 2a. Try to find the spine to trace using argmax
"""

# ╔═╡ 60e76790-a3a5-4461-a5e2-f3be7990d404
md"""
To obtain the trace, we first measure the Y-value at each X-value. we'll start with the trivial approach of using `argmax`:
"""

# ╔═╡ 7993a64d-fc7a-4b28-85a2-4bef95e6ecd1
yvals = argmax.(eachcol(img_array))

# ╔═╡ 0f7fd20f-5f27-434b-a030-2422695225e9
xvals = axes(img_data, 2)

# ╔═╡ fb28eae5-60da-4433-afd8-2e5cc7ef21e3
scatter(
    xvals, yvals;
    marker = :x,
    axis = (xlabel = "X position", ylabel = "Argmax trace data"),
)

# ╔═╡ 7e55acf3-740f-493c-99c5-75302d018c1e
md"""
There's a pretty clear line going through the center, which represents our signal, but there are also a lot of erroneous data points.

We can get rid of most of the bad data just by filtering it out using a pixel mask:
"""

# ╔═╡ c75e9173-a478-4ecc-9851-8a8550836aab
bad_pixels = @. (yvals < 400) | (yvals > 500)

# ╔═╡ 0fbb7cf4-9c7a-436c-9725-182e2502137c
let
    fig = Figure()

    ax = Axis(
        fig[1, 1];
        xlabel = "X position",
        ylabel = "Argmax trace data",
    )

    scatter!(ax, xvals, yvals; marker = :x)
    scatter!(
        ax, xvals[bad_pixels], yvals[bad_pixels];
        marker = :x,
        color = :red,
    )

    fig
end

# ╔═╡ 1551faad-1ba3-4376-b718-e8a36b37f0e0
scatter(
    xvals[.!bad_pixels], yvals[.!bad_pixels];
    marker = :x,
    axis = (xlabel = "X position", ylabel = "Argmax trace data"),
)

# ╔═╡ 00243c18-5f81-46d9-807e-c1e74e94c95c
md"""
We can be a little more precise by "zooming in" along the y-axis, so we refine the mask again to be over a narrower range:
"""

# ╔═╡ b79f8c06-48cd-42b6-887f-ed4cc31e12a1
bad_pixels_tighter = @. (yvals < 425) | (yvals > 460)

# ╔═╡ 7f8937b5-fb4e-4c18-b119-ee8bf4204a81
scatter(
    xvals[.!bad_pixels_tighter], yvals[.!bad_pixels_tighter];
    marker = :x,
    axis = (xlabel = "X position", ylabel = "Argmax trace data"),
)

# ╔═╡ 5360a684-b064-4353-9f06-2b3b44a16f8b
md"""
The stuff at x > 1100 looks bad, but there's still signal out there. We can see there is clear signal out to nearly pixel ~1400:

"""

# ╔═╡ 4ce20eff-0704-48e5-b6b9-d73512177ee4
slice_window = 425:475

# ╔═╡ d324c2ad-e9c1-4bec-aa4b-7305466263ae
img_slice = img_array[slice_window, :]

# ╔═╡ e25f71c5-af24-4a5c-a85a-33dfe657df32
function aspect!(fig)
    rowsize!(fig.layout, 1, Aspect(1, 1 / 3.0))
    resize_to_layout!(fig)
    return fig
end

# ╔═╡ 126da840-d43d-46cd-b2ff-e181587d6358
let
    fig, ax, p = image(img_slice')
    aspect!(fig)
end

# ╔═╡ f2cdbcf8-2956-4564-a7e6-3d56b2bdd9a9
md"""
## Step 2b: Use moment analysis to extract a spine to trace

We can use [moments](https://en.wikipedia.org/wiki/Moment_(mathematics)) to provide a different, possibly better, estimate of where the trace's center is. The advantage of moment analysis is that we're using all of the data to estimate the vertical position, not just the single brightest value, which is what we used above.

Note that we need to subtract off the background to avoid a bias toward the center, so we use the median of each image column as our background estimate.

The first-order moment is the intensity-weighted mean position:

```math
m_1 = \frac{\sum_i x_i f(x_i)}{\sum_i f(x_i)}
```

where ``x_i`` is the position and ``f(x_i)`` is the intensity at that position. ``f(x_i)`` must be zero in the signal-free region for ``m_1`` to return an accurate estimate of the location of the peak.
"""

# ╔═╡ f600d2f5-b890-4fef-9be5-e2a8517e96cb
background = median(img_array; dims = 1)

# ╔═╡ a665f59d-b446-4ff3-abb2-5e400b8eead4
weighted_yaxis_values = map(eachcol(img_slice .- background)) do col_weights
    mean(slice_window, weights(col_weights))
end

# ╔═╡ 6d2c22cb-2be2-44da-a51e-b97bed017888
scatter(
    xvals, weighted_yaxis_values;
    marker = :x,
    axis = (xlabel = "X position", ylabel = "Moment-1 estimated Y-value trace"),
)

# ╔═╡ 8180ac7f-a742-410f-baf2-f57ead7374c8
md"""
Overplot the "weighted" centroid locations on the data to verify they look reasonable:
"""

# ╔═╡ 87cdc911-a696-41ef-9ae0-6eee2f4496ed
begin
    # Pixel centers sit on the integer coordinates, so the image edges
    # extend half a pixel past the first and last centers
    pixel_edges(r) = CI(first(r) - 0.5, last(r) + 0.5)

    plot_slice(img, window; kwargs...) = image(
        pixel_edges(axes(img, 2)), pixel_edges(window), img';
        axis = (xlabel = "X position", ylabel = "Y position"),
        kwargs...,
    )

    plot_slice!(ax, img, window; kwargs...) = image!(
        ax, pixel_edges(axes(img, 2)), pixel_edges(window), img';
        kwargs...,
    )
end

# ╔═╡ 9f592e18-dabe-42d7-ba56-71236b65acfb
let
    fig, ax, p = plot_slice(img_slice, slice_window)

    scatter!(
        ax, xvals, weighted_yaxis_values;
        marker = :x,
        color = :lightgreen,
        alpha = 0.5,
    )

    aspect!(fig)
end

# ╔═╡ 49d0b7fd-36d2-4fd4-9481-6ed6a771f9a4
md"""

We can also compare the argmax and weighted approaches. They agree well at x < 1200, but there are simply more points from the weighted approach at x > 1200:
"""

# ╔═╡ 74a57803-7ec0-46da-9043-f5e79d6a5f49
let
    fig = Figure()

    ax = Axis(fig[1, 1]; xlabel = "X position", ylabel = "Y position")

    scatter!(
        ax, xvals, weighted_yaxis_values;
        marker = :x, alpha = 0.5, label = "Weighted"
    )

    scatter!(
        ax, xvals[.!bad_pixels_tighter], yvals[.!bad_pixels_tighter];
        marker = :cross, alpha = 0.5, label = "Argmax"
    )

    axislegend(ax; position = :lt)

    fig
end

# ╔═╡ 585a2b45-921a-49e9-ba6c-f33ceb45770a
md"""
That's a decent set of data, we'll use the moments instead of the argmax. There's still some data to flag out, though:
"""

# ╔═╡ cee07b11-9ea2-4306-bc3a-5a9b778c8b32
bad_moments = @. (weighted_yaxis_values > 460) | (weighted_yaxis_values < 430)

# ╔═╡ b139a69b-6819-4033-af56-1c0492327a8b
let
    fig = Figure()

    ax = Axis(fig[1, 1]; xlabel = "X position", ylabel = "Y position")

    scatter!(
        ax, xvals[.!bad_moments], weighted_yaxis_values[.!bad_moments];
        marker = :x, alpha = 0.5, label = "Weighted"
    )

    scatter!(
        ax, xvals[.!bad_pixels_tighter], yvals[.!bad_pixels_tighter];
        marker = :cross, alpha = 0.5, label = "Argmax"
    )

    axislegend(ax; position = :lt)

    fig
end

# ╔═╡ b89f999b-5dab-4243-ba5d-c04edef61753
md"""
## Step 3. Fit the trace profile

We want a model ``f(x)`` that gives the y-value of the centroid as a function of x.
"""

# ╔═╡ da811ced-bd73-454f-9676-e88e05b89cc7
md"""
We start by fitting a 2nd-order polynomial. This is a linear least-squares problem, so just like in the [learn.JuliaAstro > Modeling 1: Linear model fitting](https://learn.juliaastro.org/tutorials/1_models-linear/) tutorial, we can solve it by building a design matrix and using Julia's matrix division operator (`\`):
"""

# ╔═╡ 607c441c-7711-49b9-8ccd-80d00cd74c33
begin
    # Least-squares polynomial fit and its corresponding model function
    polyfit(x, y, degree) = [xi^p for xi in x, p in 0:degree] \ y
    polymodel(coeffs) = x -> evalpoly(x, coeffs)
end

# ╔═╡ 5a4c75a3-12bb-4ffd-bf47-2de82f34b6e2
trace_coeffs_deg2 = polyfit(xvals[.!bad_moments], weighted_yaxis_values[.!bad_moments], 2)

# ╔═╡ 24c2c441-7e36-49e9-8213-5e936ddd0f77
let
    fig, ax, p = scatter(
        xvals[.!bad_moments], weighted_yaxis_values[.!bad_moments];
        marker = :x, alpha = 0.5,
        axis = (xlabel = "X position", ylabel = "Y position"),
    )

    lines!(ax, xvals, polymodel(trace_coeffs_deg2).(xvals); color = :red)

    fig
end

# ╔═╡ fe8f537d-a3fe-4d8b-a460-5a833c78ba81
md"""
We plot and examine the residuals to visually inspect whether the fit is good:
"""

# ╔═╡ 70ad2443-9224-40c6-95ce-3496f78cf286
scatter(
    xvals[.!bad_moments],
    weighted_yaxis_values[.!bad_moments] .- polymodel(trace_coeffs_deg2).(xvals[.!bad_moments]);
    marker = :x,
    axis = (xlabel = "X position", ylabel = "Residual (data-model)"),
)

# ╔═╡ 762254d0-7b3a-4d5e-a688-58ede61dece3
md"""
The curvature seen at the left is a sign of a suboptimal fit. Specifically, curvature in the residual indicates that we need to use a _higher order_ model - i.e., we need more terms in the polynomial. We change `degree=2` to `degree=3`:
"""

# ╔═╡ c9d87700-0414-494e-93e4-3655c3a4bcb1
trace_coeffs_deg3 = polyfit(xvals[.!bad_moments], weighted_yaxis_values[.!bad_moments], 3)

# ╔═╡ aa8a180c-5ba0-4135-b2c4-7fe39411e414
let
    fig, ax, p = scatter(
        xvals[.!bad_moments], weighted_yaxis_values[.!bad_moments];
        marker = :x, alpha = 0.5,
        axis = (xlabel = "X position", ylabel = "Y position"),
    )
    lines!(ax, xvals, polymodel(trace_coeffs_deg3).(xvals); color = :red)
    fig
end

# ╔═╡ b572d01a-e26a-43d5-9a50-0b3b4a2b5cae
scatter(
    xvals[.!bad_moments],
    weighted_yaxis_values[.!bad_moments] .- polymodel(trace_coeffs_deg3).(xvals[.!bad_moments]);
    marker = :x,
    axis = (xlabel = "X position", ylabel = "Residual (data-model)"),
)

# ╔═╡ b6bd7a78-76fa-4db6-a40c-029811fc099f
md"""
Arguably, we should toss out the data at >1400 pixels since there's no clear signal there. We'll come back to this...

Again, we should verify the trace by overplotting on the original data:
"""

# ╔═╡ 7d4aaeca-25be-429f-a0c5-42b99528917a
let
    fig, ax, p = plot_slice(img_slice, slice_window)

    lines!(ax, xvals, polymodel(trace_coeffs_deg3).(xvals); color = :white)

    aspect!(fig)
end

# ╔═╡ 566652eb-38bb-44f6-877c-c6351f238c4f
md"""
Seeing the curve up in the model to the right -- which we do not observe in the data -- suggests we should re-fit without including the x > 1200 data at all:
"""

# ╔═╡ e01eb8b1-be5c-4949-abc8-92aa22cacc84
good_fit_region = @. !bad_moments && (xvals < 1200)

# ╔═╡ 22195b75-6ba7-476c-a660-c89a4a1ffe3d
trace_coeffs = polyfit(
    xvals[good_fit_region],
    weighted_yaxis_values[good_fit_region],
    3
)

# ╔═╡ fdc324fe-45ea-44e8-abc5-00dfbf874a65
md"""
We now have a satisfactory fit:
"""

# ╔═╡ 21f21c25-2c84-46ad-8c28-606ba3dee005
let
    fig, ax, p = plot_slice(img_slice, slice_window)

    lines!(ax, xvals, polymodel(trace_coeffs).(xvals); color = :white)

    aspect!(fig)
end

# ╔═╡ 36d4f00f-ac37-4f68-a4d8-9ed3416e6e58
let
    fig, ax, p = scatter(
        xvals[.!bad_moments], weighted_yaxis_values[.!bad_moments];
        marker = :x, alpha = 0.5,
        axis = (xlabel = "X position", ylabel = "Y position"),
    )

    lines!(ax, xvals, polymodel(trace_coeffs).(xvals); color = :red)

    fig
end

# ╔═╡ dc07bba2-44a8-4f71-85a4-50a7e2843214
let
    residuals = weighted_yaxis_values .- polymodel(trace_coeffs).(xvals)

    beyond_fit_region = @. !bad_moments && (xvals >= 1200)

    fig = Figure()

    ax = Axis(
        fig[1, 1];
        xlabel = "X position",
        ylabel = "Residual (data-model)",
    )

    scatter!(
        ax, xvals[good_fit_region], residuals[good_fit_region];
        marker = :x,
    )

    scatter!(
        ax, xvals[beyond_fit_region], residuals[beyond_fit_region];
        marker = :cross, color = :red, alpha = 0.5,
    )

    ylims!(ax, -5, 5)

    fig
end

# ╔═╡ 3e5be619-d3f7-4873-9a64-ca91a5bc400b
md"""
## Step 4. Obtain a trace profile

Now we can extract the data along that trace.

We want to take a "profile" of the trace to see how many pixels on either side of the line we should include:
"""

# ╔═╡ 61f53e5f-f305-4053-aba7-f0edf83ea54d
let
    fig, ax, p = plot_slice(img_slice, slice_window)

    trace = polymodel(trace_coeffs).(xvals)

    band!(ax, xvals, trace .- 15, trace .+ 15; color = (:lightgreen, 0.2))

    aspect!(fig)
end

# ╔═╡ 478365cb-402e-4ca8-bc08-3890d5c438e3
begin
    # Start by taking +/- 15 pixels
    npixels_to_cut = 15

    trace_center = polymodel(trace_coeffs).(xvals)

    cutouts = stack(
        img_array[(round(Int, yval) - npixels_to_cut):(round(Int, yval) + npixels_to_cut), ii]
            for (yval, ii) in zip(trace_center, xvals)
    )
end

# ╔═╡ bc95367e-2562-4bac-8a2a-c9b63eb1c42e
md"""
That last step deserves some explanation:
```julia
cutouts = stack(
    image_array[round(Int, yval) - npixels_to_cut:round(Int, yval) + npixels_to_cut, ii]
    for (yval, ii) in zip(trace_center, xvals)
)
```

- `zip(trace_center, xvals)` takes each trace y-value and each x-value and 'zips' them together, so each iteration of the generator has one x, y pair
- `image_array[round(Int, yval) - npixels_to_cut:round(Int, yval) + npixels_to_cut, ii]` is taking a single pixel along the x-direction (the second dimension, `ii`) and a range of pixels along the y-direction, i.e., `y ± n`
- `stack` then combines these vertical cutouts into a matrix with one column per x-pixel

We can see the result visually:
"""

# ╔═╡ 5de96dd1-cc54-44f7-80fa-e6e5a42d69c5
let
    fig = Figure()

    ax1 = Axis(fig[1, 1]; title = "We go from this...")

    plot_slice!(ax1, img_slice, slice_window)

    ax2 = Axis(fig[2, 1]; title = "...to this")

    plot_slice!(ax2, cutouts, -npixels_to_cut:npixels_to_cut)

    aspect!(fig)
end

# ╔═╡ c1de2f91-62e2-4a2e-aaf5-a9b5d048ef8f
md"""
Then we average along the X-direction to get the trace profile:
"""

# ╔═╡ 768e0d76-add8-4669-aea3-6c03e9562fe5
mean(cutouts .- background; dims = 1)

# ╔═╡ 39c61b6b-0500-4bbc-a735-833635251806
begin
    trace_profile_xaxis = -npixels_to_cut:npixels_to_cut

    mean_trace_profile = vec(mean(cutouts .- background, dims = 2))

    lines(
        trace_profile_xaxis, mean_trace_profile;
        axis = (;
            xlabel = "Distance from center",
            ylabel = "Average source profile",
        ),
    )
end

# ╔═╡ 136afb89-0145-4e46-ac5c-11dac62565b5
md"""
We want to fit that profile with a Gaussian for future use. This is a nonlinear fit, so we define our own Gaussian model and solve for its parameters with [Optim.jl](https://julianlsolvers.github.io/Optim.jl/stable/)'s Nelder-Mead solver — the same algorithm we used via Optimization.jl in the [learn.JuliaAstro > Modeling 1: Linear model fitting](https://learn.juliaastro.org/tutorials/1_models-linear/) tutorial:

!!! todo
    Switch to Optimization.jl (as in the Modeling 1 tutorial) once the temporary Makie pin in the Packages cell is removed; that pin restricts SciMLBase to versions incompatible with recent Optimization.jl releases.
"""

# ╔═╡ aed503aa-dac5-4385-b3c1-d419aabb03dc
begin
    gaussian(x, A, μ, σ) = A * exp(-(x - μ)^2 / (2σ^2))

    function trace_profile_objective(u)
        A, μ, σ = u
        residuals = @. mean_trace_profile - gaussian(trace_profile_xaxis, A, μ, σ)
        return sum(abs2, residuals)
    end

    # Initial guess: amplitude = profile max, centered at 0, width of 5 pixels
    guess = [maximum(mean_trace_profile), 0.0, 5.0]

    fitted_trace_profile = minimizer(optimize(trace_profile_objective, guess, NelderMead()))
end

# ╔═╡ 3dbe6187-30d3-4751-be15-5909f7093294
model_trace_profile = gaussian.(trace_profile_xaxis, fitted_trace_profile...)

# ╔═╡ 09bd9a92-b575-4ce9-938e-9588773aaf09
let
    fig = Figure()

    ax = Axis(
        fig[1, 1];
        xlabel = "Distance from center", ylabel = "Average source profile",
    )

    lines!(ax, trace_profile_xaxis, mean_trace_profile; label = "data")

    lines!(ax, trace_profile_xaxis, model_trace_profile; label = "model")

    axislegend(ax)

    fig
end

# ╔═╡ e844f47b-e418-4361-804b-f63ac16bfaad
md"""
Both the empirical trace profile `mean_trace_profile` and the modeled `model_trace_profile` can reasonably be used; the latter is more convenient to serialize (i.e., write to disk or on paper).
"""

# ╔═╡ 60a5c6f6-e701-4f16-9fe0-f55156995ed1
md"""
## Step 5. Extract the traced spectrum

We can obtain our spectrum by directly averaging the pixels along the trace:
"""

# ╔═╡ c078171c-15a7-493d-bd76-72e97ed5c3fc
begin
    average_spectrum = vec(mean(cutouts .- background, dims = 1))

    lines(average_spectrum)
end

# ╔═╡ 4cab94a9-c19a-47ba-9a8d-823e20caf9ae
md"""
Or, we can obtain our spectrum by taking the trace-weighted average:
"""

# ╔═╡ 52fbc382-a66b-4776-989f-d70663601a77
trace_avg_spectrum = [
    mean(cutout .- bg, weights(mean_trace_profile))
        for (cutout, bg) in zip(eachcol(cutouts), background)
]

# ╔═╡ fdc85139-8047-4434-b059-4572717adefa
md"""
We can also do this with the Gaussian weights:
"""

# ╔═╡ 16a47d2b-978b-4815-9ed5-580c45d14d01
gaussian_trace_avg_spectrum = [
    mean(cutout .- bg, weights(model_trace_profile))
        for (cutout, bg) in zip(eachcol(cutouts), background)
]

# ╔═╡ 2fe72675-ba6a-4d40-bec0-1deba1b5663c
let
    fig = Figure()

    ax = Axis(
        fig[1, 1];
        xlabel = "X position", ylabel = "Extracted spectrum",
    )

    lines!(ax, average_spectrum; label = "Direct average")

    lines!(ax, trace_avg_spectrum; label = "Trace-weighted average")

    lines!(
        ax, gaussian_trace_avg_spectrum;
        label = "Gaussian-model-trace-weighted average",
        color = :red, alpha = 0.5, linewidth = 0.5,
    )

    axislegend(ax)

    fig
end

# ╔═╡ 73ab28a2-45de-4744-91dd-06f5feeed894
md"""
In general, the trace-weighted average will have higher signal-to-noise, as seen here (while we haven't measured the noise, it is approximately constant across the image).

Note that the Gaussian model and the direct trace yield nearly identical results
"""

# ╔═╡ 48ac66ec-7169-4ca9-87d6-481372fb55b7
md"""
## Step 6: Repeat for another star

In this last step, we go through all the above steps again for another star (Deneb), but with less explanation.
"""

# ╔═╡ aa2995fe-045e-47ff-ab1e-99edec484196
begin
    image_data2 = load(download("https://raw.githubusercontent.com/astropy/astropy-tutorials/b2fb15153efc44c14bbc7deab36aec52f58d1367/tutorials/SpectroscopicDataReductionBasics/deneb_3s_13.63g_1.bmp"))
    image_array2 = rawview(blue.(image_data2))
end

# ╔═╡ ebdc9660-6c34-4184-8616-9d3cc22035d5
let
    fig, ax, p = image(image_array2')
    Colorbar(fig[1, 2], p)
    fig
end

# ╔═╡ 7d4136f5-3d28-4246-a18e-af7ed3a1ef28
slice_window2 = 470:520

# ╔═╡ d0c3cc31-6bd9-497e-af22-6e16625b42b9
img_slice2 = image_array2[slice_window2, :]

# ╔═╡ a30cafd3-76ff-4b57-96c5-3a2ed7bb16b6
let
    fig, ax, p = plot_slice(
        img_slice2, slice_window2;
        colormap = :viridis
    )
    aspect!(fig)
end

# ╔═╡ 0d1c84ab-8283-41aa-906f-fbec8a8d77a0
begin
    background2 = median.(eachcol(image_array2))
    weighted_yaxis_values2 = map(eachcol(img_slice2 .- background2')) do col_weights
        mean(slice_window2, weights(col_weights))
    end
    trace_coeffs2 = polyfit(xvals, weighted_yaxis_values2, 3)
    trace_center2 = polymodel(trace_coeffs2).(xvals)
end;

# ╔═╡ f6d5f6a9-d270-4c2e-a46c-bdad74c99436
let
    fig, ax, p = scatter(
        xvals, weighted_yaxis_values2;
        marker = :x, alpha = 0.5,
        axis = (xlabel = "X position", ylabel = "Y position"),
    )
    lines!(ax, xvals, trace_center2; color = :red)
    fig
end

# ╔═╡ b23d2a6c-a7d7-4ccf-afe9-c15c5d05dc75
let
    fig, ax, p = plot_slice(img_slice2, slice_window2)

    scatter!(ax, xvals, weighted_yaxis_values2; marker = :cross, color = :white, alpha = 0.25)
    lines!(ax, xvals, trace_center2; color = :red)

    aspect!(fig)
end

# ╔═╡ 885b28ab-5113-458d-9683-17dd334a852a
spectrum2 = [
    mean(
            image_array2[(round(Int, yval) - npixels_to_cut):(round(Int, yval) + npixels_to_cut), ii] .- bg,
            weights(mean_trace_profile),
        )
        for (yval, ii, bg) in zip(trace_center2, xvals, background2)
]

# ╔═╡ 77210ab4-a8b0-4c2b-aec5-cc39c63e3447
lines(spectrum2)

# ╔═╡ 33d9ca20-bbd8-4906-bf55-386d111c19c8
md"""
In the next tutorial, Spectroscopic Data Reduction 2, we'll work on the wavelength calibration.
"""

# ╔═╡ 011d0d6b-5bc4-4399-80bf-63f5f401219b
md"""
# Notebook setup 🔧
"""

# ╔═╡ 88e74f38-c943-4757-82bf-96e39d67da0c
TableOfContents()

# ╔═╡ 77442012-d0fc-4e1d-994f-e96c63f74040
function frontmatter(path)
    prefix = "#> "
    is_fm = startswith(prefix)
    block = Iterators.takewhile(is_fm, Iterators.dropwhile(!is_fm, eachline(path)))
    toml = TOML.parse(join(chopprefix.(block, prefix), "\n"))
    return toml["frontmatter"]
end

# ╔═╡ ff0a7c90-cb17-4a19-ab07-09d50c9ab98d
function keywords(kind = "note", title = "Keywords")
    nb_path = split(@__FILE__, "#==#") |> first |> string
    tags = (nb_path |> frontmatter)["tags"]
    header = "!!! $kind \"$title\""
    body = join(("`$tag`" for tag in tags), " ")
    return Markdown.parse("$header\n    $body")
end

# ╔═╡ 59b4d5c4-6100-11ef-1aa0-892e00d83607
md"""
# Spectroscopic Data Reduction Part 1: Tracing

This notebook is modified from <https://learn.astropy.org/tutorials/1_SpectroscopicTraceTutorial.html>

_Original authors: Adam Ginsburg, Kelle Cruz, Lia Corrales, Jonathan Sick, Adrian Price-Whelan_

!!! tip "Learning Goals"

    - Open a two-dimensional spectrum from an image file (bitmap)
    - Fit a spectroscopic trace

$(keywords())
"""

# ╔═╡ Cell order:
# ╟─59b4d5c4-6100-11ef-1aa0-892e00d83607
# ╟─794e3ad1-c003-4016-87b8-4a21539373f3
# ╟─0b11ec87-3d67-4aed-81d8-23f711fda033
# ╠═865165c9-ff3d-4f61-9fd1-4e682be885ed
# ╟─0a4d0e77-864d-4b67-883e-754341a220a6
# ╟─1bc961ee-6383-418f-ad59-1f32a4073a24
# ╠═07f1c80c-f8ce-4db9-b581-f8cf4dc00382
# ╟─e0ddd7dd-0db9-4219-9154-43f149eb785d
# ╠═71cb2b8b-48a6-4339-b096-3a9666e702d6
# ╠═96bbbbb3-4fd5-4a1b-8285-89dd79dc8220
# ╠═3a3290ba-35f2-47be-9597-7ca7b696bfee
# ╠═06d967b8-cc56-4edc-91ce-670a2362c4be
# ╟─ff8d460e-04c6-474f-b975-b50d04572e79
# ╠═41bd8731-4c7f-43a7-967f-50731b3133b3
# ╠═79fd7099-5b53-4b7e-8b5a-f43aea125aa7
# ╟─f22689b3-3f91-47eb-9eec-79a7f2b19810
# ╟─b2e89b4d-358a-4dd4-a71e-007d01043953
# ╟─d69c30fb-b9e1-4dab-9429-3e5daaba74fa
# ╟─60e76790-a3a5-4461-a5e2-f3be7990d404
# ╠═7993a64d-fc7a-4b28-85a2-4bef95e6ecd1
# ╠═0f7fd20f-5f27-434b-a030-2422695225e9
# ╠═fb28eae5-60da-4433-afd8-2e5cc7ef21e3
# ╟─7e55acf3-740f-493c-99c5-75302d018c1e
# ╠═c75e9173-a478-4ecc-9851-8a8550836aab
# ╠═0fbb7cf4-9c7a-436c-9725-182e2502137c
# ╠═1551faad-1ba3-4376-b718-e8a36b37f0e0
# ╟─00243c18-5f81-46d9-807e-c1e74e94c95c
# ╠═b79f8c06-48cd-42b6-887f-ed4cc31e12a1
# ╠═7f8937b5-fb4e-4c18-b119-ee8bf4204a81
# ╟─5360a684-b064-4353-9f06-2b3b44a16f8b
# ╠═4ce20eff-0704-48e5-b6b9-d73512177ee4
# ╠═d324c2ad-e9c1-4bec-aa4b-7305466263ae
# ╠═126da840-d43d-46cd-b2ff-e181587d6358
# ╟─e25f71c5-af24-4a5c-a85a-33dfe657df32
# ╟─f2cdbcf8-2956-4564-a7e6-3d56b2bdd9a9
# ╠═f600d2f5-b890-4fef-9be5-e2a8517e96cb
# ╠═a665f59d-b446-4ff3-abb2-5e400b8eead4
# ╠═6d2c22cb-2be2-44da-a51e-b97bed017888
# ╟─8180ac7f-a742-410f-baf2-f57ead7374c8
# ╠═9f592e18-dabe-42d7-ba56-71236b65acfb
# ╠═87cdc911-a696-41ef-9ae0-6eee2f4496ed
# ╟─49d0b7fd-36d2-4fd4-9481-6ed6a771f9a4
# ╠═74a57803-7ec0-46da-9043-f5e79d6a5f49
# ╟─585a2b45-921a-49e9-ba6c-f33ceb45770a
# ╠═cee07b11-9ea2-4306-bc3a-5a9b778c8b32
# ╠═b139a69b-6819-4033-af56-1c0492327a8b
# ╟─b89f999b-5dab-4243-ba5d-c04edef61753
# ╟─da811ced-bd73-454f-9676-e88e05b89cc7
# ╠═607c441c-7711-49b9-8ccd-80d00cd74c33
# ╠═5a4c75a3-12bb-4ffd-bf47-2de82f34b6e2
# ╠═24c2c441-7e36-49e9-8213-5e936ddd0f77
# ╟─fe8f537d-a3fe-4d8b-a460-5a833c78ba81
# ╠═70ad2443-9224-40c6-95ce-3496f78cf286
# ╟─762254d0-7b3a-4d5e-a688-58ede61dece3
# ╠═c9d87700-0414-494e-93e4-3655c3a4bcb1
# ╠═aa8a180c-5ba0-4135-b2c4-7fe39411e414
# ╠═b572d01a-e26a-43d5-9a50-0b3b4a2b5cae
# ╟─b6bd7a78-76fa-4db6-a40c-029811fc099f
# ╠═7d4aaeca-25be-429f-a0c5-42b99528917a
# ╟─566652eb-38bb-44f6-877c-c6351f238c4f
# ╠═e01eb8b1-be5c-4949-abc8-92aa22cacc84
# ╠═22195b75-6ba7-476c-a660-c89a4a1ffe3d
# ╟─fdc324fe-45ea-44e8-abc5-00dfbf874a65
# ╠═21f21c25-2c84-46ad-8c28-606ba3dee005
# ╠═36d4f00f-ac37-4f68-a4d8-9ed3416e6e58
# ╠═dc07bba2-44a8-4f71-85a4-50a7e2843214
# ╟─3e5be619-d3f7-4873-9a64-ca91a5bc400b
# ╠═61f53e5f-f305-4053-aba7-f0edf83ea54d
# ╠═478365cb-402e-4ca8-bc08-3890d5c438e3
# ╟─bc95367e-2562-4bac-8a2a-c9b63eb1c42e
# ╠═5de96dd1-cc54-44f7-80fa-e6e5a42d69c5
# ╟─c1de2f91-62e2-4a2e-aaf5-a9b5d048ef8f
# ╠═768e0d76-add8-4669-aea3-6c03e9562fe5
# ╠═39c61b6b-0500-4bbc-a735-833635251806
# ╟─136afb89-0145-4e46-ac5c-11dac62565b5
# ╠═aed503aa-dac5-4385-b3c1-d419aabb03dc
# ╠═3dbe6187-30d3-4751-be15-5909f7093294
# ╠═09bd9a92-b575-4ce9-938e-9588773aaf09
# ╟─e844f47b-e418-4361-804b-f63ac16bfaad
# ╟─60a5c6f6-e701-4f16-9fe0-f55156995ed1
# ╠═c078171c-15a7-493d-bd76-72e97ed5c3fc
# ╟─4cab94a9-c19a-47ba-9a8d-823e20caf9ae
# ╠═52fbc382-a66b-4776-989f-d70663601a77
# ╟─fdc85139-8047-4434-b059-4572717adefa
# ╠═16a47d2b-978b-4815-9ed5-580c45d14d01
# ╠═2fe72675-ba6a-4d40-bec0-1deba1b5663c
# ╟─73ab28a2-45de-4744-91dd-06f5feeed894
# ╟─48ac66ec-7169-4ca9-87d6-481372fb55b7
# ╠═aa2995fe-045e-47ff-ab1e-99edec484196
# ╠═ebdc9660-6c34-4184-8616-9d3cc22035d5
# ╠═7d4136f5-3d28-4246-a18e-af7ed3a1ef28
# ╠═d0c3cc31-6bd9-497e-af22-6e16625b42b9
# ╠═a30cafd3-76ff-4b57-96c5-3a2ed7bb16b6
# ╠═0d1c84ab-8283-41aa-906f-fbec8a8d77a0
# ╠═f6d5f6a9-d270-4c2e-a46c-bdad74c99436
# ╠═b23d2a6c-a7d7-4ccf-afe9-c15c5d05dc75
# ╠═885b28ab-5113-458d-9683-17dd334a852a
# ╠═77210ab4-a8b0-4c2b-aec5-cc39c63e3447
# ╟─33d9ca20-bbd8-4906-bf55-386d111c19c8
# ╟─011d0d6b-5bc4-4399-80bf-63f5f401219b
# ╠═88e74f38-c943-4757-82bf-96e39d67da0c
# ╟─ff0a7c90-cb17-4a19-ab07-09d50c9ab98d
# ╟─77442012-d0fc-4e1d-994f-e96c63f74040
# ╠═77cf408c-81de-483f-90ee-619842b67d73
