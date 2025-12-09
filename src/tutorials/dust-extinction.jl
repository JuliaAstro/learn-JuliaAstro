### A Pluto.jl notebook ###
# v0.20.21

#> [frontmatter]
#> title = "Dust Extinction"
#> layout = "layout.jlhtml"
#> date = "2025-11-25"
#> description = "Analyzing dust extinction."
#> tags = ["dust", "dust extinction", "units", "photometry", "extinction", "physics", "observational astronomy"]
#> 
#>     [[frontmatter.author]]
#>     name = "Ian Weaver"
#>     url = "https://github.com/icweaver"

using Markdown
using InteractiveUtils

# ╔═╡ b6b27fe2-c7f7-11f0-8b42-052bc2026e99
begin
	import Pkg
	Pkg.activate(Base.current_project())
    Pkg.instantiate()

	# Analysis
	using DustExtinction
	
	# Data handling
	using DataFramesMeta: DataFrame, Not, @rsubset
	using VirtualObservatory: execute, TAPService
	using FITSFiles: fits, info
	using Downloads: download
	using CodecZlib: GzipDecompressor
	
	# Units
	using DynamicQuantities: DynamicQuantities, @u_str, @us_str
	using DynamicQuantities.Constants: c as c0

	# Plots
	using Makie: Cycled, DQConversion, IntervalsBetween
	using CairoMakie: Figure, Axis, axislegend, lines, lines!, scatter!, ylims!
	using LaTeXStrings: @L_str
	using MathTeXEngine: set_texfont_family!, FontFamily
	set_texfont_family!(FontFamily("TeXGyreHeros"))
end;

# ╔═╡ 89366506-4854-4859-b109-2d6cd9456bb5
using DynamicQuantities: @register_unit

# ╔═╡ e7266de8-ba60-45b3-8988-6143e7098b72
using DynamicQuantities: ustrip

# ╔═╡ 95d05163-b6d2-4a01-a782-6da05479ee0a
using PlutoUI: TableOfContents

# ╔═╡ f18df5a0-1dda-4371-aea6-8ecbce67908c
md"""
# Dust Extinction

This notebook is modified from <https://learn.astropy.org/tutorials/color-excess.html>

!!! tip ""
	## Learning Goals
	
	* Investigate extinction curve shapes

	* Deredden spectral energy distributions and spectra

	* Calculate photometric extinction and reddening

	* Calculate synthetic photometry for a dust-reddened star

	* Convert from frequency to wavelength


!!! note ""
	## Keywords

	dust, dust extinction, units, photometry, extinction, physics, observational astronomy

!!! warning ""
	## Summary

	
"""

# ╔═╡ 2c3cabef-29fd-41cf-b3df-94f58748c22b
md"""
### Imports
"""

# ╔═╡ 555cf1bd-5b52-4d26-aaa2-033a9f4d831a
begin
	@register_unit erg exp10(-7) * u"J"
	@register_unit Å exp10(-10) * u"m"
	@register_unit Jy exp10(-26) * u"W/m^2/Hz"
end

# ╔═╡ 3680c62b-7ce8-4fbc-8cd0-249d499e990f
md"""
!!! todo
	Upstreaming to DQ.jl <https://github.com/JuliaPhysics/DynamicQuantities.jl/pull/196>
"""

# ╔═╡ f3bf7784-b2da-4ed6-9cb5-0799d81d65ae
md"""
## Example 1: Investigate Extinction Models
"""

# ╔═╡ ed46db90-9802-4615-928f-bcad4bd19804
md"""
!!! todo
	Add DQ support to DustExtinction.jl
"""

# ╔═╡ 0250218a-6670-483f-8122-045ca65d28e4
# Avoid floating point issue right at boundary
wav = let
	ϵ = eps()
	( (0.1 + ϵ) : 0.001 : 3.0 )u"μm"
end

# ╔═╡ eb1c1b16-c391-45a5-81c2-5fd3dbb12c50
let
	fig = Figure()
	
	ax = Axis(fig[1, 1];
		dim1_conversion = DQConversion(us"1/μm"),
		xlabel = "λ⁻¹",
		ylabel = "A(λ) / A(V)",
	)

	inv_wav = inv.(collect(wav))

	for model in (CCM89, F99)
		for Rv in (2.0, 3.0, 4.0)
			ext = model(; Rv)
			lines!(inv_wav, ext.(ustrip(u"Å", wav));
				label = string(nameof(model), ", ", "R = ", Rv),
			)
		end
	end

	axislegend(; position = :lt)
	
	fig
end

# ╔═╡ 943661a3-3b40-4c01-a624-79726e56b385
md"""
## Example 2: Deredden a spectrum
"""

# ╔═╡ 01195fc5-9ab9-4b9d-a4f9-a23f41bc6f2a
md"""
### MAST IUE Spectrum

!!! note
	For more on querying MAST's TAP service, see: <https://mast.stsci.edu/vo-tap/>
"""

# ╔═╡ d57bbe7f-b21a-479a-8249-a96d6aac3e95
df_spectra = execute(
	 TAPService("https://mast.stsci.edu/vo-tap/api/v0.1/caom/"),
	 """
	SELECT *
	FROM ivoa.ObsCore
	WHERE target_name = 'HD 147933'
	AND CONTAINS(
		POINT('ICRS', s_ra, s_dec),
		CIRCLE('ICRS', 246.396, -23.447, .00028)
	) = 1
	AND dataproduct_type = 'spectrum'
	 """
) |> DataFrame

# ╔═╡ 70cefdd5-4afb-4e09-b2df-75983bb40e85
md"""
Let's look at a spectrum from `lwr05639`:
"""

# ╔═╡ 285b4bc6-9689-4995-8e8e-dc6acfbd7f81
@rsubset df_spectra :obs_id == "lwr05639" && :access_format == "image/fits"

# ╔═╡ 207bfd62-2332-413d-8034-15d9927ea57e
md"""
Which we can download, uncompress, and load the fits file from the `access_url` provided in the table above. Let's do this for one of the spectra entries:
"""

# ╔═╡ dfb6b2d7-efcb-4505-9f69-d51eafb7db9f
# Download
fpath = download("http://archive.stsci.edu/pub/iue/data/lwr/05000/lwr05639.mxlo.gz")

# ╔═╡ 85b04476-2ce0-49dd-9bef-8f62c819e021
function decompress_gz_to_iobuffer(filepath::String)
    # Read the compressed file
    compressed_data = read(filepath)
    
    # Decompress directly to IOBuffer
    io = IOBuffer()
    write(io, transcode(GzipDecompressor, compressed_data))
    seekstart(io)  # Reset position to beginning
    
    return io
end

# ╔═╡ c4bc9abf-d0c8-4e5c-9c76-96b36c1f68f5
# Uncompress
io = decompress_gz_to_iobuffer(fpath)

# ╔═╡ 216027cd-7737-4056-9abc-d71c899f7568
md"""
!!! todo
	Move decompression bit to FITSFiles.jl? Does something like this already exist in Base somewhere?
"""

# ╔═╡ 03d97df3-1498-4c32-b33f-d7e8e7f26751
# Load
hdus = fits(io)

# ╔═╡ 904a9095-e2a2-4c3e-8070-719033a5c6b3
md"""
Taking a look, we see that we have the following HDUs:
"""

# ╔═╡ 466af93a-14ac-4461-ab1e-855209ce0b43
info(hdus)

# ╔═╡ 207dca91-36b1-48fa-b917-918e6a1edebd
md"""
Let's load in the table from the second HDU:
"""

# ╔═╡ 5759aefa-2e5f-4f96-afa7-897f5171082f
t = hdus[2]

# ╔═╡ 3e55d0af-a930-40da-b5ea-b822ce01140d
t.data

# ╔═╡ f97b20a0-235b-4efe-aad0-75733290dc0a
md"""
and store the wavelength and flux information into `wave_spectrum` and `flux_spectrum`, respectively:
"""

# ╔═╡ ccb3ce4f-731e-433b-9479-773aa97ac0ce
wav_spectrum, flux_spectrum = let
	flux_spectrum_raw = vec(t.data.FLUX)

	wav_spectrum_raw = range(;
		start = first(t.data.WAVELENGTH),
		step = first(t.data.DELTAW),
		length = length(flux_spectrum_raw),
	)

	# Remove spurious negative fluxes
	idxs_bad = findall(<(0), flux_spectrum_raw)

	(
		wav_spectrum_raw[Not(idxs_bad)]u"Å", 
		flux_spectrum_raw[Not(idxs_bad)]u"erg/s/cm^2/Å",
	)
end

# ╔═╡ e92c118e-7547-4fbe-8a6b-9487beb99ff9
md"""
We turn next to loading in some corresponding photometry.
"""

# ╔═╡ 37f294a1-448d-478b-a593-bfeabb1c84f6
md"""
### SIMBAD Photometry

See here for more <https://simbad.cds.unistra.fr/simbad/sim-tap/>
"""

# ╔═╡ a439386d-2b9b-44e8-9c41-cb2ec787024e
phot_simbad = execute(
	 TAPService("https://simbad.u-strasbg.fr/simbad/sim-tap/"),
	"""
	select U, B, V from allfluxes
	join ident using(oidref)
	where id = 'HD 147933'
	"""
)

# ╔═╡ 39eb85ea-9eac-4c5b-91a8-9252418dfc33
Umag, Bmag, Vmag = first(phot_simbad)

# ╔═╡ 1fac575d-9ff7-4fda-bae9-7dba874396fc
md"""
We need to convert these magnitudes to flux so that they can be plotted on the same scale as our spectrum.
"""

# ╔═╡ 844446f9-88e1-468f-bcf2-ad638aaeb844
# Central wavelengths, stored as an array for plotting later
wav_phot = wav_U, wav_B, wav_V = [0.3660, 0.4400, 0.5530] .* u"μm"

# ╔═╡ 101318a4-4bd2-4253-9e88-a814cbb6578b
# Zero-points flux in frequency space
flux_U0_nu, flux_B0_nu, flux_V0_nu = (1.81e-23, 4.26e-23, 3.64e-23) .* u"W/m^2/Hz"

# ╔═╡ df87135d-59cb-4a2f-be9f-f76c4203c645
md"""
!!! note
	More photometric definitions here: <https://ned.ipac.caltech.edu/help/photoband.lst>
"""

# ╔═╡ 36f69f11-10d8-4d42-b587-2b885f497973
to_Flam(Fnu, wav_cen) = Fnu * (c0 / wav_cen^2)

# ╔═╡ 4ee58446-1563-44c8-b31f-6f9c3722f707
# Zero-point flux in wavelength space
flux_U0, flux_B0, flux_V0 = (
	to_Flam(flux_U0_nu, wav_U),
	to_Flam(flux_B0_nu, wav_B),
	to_Flam(flux_V0_nu, wav_V),
);

# ╔═╡ 855cfb7b-5e66-415c-9d48-f66b38e835ff
to_flux(flux_0, mag) = flux_0 * exp10(-0.4 * mag)

# ╔═╡ 1c4d9d17-ab51-4ecc-9caf-c068ffc65fd2
flux_phot = [
	to_flux(flux_U0, Umag),
	to_flux(flux_B0, Bmag),
	to_flux(flux_V0, Vmag),
]

# ╔═╡ 637830ec-1b5a-4eaa-bb76-55acd3b08985
md"""
!!! note
	See [PhotometricFilters.jl](https://juliaastro.org/PhotometricFilters) for a comprehensive treatment of photometric filter curves.
"""

# ╔═╡ e59fa830-195f-45fe-bdc0-b4a2b23cd8c4
md"""
### Initial spectrum + photometry
"""

# ╔═╡ 84386f76-3972-4985-9dc4-05501bbc5f41
let
	# Spectrum
	fig, ax, p = lines(wav_spectrum, flux_spectrum;
		axis = (;
			dim1_conversion = DQConversion(us"Å"),
			dim2_conversion = DQConversion(us"erg/s/cm^2/Å"),
			xlabel = "Wavelength",
			ylabel = "Flux",
			title = "ρ Oph",
		),
		color = Cycled(2),
		label = "IEU Spectrum",
	)

	# Photometry
	scatter!(ax, wav_phot, flux_phot; color = Cycled(6), markersize = 15, label = "U, B, V")

	# Legend and axis limits
	axislegend()
	ylims!(ax, (0u"erg/s/cm^2/Å", (3.1e-10)u"erg/s/cm^2/Å"))
	fig
end

# ╔═╡ dba6f2b7-890a-4e7d-a86c-02e99f985796
md"""
### Dereddened spectrum + photometry
"""

# ╔═╡ b1207731-73be-4390-a329-30f8c1e7bc94
flux_spectrum_deredden = deredden.(
	F99(; Rv = 5), ustrip.(u"Å", wav_spectrum), flux_spectrum;
	Av = 5 * 0.5 # Aᵥ = Rᵥ * E(B - V)
) .|> us"erg/s/cm^2/Å"

# ╔═╡ 7f53c04c-5e51-4fb4-bd7c-64bc2b9d2944
flux_phot_deredden = deredden.(
	F99(; Rv = 5), ustrip.(u"Å", wav_phot), flux_phot;
	Av = 5 * 0.5 # Aᵥ = Rᵥ * E(B - V)
) .|> us"erg/s/cm^2/Å"

# ╔═╡ ff040052-e0b3-4447-811e-4f743456aed3
let
	# Spectrum
	fig, ax, p = lines(wav_spectrum, flux_spectrum_deredden;
		axis = (;
			dim1_conversion = DQConversion(us"Å"),
			dim2_conversion = DQConversion(us"erg/s/cm^2/Å"),
			yscale = log10,
			yminorticksvisible = true,
			yminorticks = IntervalsBetween(9),
			xlabel = "Wavelength",
			ylabel = "Flux",
			title = "ρ Oph",
		),
		color = Cycled(5),
		label = "IUE spectrum dereddened",
	)
	lines!(wav_spectrum, flux_spectrum;
		color = Cycled(2),
		label = "IUE spectrum",
	)

	# Photometry
	scatter!(ax, wav_phot, flux_phot_deredden;
		color = Cycled(1),
		markersize = 15,
		label = "U, B, V dereddened",
	)
	scatter!(ax, wav_phot, flux_phot;
		color = Cycled(6),
		markersize = 15,
		label = "U, B, V",
	)

	# Legend and axis limits
	axislegend()
	ylims!(ax, (1e-11u"erg/s/cm^2/Å", 5e-8u"erg/s/cm^2/Å"))
	fig
end

# ╔═╡ 8a2e5f54-ec9d-46c5-88c6-aa647e7ea036
md"""
# Notebook setup 🔧
"""

# ╔═╡ 8e118239-43e0-427c-bdf1-80c9d5fe130c
TableOfContents(; depth = 4)

# ╔═╡ Cell order:
# ╟─f18df5a0-1dda-4371-aea6-8ecbce67908c
# ╟─2c3cabef-29fd-41cf-b3df-94f58748c22b
# ╠═b6b27fe2-c7f7-11f0-8b42-052bc2026e99
# ╠═89366506-4854-4859-b109-2d6cd9456bb5
# ╠═555cf1bd-5b52-4d26-aaa2-033a9f4d831a
# ╟─3680c62b-7ce8-4fbc-8cd0-249d499e990f
# ╟─f3bf7784-b2da-4ed6-9cb5-0799d81d65ae
# ╠═e7266de8-ba60-45b3-8988-6143e7098b72
# ╟─ed46db90-9802-4615-928f-bcad4bd19804
# ╠═0250218a-6670-483f-8122-045ca65d28e4
# ╠═eb1c1b16-c391-45a5-81c2-5fd3dbb12c50
# ╟─943661a3-3b40-4c01-a624-79726e56b385
# ╟─01195fc5-9ab9-4b9d-a4f9-a23f41bc6f2a
# ╠═d57bbe7f-b21a-479a-8249-a96d6aac3e95
# ╟─70cefdd5-4afb-4e09-b2df-75983bb40e85
# ╠═285b4bc6-9689-4995-8e8e-dc6acfbd7f81
# ╟─207bfd62-2332-413d-8034-15d9927ea57e
# ╠═dfb6b2d7-efcb-4505-9f69-d51eafb7db9f
# ╠═c4bc9abf-d0c8-4e5c-9c76-96b36c1f68f5
# ╟─85b04476-2ce0-49dd-9bef-8f62c819e021
# ╟─216027cd-7737-4056-9abc-d71c899f7568
# ╠═03d97df3-1498-4c32-b33f-d7e8e7f26751
# ╟─904a9095-e2a2-4c3e-8070-719033a5c6b3
# ╠═466af93a-14ac-4461-ab1e-855209ce0b43
# ╟─207dca91-36b1-48fa-b917-918e6a1edebd
# ╠═5759aefa-2e5f-4f96-afa7-897f5171082f
# ╠═3e55d0af-a930-40da-b5ea-b822ce01140d
# ╟─f97b20a0-235b-4efe-aad0-75733290dc0a
# ╠═ccb3ce4f-731e-433b-9479-773aa97ac0ce
# ╟─e92c118e-7547-4fbe-8a6b-9487beb99ff9
# ╟─37f294a1-448d-478b-a593-bfeabb1c84f6
# ╠═a439386d-2b9b-44e8-9c41-cb2ec787024e
# ╠═39eb85ea-9eac-4c5b-91a8-9252418dfc33
# ╟─1fac575d-9ff7-4fda-bae9-7dba874396fc
# ╠═844446f9-88e1-468f-bcf2-ad638aaeb844
# ╠═101318a4-4bd2-4253-9e88-a814cbb6578b
# ╟─df87135d-59cb-4a2f-be9f-f76c4203c645
# ╠═4ee58446-1563-44c8-b31f-6f9c3722f707
# ╠═1c4d9d17-ab51-4ecc-9caf-c068ffc65fd2
# ╠═36f69f11-10d8-4d42-b587-2b885f497973
# ╟─855cfb7b-5e66-415c-9d48-f66b38e835ff
# ╟─637830ec-1b5a-4eaa-bb76-55acd3b08985
# ╟─e59fa830-195f-45fe-bdc0-b4a2b23cd8c4
# ╠═84386f76-3972-4985-9dc4-05501bbc5f41
# ╟─dba6f2b7-890a-4e7d-a86c-02e99f985796
# ╠═b1207731-73be-4390-a329-30f8c1e7bc94
# ╠═7f53c04c-5e51-4fb4-bd7c-64bc2b9d2944
# ╠═ff040052-e0b3-4447-811e-4f743456aed3
# ╟─8a2e5f54-ec9d-46c5-88c6-aa647e7ea036
# ╠═95d05163-b6d2-4a01-a782-6da05479ee0a
# ╠═8e118239-43e0-427c-bdf1-80c9d5fe130c
