sections = ["FITS files", "units", "cosmology"]

Dict(
    "main" => [uppercase(section) => collections[section].pages for section in sections],
)
