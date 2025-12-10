sections = ["FITS files", "units", "models", "cosmology", "dust"]

Dict(
    "main" => [uppercase(section) => collections[section].pages for section in sections],
)
