sections = ["FITS files", "units", "cosmology", "dust"]

Dict(
    "main" => [uppercase(section) => collections[section].pages for section in sections],
)
