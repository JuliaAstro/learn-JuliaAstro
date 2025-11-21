sections = ["FITS", "units"]

Dict(
     "main" => [uppercase(section) => collections[section].pages for section in sections],
)
