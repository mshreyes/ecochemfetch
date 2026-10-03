from chemfetch import fetch_chemicals_to_csv

pollutants = [
    "DDT",
    "Atrazine",
    "Bisphenol A",
    "Perfluorooctanoic acid",
    "Trichloroethylene"
]

fetch_chemicals_to_csv(
    names=pollutants,
    output_filepath="example_output.csv",
    delay_seconds=0.3
)