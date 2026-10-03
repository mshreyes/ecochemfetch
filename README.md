# ecochemfetch

An open-source Python tool to fetch comprehensive chemical identifiers, structural properties, CACTVS fingerprints, toxicology profiles, environmental fate metrics, regulatory codes, and mass spectra records directly into CSV files.

---

## Features

* **Zero API Keys**: Uses public scientific REST APIs (PubChem, Cactus, Wikidata, ChEBI, MassBank Europe).


* **Exhaustive Synonyms**: Fetches complete synonym and alias lists from PubChem without artificial caps.


* **881-bit Fingerprint Decoding**: Converts PubChem's Base64 `Fingerprint2D` directly into raw binary CACTVS bitstrings (`01001...`).


* **Resilient Network Handling**: Configured with exponential backoff retries for rate-limiting (`429`) and server drops (`5xx`).


* **Safe CSV Streaming**: Appends records row-by-row on disk, preventing data loss during batch queries.



---

## Data Points Captured

* **Identifiers**: CID, IUPAC Name, SMILES (Canonical & Isomeric), InChI, InChIKey, Formula, Molecular Weight, Monoisotopic Mass, Exact Mass.


* **Registries & Regulatory**: CASRN, EU EC/EINECS Number, US EPA CompTox ID (`DTXSID`), UN Dangerous Goods Transport Number, RTECS Code, ChEBI ID.


* **Classification & Roles**: ChEBI biological and chemical application roles (e.g., *herbicide, persistent organic pollutant, endocrine disruptor*).


* **GHS & Environmental Hazards**: GHS Signal word (`DANGER`, `WARNING`), environmental hazard statements (H400-series aquatic toxicity), complete GHS statement set.


* **Toxicology & Safety**: Acute toxicity doses and routes ($LD_{50}$, $LC_{50}$), IARC Carcinogenicity classifications.


* **Exposure & Ecotoxicity**: OSHA PEL, NIOSH REL, ACGIH TLV air limits, aquatic ecotoxicity benchmarks ($LC_{50}$/$EC_{50}$ for indicator species).


* **Environmental Fate & Transport**: Environmental fate summaries, biodegradation kinetics, Bioconcentration Factor (BCF), Soil Adsorption Coefficient ($K_{oc}$).


* **Mass Spectrometry**: MassBank Europe reference spectra record count, accession IDs, and instrument/mode metadata.


* **Structural Descriptors**: $XLogP$, TPSA, Complexity, Charge, H-Bond Donors/Acceptors, Rotatable Bonds, Heavy Atoms, Stereocenter counts, and 881-bit CACTVS fingerprints.



---

## Installation

Inside your project directory where `pyproject.toml` is located:

### Using `uv`

```bash
uv pip install -e .

```

### Using standard `pip`

```bash
pip install -e .

```

---

## Project Structure

```text
ecochemfetch/
├── pyproject.toml
├── README.md
└── src/
    └── ecochemfetch/
        ├── __init__.py
        └── client.py

```

---

## Usage

### 1. Batch Export to CSV

```python
from ecochemfetch import fetch_chemicals_to_csv

chemicals = [
    "Atrazine",
    "Perfluorooctanoic acid",
    "DDT",
    "Bisphenol A",
    "Chlorpyrifos",
    "Trichloroethylene"
]

fetch_chemicals_to_csv(
    names=chemicals,
    output_filepath="example_output.csv",
    delay_seconds=0.35  # Politeness interval to stay within PubChem rate limits
)

```

### 2. Single Chemical Lookup

```python
from ecochemfetch import fetch_chemical_info

info = fetch_chemical_info("Atrazine")

print("CID               :", info["cid"])
print("CAS RN            :", info["cas_rn"])
print("EPA DTXSID        :", info["epa_dtxsid"])
print("EC Number         :", info["ec_number"])
print("UN Transport      :", info["un_transport_number"])
print("GHS Signal        :", info["ghs_signal"])
print("Acute Toxicity    :", info["acute_toxicity_ld50"])
print("Soil Koc          :", info["soil_adsorption_koc"])
print("MassBank Hits     :", info["massbank_record_count"])
print("ChEBI Roles       :", info["chebi_roles"])
print("Total Synonyms    :", info["total_synonyms_count"])

```

---

| Field Name | Type | Description |
| :--- | :--- | :--- |
| `query` | `string` | Query chemical name submitted |
| `status` | `string` | `Found`, `Not Found`, or failure description |
| `cid` | `integer` | PubChem Compound ID |
| `epa_dtxsid` | `string` | US EPA CompTox Chemicals Dashboard ID |
| `cas_rn` | `string` | CAS Registry Numbers resolved from CIR/PubChem |
| `ec_number` | `string` | European Community (EINECS) number via Wikidata |
| `un_transport_number` | `string` | UN dangerous goods transport code (e.g., `UN 2761`) |
| `rtecs_id` | `string` | NIOSH RTECS toxicity registry identifier |
| `chebi_id` | `string` | EMBL-EBI ChEBI identifier |
| `chebi_roles` | `string` | Biological, environmental, and application roles |
| `iupac_name` | `string` | Systematic chemical nomenclature |
| `molecular_formula` | `string` | Molecular formula |
| `molecular_weight` | `float` | Standard molecular mass |
| `canonical_smiles` | `string` | 1D Canonical SMILES representation |
| `isomeric_smiles` | `string` | 1D SMILES including stereochemical definitions |
| `inchikey` | `string` | 27-character fixed InChIKey hash |
| `xlogp` | `float` | Octanol-water partition coefficient |
| `tpsa` | `float` | Topological polar surface area (Å²) |
| `ghs_signal` | `string` | Signal words (`DANGER`, `WARNING`) |
| `environmental_hazards` | `string` | Extracted H400-series aquatic/environmental warnings |
| `ghs_hazard_statements` | `string` | Complete list of GHS hazard statements |
| `acute_toxicity_ld50` | `string` | Acute LD50 / LC50 doses, routes, and test species |
| `carcinogenicity_iarc` | `string` | IARC carcinogenicity evaluation summary |
| `exposure_limits_air` | `string` | OSHA PEL, NIOSH REL, and ACGIH TLV exposure thresholds |
| `aquatic_ecotoxicity_benchmarks` | `string` | LC50 / EC50 ecotox benchmarks for aquatic species |
| `env_fate_summary` | `string` | Environmental partitioning and degradation summary |
| `biodegradation` | `string` | Aerobic/anaerobic microbial persistence and half-lives |
| `bioconcentration_bcf` | `string` | Aquatic bioconcentration factors (BCF) |
| `soil_adsorption_koc` | `string` | Soil organic carbon-water partition coefficient (Koc) |
| `massbank_record_count` | `integer` | Count of matching MS/MS spectra in MassBank Europe |
| `massbank_accessions` | `string` | MassBank record accession identifiers |
| `massbank_instruments` | `string` | Mass spectrometry instrument type and ionization mode |
| `cactvs_fp_bits_881` | `string` | Decoded 881-bit binary CACTVS structural fingerprint |
| `total_synonyms_count` | `integer` | Total recorded synonym count in PubChem |
| `all_synonyms` | `string` | Semicolon-delimited list of all synonyms and trade names |


---

## License

MIT

## Disclaimer

Google Gemini 3.8 Flash was used to tidy up the code and generate the README.md file.
