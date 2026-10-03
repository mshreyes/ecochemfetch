import base64
import csv
import logging
import os
import re
import sys
import time
from typing import Any, Dict, List, Optional
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry

# Allow writing large CSV cells (needed for extensive synonym blocks and MS annotations)
csv.field_size_limit(sys.maxsize)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger("ecochemfetch")

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound"
CACTUS_BASE = "https://cactus.nci.nih.gov/chemical/structure"
WIKIDATA_SPARQL_URL = "https://query.wikidata.org/sparql"
CHEBI_REST_BASE = "https://www.ebi.ac.uk/webservices/chebi/2.0/test/getCompleteEntity"
MASSBANK_EU_API = "https://massbank.eu/MassBank-api/records"

DEFAULT_TIMEOUT = 15

CAS_REGEX = re.compile(r"^\b[1-9]\d{1,6}-\d{2}-\d\b$")
DTXSID_REGEX = re.compile(r"^DTXSID\d+$")

ALL_PROPERTIES = [
    "MolecularFormula", "MolecularWeight", "CanonicalSMILES", "IsomericSMILES",
    "InChI", "InChIKey", "IUPACName", "Title", "XLogP", "ExactMass",
    "MonoisotopicMass", "TPSA", "Complexity", "Charge", "HBondDonorCount",
    "HBondAcceptorCount", "RotatableBondCount", "HeavyAtomCount",
    "IsotopeAtomCount", "AtomStereoCount", "DefinedAtomStereoCount",
    "UndefinedAtomStereoCount", "BondStereoCount", "DefinedBondStereoCount",
    "UndefinedBondStereoCount", "CovalentUnitCount", "Fingerprint2D",
]


def _create_resilient_session(retries: int = 3, backoff_factor: float = 0.5) -> requests.Session:
    session = requests.Session()
    retry_strategy = Retry(
        total=retries,
        backoff_factor=backoff_factor,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({"User-Agent": "ecochemfetch/0.6.0 (Open Chemistry Pipeline; mailto:info@example.com)"})
    return session


SESSION = _create_resilient_session()


def _safe_get_json(url: str, params: Optional[Dict[str, Any]] = None, headers: Optional[Dict[str, str]] = None) -> Optional[Any]:
    try:
        res = SESSION.get(url, params=params, headers=headers, timeout=DEFAULT_TIMEOUT)
        if res.status_code == 200:
            return res.json()
    except Exception as e:
        logger.debug(f"Request failed for {url}: {e}")
    return None


def _decode_pubchem_fingerprint(b64_str: Optional[str]) -> str:
    if not b64_str:
        return ""
    try:
        raw_bytes = base64.b64decode(b64_str)
        bits = "".join(f"{b:08b}" for b in raw_bytes)
        return bits[32:32 + 881]
    except Exception:
        return ""


def _extract_strings_from_node(node: Any) -> List[str]:
    found: List[str] = []
    if isinstance(node, dict):
        if "StringWithMarkup" in node:
            for item in node["StringWithMarkup"]:
                text = item.get("String", "").strip()
                if text:
                    found.append(text)
        elif "String" in node:
            found.append(node["String"].strip())
        for v in node.values():
            found.extend(_extract_strings_from_node(v))
    elif isinstance(node, list):
        for item in node:
            found.extend(_extract_strings_from_node(item))
    return found


def _fetch_pug_view_section(cid: int, heading: str) -> Optional[Dict[str, Any]]:
    url = f"{PUBCHEM_VIEW_BASE}/{cid}/JSON?heading={requests.utils.quote(heading)}"
    return _safe_get_json(url)


def _get_exposure_and_ecotox(cid: int) -> Dict[str, str]:
    """
    Extracts OSHA/NIOSH occupational exposure limits and aquatic ecotoxicity benchmarks.
    """
    exposure_data = _fetch_pug_view_section(cid, "Exposure Limits")
    ecotox_data = _fetch_pug_view_section(cid, "Ecotoxicity Values")

    exposure_limits = []
    ecotox_values = []

    if exposure_data:
        for line in _extract_strings_from_node(exposure_data):
            if any(term in line for term in ["PEL", "REL", "TLV", "TWA", "mg/m3", "ppm"]):
                exposure_limits.append(line)

    if ecotox_data:
        for line in _extract_strings_from_node(ecotox_data):
            if any(term in line for term in ["LC50", "EC50", "NOEC", "Daphnia", "trout"]):
                ecotox_values.append(line)

    clean_exposure = list(dict.fromkeys(exposure_limits))[:3]
    clean_ecotox = list(dict.fromkeys(ecotox_values))[:3]

    return {
        "exposure_limits_air": " | ".join(clean_exposure),
        "aquatic_ecotoxicity_benchmarks": " | ".join(clean_ecotox),
    }


def _get_toxicity_dossier(cid: int) -> Dict[str, str]:
    data = _fetch_pug_view_section(cid, "Toxicity")
    if not data:
        return {"acute_toxicity_ld50": "", "carcinogenicity_iarc": ""}

    ld50_entries = []
    iarc_entries = []

    def scan_toxicity(node: Any):
        if isinstance(node, dict):
            heading = node.get("TOCHeading", "")
            if heading == "Non-Human Toxicity Values" or "Acute" in heading:
                for line in _extract_strings_from_node(node):
                    if re.search(r"\b(LD50|LC50)\b", line, re.IGNORECASE):
                        ld50_entries.append(line)
            if "IARC" in heading or "Carcinogen" in heading:
                for line in _extract_strings_from_node(node):
                    if "Group" in line or "Carcinogenic" in line or "IARC" in line:
                        iarc_entries.append(line)
            for v in node.values():
                scan_toxicity(v)
        elif isinstance(node, list):
            for item in node:
                scan_toxicity(item)

    scan_toxicity(data)
    clean_ld50 = list(dict.fromkeys(ld50_entries))[:6]
    clean_iarc = list(dict.fromkeys(iarc_entries))[:3]

    return {
        "acute_toxicity_ld50": " | ".join(clean_ld50),
        "carcinogenicity_iarc": " | ".join(clean_iarc),
    }


def _get_environmental_fate_dossier(cid: int) -> Dict[str, str]:
    data = _fetch_pug_view_section(cid, "Environmental Fate & Transport")
    if not data:
        return {
            "env_fate_summary": "",
            "biodegradation": "",
            "bioconcentration_bcf": "",
            "soil_adsorption_koc": "",
        }

    fate_summary: List[str] = []
    biodegradation: List[str] = []
    bcf_values: List[str] = []
    koc_values: List[str] = []

    def scan_fate(node: Any):
        if isinstance(node, dict):
            heading = node.get("TOCHeading", "")
            if heading == "Environmental Fate/Exposure Summary":
                fate_summary.extend(_extract_strings_from_node(node))
            elif "Biodegradation" in heading:
                biodegradation.extend(_extract_strings_from_node(node))
            elif "Bioconcentration" in heading or heading == "Bioaccumulation":
                bcf_values.extend(_extract_strings_from_node(node))
            elif "Soil Adsorption" in heading or "Mobility in Soil" in heading:
                koc_values.extend(_extract_strings_from_node(node))
            for v in node.values():
                scan_fate(v)
        elif isinstance(node, list):
            for item in node:
                scan_fate(item)

    scan_fate(data)

    def clean_and_slice(items: List[str], max_items: int = 3) -> str:
        filtered = [it for it in items if len(it) > 10 and not it.startswith("Download")]
        return " | ".join(list(dict.fromkeys(filtered))[:max_items])

    return {
        "env_fate_summary": clean_and_slice(fate_summary, 1),
        "biodegradation": clean_and_slice(biodegradation, 2),
        "bioconcentration_bcf": clean_and_slice(bcf_values, 2),
        "soil_adsorption_koc": clean_and_slice(koc_values, 2),
    }


def _get_environmental_ghs(cid: int) -> Dict[str, str]:
    data = _fetch_pug_view_section(cid, "GHS Classification")
    if not data:
        return {"ghs_signal": "", "ghs_hazard_statements": "", "environmental_hazards": ""}

    signals = set()
    hazards = set()
    env_statements = set()

    def search_ghs(node: Any):
        if isinstance(node, dict):
            if "Information" in node:
                for info in node.get("Information", []):
                    name = info.get("Name", "")
                    val = info.get("Value", {})
                    if "Signal" in name:
                        for s in val.get("StringWithMarkup", []):
                            signals.add(s.get("String", "").strip())
                    if "Hazard Statements" in name or "GHS Hazard Statements" in name:
                        for s in val.get("StringWithMarkup", []):
                            text = s.get("String", "").strip()
                            if text:
                                hazards.add(text)
                                if re.search(r"\bH4\d{2}\b|aquatic|environment|ozone", text, re.IGNORECASE):
                                    env_statements.add(text)
            for v in node.values():
                search_ghs(v)
        elif isinstance(node, list):
            for item in node:
                search_ghs(item)

    search_ghs(data)
    return {
        "ghs_signal": ", ".join(sorted(signals)),
        "ghs_hazard_statements": " | ".join(sorted(hazards)),
        "environmental_hazards": " | ".join(sorted(env_statements)),
    }


def _get_cas_from_cactus(name: str) -> List[str]:
    url = f"{CACTUS_BASE}/{requests.utils.quote(name)}/cas"
    try:
        res = SESSION.get(url, timeout=DEFAULT_TIMEOUT)
        if res.status_code == 200 and res.text:
            lines = [line.strip() for line in res.text.splitlines() if line.strip()]
            return [c for c in lines if CAS_REGEX.match(c)]
    except Exception:
        pass
    return []


def _get_wikidata_regulatory(inchikey: str) -> Dict[str, str]:
    if not inchikey:
        return {"ec_number": "", "un_number": "", "rtecs_id": "", "chebi_id": ""}

    sparql_query = f"""
    SELECT ?ec ?un ?rtecs ?chebi WHERE {{
      ?item wdt:P235 "{inchikey}" .
      OPTIONAL {{ ?item wdt:P232 ?ec. }}
      OPTIONAL {{ ?item wdt:P695 ?un. }}
      OPTIONAL {{ ?item wdt:P657 ?rtecs. }}
      OPTIONAL {{ ?item wdt:P683 ?chebi. }}
    }} LIMIT 5
    """
    headers = {"Accept": "application/sparql-results+json"}
    data = _safe_get_json(WIKIDATA_SPARQL_URL, params={"query": sparql_query}, headers=headers)
    if not data:
        return {"ec_number": "", "un_number": "", "rtecs_id": "", "chebi_id": ""}

    bindings = data.get("results", {}).get("bindings", [])
    ecs, uns, rtecs_list, chebis = set(), set(), set(), set()
    for row in bindings:
        if "ec" in row:
            ecs.add(row["ec"]["value"])
        if "un" in row:
            uns.add(f"UN {row['un']['value']}")
        if "rtecs" in row:
            rtecs_list.add(row["rtecs"]["value"])
        if "chebi" in row:
            chebis.add(row["chebi"]["value"])

    return {
        "ec_number": "; ".join(sorted(ecs)),
        "un_number": "; ".join(sorted(uns)),
        "rtecs_id": "; ".join(sorted(rtecs_list)),
        "chebi_id": sorted(chebis)[0] if chebis else "",
    }


def _get_chebi_roles(chebi_id: str) -> str:
    if not chebi_id:
        return ""
    url = f"{CHEBI_REST_BASE}?chebiId={chebi_id}"
    try:
        res = SESSION.get(url, timeout=DEFAULT_TIMEOUT)
        if res.status_code == 200 and res.text:
            matches = re.findall(r"(.*?)", res.text)
            roles = [m.strip() for m in matches if len(m) > 3 and not m.startswith("CHEBI")]
            return " | ".join(list(dict.fromkeys(roles))[:5])
    except Exception:
        pass
    return ""


def _get_massbank_spectra(inchikey: str) -> Dict[str, Any]:
    """
    Queries MassBank Europe REST API to check for MS/MS reference spectra records.
    Returns count, top accessions, and analytical instrument metadata.
    """
    if not inchikey:
        return {"massbank_record_count": 0, "massbank_accessions": "", "massbank_instruments": ""}

    data = _safe_get_json(MASSBANK_EU_API, params={"inchi_key": inchikey})
    if not data or not isinstance(data, list):
        return {"massbank_record_count": 0, "massbank_accessions": "", "massbank_instruments": ""}

    accessions = []
    instruments = []

    for item in data:
        if isinstance(item, dict):
            acc = item.get("accession") or item.get("id")
            if acc:
                accessions.append(acc)
            # Collect instrument & mode details if provided
            inst = item.get("instrument_type") or item.get("ms_type")
            ion_mode = item.get("ion_mode")
            meta = " ".join(filter(None, [inst, ion_mode]))
            if meta:
                instruments.append(meta)

    return {
        "massbank_record_count": len(data),
        "massbank_accessions": "; ".join(accessions[:5]),
        "massbank_instruments": " | ".join(list(dict.fromkeys(instruments))[:3]),
    }


def fetch_chemical_info(name: str) -> Dict[str, Any]:
    """
    Fetches comprehensive physical, structural, toxicological,
    regulatory, exposure limits, and mass spectral profiles for a chemical.
    """
    # 1. Resolve CID
    cid_url = f"{PUBCHEM_BASE}/compound/name/{requests.utils.quote(name)}/cids/JSON"
    cid_data = _safe_get_json(cid_url)
    if not cid_data:
        return {"query": name, "status": "Not Found"}
    cids = cid_data.get("IdentifierList", {}).get("CID", [])
    if not cids:
        return {"query": name, "status": "Not Found"}
    cid = cids[0]

    # 2. General 2D Properties
    prop_list = ",".join(ALL_PROPERTIES)
    props_url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/{prop_list}/JSON"
    props_data = _safe_get_json(props_url)
    props = props_data.get("PropertyTable", {}).get("Properties", [])[0] if props_data else {}

    # 3. Exhaustive Synonyms & Identifiers
    syn_url = f"{PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON"
    syn_data = _safe_get_json(syn_url)
    all_synonyms: List[str] = []
    if syn_data:
        info = syn_data.get("InformationList", {}).get("Information", [])
        if info:
            all_synonyms = info[0].get("Synonym", [])

    cas_numbers = _get_cas_from_cactus(name)
    if not cas_numbers:
        cas_numbers = list(dict.fromkeys(s for s in all_synonyms if CAS_REGEX.match(s)))

    dtxsid_list = list(dict.fromkeys(s for s in all_synonyms if DTXSID_REGEX.match(s)))

    # 4. Regulatory, Toxicity & Environmental Dossiers
    env_ghs = _get_environmental_ghs(cid)
    tox_data = _get_toxicity_dossier(cid)
    fate_data = _get_environmental_fate_dossier(cid)
    exposure_ecotox = _get_exposure_and_ecotox(cid)

    # 5. Wikidata & ChEBI Roles
    inchikey = props.get("InChIKey", "")
    wiki_data = _get_wikidata_regulatory(inchikey)
    chebi_roles = _get_chebi_roles(wiki_data.get("chebi_id", ""))

    # 6. MassBank Europe MS/MS Spectra Lookup
    mb_data = _get_massbank_spectra(inchikey)

    raw_fp = props.get("Fingerprint2D", "")

    return {
        "query": name,
        "status": "Found",
        "cid": cid,
        # Environmental & Regulatory Identifiers
        "epa_dtxsid": "; ".join(dtxsid_list),
        "cas_rn": "; ".join(cas_numbers),
        "ec_number": wiki_data["ec_number"],
        "un_transport_number": wiki_data["un_number"],
        "rtecs_id": wiki_data["rtecs_id"],
        "chebi_id": wiki_data["chebi_id"],
        "chebi_roles": chebi_roles,
        # Structure & Basic Properties
        "iupac_name": props.get("IUPACName", ""),
        "molecular_formula": props.get("MolecularFormula", ""),
        "molecular_weight": props.get("MolecularWeight", ""),
        "canonical_smiles": props.get("CanonicalSMILES", ""),
        "isomeric_smiles": props.get("IsomericSMILES", ""),
        "inchikey": inchikey,
        "xlogp": props.get("XLogP", ""),
        "tpsa": props.get("TPSA", ""),
        # Environmental Hazard & GHS
        "ghs_signal": env_ghs["ghs_signal"],
        "environmental_hazards": env_ghs["environmental_hazards"],
        "ghs_hazard_statements": env_ghs["ghs_hazard_statements"],
        # Toxicity & Exposure
        "acute_toxicity_ld50": tox_data["acute_toxicity_ld50"],
        "carcinogenicity_iarc": tox_data["carcinogenicity_iarc"],
        "exposure_limits_air": exposure_ecotox["exposure_limits_air"],
        "aquatic_ecotoxicity_benchmarks": exposure_ecotox["aquatic_ecotoxicity_benchmarks"],
        # Environmental Fate & Transport
        "env_fate_summary": fate_data["env_fate_summary"],
        "biodegradation": fate_data["biodegradation"],
        "bioconcentration_bcf": fate_data["bioconcentration_bcf"],
        "soil_adsorption_koc": fate_data["soil_adsorption_koc"],
        # Mass Spectrometry (MassBank EU)
        "massbank_record_count": mb_data["massbank_record_count"],
        "massbank_accessions": mb_data["massbank_accessions"],
        "massbank_instruments": mb_data["massbank_instruments"],
        # Fingerprint & Exhaustive Synonyms
        "cactvs_fp_bits_881": _decode_pubchem_fingerprint(raw_fp),
        "total_synonyms_count": len(all_synonyms),
        "all_synonyms": "; ".join(all_synonyms),
    }


def fetch_chemicals_to_csv(
    names: List[str],
    output_filepath: str = "environmental_toxicology.csv",
    delay_seconds: float = 0.35,
) -> None:
    """
    Streams complete chemical, toxicological, spectral, and fate records to a CSV.
    """
    all_headers = [
        "query", "status", "cid", "epa_dtxsid", "cas_rn", "ec_number",
        "un_transport_number", "rtecs_id", "chebi_id", "chebi_roles",
        "iupac_name", "molecular_formula", "molecular_weight", "canonical_smiles",
        "isomeric_smiles", "inchikey", "xlogp", "tpsa", "ghs_signal",
        "environmental_hazards", "ghs_hazard_statements", "acute_toxicity_ld50",
        "carcinogenicity_iarc", "exposure_limits_air", "aquatic_ecotoxicity_benchmarks",
        "env_fate_summary", "biodegradation", "bioconcentration_bcf",
        "soil_adsorption_koc", "massbank_record_count", "massbank_accessions",
        "massbank_instruments", "cactvs_fp_bits_881", "total_synonyms_count", "all_synonyms"
    ]

    write_header = not os.path.exists(output_filepath)

    with open(output_filepath, mode="a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=all_headers, restval="")
        if write_header:
            writer.writeheader()

        total = len(names)
        for idx, name in enumerate(names, start=1):
            clean_name = name.strip()
            if not clean_name:
                continue

            print(f"[{idx}/{total}] Fetching environmental dossier: {clean_name}...")
            try:
                record = fetch_chemical_info(clean_name)
            except Exception as e:
                logger.error(f"Error fetching '{clean_name}': {e}")
                record = {"query": clean_name, "status": f"Error: {e}"}

            writer.writerow(record)
            f.flush()

            if delay_seconds > 0:
                time.sleep(delay_seconds)

    print(f"\nFinished. Comprehensive records saved to '{output_filepath}'.")