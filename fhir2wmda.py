#!/usr/bin/env python3
"""
fhir_to_wmda_post.py

Convert a FHIR Bundle JSON into a WMDA Search & Match createPatient payload
and optionally POST it to the WMDA API.

What this implementation supports from the provided Fhir_Input.json:
- Patient.identifier[use=usual].value -> patientId
- Patient.birthDate -> dateOfBirth
- Patient.gender -> sex (male->M, female->F, other->O, unknown->U)
- Genotype Observations (LOINC 84413-4) -> hla.<locus>.field1 / field2
  for loci such as HLA-A, HLA-B, HLA-C, DPA1, DPB1, DQA1, DQB1, DRB1/3/4/5

Optional WMDA fields can be passed on the CLI:
- --legal-terms
- --ethnicity
- --pool-country-code
- --transplant-centre-id
- --abo
- --rhesus
- --weight

Examples:
    python fhir_to_wmda_post.py Fhir_Input.json --pretty
    python fhir_to_wmda_post.py Fhir_Input.json --pretty --output wmda_payload.json
    python fhir_to_wmda_post.py Fhir_Input.json --post --base-url https://api.example.org \
        --token YOUR_TOKEN --legal-terms

Optional dependency for HTTP POST support:
    pip install requests
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any


SEX_MAP = {
    "male": "M",
    "female": "F",
    "other": "O",
    "unknown": "U",
}

GENE_TO_LOCUS = {
    "HLA-A": "a",
    "HLA-B": "b",
    "HLA-C": "c",
    "HLA-DPA1": "dpa1",
    "HLA-DPB1": "dpb1",
    "HLA-DQA1": "dqa1",
    "HLA-DQB1": "dqb1",
    "HLA-DRB1": "drb1",
    "HLA-DRB3": "drb3",
    "HLA-DRB4": "drb4",
    "HLA-DRB5": "drb5",
}


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def ensure_list(value: Any) -> list[Any]:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def first(seq, default=None):
    for item in seq:
        return item
    return default


def build_resource_index(bundle: dict[str, Any]) -> list[dict[str, Any]]:
    if bundle.get("resourceType") == "Bundle":
        return [entry["resource"] for entry in bundle.get("entry", []) if "resource" in entry]
    return [bundle]


def find_first_resource(resources: list[dict[str, Any]], resource_type: str) -> dict[str, Any] | None:
    return first(r for r in resources if r.get("resourceType") == resource_type)


def find_resources(resources: list[dict[str, Any]], resource_type: str) -> list[dict[str, Any]]:
    return [r for r in resources if r.get("resourceType") == resource_type]


def get_patient_id(patient: dict[str, Any]) -> str | None:
    identifiers = ensure_list(patient.get("identifier"))
    usual = first(i for i in identifiers if i.get("use") == "usual" and i.get("value"))
    if usual:
        return usual.get("value")
    first_non_empty = first(i for i in identifiers if i.get("value"))
    return first_non_empty.get("value") if first_non_empty else None


def map_sex(patient: dict[str, Any]) -> str | None:
    return SEX_MAP.get(patient.get("gender"))


def get_gene_name(observation: dict[str, Any]) -> str | None:
    for component in ensure_list(observation.get("component")):
        code = component.get("code", {})
        codings = ensure_list(code.get("coding"))
        if any(c.get("system") == "http://loinc.org" and c.get("code") == "48018-6" for c in codings):
            value = component.get("valueCodeableConcept", {})
            for coding in ensure_list(value.get("coding")):
                display = coding.get("display")
                code_val = coding.get("code")
                if display:
                    return display
                if code_val:
                    return code_val
    return None


def get_genotype_gl_string(observation: dict[str, Any]) -> str | None:
    value_cc = observation.get("valueCodeableConcept", {})
    for coding in ensure_list(value_cc.get("coding")):
        code = coding.get("code")
        if code:
            return code
    return None


def normalize_allele(raw: str) -> str:
    """
    Examples:
      HLA-A:01:01G     -> 01:01
      HLA-A*01:02      -> 01:02
      HLA-B*15:01:01G  -> 15:01
      HLA-C*03:04:01G  -> 03:04
    """
    m = re.search(r'(\d{2,3}):(\d{2,3})', raw)
    if not m:
        return raw
    return f"{m.group(1)}:{m.group(2)}"


def parse_gl_string(gl_string: str) -> list[str]:
    """
    Input:
      hla#3.23.0#HLA-A:01:01G+HLA-A*01:02
    Output:
      ["01:01", "01:02"]
    """
    if "#" in gl_string:
        gl_string = gl_string.split("#")[-1]
    parts = [p.strip() for p in gl_string.split("+") if p.strip()]
    return [normalize_allele(p) for p in parts]


def is_genotype_observation(observation: dict[str, Any]) -> bool:
    code = observation.get("code", {})
    for coding in ensure_list(code.get("coding")):
        if coding.get("system") == "http://loinc.org" and coding.get("code") == "84413-4":
            return True
    return False


def extract_hla(observations: list[dict[str, Any]]) -> dict[str, dict[str, str]]:
    hla: dict[str, dict[str, str]] = {}

    for obs in observations:
        if not is_genotype_observation(obs):
            continue

        gene_name = get_gene_name(obs)
        gl_string = get_genotype_gl_string(obs)
        if not gene_name or not gl_string:
            continue

        locus = GENE_TO_LOCUS.get(gene_name)
        if not locus:
            continue

        alleles = parse_gl_string(gl_string)
        if not alleles:
            continue

        entry: dict[str, str] = {}
        if len(alleles) >= 1:
            entry["field1"] = alleles[0]
        if len(alleles) >= 2:
            entry["field2"] = alleles[1]
        hla[locus] = entry

    return hla


def compact(value: Any) -> Any:
    if isinstance(value, dict):
        out = {k: compact(v) for k, v in value.items()}
        return {k: v for k, v in out.items() if v is not None and v != {} and v != []}
    if isinstance(value, list):
        out = [compact(v) for v in value]
        return [v for v in out if v is not None and v != {} and v != []]
    return value


def build_payload(
    bundle: dict[str, Any],
    *,
    legal_terms: bool | None = None,
    ethnicity: str | None = None,
    pool_country_code: str | None = None,
    transplant_centre_id: str | None = None,
    abo: str | None = None,
    rhesus: str | None = None,
    weight: float | None = None,
) -> dict[str, Any]:
    resources = build_resource_index(bundle)
    patient = find_first_resource(resources, "Patient")
    if not patient:
        raise ValueError("No Patient resource found in input")

    observations = find_resources(resources, "Observation")

    payload: dict[str, Any] = {
        "patientId": get_patient_id(patient),
        "hla": extract_hla(observations),
        "dateOfBirth": patient.get("birthDate"),
        "sex": map_sex(patient),
        "ethnicity": ethnicity,
        "poolCountryCode": pool_country_code,
        "transplantCentreId": transplant_centre_id,
        "abo": abo,
        "rhesus": rhesus,
        "weight": weight,
        "legalTerms": legal_terms,
    }

    return compact(payload)


def post_to_wmda(payload: dict[str, Any], base_url: str, token: str | None, timeout: int = 30):
    try:
        import requests
    except ImportError as exc:
        raise SystemExit(
            "The optional dependency 'requests' is required for --post. "
            "Install it with: pip install requests"
        ) from exc

    url = base_url.rstrip("/") + "/patients/createPatient"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
    }
    if token:
        headers["Authorization"] = f"Bearer {token}"

    response = requests.post(url, json=payload, headers=headers, timeout=timeout)
    response.raise_for_status()
    return response


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert FHIR Bundle JSON to WMDA createPatient payload and optionally POST it")
    parser.add_argument("input", type=Path, help="Path to FHIR input JSON")
    parser.add_argument("-o", "--output", type=Path, help="Write WMDA payload to file")
    parser.add_argument("--pretty", action="store_true", help="Pretty-print JSON payload")
    parser.add_argument("--legal-terms", action="store_true", help="Set legalTerms=true")
    parser.add_argument("--ethnicity", help="Set ethnicity")
    parser.add_argument("--pool-country-code", help="Set poolCountryCode")
    parser.add_argument("--transplant-centre-id", help="Set transplantCentreId")
    parser.add_argument("--abo", help="Set ABO value")
    parser.add_argument("--rhesus", help="Set rhesus value")
    parser.add_argument("--weight", type=float, help="Set weight")

    parser.add_argument("--post", action="store_true", help="POST the generated payload to WMDA")
    parser.add_argument("--base-url", help="WMDA API base URL, e.g. https://host/api/v1")
    parser.add_argument("--token", help="Bearer token for Authorization header")
    parser.add_argument("--timeout", type=int, default=30, help="HTTP timeout in seconds")
    parser.add_argument("--save-response", type=Path, help="Write HTTP response body to file")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    bundle = load_json(args.input)

    payload = build_payload(
        bundle,
        legal_terms=True if args.legal_terms else None,
        ethnicity=args.ethnicity,
        pool_country_code=args.pool_country_code,
        transplant_centre_id=args.transplant_centre_id,
        abo=args.abo,
        rhesus=args.rhesus,
        weight=args.weight,
    )

    text = json.dumps(payload, indent=2 if args.pretty else None, ensure_ascii=False)

    if args.output:
        args.output.write_text(text + ("\n" if not text.endswith("\n") else ""), encoding="utf-8")
    else:
        print(text)

    if args.post:
        if not args.base_url:
            raise SystemExit("--post requires --base-url")

        response = post_to_wmda(payload, args.base_url, args.token, timeout=args.timeout)

        print("\nHTTP", response.status_code, file=sys.stderr)
        try:
            response_text = json.dumps(response.json(), indent=2, ensure_ascii=False)
        except ValueError:
            response_text = response.text

        if args.save_response:
            args.save_response.write_text(response_text + ("\n" if not response_text.endswith("\n") else ""), encoding="utf-8")
        else:
            print(response_text)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
