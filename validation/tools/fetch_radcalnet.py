"""Fetch the RadCalNet records a validation case used and recompute its reference side.

Uses spectrAccess (https://github.com/SpectraWorksOSS/SpectrAccess, ``pip install
spectraccess``) to download the exact RadCalNet daily files named in a case's
``case.json``, select the record at the case's observation time, and write the
spectrum to CSV. With ``--srf`` it also band-integrates the spectrum the way the
case CSVs were made, so the ``radcalnet_*`` columns can be reproduced and diffed.

RadCalNet needs a free portal account. Store it once with
``python -c "import spectraccess; spectraccess.login('radcalnet')"``
(OS keyring); this script never reads credentials from files or arguments.

Examples::

    python validation/tools/fetch_radcalnet.py validation/cases/landsat9_rvus_2025-06-16
    python validation/tools/fetch_radcalnet.py validation/cases/landsat9_rvus_2025-06-16 \
        --srf oli2_rsr.csv

``--srf`` takes a CSV with columns ``band,wavelength_nm,response`` (for example
the USGS Landsat relative spectral response tables, one row per band and
wavelength). Band value = SRF-weighted mean of the RadCalNet spectrum on its
10 nm grid. Band uncertainty = SRF-weighted mean of the per-wavelength
uncertainty, i.e. errors treated as fully correlated across wavelengths within
the band (an upper bound). RadCalNet states no coverage factor; the case files
read it as k=1, which is our declared assumption, not RadCalNet's statement.
``--local`` parses already-downloaded RadCalNet files instead of fetching.
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def _records(case: dict) -> list[dict]:
    rc = case["radcalnet"]
    out = list(rc.get("records", []))
    for key in ("surface_record", "toa_record"):
        if key in rc:
            out.append(rc[key])
    return out


def _load(record: dict, local_dir: Path | None, out_dir: Path) -> tuple[str, pd.DataFrame]:
    from spectraccess.connectors.radcalnet import RadCalNetConnector

    site, filename = record["record_file"].split("/", 1)
    connector = RadCalNetConnector()
    if local_dir is not None:
        raw = (local_dir / filename).read_bytes()
    else:
        targets = [t for t in connector.discover(site, kind=None) if t.filename == filename]
        if not targets:
            raise SystemExit(f"{filename}: not listed by RadCalNet for site {site}")
        raw = Path(connector.fetch(targets[0], dest=out_dir / filename)).read_bytes()
    return filename, connector.parse(raw, source_file=filename)


def _at_time(frame: pd.DataFrame, observed_at: str) -> pd.DataFrame:
    when = pd.Timestamp(datetime.fromisoformat(observed_at))
    stamps = pd.to_datetime(frame["timestamp"], utc=True)
    rows = frame[stamps == when]
    if rows.empty:
        raise SystemExit(f"no RadCalNet record at {observed_at}")
    return rows.sort_values("wavelength_nm")


def _band_integrate(spectrum: pd.DataFrame, quantity: str, srf: pd.DataFrame) -> list[dict]:
    wl = spectrum["wavelength_nm"].to_numpy(float)
    val = spectrum[quantity].to_numpy(float)
    unc = spectrum[f"{quantity}_unc"].to_numpy(float)
    rows = []
    for band, group in srf.groupby("band", sort=False):
        g = group.sort_values("wavelength_nm")
        w = np.interp(wl, g["wavelength_nm"].to_numpy(float), g["response"].to_numpy(float),
                      left=0.0, right=0.0)
        ok_v = np.isfinite(val) & (w > 0)
        ok_u = ok_v & np.isfinite(unc)
        rows.append({
            "band": band,
            quantity: float(np.sum(w[ok_v] * val[ok_v]) / np.sum(w[ok_v])) if ok_v.any() else None,
            f"{quantity}_unc_k1": float(np.sum(w[ok_u] * unc[ok_u]) / np.sum(w[ok_u])) if ok_u.any() else None,
            "srf_weight_with_valid_value": float(np.sum(w[ok_v]) / np.sum(w)) if np.sum(w) > 0 else 0.0,
        })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("case_dir", type=Path)
    ap.add_argument("--srf", type=Path, help="CSV with band,wavelength_nm,response")
    ap.add_argument("--local", type=Path, help="directory holding already-downloaded RadCalNet files")
    ap.add_argument("--out", type=Path, help="output directory (default: <case_dir>/radcalnet)")
    args = ap.parse_args()

    case = json.loads((args.case_dir / "case.json").read_text())
    out_dir = args.out or args.case_dir / "radcalnet"
    out_dir.mkdir(parents=True, exist_ok=True)
    srf = pd.read_csv(args.srf) if args.srf else None

    for record in _records(case):
        filename, frame = _load(record, args.local, out_dir)
        quantity = "toa_reflectance" if filename.endswith(".output") else "surface_reflectance"
        spectrum = _at_time(frame, record["observed_at_utc"])
        stem = Path(filename).name.replace(".", "_")
        spectrum.to_csv(out_dir / f"{stem}_spectrum.csv", index=False)
        print(f"{filename} @ {record['observed_at_utc']}: {len(spectrum)} wavelengths")
        if srf is not None:
            bands = _band_integrate(spectrum, quantity, srf)
            with (out_dir / f"{stem}_bands.csv").open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(bands[0]), lineterminator="\n")
                w.writeheader()
                w.writerows(bands)
            for b in bands:
                print(f"  {b['band']}: {b[quantity]} +/- {b[f'{quantity}_unc_k1']}")


if __name__ == "__main__":
    main()
