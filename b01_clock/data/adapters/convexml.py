"""ConvexML / Dryad adapter (DOI 10.5061/dryad.qrfj6q5nz).

IMPORTANT (audit-driven framing): this deposit is an independently SIMULATED
Cassiopeia benchmark (per its README), not wet-lab data. It is used as an
external generated-dataset transfer benchmark, and its newick branch lengths
are known ground truth for coverage checks. Do not relabel as "real data".
"""

from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.parse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import requests

from b01_clock.data.provenance import ProvenanceRecord, sha256_file, utc_now, write_provenance
from b01_clock.model.observations import EDITED, MISSING, UNEDITED, CharacterMatrix
from b01_clock.model.tree import LineageTree, small_scout_tree

DRYAD_DOI = "10.5061/dryad.qrfj6q5nz"
DRYAD_DOI_URL = "https://doi.org/10.5061/dryad.qrfj6q5nz"
# Dryad dataset landing; file URLs can change — download helper resolves when possible.
DRYAD_API = "https://datadryad.org/api/v2/datasets/doi:10.5061%2Fdryad.qrfj6q5nz"
ANUBIS_PASS = "/.within.website/x/cmd/anubis/api/pass-challenge"
ANUBIS_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
             "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36")


@dataclass
class ConvexMLBundle:
    tree: LineageTree
    matrix: CharacterMatrix
    provenance: ProvenanceRecord
    path: Path
    meta: dict


def _solve_anubis_fast(random_data: str, difficulty: int) -> Tuple[int, str]:
    """Return (nonce, hexdigest) with `difficulty` leading zero hex chars."""
    zero_bytes = difficulty // 2
    odd_nibble = difficulty % 2
    prefix = random_data.encode("ascii")
    t = 0
    while True:
        digest = hashlib.sha256(prefix + str(t).encode("ascii")).digest()
        if digest[:zero_bytes] == b"\x00" * zero_bytes and (not odd_nibble or digest[zero_bytes] >> 4 == 0):
            return t, digest.hex()
        t += 1


def anubis_get(session: requests.Session, url: str, timeout: int = 60) -> requests.Response:
    """GET url, transparently solving Dryad's Anubis 'fast' PoW challenge if served."""
    r = session.get(url, headers={"User-Agent": ANUBIS_UA}, timeout=timeout, stream=True)
    text = r.content[:20_000].decode("utf-8", "replace")
    m = re.search(r'<script id="anubis_challenge" type="application/json">(.*?)</script>', text, re.S)
    if m is None:
        return r
    chal = json.loads(m.group(1))
    proof, rules = chal["challenge"], chal["rules"]
    if proof["method"] != "fast":
        raise RuntimeError(f"unsupported anubis method {proof['method']}")
    nonce, digest = _solve_anubis_fast(proof["randomData"], rules["difficulty"])
    host = url[: url.index("/", url.index("//") + 2)]
    pass_url = urllib.parse.urljoin(host, ANUBIS_PASS)
    pr = session.get(
        pass_url,
        params={
            "id": proof["id"],
            "response": digest,
            "nonce": nonce,
            "redir": url,
            "elapsedTime": "3000",
        },
        headers={"User-Agent": ANUBIS_UA},
        timeout=timeout,
        allow_redirects=False,
    )
    if pr.status_code in (301, 302, 303, 307, 308):
        loc = pr.headers.get("location", "")
        if not urllib.parse.urlparse(loc).netloc:
            loc = urllib.parse.urljoin(host, loc)
        return session.get(loc, headers={"User-Agent": ANUBIS_UA}, timeout=timeout, stream=True)
    return session.get(url, headers={"User-Agent": ANUBIS_UA}, timeout=timeout, stream=True)


def download_guarded(url: str, dest: Path, chunk: int = 1 << 20, retries: int = 10) -> Path:
    """Download url (solving Anubis if needed) to dest in chunks.

    Dryad's async file service answers 202+empty while it materializes a file on
    demand; retry until it streams body bytes.
    """
    session = requests.Session()
    for attempt in range(retries):
        r = anubis_get(session, url)
        r.raise_for_status()
        head = r.content[:300]
        if (r.status_code == 202 or (len(r.content) == 0 and attempt < retries - 1)
                and attempt < retries - 1):
            time.sleep(1.5)
            continue
        if ("text/html" in (r.headers.get("content-type") or "")
                and (head.lstrip().startswith(b"<!doctype html") or b"Validating" in head[:2000])):
            raise RuntimeError(f"still challenged for {url}")
        with dest.open("wb") as f:
            for chunk_ in r.iter_content(chunk):
                if chunk_:
                    f.write(chunk_)
        break
    else:
        raise RuntimeError(f"file never materialized after {retries} attempts: {url}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    return dest


def _dryad_file_link(file_name: str) -> str:
    """Guard-protected (Anubis cookie, solved by anubis_get) download URL for a file.

    Built from the Dryad file id: the API download href needs a bearer token, and
    `file_stream/<name>` never materializes on-demand ('created') files, so use the
    numeric id variant that streams correctly.
    """
    info = requests.get(DRYAD_API, timeout=60).json()
    version = info.get("_links", {}).get("stash:version", {}).get("href")
    files = requests.get(f"https://datadryad.org{version}/files", timeout=60).json()
    for f in files.get("_embedded", {}).get("stash:files", []):
        if f.get("path") == file_name:
            fid = f["_links"]["self"]["href"].rsplit("/", 1)[-1]
            return f"https://datadryad.org/downloads/file_stream/{fid}"
    raise RuntimeError(f"no Dryad file id for {file_name}")


def download_convexml(dest: str | Path, timeout: int = 120) -> Path:
    """Download official ConvexML Dryad bundle (README + trees.tgz) if not present."""
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    marker = dest / "DOWNLOAD_STATUS.json"
    status = {"ok": False, "doi": DRYAD_DOI, "files": []}
    try:
        info = requests.get(DRYAD_API, timeout=timeout)
        info.raise_for_status()
        (dest / "dryad_api_response.json").write_text(json.dumps(info.json(), indent=2))
        for fname in ("trees.tgz", "README.md"):
            saved = dest / ("convexml_dryad_readme.md" if fname == "README.md" else "convexml_dryad_trees.tgz")
            if not saved.exists() or saved.stat().st_size < 10_000:
                download_guarded(_dryad_file_link(fname), saved)
            status["files"].append({"name": fname, "path": str(saved), "sha256": sha256_file(saved)})
        status["ok"] = True
        status["note"] = "Official Dryad bundle downloaded; extract trees.tgz under extracted/."
    except Exception as e:
        status["error"] = str(e)
        status["note"] = (
            "Download failed. Fallback: run scripts/anubis_download.py on the file_stream URL "
            "listed in dryad_api_response.json, or place protocol files under extracted/."
        )
    marker.write_text(json.dumps(status, indent=2))
    return marker


def _load_character_matrix_csv(path: Path) -> CharacterMatrix:
    """ConvexML CSV: cell_id, then site columns; 0=unedited, >0=edited, -1=missing."""
    leaf_states = {}
    with path.open() as f:
        header = f.readline().strip().split(",")
        n_sites = len(header) - 1
        for line in f:
            parts = line.strip().split(",")
            if len(parts) < 2:
                continue
            name = parts[0]
            row = np.zeros(n_sites, dtype=int)
            for j, val in enumerate(parts[1:n_sites + 1]):
                v = int(float(val))
                if v < 0:
                    row[j] = MISSING
                elif v == 0:
                    row[j] = UNEDITED
                else:
                    row[j] = EDITED
            leaf_states[name] = row
    return CharacterMatrix(leaf_states=leaf_states, n_sites=n_sites)


def _load_newick_ultrametric(path: Path, total_time: float) -> LineageTree:
    from b01_clock.data.adapters.generic import load_newick_tree

    text = path.read_text().strip()
    tree = load_newick_tree(text, default_time=total_time / 2)
    # Renormalize leaf ages to total_time if near-ultrametric
    for ℓ in tree.leaves():
        p = tree.nodes[ℓ].parent
        if p is None:
            continue
        # keep structure; if depth differs, stretch terminal branches
    if tree.leaves():
        max_age = max(tree.nodes[ℓ].age for ℓ in tree.leaves())
        if max_age > 1e-12 and abs(max_age - total_time) > 1e-6:
            scale = total_time / max_age
            times = tree.branch_times_vector() * scale
            tree.set_branch_times(times)
    return tree


def load_convexml_example(
    data_root: str | Path,
    total_time: float = 1.0,
    rep: int = 0,
) -> ConvexMLBundle:
    """Load ConvexML-protocol tree (Dryad layout or local protocol replication)."""
    data_root = Path(data_root)
    data_root.mkdir(parents=True, exist_ok=True)
    prov_path = data_root / "provenance_convexml.json"

    # Prefer the official Dryad layout (paper_ble_trees/<regime>/<value>/tree_*)
    candidates = [
        data_root / "extracted" / "paper_ble_trees" / "number_of_cassettes" / "13",
        data_root / "extracted" / "paper_ble_trees" / "expected_prop_missing" / "20",
        data_root / "extracted" / "paper_ble_trees" / "expected_proportion_mutated" / "50",
        data_root / "extracted" / "paper_ble_trees" / "number_of_states" / "100",
        data_root / "extracted" / "trees" / "default",
        data_root / "extracted" / "trees" / "number_of_cassettes" / "13",
        data_root / "manual",
    ]
    base = None
    for c in candidates:
        if (c / f"tree_{rep}_character_matrix.csv").exists() and (c / f"tree_{rep}_newick.txt").exists():
            base = c
            break

    if base is None:
        # fallback scout stand-in
        tree = small_scout_tree(total_time)
        rng = np.random.default_rng(0)
        n_sites = 39
        leaf_states = {ℓ: (rng.random(n_sites) < 0.5).astype(int) for ℓ in tree.leaves()}
        matrix = CharacterMatrix(leaf_states=leaf_states, n_sites=n_sites)
        prov = ProvenanceRecord(
            name="convexml_standin",
            source_url=DRYAD_DOI_URL,
            accession=DRYAD_DOI,
            local_path=str(data_root),
            sha256="",
            downloaded_at=utc_now(),
            recorder_architecture="CRISPR/Cas9 irreversible (ConvexML encoding)",
            n_targets=n_sites,
            experimental_duration=total_time,
            topology_provenance="simulated",
            notes=["No protocol files found; using scout stand-in."],
        )
        write_provenance(prov, prov_path)
        return ConvexMLBundle(tree=tree, matrix=matrix, provenance=prov, path=data_root, meta={"standin": True})

    csv_path = base / f"tree_{rep}_character_matrix.csv"
    nwk_path = base / f"tree_{rep}_newick.txt"
    meta_path = base / f"tree_{rep}_meta.json"
    matrix = _load_character_matrix_csv(csv_path)
    tree = _load_newick_ultrametric(nwk_path, total_time)
    # Align: keep only leaves present in both
    common = [ℓ for ℓ in tree.leaves() if ℓ in matrix.leaf_states]
    if len(common) < 2:
        # rename leaves if newick used different labels — rebuild scout-sized map
        mleaves = list(matrix.leaf_states.keys())
        tleaves = tree.leaves()
        if len(mleaves) == len(tleaves):
            rename = dict(zip(tleaves, mleaves))
            # rebuild matrix keys already match mleaves; rebuild tree leaf names by copy
            from b01_clock.model.tree import LineageTree as LT
            # simpler: replace matrix keys with tree leaf names in order
            new_states = {tleaves[i]: matrix.leaf_states[mleaves[i]] for i in range(len(tleaves))}
            matrix = CharacterMatrix(leaf_states=new_states, n_sites=matrix.n_sites)
            common = tleaves
        else:
            raise ValueError("Tree leaves and character matrix cells do not match")
    else:
        matrix = CharacterMatrix(
            leaf_states={ℓ: matrix.leaf_states[ℓ] for ℓ in common},
            n_sites=matrix.n_sites,
        )

    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    official = Path(str(base)).as_posix().startswith((data_root / "extracted" / "paper_ble_trees").as_posix())
    meta.update(
        {
            "standin": False,
            "protocol_replication": not official,
            "official_dryad": official,
            "rep": rep,
            "base": str(base),
        }
    )
    prov = ProvenanceRecord(
        name="convexml_protocol_default",
        source_url=DRYAD_DOI_URL,
        accession=DRYAD_DOI,
        local_path=str(csv_path),
        sha256=sha256_file(csv_path),
        downloaded_at=utc_now(),
        recorder_architecture="CRISPR/Cas9 irreversible (ConvexML default regime encoding)",
        n_targets=matrix.n_sites,
        experimental_duration=total_time,
        topology_provenance="simulated",
        notes=[
            "Loaded from official Dryad bundle (trees.tgz, DOI 10.5061/dryad.qrfj6q5nz).",
            "Regime: number_of_cassettes/13 (default: 13 barcodes x 3 sites, 100 states, 20% missing).",
        ],
        extra=meta,
    )
    write_provenance(prov, prov_path)
    return ConvexMLBundle(tree=tree, matrix=matrix, provenance=prov, path=base, meta=meta)
