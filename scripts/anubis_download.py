"""CLI for the Anubis-aware Dryad downloader (b01_clock.data.adapters.convexml)."""

import sys
from pathlib import Path

from b01_clock.data.adapters.convexml import download_guarded

if __name__ == "__main__":
    _url, _dest = sys.argv[1], Path(sys.argv[2])
    download_guarded(_url, _dest)
    print(f"saved {_dest} ({_dest.stat().st_size} bytes)")