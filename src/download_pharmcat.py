#!/usr/bin/env python3
"""
Download the PharmCAT JAR from GitHub releases into <project_root>/lib/.

Usage:
    python src/download_pharmcat.py
    python src/download_pharmcat.py --version 3.2.0
"""

import argparse
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

GITHUB_API = "https://api.github.com/repos/PharmGKB/PharmCAT/releases"
LIB_DIR = Path(__file__).resolve().parent.parent / "lib"


def _request(url: str) -> dict:
    req = urllib.request.Request(
        url, headers={"Accept": "application/vnd.github.v3+json"},
    )
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read().decode())


def get_latest_release_url() -> tuple[str, str]:
    data = _request(f"{GITHUB_API}/latest")
    version = data["tag_name"].lstrip("v")
    return _pick_jar_asset(data, version), version


def get_versioned_release_url(version: str) -> str:
    data = _request(f"{GITHUB_API}/tags/v{version}")
    return _pick_jar_asset(data, version)


def _pick_jar_asset(release_data: dict, version: str) -> str:
    for asset in release_data.get("assets", []) or []:
        name = asset["name"]
        if name.endswith("-all.jar") or (
            name.endswith(".jar") and "pharmcat" in name.lower()
        ):
            return asset["browser_download_url"]
    raise RuntimeError(f"No JAR asset found in release {version}")


def download_file(url: str, dest: Path) -> None:
    print(f"Downloading: {url}")
    print(f"Destination: {dest}")
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as resp:
        total = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        chunk_size = 256 * 1024  # 256 KB
        with open(dest, "wb") as f:
            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if total:
                    pct = downloaded / total * 100
                    print(
                        f"\r  {downloaded / 1024 / 1024:.1f} / "
                        f"{total / 1024 / 1024:.1f} MB ({pct:.0f}%)",
                        end="", flush=True,
                    )
    print()


def main() -> int:
    p = argparse.ArgumentParser(description="Download PharmCAT JAR.")
    p.add_argument("--version", "-V", default=None,
                   help="PharmCAT version (default: latest).")
    args = p.parse_args()

    LIB_DIR.mkdir(parents=True, exist_ok=True)

    try:
        if args.version:
            url = get_versioned_release_url(args.version)
            version = args.version
        else:
            url, version = get_latest_release_url()
        filename = url.rsplit("/", 1)[-1]
        dest = LIB_DIR / filename
        if dest.is_file():
            print(f"PharmCAT v{version} already downloaded: {dest}")
            return 0
        download_file(url, dest)
        print(f"PharmCAT v{version} downloaded to: {dest}")
        return 0
    except urllib.error.HTTPError as e:
        print(f"HTTP error: {e.code} {e.reason}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
