#!/usr/bin/env python3
"""
Download the latest PharmCAT JAR file from GitHub releases.

Usage:
    python download_pharmcat.py
    python download_pharmcat.py --version 2.15.4
"""

import argparse
import json
import sys
import urllib.request
import urllib.error
from pathlib import Path

GITHUB_API = "https://api.github.com/repos/PharmGKB/PharmCAT/releases"
LIB_DIR = Path(__file__).parent / "lib"


def get_latest_release_url() -> tuple[str, str]:
    """Get the download URL and version of the latest PharmCAT release."""
    req = urllib.request.Request(
        f"{GITHUB_API}/latest",
        headers={"Accept": "application/vnd.github.v3+json"},
    )
    with urllib.request.urlopen(req) as resp:
        data = json.loads(resp.read().decode())

    version = data["tag_name"].lstrip("v")
    for asset in data.get("assets", []):
        name = asset["name"]
        if name.endswith("-all.jar") or (name.endswith(".jar") and "pharmcat" in name.lower()):
            return asset["browser_download_url"], version

    raise RuntimeError(f"No JAR asset found in release {version}")


def get_versioned_release_url(version: str) -> str:
    """Get download URL for a specific PharmCAT version."""
    req = urllib.request.Request(
        f"{GITHUB_API}/tags/v{version}",
        headers={"Accept": "application/vnd.github.v3+json"},
    )
    with urllib.request.urlopen(req) as resp:
        data = json.loads(resp.read().decode())

    for asset in data.get("assets", []):
        name = asset["name"]
        if name.endswith("-all.jar") or (name.endswith(".jar") and "pharmcat" in name.lower()):
            return asset["browser_download_url"]

    raise RuntimeError(f"No JAR asset found in release v{version}")


def download_file(url: str, dest: Path) -> None:
    """Download a file with progress reporting."""
    print(f"Downloading: {url}")
    print(f"Destination: {dest}")

    req = urllib.request.Request(url)
    with urllib.request.urlopen(req) as resp:
        total = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        chunk_size = 1024 * 256  # 256 KB

        with open(dest, "wb") as f:
            while True:
                chunk = resp.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                if total > 0:
                    pct = downloaded / total * 100
                    mb = downloaded / (1024 * 1024)
                    total_mb = total / (1024 * 1024)
                    print(f"\r  {mb:.1f} / {total_mb:.1f} MB ({pct:.0f}%)", end="", flush=True)

    print()  # newline after progress


def main() -> int:
    parser = argparse.ArgumentParser(description="Download PharmCAT JAR file.")
    parser.add_argument(
        "--version", "-V",
        default=None,
        help="PharmCAT version to download (default: latest)",
    )
    args = parser.parse_args()

    LIB_DIR.mkdir(parents=True, exist_ok=True)

    try:
        if args.version:
            url = get_versioned_release_url(args.version)
            version = args.version
        else:
            url, version = get_latest_release_url()

        filename = url.split("/")[-1]
        dest = LIB_DIR / filename

        if dest.exists():
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
