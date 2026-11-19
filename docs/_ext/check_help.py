"""Check the MATLAB help-text converter on every documented function.

Runs a full Sphinx build into a temporary folder with the converter's dump
switched on, then:

* prints the reST produced for each documented object (``--show all``,
  ``--show NAME``; default: none), with the raw help text if ``--raw``;
* lists objects shown verbatim because parsing failed;
* lists docutils/Sphinx warnings raised inside docstrings;
* scans the HTML for markup that leaked through as text (":func:",
  ":param", ":returns:", ".. code-block", "``");
* counts objects rendered with Parameters / Returns lists.

Usage (from the repository root)::

    ~/.venvs/matlab-docs/bin/python docs/_ext/check_help.py [--show all|NAME] [--raw]

Exit status is 1 if any check fails.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import tempfile
from pathlib import Path

from sphinx.cmd.build import build_main

DOCS = Path(__file__).resolve().parent.parent


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--show", default="", help="'all' or an object name to print")
    ap.add_argument("--raw", action="store_true", help="also print the raw help text")
    args = ap.parse_args()

    tmp = Path(tempfile.mkdtemp(prefix="matlab_help_"))
    dump, warn = tmp / "dump.json", tmp / "warnings.txt"
    os.environ["MATLAB_HELP_DUMP"] = str(dump)
    code = build_main(["-E", "-q", "-b", "html", "-w", str(warn), str(DOCS), str(tmp / "html")])
    if code != 0:
        print(f"Sphinx build failed (exit {code})")
        return 1
    entries = json.loads(dump.read_text())

    for e in entries:
        if args.show == "all" or (args.show and e["name"].split(".")[-1] == args.show.split(".")[-1]):
            print("=" * 78)
            print(f"{e['what']} {e['name']}" + ("   [VERBATIM FALLBACK]" if e["fallback"] else ""))
            if args.raw:
                print("-" * 30 + " raw " + "-" * 43)
                print("\n".join(e["raw"]))
            print("-" * 30 + " reST " + "-" * 42)
            print("\n".join(e["rst"]))

    failed = False
    fallbacks = [e["name"] for e in entries if e["fallback"]]
    warnings = [l for l in warn.read_text().splitlines() if l.strip()] if warn.exists() else []
    doc_warn = [l for l in warnings if "docstring of" in l or "matlab_help" in l]
    other_warn = [l for l in warnings if l not in doc_warn and not l.startswith(" ")]

    leaks = []
    for page in sorted((tmp / "html" / "api").glob("*.html")):
        html = page.read_text()
        body = re.sub(r"<(script|style)\b.*?</\1>", "", html, flags=re.S)
        text = re.sub(r"<[^>]+>", " ", body)
        for pat in (":mat:func:", ":func:", ":param ", ":type ", ":returns:", ".. code-block", ".. rubric", "``"):
            n = text.count(pat)
            if n:
                leaks.append(f"{page.name}: {n} x {pat!r}")

    with_params = sum(1 for e in entries if any(l.startswith(":param ") for l in e["rst"]))
    with_returns = sum(1 for e in entries if any(l.startswith(":returns:") for l in e["rst"]))

    print("=" * 78)
    print(f"documented objects with help text : {len(entries)}")
    print(f"  with Parameters                 : {with_params}")
    print(f"  with Returns                    : {with_returns}")
    print(f"verbatim fallbacks                : {len(fallbacks)}")
    for n in fallbacks:
        print(f"    {n}")
    print(f"warnings inside docstrings        : {len(doc_warn)}")
    for l in doc_warn:
        print(f"    {l}")
    print(f"other build warnings              : {len(other_warn)}")
    for l in other_warn:
        print(f"    {l}")
    print(f"markup leaked into HTML           : {len(leaks)}")
    for l in leaks:
        print(f"    {l}")
    failed = bool(fallbacks or doc_warn or leaks)
    print("RESULT:", "FAIL" if failed else "OK")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
