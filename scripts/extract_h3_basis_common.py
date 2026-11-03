#!/usr/bin/env python3
"""Replace the duplicated H3 /indy/ COMMON block with the bundled core API."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


INDEX_ROUTINE = re.compile(
    r"\n!> @brief Evaluate indexa3[^\n]*\n"
    r"subroutine indexa3\( mbig, mtop \).*?"
    r"\nend subroutine indexa3\n",
    flags=re.DOTALL,
)


def replace_once(source: str, old: str, new: str, description: str) -> str:
    """Replace one required source fragment or fail with a useful message."""

    count = source.count(old)
    if count != 1:
        raise ValueError("{}: expected one match, found {}".format(
            description, count
        ))
    return source.replace(old, new, 1)


def migrate_source(source: str) -> str:
    """Return one H3 source using the private bundled basis table."""

    use_match = re.search(
        r"^(\s*use\s+[a-z0-9_]+data,\s*only:[^\n]+)$",
        source,
        flags=re.IGNORECASE | re.MULTILINE,
    )
    if use_match is None:
        raise ValueError("fit-data import not found")
    source = source[:use_match.end()] + (
        "\n    use chempotpy_h3_basis, only: get_h3_basis_term, "
        "initialize_h3_basis"
    ) + source[use_match.end():]

    source = replace_once(
        source,
        "    integer :: ii, iind, ipar, it, jind, kind, maxi, maxi2, "
        "maxi2a, maxi2b\n"
        "    save\n"
        "    integer, parameter :: n2 = 5000\n"
        "    common /indy/iind(n2),jind(n2),kind(n2),ipar(n2),maxi,maxi2\n",
        "    integer :: coefficient_index, first_power, ii, it, maxi2a, "
        "maxi2b, second_power, third_power\n"
        "    save\n",
        "COMMON declaration",
    )
    source = replace_once(
        source,
        "    maxi2 = 0\n"
        "    maxi  = 0\n"
        "    call indexa3(12,12)\n"
        "    maxi2a = maxi2\n\n"
        "    maxi = 0\n"
        "    call indexa3(11,11)\n"
        "    maxi2b = maxi2-maxi2a\n",
        "    call initialize_h3_basis( 12, 11, maxi2a, maxi2b )\n",
        "basis initialization",
    )
    source = replace_once(
        source,
        "    do ii=1,maxi2a\n"
        "        rprod = rho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))\n"
        "        sum1  = sum1+xlin1(ipar(ii))*rprod\n"
        "        sum2  = sum2+xlin2(ipar(ii))*rprod\n"
        "    enddo\n\n"
        "    do ii=maxi2a+1,maxi2a+maxi2b\n"
        "        rprod = rho1t(kind(ii))*rho2t(jind(ii))*rho3t(iind(ii))\n"
        "        sum3  = sum3+xlin3(ipar(ii))*rprod\n"
        "        sum4  = sum4+xlin4(ipar(ii))*rprod\n"
        "    enddo\n",
        "    do ii=1,maxi2a\n"
        "        call get_h3_basis_term( ii, first_power, second_power, "
        "third_power, coefficient_index )\n"
        "        rprod = rho1t(third_power)*rho2t(second_power)*"
        "rho3t(first_power)\n"
        "        sum1  = sum1+xlin1(coefficient_index)*rprod\n"
        "        sum2  = sum2+xlin2(coefficient_index)*rprod\n"
        "    enddo\n\n"
        "    do ii=maxi2a+1,maxi2a+maxi2b\n"
        "        call get_h3_basis_term( ii, first_power, second_power, "
        "third_power, coefficient_index )\n"
        "        rprod = rho1t(third_power)*rho2t(second_power)*"
        "rho3t(first_power)\n"
        "        sum3  = sum3+xlin3(coefficient_index)*rprod\n"
        "        sum4  = sum4+xlin4(coefficient_index)*rprod\n"
        "    enddo\n",
        "basis evaluation",
    )
    source, count = INDEX_ROUTINE.subn("\n", source, count=1)
    if count != 1:
        raise ValueError("indexa3 routine not found")
    return source


def main() -> int:
    """Migrate requested H3 sources and report the number changed."""

    parser = argparse.ArgumentParser()
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    changed = 0
    for path in arguments.sources:
        original = path.read_text(errors="strict")
        migrated = migrate_source(original)
        changed += migrated != original
        if arguments.write and migrated != original:
            path.write_text(migrated)
    print("fortran_sources_changed: {}".format(changed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
