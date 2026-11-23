#!/usr/bin/env python3
"""Create a modern CHEMPOTPY PES or DPEM source scaffold."""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Mapping, Optional, Sequence, Tuple

try:
    from scripts.format_fortran import format_source_text
except ModuleNotFoundError:
    from format_fortran import format_source_text


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = REPOSITORY_ROOT / "chempotpy"
TEMPLATE_ROOT = REPOSITORY_ROOT / "templates"
NAME_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_+\-]*$")
PLACEHOLDER_PATTERN = re.compile(r"{{([A-Z0-9_]+)}}")


def validated_name(value: str, label: str) -> str:
    """Return a safe system or surface path component."""

    normalized = value[:-4] if value.lower().endswith(".f90") else value
    if not NAME_PATTERN.fullmatch(normalized):
        raise ValueError(
            "{} must start with a letter and contain only letters, digits, "
            "underscores, plus signs, or hyphens".format(label)
        )
    return normalized


def module_identifier(system: str, surface: str) -> str:
    """Return a valid, deterministic Fortran module-name prefix."""

    combined = "chempotpy_{}_{}".format(system, surface).lower()
    return re.sub(r"[^a-z0-9_]", "_", combined)


def parsed_atom_symbols(value: Optional[str], natoms: int) -> Tuple[str, ...]:
    """Parse and validate a comma-separated canonical atom ordering."""

    if value is None:
        return ()
    symbols = tuple(item.strip() for item in value.split(",") if item.strip())
    if len(symbols) != natoms:
        raise ValueError(
            "--atoms must provide exactly {} comma-separated symbols".format(
                natoms
            )
        )
    if any(re.fullmatch(r"[A-Z][a-z]?", symbol) is None for symbol in symbols):
        raise ValueError("--atoms contains an invalid element symbol")
    return symbols


def render_template(path: Path, substitutions: Mapping[str, str]) -> str:
    """Render one repository template and reject unresolved placeholders."""

    rendered = path.read_text()
    for name, value in substitutions.items():
        rendered = rendered.replace("{{{{{}}}}}".format(name), value)
    unresolved = sorted(set(PLACEHOLDER_PATTERN.findall(rendered)))
    if unresolved:
        raise ValueError(
            "template has unresolved placeholders: {}".format(
                ", ".join(unresolved)
            )
        )
    return rendered


def scaffold_surface(
    system: str,
    surface: str,
    natoms: int,
    nstates: int,
    entry_point: str,
    functional_form: str,
    common_name: str,
    reference: str,
    derivatives: int = 0,
    bodies: Optional[str] = None,
    atoms: Optional[str] = None,
    package_root: Path = PACKAGE_ROOT,
    template_root: Path = TEMPLATE_ROOT,
) -> Tuple[Path, ...]:
    """Write a non-overwriting PES scaffold and return the created paths."""

    system = validated_name(system, "system")
    surface = validated_name(surface, "surface")
    if natoms < 1 or nstates < 1:
        raise ValueError("natoms and nstates must be positive")
    if derivatives < 0:
        raise ValueError("derivatives must be non-negative")
    normalized_entry = entry_point.lower()
    if normalized_entry not in {"pes", "dpem"}:
        raise ValueError("entry_point must be 'pes' or 'dpem'")
    if not functional_form.strip() or not reference.strip():
        raise ValueError("functional form and reference must not be empty")

    system_directory = package_root / system
    initialization = system_directory / "__init__.py"
    source = system_directory / (surface + ".f90")
    if source.exists():
        raise FileExistsError("refusing to overwrite {}".format(source))

    atom_symbols = parsed_atom_symbols(atoms, natoms)
    create_system = not initialization.exists()
    if create_system and not atom_symbols:
        raise ValueError(
            "--atoms is required when creating a new system directory"
        )

    substitutions = {
        "SYSTEM": system,
        "SURFACE": surface,
        "MODULE_NAME": module_identifier(system, surface),
        "FUNCTIONAL_FORM": functional_form.strip(),
        "COMMON_NAME": common_name.strip() or surface,
        "REFERENCE": reference.strip(),
        "DERIVATIVES": str(derivatives),
        "BODIES": str(bodies if bodies is not None else natoms),
        "NATOMS": str(natoms),
        "NSTATES": str(nstates),
    }
    template = template_root / (normalized_entry + "_module.f90.in")
    source_text = format_source_text(render_template(template, substitutions))

    system_directory.mkdir(parents=True, exist_ok=True)
    created = []
    if create_system:
        initialization_text = render_template(
            template_root / "system_init.py.in",
            {
                "SYSTEM": system,
                "SURFACE": surface,
                "ATOM_SYMBOLS": repr(atom_symbols),
            },
        )
        initialization.write_text(initialization_text)
        created.append(initialization)
    source.write_text(source_text)
    created.append(source)
    return tuple(created)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Parse command-line arguments and create one PES scaffold."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("system")
    parser.add_argument("surface")
    parser.add_argument("--natoms", type=int, required=True)
    parser.add_argument("--nstates", type=int, required=True)
    parser.add_argument(
        "--entry-point", choices=("pes", "dpem"), default="pes"
    )
    parser.add_argument("--functional-form", required=True)
    parser.add_argument("--common-name", default="")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--derivatives", type=int, default=0)
    parser.add_argument(
        "--bodies",
        help="body count or 'variable'; defaults to natoms",
    )
    parser.add_argument(
        "--atoms",
        help="canonical comma-separated atom order; required for a new system",
    )
    args = parser.parse_args(argv)
    try:
        created = scaffold_surface(
            system=args.system,
            surface=args.surface,
            natoms=args.natoms,
            nstates=args.nstates,
            entry_point=args.entry_point,
            functional_form=args.functional_form,
            common_name=args.common_name,
            reference=args.reference,
            derivatives=args.derivatives,
            bodies=args.bodies,
            atoms=args.atoms,
        )
    except (FileExistsError, OSError, ValueError) as error:
        parser.error(str(error))
    for path in created:
        print(path.relative_to(REPOSITORY_ROOT))
    print(
        "Next: replace CHEMPOTPY_IMPLEMENTATION_REQUIRED, add independently "
        "verified references, and follow CONTRIBUTING.md."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
