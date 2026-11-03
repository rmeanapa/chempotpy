.PHONY: all audit build check check-binaries check-fortran check-style clean test test-surfaces

PYTHON ?= python3
BUILD_FLAGS ?=
BUILD_OUTPUT ?= chempotpy

audit:
	$(PYTHON) scripts/audit_pes.py

check-fortran:
	$(PYTHON) scripts/check_fortran.py

check-binaries:
	$(PYTHON) scripts/check_binaries.py --output-root $(BUILD_OUTPUT)

check-style:
	$(PYTHON) scripts/check_fortran_style.py

build:
	$(PYTHON) scripts/build_pes.py $(BUILD_FLAGS)

all: build

check: audit test check-fortran check-style

test:
	$(PYTHON) -m unittest discover -s tests -v

test-surfaces:
	$(PYTHON) -m unittest tests.test_all_surfaces -v

clean:
	$(PYTHON) scripts/clean_generated.py
