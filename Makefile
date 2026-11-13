.PHONY: all audit build check check-binaries check-conda check-embedded-parameters check-fortran check-generated check-modernization check-modernization-source check-numerical-regression check-python-extensions check-style clean test test-surfaces

PYTHON ?= python3
BUILD_FLAGS ?=
BUILD_OUTPUT ?= chempotpy
NUMERICAL_BASELINE ?= tests/data/all_surface_references.json

audit:
	$(PYTHON) scripts/audit_pes.py

check-fortran:
	$(PYTHON) scripts/check_fortran.py

check-binaries:
	$(PYTHON) scripts/check_binaries.py --output-root $(BUILD_OUTPUT)

check-python-extensions:
	$(PYTHON) scripts/check_python_extensions.py --output-root $(BUILD_OUTPUT)

check-conda: build check-binaries check-python-extensions

check-style:
	$(PYTHON) scripts/check_fortran_style.py

check-modernization-source:
	$(PYTHON) scripts/check_modernization.py

check-embedded-parameters:
	$(PYTHON) scripts/generate_embedded_network_parameters.py --check
	$(PYTHON) scripts/generate_embedded_pes_parameters.py --check

check-generated: check-embedded-parameters
	$(PYTHON) scripts/generate_physical_constants.py --check

build:
	$(PYTHON) scripts/build_pes.py $(BUILD_FLAGS)

all: build

check: audit check-modernization-source check-generated test check-fortran check-style

check-modernization: check

test:
	$(PYTHON) -m unittest discover -s tests -v

test-surfaces:
	$(PYTHON) -m unittest tests.test_all_surfaces -v

check-numerical-regression:
	$(PYTHON) scripts/check_numerical_regression.py \
		--baseline $(NUMERICAL_BASELINE) --jobs 4

clean:
	$(PYTHON) scripts/clean_generated.py
