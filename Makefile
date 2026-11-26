.PHONY: all audit build check check-binaries check-conda check-embedded-parameters check-fortran check-fortran-warnings check-generated check-modernization check-modernization-source check-new-pes check-numerical-regression check-python-extensions check-style clean test test-surfaces

CONDA_ENV ?= chempotpy
PYTHON ?= $(shell \
	if python3 -c "import numpy" >/dev/null 2>&1 && \
		command -v meson >/dev/null 2>&1 && \
		command -v ninja >/dev/null 2>&1 && \
		command -v "$${FC:-gfortran}" >/dev/null 2>&1; then \
		printf '%s' python3; \
	elif command -v conda >/dev/null 2>&1 && \
		conda run -n $(CONDA_ENV) sh -c \
		'python -c "import numpy" >/dev/null 2>&1 && \
		command -v meson >/dev/null 2>&1 && \
		command -v ninja >/dev/null 2>&1 && \
		command -v "$${FC:-gfortran}" >/dev/null 2>&1'; then \
		printf '%s' 'conda run -n $(CONDA_ENV) python'; \
	else \
		printf '%s' python3; \
	fi)
BUILD_FLAGS ?=
BUILD_OUTPUT ?= chempotpy
NUMERICAL_BASELINE ?= tests/data/all_surface_references.json
BASE_REF ?= HEAD^

audit:
	$(PYTHON) scripts/audit_pes.py

check-fortran:
	$(PYTHON) scripts/check_fortran.py

check-fortran-warnings:
	$(PYTHON) scripts/check_fortran.py --include-support --warnings-as-errors

check-binaries:
	$(PYTHON) scripts/check_binaries.py --output-root $(BUILD_OUTPUT)

check-python-extensions:
	$(PYTHON) scripts/check_python_extensions.py --output-root $(BUILD_OUTPUT)

check-conda: build check-binaries check-python-extensions

check-style:
	$(PYTHON) scripts/check_fortran_style.py

check-modernization-source:
	$(PYTHON) scripts/check_modernization.py

check-new-pes:
	$(PYTHON) scripts/check_new_pes.py --base-ref $(BASE_REF)

check-embedded-parameters:
	$(PYTHON) scripts/generate_embedded_network_parameters.py --check
	$(PYTHON) scripts/generate_embedded_pes_parameters.py --check
	$(PYTHON) scripts/generate_hcl_h2o_parameters.py --check
	$(PYTHON) scripts/generate_h3o2m_sgdml_parameters.py --check
	$(PYTHON) scripts/generate_nh3_surfgen_parameters.py --check

check-generated: check-embedded-parameters
	$(PYTHON) scripts/generate_constants.py --check
	$(PYTHON) scripts/generate_provider_catalog.py --check
	$(PYTHON) scripts/generate_readme_catalog.py --check

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
