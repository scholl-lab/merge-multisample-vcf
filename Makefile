# =============================================================================
# Makefile — lint, format, and typecheck for merge-multisample-vcf
# =============================================================================
# Install dev tools:  pip install ruff mypy snakefmt shellcheck-py
#
# Quick reference:
#   make lint        — check everything (no changes)
#   make format      — auto-fix formatting
#   make fix         — auto-fix lint issues + format
#   make typecheck   — run mypy
#   make test        — run all tests
#   make clean       — remove caches
# =============================================================================

SHELL := /bin/bash
.DEFAULT_GOAL := help

# --- File sets ---------------------------------------------------------------
PY_FILES   := scripts/generate_config.py workflow/rules/helpers.py
SMK_FILES  := workflow/rules/*.smk workflow/Snakefile
SH_FILES   := scripts/run_snakemake.sh

# =============================================================================
# Composite targets
# =============================================================================

.PHONY: lint
lint: lint-py lint-smk lint-sh typecheck  ## Run all checks (no modifications)

.PHONY: format
format: format-py format-smk format-sh  ## Auto-format all files

.PHONY: fix
fix: fix-py format-smk format-sh  ## Auto-fix lint issues + format

.PHONY: check
check: lint  ## Alias for lint

# =============================================================================
# Python  (ruff + mypy)
# =============================================================================

.PHONY: lint-py
lint-py:  ## Check Python: ruff lint + format check
	ruff check $(PY_FILES)
	ruff format --check $(PY_FILES)

.PHONY: format-py
format-py:  ## Format Python files
	ruff format $(PY_FILES)

.PHONY: fix-py
fix-py:  ## Auto-fix Python lint issues + format
	ruff check --fix $(PY_FILES)
	ruff format $(PY_FILES)

.PHONY: typecheck
typecheck:  ## Run mypy type checking
	mypy $(PY_FILES)

# =============================================================================
# Snakemake  (snakefmt)
# =============================================================================

.PHONY: lint-smk
lint-smk:  ## Check Snakemake formatting
	snakefmt --check $(SMK_FILES)

.PHONY: format-smk
format-smk:  ## Format Snakemake files
	snakefmt $(SMK_FILES)

# =============================================================================
# Shell  (shellcheck)
# =============================================================================

.PHONY: lint-sh
lint-sh:  ## Check shell scripts (shellcheck, warning level)
	shellcheck -S warning $(SH_FILES)

.PHONY: format-sh
format-sh:  ## Fix shell script line endings (LF)
	@for f in $(SH_FILES); do \
		sed -i 's/\r$$//' "$$f" && echo "  LF: $$f"; \
	done

# =============================================================================
# Tests  (pytest)
# =============================================================================

.PHONY: test
test:  ## Run all tests
	pytest

.PHONY: test-unit
test-unit:  ## Run unit tests (skip dry-run tests)
	pytest -m "not dryrun"

.PHONY: test-dryrun
test-dryrun:  ## Run dry-run tests only
	pytest -m dryrun

.PHONY: test-cov
test-cov:  ## Run tests with coverage report
	pytest --cov=scripts --cov=workflow/rules --cov-report=term-missing

# =============================================================================
# Utilities
# =============================================================================

.PHONY: install-dev
install-dev:  ## Install dev linting tools
	pip install ruff mypy snakefmt shellcheck-py pytest pytest-cov

.PHONY: clean
clean:  ## Remove tool caches
	rm -rf .mypy_cache .ruff_cache __pycache__ scripts/__pycache__ workflow/rules/__pycache__

.PHONY: help
help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'
