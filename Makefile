APPNAME=$(shell grep -m 1 name pyproject.toml|cut -f2 -d'"')
TARGETS=build clean dependencies deploy editable_install install quickstart test uninstall
VERSION=$(shell grep version pyproject.toml|cut -f2 -d'"')

all:
	@echo "Try one of: ${TARGETS}"

build: clean dependencies
	poetry build

bump-version:
	git-cliff --bumped-version | sed 's/^v//' > version.txt
	python bump-version.py
	git-cliff --bump > docs/changelog.md

clean:
	find . -name '*.pyc' -delete
	find . -type d -name '__pycache__' -exec rm -rf {} +
	rm -rf dist build

dependencies:
	poetry install --without dev --no-root

dependencies_dev:
	poetry install --only dev --no-root

editable_install:
	pip install --editable .

install: build
	pip install dist/*.whl

pre-commit:
	if [ ! -f .git/hooks/pre-commit ]; then pre-commit install; fi
	pre-commit run --all-files

tag:
	git tag v${VERSION}

test: unit-test functional_test_00 functional_test_01 functional_test_02
	@echo "End-to-End tests"


unit-test:
	@echo "Unit Testing"
	pytest --cov=src/gwaspipe/ tests

functional_test_00:
	@echo "Functional test 00"
	gwaspipe \
	  -c examples/config_sumstats_harmonization.yml \
	  -i examples/input_data.tsv.gz \
	  -f regenie \
	  -o results

functional_test_01:
	@echo "Functional test 01"
	gwaspipe \
	  -c examples/config_sumstats_harmonization.yml \
	  -i examples/input_data_01.csv.gz \
	  -f gtex \
	  -s ',' \
	  -o results

functional_test_02:
	@echo "Functional test 02"
	gwaspipe \
	  -c examples/config_sumstats_harmonization_check_ambiguous_snps.yml \
	  -i examples/input_data_02.csv.gz \
	  -f regenie \
	  -s ' ' \
	  -o results

uninstall:
	pip uninstall -y ${APPNAME}
