APPNAME=$(shell grep -m 1 name pyproject.toml|cut -f2 -d'"')
TARGETS=build clean dependencies deploy editable_install install test uninstall
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

test:
	@echo "Testing"
	pytest --cov=src/gwaspipe/ tests
	python -m unittest discover -s tests

uninstall:
	pip uninstall -y ${APPNAME}
