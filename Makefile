TARGETS=build clean dependencies

all:
	@echo "Try one of: ${TARGETS}"

clean:
	find . -name '*.pyc' -delete
	find . -type d -name '__pycache__' -exec rm -rf {} +
	rm -rf dist build

dependencies:
	conda env update --file environment.yml

dev-dependencies: dependencies
	conda env update --file environment_dev.yml

pre-commit:
	if [ ! -f .git/hooks/pre-commit ]; then pre-commit install; fi
	pre-commit run --all-files

test:
	@echo "Testing"
	python -m unittest discover -s tests
