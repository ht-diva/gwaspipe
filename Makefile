TARGETS=build clean dependencies

all:
	@echo "Try one of: ${TARGETS}"

clean:
	find . -name '*.pyc' -delete
	find . -type d -name '__pycache__' -exec rm -rf {} +
	rm -rf dist build

dependencies:
	conda env update --file environment.yml
