PYTHON ?= python
CYTHON ?= cython
NOSETESTS ?= nosetests
CTAGS ?= ctags

all: clean inplace test

inplace:
	$(PYTHON) setup.py build_ext -i


test: test-code

test-code: inplace
	$(NOSETESTS) -s -v iced

test-coverage:
	rm -rf coverage .coverage
	$(NOSETESTS) -s -v --with-coverage iced --cover-package iced

clean-ctags:
	rm -f tags

clean: clean-ctags
	$(PYTHON) setup.py clean
	rm -rf dist
	rm -rf build

trailing-spaces:
	find iced -name "*.py" -exec perl -pi -e 's/[ \t]*$$//' {} \;

cython:
	find iced -name "*.pyx" -exec $(CYTHON) {} \;
