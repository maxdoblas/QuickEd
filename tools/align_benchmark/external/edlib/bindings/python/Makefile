default: build

.PHONY:
edlib: $(shell find ../../edlib)
	# create a clean (maybe updated) copy of edlib src
	rm -rf edlib && cp -r ../../edlib .

pyedlib.bycython.cpp: edlib.pyx cedlib.pxd
	python -m pip install cython
	cython --cplus edlib.pyx -o edlib.bycython.cpp

# To build package, README.rst is needed, because it goes into long description of package,
# which is what is visible on PyPI.
# However, to generate README.rst from README-tmpl.rst, built package is needed (for `import edlib` in cog)!
# Therefore, we first build package without README.rst, use it to generate README.rst,
# and then finally build package again but with README.rst.

BUILD_SOURCE_FILES=edlib pyedlib.bycython.cpp setup.py

buildWithoutREADME.rst: ${BUILD_SOURCE_FILES}
	EDLIB_OMIT_README_RST=1 python setup.py build_ext -i

README.rst: buildWithoutREADME.rst README-tmpl.rst
	python -m pip install cogapp
	cog -d -o README.rst README-tmpl.rst	

BUILD_FILES=${BUILD_SOURCE_FILES} README.rst

build: ${BUILD_FILES}
	python setup.py build_ext -i

sdist: edlib pyedlib.bycython.cpp setup.py README.rst MANIFEST.in
	python setup.py sdist

publish: clean sdist
	twine upload dist/*

clean:
	rm -rf edlib dist edlib.egg-info build
	rm -f edlib.c *.bycython.* edlib.*.so
	rm -f README.rst
