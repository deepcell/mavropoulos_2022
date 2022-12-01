# Minimal makefile for Sphinx documentation

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = pipeline
SOURCEDIR     = docs
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Delete current documentation for a fresh build, build API autodoc stub files, and then build html docs.
all:
	rm -rf ${BUILDDIR} ${SOURCEDIR}/dev_manual/api
	make api
	make html

# Build the API autodoc stub files needed for the developer manual. Options are as follows:
#   --no-toc: don't create a table of contents file; the index.rst file in the dev_manual suffices
#   --output-dir: directory to place all output
#   --separate: put documentation for each module on its own page
api:
	sphinx-apidoc --no-toc --separate --output-dir ${SOURCEDIR}/dev_manual/api ${SPHINXPROJ}

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
