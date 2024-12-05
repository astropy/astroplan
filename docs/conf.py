# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# Astropy documentation build configuration file.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# Note that not all possible configuration values are present in this file.
#
# All configuration values have a default. Some values are defined in
# the global Astropy configuration which is loaded here before anything else.
# See astropy.sphinx.conf for which values are set there.

from configparser import ConfigParser
import sys
import datetime
from importlib import metadata

try:
    from sphinx_astropy.conf.v2 import *  # noqa
except ImportError:
    print('ERROR: the documentation requires the sphinx-astropy package to '
          'be installed')
    sys.exit(1)

# Get configuration information from setup.cfg
conf = ConfigParser()
conf.read([os.path.join(os.path.dirname(__file__), '..', 'setup.cfg')])
setup_cfg = dict(conf.items('metadata'))


# -- General configuration ----------------------------------------------------
# By default, highlight as Python 3.
highlight_language = 'python3'

# If your documentation needs a minimal Sphinx version, state it here.
# needs_sphinx = '1.7'

# Extend astropy intersphinx_mapping with packages we use here
intersphinx_mapping['astroquery'] = ('http://astroquery.readthedocs.io/en/latest/', None)

# Exclude astropy intersphinx_mapping for unused packages
del intersphinx_mapping['h5py']  # noqa

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns.append('_templates')  # noqa
# Exclude template PSF block specification documentation
exclude_patterns.append('psf_spec/*')  # noqa

plot_formats = ['png', 'hires.png', 'pdf', 'svg']

# This is added to the end of RST files - a good place to put
# substitutions to be used globally.
rst_epilog = """
.. _Astropy: https://www.astropy.org/
"""

# -- Project information ------------------------------------------------------

# This does not *have* to match the package name, but typically does
project = setup_cfg['name']
author = setup_cfg['author']
copyright = '{0}, {1}'.format(
    datetime.datetime.now().year, setup_cfg['author'])

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

# The full version, including alpha/beta/rc tags.
release = metadata.version(project)
# The short X.Y version.
version = ".".join(release.split(".")[:2])

# Only include dev docs in dev version.
dev = "dev" in release


# -- Options for HTML output ---------------------------------------------------

# A NOTE ON HTML THEMES
# The global astropy configuration uses a custom theme, 'bootstrap-astropy',
# which is installed along with astropy. A different theme can be used or
# the options for this theme can be modified by overriding some of the
# variables set in the global configuration. The variables set in the
# global configuration are listed below, commented out.

# Add any paths that contain custom themes here, relative to this directory.
# To use a different custom theme, add the directory containing the theme.
#html_theme_path = []

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes. To override the custom theme, set this to the
# name of a builtin theme or the name of a custom theme in html_theme_path.
# html_theme = None

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#html_favicon = ''

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
#html_last_updated_fmt = ''

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = '{0} v{1}'.format(project, release)

# Output file base name for HTML help builder.
htmlhelp_basename = project + 'doc'

# Prefixes that are ignored for sorting the Python module index
modindex_common_prefix = ["astroplan."]

# -- Options for LaTeX output --------------------------------------------------

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [('index', project + '.tex', project + u' Documentation',
                    author, 'manual')]


# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [('index', project.lower(), project + u' Documentation',
              [author], 1)]


# -- Options for the edit_on_github extension ----------------------------------

if setup_cfg.get('edit_on_github').lower() == 'true':

    extensions += ['sphinx_astropy.ext.edit_on_github']

    edit_on_github_project = setup_cfg['github_project']
    edit_on_github_branch = "main"

    edit_on_github_source_root = ""
    edit_on_github_doc_root = "docs"

# Make appropriate substitutions to mock internet querying methods
# within the tests.
# Currently this is not needed because of the content of the tests, but we leave
# it here in case it's needed again in the future.  BUT beware of:
# https://github.com/astropy/astroplan/issues/96
#from astroplan.utils import _mock_remote_data
#_mock_remote_data()

# Add additional Sphinx extensions:
extensions += [
    'matplotlib.sphinxext.plot_directive',
    'sphinx.ext.graphviz'
]

# -- Resolving issue number to links in changelog -----------------------------
github_issues_url = 'https://github.com/astropy/astroplan/issues/'

# -- Options for linkcheck output -------------------------------------------
linkcheck_retry = 5
linkcheck_ignore = [
    r'https://github\.com/astropy/astroplan/(?:issues|pull)/\d+',
]
linkcheck_timeout = 180
linkcheck_anchors = False

# -- Turn on nitpicky mode for sphinx (to warn about references not found) ----
nitpicky = True

html_copy_source = False

html_theme_options.update(  # noqa: F405
    {
        "github_url": "https://github.com/astropy/astroplan",
        "external_links": [
            {"name": "astropy docs", "url": "https://docs.astropy.org/en/stable/"},
        ],
        "use_edit_page_button": True,
    }
)

html_context = {
    "default_mode": "light",
    "to_be_indexed": ["stable", "latest"],
    "is_development": dev,
    "github_user": "astropy",
    "github_repo": "astroplan",
    "github_version": "main",
    "doc_path": "docs",
}

# Add any extra paths that contain custom files (such as robots.txt or
# .htaccess) here, relative to this directory. These files are copied
# directly to the root of the documentation.
html_extra_path = ["robots.txt"]

#
# Some warnings are impossible to suppress, and you can list specific references
# that should be ignored in a nitpick-exceptions file which should be inside
# the docs/ directory. The format of the file should be:
#
# <type> <class>
#
# for example:
#
# py:class astropy.io.votable.tree.Element
# py:class astropy.io.votable.tree.SimpleElement
# py:class astropy.io.votable.tree.SimpleElementWithContent
#
# Uncomment the following lines to enable the exceptions:
#
# for line in open('nitpick-exceptions'):
#     if line.strip() == "" or line.startswith("#"):
#         continue
#     dtype, target = line.split(None, 1)
#     target = target.strip()
#     nitpick_ignore.append((dtype, str(target)))
