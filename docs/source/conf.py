# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'DOMINO'
copyright = '2025, Pan Jia'
author = 'Pan Jia'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [ 
    'myst_parser',
    # 'myst_nb',
    "sphinxcontrib.mermaid",
    # 'sphinx-prompt',
    "sphinx_copybutton",
]


source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}


myst_enable_extensions = [
    "tasklist",
    "deflist",
    "dollarmath",
    "colon_fence",
    "linkify"
]

templates_path = ['_templates']
exclude_patterns = []


html_theme_options = {
    'navigation_depth': 3,        
    'collapse_navigation': False,  
    'sticky_navigation': True,    
    'includehidden': True,        
    'titles_only': False          
}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']

html_css_files = [
    'custom.css',
]

html_context = {
    'display_github': True,  
    'github_user': 'programmerjia',  
    'github_repo': 'DOMINO-Tutorials',      
    'github_version': 'main',          
    'conf_py_path': '/docs/source/',   
    'edit_page_url_template': 'https://github.com/{github_user}/{github_repo}/edit/{github_version}/docs/source/{path}',
}
