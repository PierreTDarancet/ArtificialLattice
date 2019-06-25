from setuptools import setup, find_packages

setup(
# -----------
name='TopoQuest',
version = 'v0.1',
url = 'https://www.anl.gov/cnm',
description = 'TopoQuest',
# -----------
author = 'ANL_NST_TMG',
author_email = 'ssrinivasan@anl.gov',
# -----------
install_require = ['numpy','scipy'],
packages = find_packages('src'),
package_dir={'':'src'},
# -----------
# -----------
#test_suite='tests',
)
