from setuptools import setup

setup(name='wraparound',
      version='1.0.0',
      description=('Wraparound: finalize Kwant systems with '
                   'multiple translational symmetries'),
      url='https://gitlab.kwant-project.org/cwg/wraparound',
      author='Christoph Groth',
      license='BSD 2-clause',
      py_modules=['wraparound'],
      install_requires=['kwant', 'numpy', 'tinyarray'],
      zip_safe=False)
