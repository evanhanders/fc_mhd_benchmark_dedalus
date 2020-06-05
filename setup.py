from setuptools import setup

setup(name='fce_dedalus',
      version='0.0.1',
      description='Tools and scripts for studying fully compressible systems in dedalus',
      url='https://github.com/evanhanders/fce_dedalus',
      author='Evan Anders',
      author_email='evan.anders@colorado.edu',
      license='GPL-3.0',
      packages=['fce_dedalus', 'fce_dedalus/logic', 'fce_dedalus/atmospheres', 'fce_dedalus/domains'],
      zip_safe=False)
