# -*- coding: utf-8 -*-

from distutils.core import setup
import shutil
shutil.copy('README.md', 'channelutil/README.md')

setup(name='channelutil',
      version='0.6',
      description='Python package to convert between asymptotic wavenumbers and energies for multichannel, inelastic collisions.',
      author="Peter Bingham",
      author_email="petersbingham@hotmail.co.uk",
      packages=['channelutil'],
      package_data={'channelutil': ['tests/*', 'README.md']}
     )
