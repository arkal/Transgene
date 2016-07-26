from __future__ import print_function

import sys

from setuptools import setup
from setuptools.command.test import test as TestCommand
from version import version

# Set up a test class
# noinspection PyAttributeOutsideInit
class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        # Sanitize command line arguments to avoid confusing Toil code attempting to parse them
        sys.argv[1:] = []
        err_number = pytest.main(self.pytest_args)
        sys.exit(err_number)


setup(name='transgene',
      version=version,
      description='Translation of genomic events to proteomic space',
      url='http://github.com/arkal/transgene',
      author='Arjun Arkal Rao',
      author_email='aarao@ucsc.edu',
      license='Apache',
      install_requires=[
          'pysam>=0.9.1',
          'swalign>=0.3.3',
      ],
      tests_require=[
          'pytest==2.8.3'],
      test_suite='transgene',
      entry_points={
          'console_scripts': [
              'transgene = transgene:run_transgene']},
      cmdclass={'test': PyTest},
      zip_safe=False)
