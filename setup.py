from setuptools import setup

setup(name='numpy_alignments',
      version='0.0.1',
      description='Numpy Alignments',
      url='http://github.com/ivargr/simple_read_mutator',
      author='Ivar Grytten',
      author_email='',
      license='MIT',
      packages=['numpy_alignments'],
      zip_safe=False,
      install_requires=['numpy'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points = {
            'console_scripts': ['numpy_alignments=numpy_alignments.command_line_interface:main'],
      })

"""
To update package:
#Update version number manually in this file

sudo python3 setup.py sdist
sudo python3 setup.py bdist_wheel
twine upload dist/simple_read_mutatot-X.tar.gz
"""
