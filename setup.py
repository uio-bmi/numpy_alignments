from setuptools import setup

setup(name='numpy_alignments',
      version='0.0.9',
      description='Numpy Alignments',
      url='http://github.com/ivargr/numpy_alignments',
      author='Ivar Grytten',
      author_email='',
      license='MIT',
      packages=['numpy_alignments'],
      zip_safe=False,
      install_requires=['numpy', 'tqdm', 'plotly', 'graph_read_simulator', 'bionumpy>=0.2.20'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points = {
            'console_scripts': ['numpy_alignments=numpy_alignments.command_line_interface:main'],
      })

"""
To update package:
#Update version number manually in this file

rm -rf dist
python3 setup.py sdist
python3 setup.py bdist_wheel
twine upload dist/*.tar.gz
"""
