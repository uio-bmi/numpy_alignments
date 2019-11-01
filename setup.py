from setuptools import setup

setup(name='numpy_alignments',
      version='0.0.2',
      description='Numpy Alignments',
      url='http://github.com/ivargr/numpy_alignments',
      author='Ivar Grytten',
      author_email='',
      license='MIT',
      packages=['numpy_alignments'],
      zip_safe=False,
      install_requires=['numpy', 'tqdm', 'plotly'],
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
twine upload dist/numpy_alignments-X.tar.gz
"""
