from setuptools import setup

def readme():
	with open('README.md') as f:
		return f.read()

setup(name='geode',
      version='0.1',
      description='A Python implementation of the R package GeoDE',
      long_description=readme(),
      url='https://github.com/wangz10/geode',
      author='Zichen Wang',
      author_email='zichen.wang@mssm.edu',
      license='GPL',
      packages=['geode'],
      install_requires=['numpy', 'scipy', 'scikit-learn'],
      test_suite="tests",
      include_package_data=True,
      zip_safe=False)
