try:
    from setuptools import setup
    from setuptools import find_packages
except ImportError:
    from distutils.core import setup

setup(
    name='rtPrimerDesign',
    version='0.0.1',
    description='Primer Design for real time PCR',
    url='https://github.com/tomharrop/rtPrimerDesign',
    author='Tom Harrop',
    author_email='twharrop@gmail.com',
    license='GPL-3',
    packages=find_packages(),
    install_requires=[
        'beautifulsoup4>=4.5.1',
        'joblib>=0.9.4',
        'lxml>=3.7.1',
        'progressbar2>=3.12.0',
        'requests>=2.12.4'],
    zip_safe=False
)
