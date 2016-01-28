try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Primer Design for real time PCR',
    'author': 'Tom Harrop',
    'url': 'https://github.com/tomharrop/rtPrimerDesign',
    'author_email': '-----',
    'version': '0.0.1',
    'packages': ['rtPrimerDesign'],
    'install_requires': ['csv', 'urllib', 'bs4', 'requests', 'progressbar']
    'name': 'rtPrimerDesign'
}

setup(**config)
