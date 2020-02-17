import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

long_description = read('README.md')

version_py = os.path.join('crispy_service', 'version.py')
version = read(version_py).strip().split('=')[-1].replace("'", "").strip()

install_requires = [
    "biopython",
    "crispy-models",
    "nearmiss",
    "redis",
    "requests",
]

tests_require = [
    "mockredispy-kblin",
    "pytest",
    "pytest-cov",
]

setup(name='crispy-service',
    version=version,
    install_requires=install_requires,
    tests_require=tests_require,
    author='Kai Blin',
    author_email='packaging@secondarymetabolites.org',
    description='Job runners for CRISPy web service',
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=['crispy_service'],
    url='https://github.com/secondarymetabolites/crispy-service/',
    license='GNU Affero General Public License',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
        'Operating System :: OS Independent',
    ],
    extras_require={
        'testing': tests_require,
    },
    entry_points={
        'console_scripts': [
            'crispy-prepare=crispy_service.prepare:main',
            'crispy-scan=crispy_service.scan:main',
            'crispy-standalone=crispy_service.standalone:main',
        ],
    },
)

