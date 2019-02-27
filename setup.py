from setuptools import setup

requirements = [
    'click',
    'biopython',
    'numpy',
    'pandas']

setup(
    name='primacy',
    version='1.0.0',
    description="Multiplex PCR primer optimization",
    author="Tara Furstenau",
    author_email='tara.furstenau@nau.edu',
    url='https://github.com/FofanovLab/Primacy',
    packages=['primacy'],
    entry_points={
        'console_scripts': [
            'primacy_cli=primacy.main:main',
            'primacy=primacy.main:gui'
        ]
    },
    install_requires=requirements,
    keywords='bio primer',
    classifiers=[
        'Programming Language :: Python :: 3.6',
    ]
)