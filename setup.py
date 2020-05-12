from setuptools import setup

requirements = [
    "click",
    "pandas",
    #"rdkit",
    # "molvs",
    "sphinx"
]

setup(
    name="ChemicalCuration",
    version="0.0.1",
    description="Library providing functions for curating chemical data according to the Trust But Verify process.",

    packages=['chemical_curation'],
    
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "tbv = chemical_curation.tbv:cli"
        ]
    }
)
