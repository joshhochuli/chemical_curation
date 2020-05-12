from setuptools import setup, find_packages

requirements = [
    "click",
    "pandas",
    #"rdkit",
    # "molvs",
    "sphinx"
]

setup(
    name="chemical_curation",
    version="0.0.1",
    description="Library providing functions for curating chemical data according to the Trust But Verify process.",

    packages=find_packages('src'),
    package_dir={'':'src'},
    
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "tbv = chemical_curation.tbv:cli"
        ]
    }
)
