from setuptools import setup, find_packages()

setup(
    name="ChemicalCuration",
    version="0.0.1",
    description="Library providing functions for curating chemical data according to the Trust But Verify process."

    packages=find_packages(),
    include_package_data=True,
    install_requires=["click",
                      "sphinx",
                      "rdkit",
                      "molvs",
                      "pandas",
                      "unittest"
                      ],

    entry_points="""
    [console_scripts]
    chemical_curation=tbv:main
    """ 
)
