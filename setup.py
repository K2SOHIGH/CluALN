from setuptools import setup, find_packages

setup(
    name='snakeclualn',
    version='0.1',
    description="""
        cluster sequences with mmseqs and perform MSA per cluster before MSA merging. 
    """,
    url='',
    author='Maxime Millet',
    author_email='maxime.luc.millet@gmail.com',
    license='MIT',
    packages=find_packages(),
    # install_requires=[
    #     "snakemake==7.8.5",
    #     "pyyaml==6.0",
    #     "pandas"        
    # ],    
    include_package_data=True,
    entry_points = {
        'console_scripts': ['clualn = snakeclualn.main:main'],
    },
    zip_safe=False
)