from setuptools import setup, find_packages

setup(
    name='snakeclualn',
    version='0.1',
    description="""
        cluster sequences with mmseqs2 and more. 
    """,
    url='',
    author='Maxime Millet',
    author_email='maxime.luc.millet@gmail.com',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    entry_points = {
        'console_scripts': ['clualn = snakeclualn.main:main']        
    },
    zip_safe=False
)