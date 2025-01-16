#from setuptools import setup


#setup(
#    use_scm_version={
#        "write_to": "xapres/_version.py",
#        "write_to_template": '__version__ = "{version}"',
#        "tag_regex": r"^(?P<prefix>v)?(?P<version>[^\+]+)(?P<suffix>.*)?$",
#    }
#)

from setuptools import setup, find_packages

setup(
    name='xapres',
    version='0.14.6',
    author='ldeo_glaciology',
    author_email='j.kingslake@columbia.edu',
    description='A package for processing data from the Autonomous phase-sensitive Radio-Echo Sounder (ApRES) using xarray.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/ldeo-glaciology/xapres',
    license='MIT',
    license_files=('LICENSE.txt',),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Topic :: Scientific/Engineering',
        'Intended Audience :: Science/Research',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'License :: OSI Approved :: MIT License',
    ],
    install_requires=[
        'numpy',
        'requests',
        'xarray',
        'gcsfs',
        'fsspec',
        'pandas',
        'tqdm',
    ],
    python_requires='>=3.6',
    packages=find_packages(exclude=['docs', 'data']),
    zip_safe=False,
)
