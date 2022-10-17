from setuptools import setup


setup(
    name='tem-simulator-scripts',
    use_scm_version={'write_to': 'temsimscripts/__init__.py'},  # True,
    setup_requires=['setuptools_scm'],
    python_requires='>=3.7.0, <3.10',
    packages=[
        'temsimscripts'
    ],
    url='https://github.com/MPI-Dortmund/tem-simulator-scripts',
    license='Mozilla Public License Version 2.0',
    author='Thorsten Wagner, Markus Stabrin, Gavin Rice',
    install_requires=[
        "numpy >= 1.20.0",
        "mrcfile",
        "scipy",
        "tqdm",
        "pandas",
        'biotite',
        "scikit-image"
    ],
    author_email='thorsten.wagner@mpi-dortmund.mpg.de',
    description='Support scripts for the TEM Simulator (V1.3)',
    long_description='',
    long_description_content_type='text/markdown',
    entry_points={
        'console_scripts': [
            'tsimscripts_gen_input.py = temsimscripts.generate_input_file:_main_',
            'tsimscripts_gen_filaments.py = temsimscripts.generate_filaments:_main_',
            'tsimscripts_gen_raw_tilt.py = temsimscripts.generate_raw_tilt:_main_',
            'tsimscripts_gen_trans_file.py = temsimscripts.generate_trans_file:_main_',
            'tsimscripts_gen_map.py = temsimscripts.generate_map:_main_',
            'tsimscripts_gen_coords.py = temsimscripts.generate_coordinates:_main_',
            'tsimscripts_extract.py = temsimscripts.extract_subvolumes:_main_',
        ]},
    scripts = [
        'pipeline/tsimscripts_pipe.sh'
    ]
)
