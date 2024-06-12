from setuptools import setup, find_packages

setup(
    name='gr_particle_path_tool',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'sympy',
        'matplotlib',
        'jupyter',
        'ipywidgets',
        'scipy'
    ],
    entry_points={
        'console_scripts': [
            # Add any console scripts here if needed
        ],
    },
    include_package_data=True,
    description='Open-source research tool for visualizing particle paths around a black hole',
    author='Your Name',
    author_email='your_email@example.com',
    url='https://github.com/your_username/GR_Particle_Path_Research_Tool',
    license='MIT',
)

