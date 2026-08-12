from setuptools import setup, find_packages

setup(
    name="cyclospora_pyeuk",
    version="2.0.0",
    description="Modern CDC Cyclospora cayetanensis MLST typing, Eukaryotyping ensemble distance engine, and outbreak cluster finder",
    author="BRC-Analytics / CDC Modernization Team",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "pandas>=1.3.0",
        "scikit-learn>=1.0.0",
        "numba>=0.53.0",
    ],
    entry_points={
        "console_scripts": [
            "cyclospora-typing=cyclospora_pyeuk.cli:main",
        ],
    },
)
