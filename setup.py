from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="QMzyme.Biopython.kdtrees", ["QMzyme/Biopython/kdtrees.c"]
        )
    ]
)
