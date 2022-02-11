import setuptools

setuptools.setup(
    name="chemicaldiagram",
    version="0.0.6",
    author="Q Ai",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=["chemicaldiagram"],
    package_data={'chemicaldiagram': ['*.json']},
    include_package_data=True,
    python_requires=">=3.7",
)
