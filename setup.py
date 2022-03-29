import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PyDCF",
    version="1.0.1",
    author="Jinsoo Park",
    author_email="jp7dec23@gmail.com",
    description="The Davis Chandrasekhar Fermi method written in python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/parkji30/PyDCF",
    project_urls={
        "Bug Tracker": "https://github.com/parkji30/PyDCF/issues/",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)