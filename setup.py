import setuptools

with open("README.md", "r") as fh:
  long_description = fh.read()

setuptools.setup(
  name="SpatialHeterogeneityAnalysis",
  version="1.0.0",
  author="Devin Jun",
  author_email="1402795334@qq.com",
  description="SHA: a Toolkit for Spatial Heterogeneity Analysis of Spatially Resolved Transcriptomes",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/Devin-jun/SHA-SpatialHeterogeneityAnalysis",
  packages=setuptools.find_packages(),
  classifiers=[
  "Programming Language :: Python :: 3",
  "License :: OSI Approved :: MIT License",
  "Operating System :: OS Independent",
  ],
)