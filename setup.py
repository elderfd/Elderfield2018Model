#!/usr/bin/env python

from distutils.core import setup

setup(
    name = "Elderfield 2017 Fungicide Model",
    version = "1.0",
    description = "Model of fungicide resistance evolution in either wheat septoria leaf blotch or grapevine powdery mildew",
    author = "James Elderfield",
    author_email = "jad.elderfield@gmail.com",
    license = "MIT",
    url = "",  # TODO: Point to bioarxiv version of paper
    packages = ["fungicide_model"],
    requires = ["numpy", "pandas", "scipy"]
)