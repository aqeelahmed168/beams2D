## Building Documentation using Sphinx

It is not as easy as it looks...


Some insights from the fluidSim packagae
https://github.com/fluiddyn/fluidsim

After the quickstart with Sphinx, try setting up
following this post
https://stackoverflow.com/a/62613202/12615865


Still need some customizations in the `index.rst` file, combine that with the
input in the `__init__.py` files to have specific control.


Look into the `config.py` and the `index.rst` file for further reference.


If you want to use the sphinx build using pip3

```shell
cd docs
python3 -m sphinx.cmd.build -b html . _build
```

To include the math (also need to install the basictex and other packages)
Note: Tested with Sphinx v4.5.0

```shell
cd docs
python3 -m sphinx.cmd.build -b html -D imgmath_latex=/Library/TeX/texbin/latex -D imgmath_dvisvgm=/Library/TeX/texbin/dvisvgm . _build/html
```