## Installing the package
```
python3 setup.py install --user
```

Or install using pip (in the main directory)
```
pip3 install -e . --user
```

Check the installation
```
pip3 list | grep -F beams2D 
```


## Install Using Virtual Environment (recommended for modifications)
Go to the main beams2D dir (use --system-site-packages to re-use numpy, etc)
```
python3 -m venv venv-beams2d --system-site-packages
source venv-beams2d/bin/activate
pip3 install -e . --user
```

For working with `JupyterLab` use
```
pip3 install ipykernel
python3 -m ipykernel install --user --name=venv-beams2d
```

Deactive the env (if not required anymore)
```
deactivate
```

### Quick start

```shell
cd tutorials/simpleBeam
python3 -m jupyterlab simpleBeam.ipynb
```