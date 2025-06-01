Autodoc Thing
https://github.com/RadioAstronomySoftwareGroup/pyuvsim/blob/main/docs/conf.py

Okay, so, as far as I remember, this is the rundown:
1. When doing sphinx-quickstart, hit 'n' instead of 'y'.
2. After running through the quickstart, you should still be in the docs folder and be able to `nano conf.py` . When editing conf.py, immediatelly add these three lines: 
```
import sys
import os
sys.path.insert(0, os.path.abspath('../src/'))
```
With the '../src/' depending on your package directory structure. src should contain the directory of your package. So, in my case, the full directory structure looks like this:
```
reRoute_Dynamics_Repo/
├─ docs/
|  ├─ conf.py
|  ├─ index.rst
|  ├─ modules.rst
|  ├─ etc...
├─ src/
|  ├─ reRoute_Dynamics/
|  |  ├─ __init__.py
|  |  ├─ module1.py
|  |  ├─ etc...
├─ etc/
```
The next edit to make to conf.py is to fill extensions  like this:
```
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx_rtd_theme', 'sphinx.ext.viewcode']
```
and then change html_theme to be html_theme = `sphinx_rtd_theme` .
3. Next, run the sphinx-apidoc -o docs/ src/package_name from right outside the docs directory.
4. Now that all that's done and dusted. we need to double check modules.rst and index.rst.
modules.rst should have 1 line in the toctree after the maxdepth, and that should be your package name. In my case, it looks like this:

  GNU nano 6.2                                                                                        modules.rst                                                                                                  
reRoute_Dynamics_Core
=====================

.. toctree::
   :maxdepth: 4

   reRoute_Dynamics
And then index.rst should be directing to your modules.rst file, like so:
^^ headers and stuff up here ^^


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules
5. Now, running make clean html while in the docs/ directory should yeild a workable result. It may spit out errors at you if your docstrings are giving it hell, but that's more easily fixable.

Quick addendum!
When trying to get ReadTheDocs to work for viewing online once it's succesfully rendered the HTML, you need to do the following steps:

    If you're sure you've got all the packages you want and need in your current python environment, run pip list --format=freeze > requirements.txt  to create a list of requirements. Some of the included packages may give you trouble, so you may need to iteratively remove lines or adjust the requirements.txt as nessecary to get ReadTheDocs to play nice
    add a .readthedocs.yaml file to your base repository directory, as per what the walkthrough requests. This is going to have adjustments depending on your os and python version, so make sure you double check those. Uncomment the 

python:
  install:
   - requirements: docs/requirements.txt

and double check that the tabular alignment is correct. Also, double check that the sphinx build configuration leads to docs/conf.py assuming you used the last set of instructions for sphinx
3. Push your repo! If you've connected it to the ReadTheDocs website already (which is as simple as logging in and hitting the big green button, really) it will automatically build the pages from the previously constructed HTML. Be prepared to push like 5 times for small dumb adjustments - the error reporting for ReadTheDocs is pretty good, so it should be easy to figure out what went wrong