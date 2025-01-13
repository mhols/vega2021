Examples
=============

Installation/Usage:
*******************
Donwload the sources from git (). Create a virtual environment for
nextra. 

.. code-block:: console
    
    $ git clone https://github.nextra

    $ cd nextra
    $ python -m venv venv-nextra
    $ source venv-nexra/bin/activate

Next install the required python packages into the virtaul environment that you 
just have activated.

.. code-block:: console
    
    $ pip install -r requirements.txt

In order to allow to acces the module you sould set the environemnt variable 
PYTHONPATH to the src folder inside the nextra package

.. code-block:: console
    
    $ export PYTHONPATH=/your/path/to/nextra/src

The reference files and line catalogs that are delivered with the package need 
to be downloaded independently.


You may verify that everithing is well installed by running a couple of test cases

.. code-block:: console
    
    $ python ./src/tests/test_vega.py


A sample session
****************

Open an interactive python shell, like ipython. 

.. code-block:: python

    In[1] from nextra import *
 