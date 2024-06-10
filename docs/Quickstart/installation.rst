Installation
================================================

Requirements
---------------
*These libraries are installed during pip installation.*

*  Python >= 3.11 (older versions might work too, but testing is done with 3.11)
*  NumPy
*  `aqme <https://aqme.readthedocs.io/en/latest/>`_
*  `MDAnalysis <https://www.mdanalysis.org>`_

Source Code (recommended)
--------------------------
The source code is available from https://github.com/hklem/QMzyme.

.. code-block:: bash
    
    git clone https://github.com/hklem/QMzyme.git
    cd QMzyme

    # The next 2 commands are optional, but recommended
    conda create -n qmzyme "python==3.11"
    conda activate qmzyme

    pip install -e . # the -e signifies developer install 

    # If you want to run tests you can install the required dependecies by executing:
    pip install 'QMzyme[test]'


From PyPi
-----------
The data files used in the Tutorials and Cookbook documentation 
section may not be included if you do the pip install. You can 
download the data files as you wish directly from the 
`GitHub repo <https://github.com/hklem/QMzyme/tree/main/QMzyme/data>`_.

.. code-block:: bash

    pip install QMzyme #NOTE the PyPi distribution still under development. User interface may break
