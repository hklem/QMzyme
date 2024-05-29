Installation
================================================

Requirements
---------------
*These libraries are installed during pip installation.*

*  Python >= 3.11 (older versions might work too, but testing is done with 3.11)
*  NumPy
*  `aqme <https://aqme.readthedocs.io/en/latest/>`_
*  `MDAnalysis <https://www.mdanalysis.org>`_

Source Code
-------------
The source code is available from https://github.com/hklem/QMzyme.

.. code-block:: bash
    
    git clone https://github.com/hklem/QMzyme/QMzyme.git
    cd QMzyme

    # The next 2 commands are optional, but recommended
    conda create -n qmzyme "python==3.11"
    conda activate qmzyme

    pip install -e . #-e signifies developer install 


From PyPi
-----------

.. code-block:: bash

    pip install QMzyme #NOTE the PyPi distribution still under development. User interface may break
