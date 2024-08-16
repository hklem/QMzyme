.. |logo| image:: Images/QMzyme_logo.png
   :width: 300

|logo|

:QMzyme version: |release|
..
   :Last updated: |today|

.. note::
   QMzyme is in a developmental stage. 
   Please note the user interface may change! 
   The first API version with guaranteed stability 
   will be released as version 1.0.0 (QMzyme==1.0.0) on PyPi (pip).

=============================================
Introduction
=============================================


**QMzyme** is a Python toolkit to facilitate (quantum mechanical) QM-based enzyme 
calculations. The :class:`~QMzyme.GenerateModel.GenerateModel` class guides the process of generating 
QM-calculation ready truncated or partitioned structures. The code 
framework can accept all input files that `MDAnalysis <https://userguide.mdanalysis.org/stable/index.html>`_ accepts to create an
MDAnalysis `Universe <https://userguide.mdanalysis.org/stable/universe.html>`_ object. The QMzyme framework works with MDAnalysis modules
to create more dynamic QMzyme data structures: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`, 
:class:`~QMzyme.QMzymeRegion.QMzymeResidue`, :class:`~QMzyme.QMzymeRegion.QMzymeRegion`, and
:class:`~QMzyme.QMzymeModel.QMzymeModel`. QMzymeModel is an abstraction of a molecular system, such as a real enzyme, that 
comprises at least one QMzymeRegion. Its utility comes from the ability to perform calculations on it. 
Calculation methods (think Hamiltonians) are assigned at the QMzymeRegion level. QMzyme, as its namesake suggests,
takes a QM-centric perspective. Therefore, the CalculateModel and Writers modules are designed for compatibility with QM-focused 
software, rather than MD-focused software, even though many MD software support QM program interfacing. The calculation
results will (ideally) be validatable through comparison to experiment, and (hopefully) provide new chemical or methodological insights. 

If you have ideas or suggestions on how to improve QMzyme please do not hesitate to engage on the `QMzyme Ideas GitHub Project space <https://github.com/users/hklem/projects/11/views/1>`_, or contribute directly by `forking the repository and submitting a pull request <https://qmzyme.readthedocs.io/en/latest/Contributing/general.html>`_. Questions can be directed via email to 
heidi.klem{AT}nist{dot}gov.


.. toctree::
   :maxdepth: 2

   Quickstart/installation

.. toctree::
   :maxdepth: 2
   :caption: Examples/Tutorials

   Examples/index

.. toctree::
   :maxdepth: 2
   :caption: Code Documentation

   API/index

.. toctree::
   :maxdepth: 2
   :caption: Contributing to QMzyme!

   Contributing/index


.. Hide the contents from the front page because they are already in
.. the side bar in the Alabaster sphinx style; requires Alabaster
.. config sidebar_includehidden=True (default)

.. Contents
.. ========
.. toctree::
   :maxdepth: 4
   :caption: Documentation
   :numbered:		
   :hidden:
   

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
