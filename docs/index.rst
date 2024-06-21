.. |logo| image:: Images/QMzyme_logo.png
   :width: 300

|logo|

:QMzyme version: |release|
..
   :Last updated: |today|

.. note::
   QMzyme is currently under-development. 
   Please note the user interface may change! 
   The first stable API version will be released 
   as version 1.0.0 (QMzyme==1.0.0) on PyPi.*

=============================================
Introduction
=============================================


**QMzyme** is a Python toolkit to facilitate (quantum mechanical) QM-based enzyme 
calculations. The :class:`~QMzyme.GenerateModel.GenerateModel` class guides the process of generating 
QM-calculation ready truncated or partitioned models. The code 
framework can accept any input files `MDAnalysis <https://userguide.mdanalysis.org/stable/index.html>`_ would accept to formulate an
MDAnalysis `Universe <https://userguide.mdanalysis.org/stable/universe.html>`_ object. The QMzyme framework works with MDAnalysis modules
to create more dynamic QMzyme objects: :class:`~QMzyme.QMzymeAtom.QMzymeAtom`, 
:class:`~QMzyme.QMzymeRegion.QMzymeResidue`, :class:`~QMzyme.QMzymeRegion.QMzymeRegion`, and
:class:`~QMzyme.QMzymeModel.QMzymeModel`. QMzymeModel is an abstraction of a molecular system, such as a real enzyme. 
Its utility comes from the ability to perform calculations on it. The calculation
results will (ideally) be validatable through comparison to experiment, and 
(hopefully) provide new chemical or methodological insights. 

.. toctree::
   :maxdepth: 2
   :caption: Table of contents

   Quickstart/installation
   Examples/index
   API/index
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
