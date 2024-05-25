================================================
QMzyme
================================================

.. contents::
   :local:


Introduction
================================================

.. introduction-start

QMzyme is a Python toolkit to facilitate (quantum mechanical) QM-based enzyme 
calculations. The GenerateModel module guides the process of generating 
QM-calculation ready truncated or partitioned enzyme models. The code 
framework can accept any input files MDAnalysis would accept to formulate an
MDAnalysis Universe object. The QMzyme framework works with MDAnalysis modules
to create more dynamic QMzyme objects: QMzymeAtom, QMzymeResidue, QMzymeRegion and
QMzymeModel. The QMzymeModel is an abstraction of a molecular system, such as a real enzyme. 
Its utility comes from the ability to perform calculations on it. The calculation
results will (ideally) be validatable through comparison to experiment, and 
(hopefully) provide new chemical or methodological insights.  

.. introduction-end


.. main_modules-start
Main Modules
===============

GenerateModel
----------------

Module where the user introduces the starting structure and makes decisions on
how the QMzymeModel should be constructed by defining one or more regions. The
GenerateModel class allows a user to define/select regions. Calculation details
are assigned at the QMzymeRegion level.


CalculateModel
----------------

Module where the user defines how a region will be treated in a calculation.
The CalculateModel classes allow the user to declare calculation-level attributes
for a QMzymeRegion instance.



.. main_modules-end

.. data_structures-start
Data Structures
===================

QMzymeModel
-------------

Composed of QMzymeRegion objects.

QMzymeRegion
-------------

Composed of QMzymeAtom and QMzymeResidue objects.

QMzymeResidue
---------------

Composed of QMzymeAtom objects.

QMzymeAtom
------------

The foundational unit that all other QMzyme objects are built from.

.. data_structures-end

.. installation-start
Installation
================================================

.. code-block:: bash

    pip install QMzyme (distribution under development.)

.. installation-end 

.. requirements-start
Requirements
================================================

*  Python >= 3.11
*  NumPy
*  `aqme <https://aqme.readthedocs.io/en/latest/>`_
*  MDAnalysis

*These libraries are installed during the pip installation.*

.. requirements-end


