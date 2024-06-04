.. _calculation_methods:

Adding a Calculation Type
============================

In QMzyme, :class:`~QMzyme.CalculateModel` contains all methods used to 
define a calculation type, its details, and assign it to a region. When a calculation
class (Ex., :class:`~QMzyme.CalculateModel.QM_Method`) is assigned to a region using 
its parent class (:class:`~QMzyme.CalculateModel.CalculationBase`) method, `.assign_to_region()` 
and the calculation type is then recorded in :class:`~QMzyme.CalculateModel.CalculateModel`.
This is required to enable :class:`~QMzyme.GenerateModel.GenerateModel` to automate writing of 
calculation input files, because :class:`~QMzyme.Writers.Writer` will be able to read
the CalculateModel data and decide what writer class to call (assuming there is 
a writer class that corresponds to your new calculation type!). See :ref:`writers` to learn more.

See :ref:`general` for general guidance on contributing to QMzyme to get started.
