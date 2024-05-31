.. _writers:

Adding a Writer
================

:class:`~QMzyme.Writers.Writer` is called directly by 
:class:`~QMzyme.GenerateModel.GenerateModel` to automate writing of 
calculation input files. If you would like to contribute by expanding the 
suite of supported input, you can write your own concrete writer class! 

See :ref:`general` for general guidance on contributing to QMzyme to get started.

If you just wish for your writer to be used directly (i.e., someone would
import your class from QMzyme.Writers and call it, providing all arguments),
then all you have to do is write the class and add appropriate tests. 

However, in order for your writer class to be automatically detectable by GenerateModel 
you will need to do some steps in addition to writing your class and creating tests:
    - register your class to :class:`~QMzyme.Writers.WritersRegistry`.
    - add code to :class:`~QMzyme.GenerateModel.GenerateModel`, see :ref:`calculation_methods`.

