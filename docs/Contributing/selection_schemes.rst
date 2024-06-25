Adding a Selection Scheme 
===========================

:class:`~QMzyme.SelectionSchemes.SelectionScheme` is an `abstract base class (ABC)
<https://docs.python.org/3/library/abc.html>`_ used to prescribe concrete selection 
scheme classes. Please note, if you submit a GitHub Pull Request to include a new 
selection scheme you will need to have also added appropriate tests in the QMzyme/tests 
directory. The docstring of the `SelectionScheme` base class describes how you can create 
your own subclass:

.. autoclass:: QMzyme.SelectionSchemes.SelectionScheme
  :members:
