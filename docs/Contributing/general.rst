.. _general:

General Guidance
====================

Fork repository
------------------

Visit https://github.com/hklem/QMzyme. On the top right of the page, click Fork 
> create a new fork (fill in details) > Create fork.

Make changes
--------------

Make sure to include docstrings and comments. The easier the code is to read, 
the more likely the changes are to be accepted.

Some modules, such as SelectionSchemes.py and TruncationSchemes.py are 
accompanied by abstract classes under the same name that prescribe how 
additional schemes should be constructed. To better understand abstract 
classes and methods, see https://docs.python.org/3/library/abc.html.

Viewing Documentation Locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also add onto the read the docs webpage in the QMzyme/docs directory.
To build the documentaion pages locally go to the QMzyme/docs dirctory and run:

.. code-block:: bash

    pip install -r requirements.txt # to install sphinx and other docs related packages
    brew install pandoc # pandoc is required to configure the jupyter notebook files used as tutorials. If you are not on MacOS see pandoc webpage for other installation methods: https://pandoc.org/installing.html
    make html
    open _build/html/index.html # a web browswer should open

There may be some WARNING messages. For the most part, they are probably okay to ignore (i.e., 'WARNING: document isn't included in any toctree')

Test Changes
------------------

If you are adding or changing code functionality, you will also be asked 
to implement tests before the pull request is accepted and merged. The tests 
are located in the /tests/ directory. If you modified an existing module, 
find the corresponding test_MODULENAME.py file and add one or more tests there.

QMzyme uses the pytest package for its testing suite. If you are unfamiliar with pytest 
you can learn more on their [documentation webpage](https://docs.pytest.org/en/8.2.x/).

Running Tests Locally
~~~~~~~~~~~~~~~~~~~~~~~
In the main QMzyme directory run:

.. code-block:: bash

    pip install -e .[test] # installs pytest and pytest-cov
    # To run only a single test file:
    pytest tests/test_FILENAME.py -vv --color yes
    # To run all tests
    pytest -vv --color yes

If any tests failed review the verbose output to understand why and fix the code accordingly.


Create pull request
------------------

See https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork.
Pull requests will only be accepted if all tests pass.
