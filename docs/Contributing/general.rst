.. _general:

General Guidance
====================

Fork repository
------------------

Visit https://github.com/hklem/QMzyme. On the top right of the page, click Fork 
> create a new fork (fill in details) > Create fork.

Make changes
------------------

Make sure to include docstrings and comments. The easier the code is to read, 
the more likely the changes are to be accepted.

Some modules, such as SelectionSchemes.py and TruncationSchemes.py are 
accompanied by abstract classes under the same name that prescribe how 
additional schemes should be constructed. To better understand abstract 
classes and methods, see https://docs.python.org/3/library/abc.html.

Test changes
------------------

If you are adding or changing code functionality, you will also be asked 
to implement tests before the pull request is accepted and merged. The tests 
are located in the /tests/ directory. If you modified an existing module, 
find the corresponding test_MODULENAME.py file and add one or more tests there.

Create pull request
------------------

See https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork.

