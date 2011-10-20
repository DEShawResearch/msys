.. Msys documentation master file, created by
   sphinx-quickstart on Wed Oct 12 12:31:22 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Msys library
====================================

Msys is a library for manipulating molecular structures and their
associated forcefields.  If you are already familiar with the Ent
library, you'll see that the two have much in common.

Msys is written in C++, with bindings for Python.  The direct Python
bindings are generated using boost::python, and are intended to be
used only by developers, or those familiar with the C++ interface.
The rest of this guide describes the high level Python interface,
which is based on the Python bindings.  

Contents:

.. toctree::
   :maxdepth: 2

   overview
   paramtable
   termtable
   system

Command line tools
==================

Msys is packaged with a set of command line tools that wrap functionality
present in the Python interface.

.. toctree::

   tools


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

