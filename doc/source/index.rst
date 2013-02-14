.. Msys documentation master file, created by
   sphinx-quickstart on Wed Oct 12 12:31:22 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Msys library
====================================

Msys is a library for manipulating molecular structures and their
associated forcefields.  Among its many features:

* Reads and writes DMS and MAE file formats;

* Implements the VMD atom selection language;

* Provides a high-performance C++ interface as well as a friendlier Python
  interface;

* Contains a collection of command line and scripting tools for examining,
  combining, and processing DMS files.


Msys is written in C++, with bindings for Python.  The direct Python
bindings are generated using boost::python, and are intended to be
used only by developers, or those familiar with the C++ interface.
The rest of this guide describes the high level Python interface,
which is based on the Python bindings.  

Contents:

.. toctree::
   :maxdepth: 2

   overview
   selections
   paramtable
   termtable
   system
   nonbonded

Command line tools
==================

Msys is packaged with a set of command line tools that wrap functionality
present in the Python interface.

.. toctree::

   tools


DMS files
=========

The preferred file format for Msys structures is the DMS file.

.. toctree::

   dmsversion
   dms


MAE files
=========

Msys tries to support MAE files to the extent possible.

.. toctree::

   mae


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

