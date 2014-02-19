.. Msys documentation master file, created by
   sphinx-quickstart on Wed Oct 12 12:31:22 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The Msys library
================

Msys is a library for inspecting and manipulating chemical structures of
the sort used in molecular simulations.  Its main features include:

 * A hierarchical organization of structure information, as opposed to the
   flat tables typically used in molecular file formats;

 * A powerful VMD-like atom selection language;

 * Representation of molecular forcefields which supports sharing of
   parameters by many groups of terms;

 * Conversion to and from many chemical file formats, with forcefield
   and chemical information preserved where possible;

 * Command line tools for scripting common tasks;

 * C++ and Python interfaces.

In the rest of this introduction, we go into more detail about the
Msys data model and how to perform common tasks.  The descriptions
will be of the interface as seen from the Python layer.  

Contents:

.. toctree::
   :maxdepth: 2

   overview
   selections
   python
   paramtable
   termtable
   system
   nonbonded
   pfx


Command line tools
==================

Msys is packaged with a set of command line tools that wrap functionality
present in the Python interface.

.. toctree::

   tools


Recipes
=======

Here are some illustrative Python scripts for situations when a command line 
tool isn't available.

.. toctree::

   recipes


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

