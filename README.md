﻿# abaqus-py

Unsorted scripts used with the Abaqus FE software.

Most functions are relying on abaqus python modules and must therefore be run in abaqus python interpreter. 

Warning: Most stuff is quickly thrown together and is not really tested. Therefore, use with care.

## Usage

Simplest way us perhaps to load the module of interest using Python's ``imp`` module (See https://docs.python.org/2/library/imp.html).

Example:
```
import imp
odbToCSV = imp.load_source('.', path_name)
odbToCSV.field_to_csv(...)
```
where path_name is the path to the module source file.
