
MATWRITE
========

MATWRITE is a Stata plugin that allows the user to export data to the MATLAB
.mat file format from within Stata. MATWRITE allows you to export all (or 
selected) variables as column vectors, and all (or selected) matrices as
matrices. It is written by Andrew Shephard (<asheph@econ.upenn.edu>).

Installation
------------

If you have web-aware Stata, the recommended way of installing is to type `ssc 
install matwrite` from within Stata. Note that this may not be the most up to 
date version.

MATWRITE uses the [Stata Plugin Interface](http://www.stata.com/plugins/). This
means that the plugin file is platform specific. To compile the plugin file on a
Unix-like operating system type:

```g++ -shared -fPIC -DSYSTEM=OPUNIX stplugin.c matwrite.cpp -o matwrite.plugin```

see the [Stata Plugin Interface](http://www.stata.com/plugins/) documentation for
information on compiling for Windows or Macintosh systems. After compiling
the plugin file (*matwrite.plugin*) copy *matwrite.plugin*, *matwrite.ado* and 
*matwrite.hlp* to your personal ado file directory.

Usage
-----

Type `help matwrite` within Stata for command syntax.


License
-------

MATWRITE is free software: you can redistribute it and/or modify it under the 
terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

MATWRITE is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
MATWRITE.  If not, see <http://www.gnu.org/licenses/>.

