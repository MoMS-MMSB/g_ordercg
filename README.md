# g_ordercg: Calculate order parameter on coarse grained models
![maintenance-status](https://img.shields.io/badge/maintenance-deprecated-red.svg)

Deuterium order parameter, as calculated by ``g_order`` and measured
experimentally by NMR spectroscopy, cannot be calculated for coarse grained models
since not all carbons are explicitly represented. In such models the
ordering of molecules can be estimated using a second-rank order parameter
defined as:
$$
    S = 1/2 * (3 * <cos^2(θ)> - 1)
$$
where θ is the angle between a predefined axis (i.e. the membrane
normal) and the bond between two successive beads. The brakets < > represent an
ensemble and time averaging. According to this formula, the order parameter will be 1 if
the molecules are perfectly aligned with the axis, -0.5 if they are perfectly
normal to it and 0 if the θ angle is distributed randomly.

The ``g_ordercg`` program calculates the second rank order parameter from
`GROMACS <http://www.gromacs.org>`_ trajectories. Several output are available:

* the order parameter for each bond as a function of time;
* the order parameter profile as a function of an axis coordinate;
* the order parameter profile as a function of the distance from a group of atoms;
* the order parameter on a grid (i.e., as a function of two axes).

This program has been written by Jonathan Barnoud <jonathan.barnoud@inserm.fr>.
It is based on GROMACS source code and is distributed under the terms
of the General Public Licence version 2 or greater. See the LICENCE file for
more information.

## Installation

``g_ordercg`` depends on the GROMACS software package, which needs to be
installed.  Only versions 4.5.x have been tested, but ``g_ordercg`` might be
compatible with other versions of GROMACS.

To install ``g_ordercg``, GROMACS needs to be loaded. You can load
it using:

    source /path_to_gromacs/bin/GMXRC

Go into the source directory of the program, then run ``make``. The
``g_ordercg`` executable should be created.  Make sure that
this executable is in the research path of your shell.

## Usage
Here we assume that ``g_ordercg`` is in the research path of your shell. To get
some help just run ``g_ordercg -h``. All available options will be listed.

A classical use would be:

    g_ordercg -f traj.xtc -s topol.tpr -n index.ndx -o order.xvg

GROMACS needs to be loaded for ``g_ordercg`` to work.

### Required arguments
Three arguments are required for any use of the program. They are:

* ``-f``: the path to the trajectory to read;
* ``-s``: the path to the topology (tpr file), except if you use the ``-t``
  option (see the `Advanced input`_ section bellow) the bonds will be read from
  this file;
* ``-n``: the path to an index file that describe the group of atoms you are
  interested in;

The ``-d`` option define the reference axis (i.e. the axis normal to the
membrane). The axis is set at Z by default.

### Output control
The following arguments control the output. You can get either one or several of
the possible outputs but you need to select at least one of them.

#### Order parameter as a function of time


The ``-o`` option can be used to get the order parameter for each bond as a
function of time. The option can also be used to define the output file name,
the default file name if ``order.xvg``.

#### Order parameter profile along an axis
Use ``-op`` to get the order parameter profile as a function of the coordinate
along an axis. The default output file is ``order_profile.xvg``.

By default the profile is plotted against the reference axis defined with the
``-d`` option, you can plot the profile against a different axis using the
``-dp`` option. The value of this option can be x, y or z, it is d by default
that means that the reference axis is used.

The number of bonds used per bin is written in the file
``profile_sampling.xvg``. You can change the name of this file using the
``-osp`` option.

The number of bins can be chosen using the ``-sl`` option. Be aware that the
number of bins is constant but not their size, be careful if your box size
fluctuates a lot.

#### Order parameter profile as a function of the distance to a group of atoms

The ``-od`` option activates the output of the order parameter profile as a
function of the distance to the center of mass of a group of atom. The default
output file name is ``order_distprof.xvg`` but you can change it with ``-od``.

The number of bonds used in each bins is recorded in ``distprof_sampling.xvg``.
The name of this file can be changed thanks to the ``-osd`` option.

The calculated distance is the distance between the middle of a bond and the
center of mass of the reference group of atoms. Distance calculation is done,
by default, in two dimensions ignoring the dimension of the reference axis. In
membrane systems this corresponds to calculate distances only in the membrane
plan. This behavior can be changed using the ``-2D`` option. The authorized
values are x, y and z to ignore one of the unit axis, but also d that ignores
the axis defined by the ``-d`` option and n to not ignore any axis and do the
distance calculation in 3D.

As for the previous output, the number of bins can be set with the ``-sl``
option.

#### Order parameter landscape
To get the order parameter landscape on a membrane plan, use the ``-og`` option.
The ``order_grid.dat`` and ``grid_sampling.dat`` files will be written. These
files are the order parameter grid and the sampling record for each cell
respectively. File names can be changed using ``-og`` for the data and ``-osg``
for the sampling output.

The plane of the landscape is the plane normal to the reference axis defined by
the ``-d`` option. You can control the number of cells in each direction using
``-sl`` and ``-sl2``. If the value of the ``-sl2`` option is negative (default),
then the value of ``-sl`` is used for both axes. To know what dimensions are
controlled by these two options, see the table below:
```
====== ======= ========
``-d`` ``-sl`` ``-sl2``
====== ======= ========
x      y       z
y      x       z
z      x       y
====== ======= ========
```
To convert the output file into an image see the `Generate pictures from
landscapes`_ section bellow.

### Advanced input

By default, the program read the list of bonds per molecule from the topology
input (``-s`` argument). The list of bonds that is then used is written in the
output files for reference.

You may not want to calculate the order parameter using all the bonds of the
molecules of interest or you might want to use vectors that are not considered
as bonds in the topology. The former case can be useful, for instance, if you
want to ignore lipid polar heads from the calculation, the latter case can be
use with some polymers. You can pass to the program an index file describing
the topology of your molecule using the ``-t`` option. The index file needs to
contain a group that specifies for each bond the internal index of the
bead/pseudo-atom at its extremity. This index group looks similar to the two
first columns of the ``[ bonds ]`` section as written in an itp file. Here is
an example for the POPC molecule in the `Martini <http://www.cgmartini.nl>`_
force field: ::
```
    [ full_POPC ]
    1 2 	
    2 3 	
    3 4 	
    4 5 	
    5 6 	
    6 7 	
    7 8 	
    3 9 	
    9 10 	
    10 11 
    11 12 
    12 13 

    [ POPC_tails ]
    5 6 	
    6 7 	
    7 8 	
    9 10 	
    10 11 
    11 12 
    12 13
```
If you are interested in itermolecular vectors you can use the ``-vectors``
option. This change drastically the behavior of the program: the reference index
group is read as a list of vectors with each pair of indices considered as the
extremity beads/pseudo-atoms of the vectors. This vectors are used for the order
parameter calculation instead of those described in the topology. Be aware that
this option is less tested than the others.

If ``-novectors`` is used (default) the program reads the reference group and
build a non-redundant list of molecules of interest from it. If you do not want
the program to behave this way, you can use the ``-first`` option, then the
program will read the indices of the reference group and consider them as the
first atoms of the molecules. The atom numbered 1 in the molecule connectivity
as read in the topology or from the ``-t`` option will refer to this atoms. Use
this option with caution.

GROMACS usual options to control trajectory reading apply here. You can use
``-b``, ``-e`` and ``-dt`` with a time in picoseconds to control the beginning
time, the end time and the time step to use for the analysis.

### Generate pictures from landscapes
The landscape output is a text file describing the order parameter values on a
grid. The file format is not XPM like most grid outputs produced by GROMACS
tools so the ``xpm2ps`` utility can not be used to produce usable pictures. The
``dispgrid`` python script aims to exploit the data and to produce pictures from
them.

You need python 2.7, with the numpy and matplotlib modules to run dispgrid.

Basic usage of dispgrid is:

    dispgrid input.dat output.png

See the help available by typing ``dispgrid -h`` for more features.
