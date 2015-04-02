Gray-Box Optimization using the Parameter-less Population Pyramid
==

To compile you will need C++11.  We use gcc version 4.8.2 for our complilation.

Our build system uses Makefiles to build.  You can compile the release version
by changing directory to Release and calling "make GBO".

All of the source code is available in the 'src' directory.

To run an experiment, call the executable with command line arguments for configuration.
This will run the default test configuration:

Release/GBO config/default.cfg

The command line accepts any number of configuration files and configuration flags given
in the form "-key value" where key is the name of the option to change and value is the
value to change it to.  Arguments override in last given order.  For example:

Release/GBO config/default.cfg -problem NearestNeighborNKQ -length 1000 -verbosity 2

This line will use the default configuration, setting the problem to Nearest Neighbor NKQ
with a genome size of 1000, and increasing the amount of information output to the screen.

Also, see main.cpp for implementation related details.

The data used in the paper is available as part of the
[GBO-using-P3.tar.gz](https://drive.google.com/open?id=0BwZEhB6mzmkadm5qbUhFZTdXU2M&authuser=0) file.
All of the graphing and statistical processing scripts are contained in data_parsing/.

