# piecewiseValues.r
# -----------------
# Begin: 2008/11/11
# Last revision: $Date: 2008-11-25 23:24:25 $ $Author: areeves $
# Version: $Revision: 1.1 $
# Project: Random number distributions for Delphi: verification code
# Website: http://www.naadsm.org
# Author: Aaron Reeves <Aaron.Reeves@colostate.edu>
# --------------------------------------------------
# Copyright (C) 2008 Animal Population Health Institute, Colorado State University
# 
# This program is free software; you can redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.


# The values in this table were generated with TPdfPiecewise.rand()
# from a function defined by the following points:
#
#  x, y
#  0, 0
# 10, 0.038468
# 20, 0.00097
# 30, 0.00097
# 45, 0.039567
# 60, 0

dat <- read.table( "piecewiseValues.txt" );


# The histogram generated from the dataset can be used as a quick visual check that the
# values produced by the rand() function follow the distribution.
#
# The Microsoft Word document in this folder contains graphical examples.

hist( dat[,1], freq = FALSE );

