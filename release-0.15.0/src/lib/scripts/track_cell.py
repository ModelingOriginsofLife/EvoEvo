
#!/usr/bin/env python
# coding: utf-8

#***************************************************************************
# Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon
# E-mail: charles.rocabert@inria.fr
# Web: http://www.evoevo.eu/
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#***************************************************************************

import sys
import os
import time as timelib
from matplotlib.pyplot import *

### Read track.txt data ###
def read_data():
	g_data   = {}
	i_data   = {}
	in_data  = {}
	out_data = {}
	score    = []
	time     = []

	#--------------------#
	# GET VARIABLE NAMES #
	#--------------------#
	f = open("track.txt", "r")
	l = f.readline()
	variables = l.strip("\n").split(" ")
	for var in variables:
		if var.startswith("g"):
			g_data[var] = []
		elif var.startswith("i") and not var.startswith("in"):
			i_data[var] = []
		elif var.startswith("in"):
			in_data[var] = []
		elif var.startswith("out"):
			out_data[var] = []

	#----------#
	# GET DATA #
	#----------#
	l = f.readline()
	while l:
		l = l.strip("\n")
		l = l.split(" ")
		if len(l) == len(variables):
			for i in range(len(variables)):
				if variables[i].startswith("g"):
					g_data[variables[i]].append(float(l[i]))
				elif variables[i].startswith("i") and not variables[i].startswith("in"):
					i_data[variables[i]].append(float(l[i]))
				elif variables[i].startswith("in"):
					in_data[variables[i]].append(float(l[i]))
				elif variables[i].startswith("out"):
					out_data[variables[i]].append(float(l[i]))
				elif variables[i].startswith("score"):
					score.append(float(l[i]))
				elif variables[i].startswith("t"):
					time.append(float(l[i]))
			l = f.readline()
		else:
			break
	f.close()

	#-----------#
	# SORT DATA #
	#-----------#
	for var in g_data.keys():
		zero = True
		for i in range(len(g_data[var])):
			if g_data[var][i] > 0.0:
				zero = False
		if zero:
			del g_data[var]

	for var in i_data.keys():
		zero = True
		for i in range(len(i_data[var])):
			if i_data[var][i] > 0.0:
				zero = False
		if zero:
			del i_data[var]

	for var in in_data.keys():
		zero = True
		for i in range(len(in_data[var])):
			if in_data[var][i] > 0.0:
				zero = False
		if zero:
			del in_data[var]

	for var in out_data.keys():
		zero = True
		for i in range(len(out_data[var])):
			if out_data[var][i] > 0.0:
				zero = False
		if zero:
			del out_data[var]

	#-------------#
	# RETURN DATA #
	#-------------#
	return g_data, i_data, in_data, out_data, score, time


#############
#   MAIN    #
#############

ion()

while 1:
	
	f  = open("last_simulation_timestep.txt", "r")
	l  = f.readline()
	ct = l.strip("\n")
	f.close()

	g_data, i_data, in_data, out_data, score, time = read_data()

	clf()

	figure(1)

	subplot(2,2,1)
	for var in g_data.keys():
		plot(g_data[var])
	title("Produced enzymes ("+str(ct)+")")

	subplot(2,2,2)
	for var in i_data.keys():
		plot(i_data[var])
	title("Inherited enzymes ("+str(ct)+")")

	subplot(2,2,3)
	for var in in_data.keys():
		plot(in_data[var])
	title("Internal metabolites ("+str(ct)+")")

	subplot(2,2,4)
	for var in out_data.keys():
		plot(out_data[var])
	title("External metabolites ("+str(ct)+")")

	draw()

	timelib.sleep(1)

