
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

### PRINT USAGE ####
# param  void
def print_usage():
  print ""
  print "=== WRITE JAVASCRIPT FILES FROM STATISTICS ==="
  print "Usage: python statistics_to_js.py [parameters]"
  print "Parameters are:"
  print "-h, --help:"
  print "    Print this help, then exit."
  print "-p, --path:"
  print "    Give the path of the simulation folder."
  print ""

### READ COMMAND LINE ARGUMENTS ###
# param  string-array args : list of command line arguments
def read_args( args ):
  if len(args) < 2:
    print("Lack of parameters, see help (-h --help)")
    sys.exit()
  else:
    PATH    = ""
    for i in range(len(args)):
      if args[i] == "-h" or args[i] == "--help":
        print_usage()
        sys.exit()
      elif args[i] == "-p" or args[i] == "--path":
        if i+1 >= len(args):
          print("Lack of parameters, see help (-h --help)")
          sys.exit()
        else:
          PATH = args[i+1]
    if PATH == "":
      print("Lack of parameters, see help (-h --help)")
      sys.exit()
    else:
      return PATH

### GET THE LIST OF VARIABLES IN A DICTIONARY ###
# param  string filename : statistical file to parse
def get_variable_names( filename ):
  f = open(filename, "r")
  l = f.readline()
  l = l.strip("\n")
  l = l.split(" ")
  dico = {}
  for i in range(len(l)):
    dico[l[i]] = i
  f.close()
  return dico

### GET THE VARIABLE EVOLUTION ###
# param  string variable_name  : variable name (key of the variables dictionary)
# param  string filename       : statistical file to parse
# param  dict   variable_dict  : variables dictionary (extracted from the same file)
# param  int    vector_length  : length of the data vector
# param  bool   data_cumulated : extract cumulated variables
# param  bool   data_reversed  : reverse data vectors
def get_variable(
  variable_name,
  filename,
  variable_dict,
  vector_length,
  data_cumulated,
  data_reversed ):
  #------------------------------------#
  # 1) get the vector length           #
  #------------------------------------#
  f = open(filename, "r")
  l = f.readline()
  l = f.readline()
  COUNT = 0
  while l:
    COUNT += 1
    l = f.readline()
  f.close()

  #------------------------------------#
  # 2) compute the mean period         #
  #------------------------------------#
  PERIOD = 1
  if COUNT > vector_length:
    PERIOD = COUNT/vector_length

  #------------------------------------#
  # 3) get the data vector             #
  #------------------------------------#
  var   = []
  count = 0
  mean  = 0.0
  line_counter = 0
  f = open(filename, "r")
  l = f.readline()
  line_counter += 1
  l = l.strip("\n")
  l = l.split(" ")
  L = len(l)
  l = f.readline()
  while l:
    line_counter += 1
    l = l.strip("\n")
    l = l.split(" ")
    if len(l) == L:
      if count < PERIOD:
        try:
          mean += float(l[variable_dict[variable_name]])
        except:
          print "Incorrect line structure in file '"+filename+"' at line "+str(line_counter)+". Exit."
          print variable_dict
          print variable_name
          return []
          #sys.exit()
        count += 1
      else:
        mean /= PERIOD
        var.append(mean)
        count = 0
        mean  = 0.0
      l = f.readline()
    else:
      break
  f.close()

  #------------------------------------#
  # 4) revert data if asked            #
  #------------------------------------#
  if data_reversed:
    reverted_data = []
    for i in range(len(var)):
      reverted_data.append(var[(len(var)-i-1)])
    var = reverted_data

  #------------------------------------#
  # 5) compute cumulated data if asked #
  #------------------------------------#
  if data_cumulated:
    cum = 0.0
    cumulated_data = []
    for i in range(len(var)):
      cum += var[i]
      cumulated_data.append(cum)
    var = cumulated_data

  #------------------------------------#
  # 6) return data vector              #
  #------------------------------------#
  return var

### WRITE DYGRAPH JS FILE WITH SPECIFIED DATA ###
# param  string-array variables_names : list of variables to extract
# param  dict         variable_dict   : variables dictionary
# param  string       js_filename     : javascript output file name
# param  string       div_id          : html div identifier
# param  string       title           : figure title
# param  string       xlab            : figure x label
# param  string       ylab            : figure y label
# param  string-array dygraph_options : dygraph object options
# param  string       dygraph_name    : dygraph object name
# param  int          vector_length   : length of the data vector
# param  bool         data_cumulated  : extract cumulated variables
# param  bool         data_reversed   : reverse data vectors
def write_dygraph_js(
  variables_names,
  variable_dict,
  js_filename,
  filename,
  div_id,
  title,
  xlab,
  ylab,
  dygraph_options,
  dygraph_name,
  vector_length,
  data_cumulated,
  data_reversed ):
  #--------------------------#
  # 1) get data              #
  #--------------------------#
  MINL = 2000000000
  data = []
  for variable_name in variables_names:
    if variable_name == "t" or variable_name == "generations":
      var = get_variable(variable_name, filename, variable_dict, vector_length, False, data_reversed)
    else:
      var = get_variable(variable_name, filename, variable_dict, vector_length, data_cumulated, data_reversed)
    if MINL > len(var):
      MINL = len(var)
    data.append(var)

  f = open(js_filename, "w")

  #--------------------------#
  # 2) write header          #
  #--------------------------#
  f.write("/****************************************************************************\n")
  f.write(" * Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon\n")
  f.write(" * E-mail: charles.rocabert@inria.fr\n")
  f.write(" * Web: http://www.evoevo.eu/\n")
  f.write(" *\n")
  f.write(" * This program is free software: you can redistribute it and/or modify\n")
  f.write(" * it under the terms of the GNU General Public License as published by\n")
  f.write(" * the Free Software Foundation, either version 3 of the License, or\n")
  f.write(" * (at your option) any later version.\n")
  f.write(" *\n")
  f.write(" * This program is distributed in the hope that it will be useful,\n")
  f.write(" * but WITHOUT ANY WARRANTY; without even the implied warranty of\n")
  f.write(" * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n")
  f.write(" * GNU General Public License for more details.\n")
  f.write(" *\n")
  f.write(" * You should have received a copy of the GNU General Public License\n")
  f.write(" * along with this program.  If not, see <http://www.gnu.org/licenses/>.\n")
  f.write(" ****************************************************************************/\n\n")

  #--------------------------#
  # 3) write data in js file #
  #--------------------------#
  f.write("var data = [];\n")
  for i in range(MINL):
    line = "data.push(["
    for j in range(len(variables_names)):
      line += str(data[j][i])+", "
    line = line.strip(", ")
    line += "]);"
    f.write(line)
  f.write("\n")

  #--------------------------#
  # 4) write dygraph class   #
  #--------------------------#
  f.write(dygraph_name+" = new Dygraph(\n")
  f.write("document.getElementById(\""+div_id+"\"),\n")
  f.write("data,\n")
  f.write("{\n")
  f.write("title: '"+title+"',\n")
  f.write("xlabel: '"+xlab+"',\n")
  f.write("ylabel: '"+ylab+"',\n")
  line = "labels: ["
  for variable_name in variables_names:
    line += "'"+variable_name +"', "
  line = line.strip(", ")
  line += "],\n"
  f.write(line)
  #----------------#
  # CUSTOM OPTIONS #
  #----------------#
  for option in dygraph_options:
    f.write(option+",\n")
  #----------------#
  # COMMON OPTIONS #
  #----------------#
  f.write("legend: 'always',\n")
  #f.write("rollPeriod: "+str(period)+",\n")
  #f.write("rollPeriod: 1,\n")
  #f.write("showRoller: true,\n")
  f.write("labelsDivStyles: { 'textAlign': 'right' },\n")
  f.write("animatedZooms: true}\n")
  f.write(");\n\n");
  f.close()


############
#   MAIN   #
############
if __name__ == '__main__':
  
  PATH = read_args(sys.argv)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 1) DEFINE OPTIONS                                     #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  stacked_option = ["stackedGraph: true"]
  filled_option  = ["fillGraph: true"]
  vector_length  = 500
  
  # POPULATION FIGURES ---------------------------------------------------------------------#
  
  filename       = PATH+"statistics/phenotype_mean.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = False
  
  # POPULATION SIZE #
  write_dygraph_js(["t", "population_size"], dico, PATH+"viewer/src/js/population_size.js", filename, "population_size_div", "Population size", "Simulation time", "Population size", stacked_option, "population_1", vector_length, data_cumulated, data_reversed)
  
  # GROWTH RATE #
  write_dygraph_js(["t", "growth_rate"], dico, PATH+"viewer/src/js/growth_rate.js", filename, "growth_rate_div", "Population growth rate", "Simulation time", "Growth rate", stacked_option, "population_2", vector_length, data_cumulated, data_reversed)
  

  # PHENOTYPE FIGURES ----------------------------------------------------------------------#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 2) GENERATE MEAN PHENOTYPE FIGURES                    #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  filename       = PATH+"statistics/phenotype_mean.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = False

  counter = 1

  # PROTEIN AMOUNTS #
  write_dygraph_js(["t", "inherited_TF_amount", "inherited_E_amount", "TF_amount", "E_amount"], dico, PATH+"viewer/src/js/mean_protein_amounts.js", filename, "mean_protein_amounts_div", "Protein amounts", "Simulation time", "Concentration", stacked_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # METABOLIC AMOUNTS #
  write_dygraph_js(["t", "inherited_metabolic_amount", "metabolic_amount", "cumulated_error"], dico, PATH+"viewer/src/js/mean_metabolic_amounts.js", filename, "mean_metabolic_amounts_div", "Metabolic amounts", "Simulation time", "Concentration", stacked_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SCORE #
  write_dygraph_js(["t", "score"], dico, PATH+"viewer/src/js/mean_score.js", filename, "mean_score_div", "Score", "Simulation time", "Score", filled_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # ENERGY #
  write_dygraph_js(["t", "energy"], dico, PATH+"viewer/src/js/mean_energy.js", filename, "mean_energy_div", "Energy", "Simulation time", "Energy", filled_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # LIFESPAN #
  write_dygraph_js(["t", "lifespan"], dico, PATH+"viewer/src/js/mean_lifespan.js", filename, "mean_lifespan_div", "Lifespan", "Simulation time", "Lifespan (simulation steps)", filled_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF DIVISIONS #
  write_dygraph_js(["t", "number_of_divisions"], dico, PATH+"viewer/src/js/mean_nb_divisions.js", filename, "mean_nb_divisions_div", "Number of divisions", "Simulation time", "Number of divisions", filled_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # TOXICITY ACCUMULATION #
  write_dygraph_js(["t", "toxicity"], dico, PATH+"viewer/src/js/mean_toxicity.js", filename, "mean_toxicity_div", "Toxicity accumulation", "Simulation time", "Toxicity level", filled_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S EXCHANGES #
  write_dygraph_js(["t", "metabolic_uptake", "metabolic_release"], dico, PATH+"viewer/src/js/mean_exchanges.js", filename, "mean_exchanges_div", "Exchanges with environment", "Simulation time", "Concentrations", stacked_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S METABOLIC GROWTH RATE #
  write_dygraph_js(["t", "metabolic_growth_rate", "Dmetabolic_growth_rate"], dico, PATH+"viewer/src/js/mean_metabolic_growth_rate.js", filename, "mean_metabolic_growth_rate_div", "Metabolic growth rate", "Simulation time", "Growth rate", filled_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S GRN NODES AND EDGES #
  write_dygraph_js(["t", "grn_nb_nodes", "grn_nb_edges"], dico, PATH+"viewer/src/js/mean_grn_nodes_edges.js", filename, "mean_grn_nodes_edges_div", "Genetic regulation network", "Simulation time", "Number of elements", filled_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S METABOLIC NODES AND EDGES #
  write_dygraph_js(["t", "metabolic_nb_nodes", "metabolic_nb_edges"], dico, PATH+"viewer/src/js/mean_metabolic_nodes_edges.js", filename, "mean_metabolic_nodes_edges_div", "Metabolic network", "Simulation time", "Number of elements", filled_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # REGULATION REDUNDANCY #
  write_dygraph_js(["t", "regulation_redundancy"], dico, PATH+"viewer/src/js/mean_regulation_redundancy.js", filename, "mean_regulation_redundancy_div", "Regulation redundancy", "Simulation time", "Number of redundant tuples", filled_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # METABOLIC REDUNDANCY #
  write_dygraph_js(["t", "metabolic_redundancy"], dico, PATH+"viewer/src/js/mean_metabolic_redundancy.js", filename, "mean_metabolic_redundancy_div", "Metabolic redundancy", "Simulation time", "Number of redundant tuples", filled_option, "mean_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 3) GENERATE BEST PHENOTYPE FIGURES                    #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  filename       = PATH+"statistics/phenotype_best.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = False

  counter = 1

  # PROTEIN AMOUNTS #
  write_dygraph_js(["t", "inherited_TF_amount", "inherited_E_amount", "TF_amount", "E_amount"], dico, PATH+"viewer/src/js/best_protein_amounts.js", filename, "best_protein_amounts_div", "Protein amounts", "Simulation time", "Concentration", stacked_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # METABOLIC AMOUNTS #
  write_dygraph_js(["t", "inherited_metabolic_amount", "metabolic_amount", "cumulated_error"], dico, PATH+"viewer/src/js/best_metabolic_amounts.js", filename, "best_metabolic_amounts_div", "Metabolic amounts", "Simulation time", "Concentration", stacked_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SCORE #
  write_dygraph_js(["t", "score"], dico, PATH+"viewer/src/js/best_score.js", filename, "best_score_div", "Score", "Simulation time", "Score", filled_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # ENERGY #
  write_dygraph_js(["t", "energy"], dico, PATH+"viewer/src/js/best_energy.js", filename, "best_energy_div", "Energy", "Simulation time", "Energy", filled_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # LIFESPAN #
  write_dygraph_js(["t", "lifespan"], dico, PATH+"viewer/src/js/best_lifespan.js", filename, "best_lifespan_div", "Lifespan", "Simulation time", "Lifespan (simulation steps)", filled_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF DIVISIONS #
  write_dygraph_js(["t", "number_of_divisions"], dico, PATH+"viewer/src/js/best_nb_divisions.js", filename, "best_nb_divisions_div", "Number of divisions", "Simulation time", "Number of divisions", filled_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # TOXICITY ACCUMULATION #
  write_dygraph_js(["t", "toxicity"], dico, PATH+"viewer/src/js/best_toxicity.js", filename, "best_toxicity_div", "Toxicity accumulation", "Simulation time", "Toxicity level", filled_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S EXCHANGES #
  write_dygraph_js(["t", "metabolic_uptake", "metabolic_release"], dico, PATH+"viewer/src/js/best_exchanges.js", filename, "best_exchanges_div", "Exchanges with environment", "Simulation time", "Concentrations", stacked_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S METABOLIC GROWTH RATE #
  write_dygraph_js(["t", "metabolic_growth_rate", "Dmetabolic_growth_rate"], dico, PATH+"viewer/src/js/best_metabolic_growth_rate.js", filename, "best_metabolic_growth_rate_div", "Metabolic growth rate", "Simulation time", "Growth rate", filled_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S GRN NODES AND EDGES #
  write_dygraph_js(["t", "grn_nb_nodes", "grn_nb_edges"], dico, PATH+"viewer/src/js/best_grn_nodes_edges.js", filename, "best_grn_nodes_edges_div", "Genetic regulation network", "Simulation time", "Number of elements", filled_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S METABOLIC NODES AND EDGES #
  write_dygraph_js(["t", "metabolic_nb_nodes", "metabolic_nb_edges"], dico, PATH+"viewer/src/js/best_metabolic_nodes_edges.js", filename, "best_metabolic_nodes_edges_div", "Metabolic network", "Simulation time", "Number of elements", filled_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # REGULATION REDUNDANCY #
  write_dygraph_js(["t", "regulation_redundancy"], dico, PATH+"viewer/src/js/best_regulation_redundancy.js", filename, "best_regulation_redundancy_div", "Regulation redundancy", "Simulation time", "Number of redundant tuples", filled_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # METABOLIC REDUNDANCY #
  write_dygraph_js(["t", "metabolic_redundancy"], dico, PATH+"viewer/src/js/best_metabolic_redundancy.js", filename, "best_metabolic_redundancy_div", "Metabolic redundancy", "Simulation time", "Number of redundant tuples", filled_option, "best_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # GENOME STRUCTURE FIGURES ----------------------------------------------------------------------#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 4) GENERATE MEAN GENOME STRUCTURE FIGURES             #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  filename       = PATH+"statistics/genome_structure_mean.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = False

  counter = 1

  # GENOME SIZE #
  write_dygraph_js(["t", "genome_size", "functional_size"], dico, PATH+"viewer/src/js/mean_genome_size.js", filename, "mean_genome_size_div", "Genome size", "Simulation time", "Number of tuples", filled_option, "mean_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF PEARLS #
  write_dygraph_js(["t", "nb_NC", "nb_E", "nb_TF", "nb_BS", "nb_P"], dico, PATH+"viewer/src/js/mean_nb_pearls.js", filename, "mean_nb_pearls_div", "Genome composition", "Simulation time", "Number of tuples by type", stacked_option, "mean_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # ENZYME TUPLE CLASSES #
  write_dygraph_js(["t", "nb_inner_enzymes", "nb_inflow_pumps", "nb_outflow_pumps"], dico, PATH+"viewer/src/js/mean_E_composition.js", filename, "mean_E_composition_div", "Enzyme tuple classes", "Simulation time", "Number of tuples by class", stacked_option, "mean_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF ENHANCERS OR OPERATORS #
  write_dygraph_js(["t", "nb_functional_regions", "nb_enhancers", "nb_operators"], dico, PATH+"viewer/src/js/mean_nb_enhancers_operators.js", filename, "mean_nb_enhancers_operators_div", "Number of enhancers and operators", "Simulation time", "Number of sites", filled_option, "mean_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SIZE OF ENHANCERS OR OPERATORS #
  write_dygraph_js(["t", "enhancer_size", "operator_size"], dico, PATH+"viewer/src/js/mean_size_enhancers_operators.js", filename, "mean_size_enhancers_operators_div", "Size of enhancers and operators", "Simulation time", "Size", filled_option, "mean_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # TYPES OF REGIONS #
  write_dygraph_js(["t", "nb_E_regions", "nb_TF_regions", "nb_mixed_regions"], dico, PATH+"viewer/src/js/mean_region_types.js", filename, "mean_region_types_div", "Types of functional regions", "Simulation time", "Number by type", stacked_option, "mean_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SIZE OF REGIONS BY TYPE #
  write_dygraph_js(["t", "E_region_size", "TF_region_size", "mixed_region_size"], dico, PATH+"viewer/src/js/mean_region_size.js", filename, "mean_region_size_div", "Size of functional regions", "Simulation time", "Mean size by type", filled_option, "mean_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SIZE OF OPERONS BY TYPE #
  write_dygraph_js(["t", "E_operon_size", "TF_operon_size", "mixed_operon_size"], dico, PATH+"viewer/src/js/mean_operon_size.js", filename, "mean_operon_size_div", "Size of operons", "Simulation time", "Mean size by type", filled_option, "mean_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 5) GENERATE BEST GENOME STRUCTURE FIGURES             #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  filename       = PATH+"statistics/genome_structure_best.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = False

  counter = 1

  # GENOME SIZE #
  write_dygraph_js(["t", "genome_size", "functional_size"], dico, PATH+"viewer/src/js/best_genome_size.js", filename, "best_genome_size_div", "Genome size", "Simulation time", "Number of tuples", filled_option, "best_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF PEARLS #
  write_dygraph_js(["t", "nb_NC", "nb_E", "nb_TF", "nb_BS", "nb_P"], dico, PATH+"viewer/src/js/best_nb_pearls.js", filename, "best_nb_pearls_div", "Genome composition", "Simulation time", "Number of tuples by type", stacked_option, "best_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # ENZYME TUPLE CLASSES #
  write_dygraph_js(["t", "nb_inner_enzymes", "nb_inflow_pumps", "nb_outflow_pumps"], dico, PATH+"viewer/src/js/best_E_composition.js", filename, "best_E_composition_div", "Enzyme tuple classes", "Simulation time", "Number of tuples by class", stacked_option, "best_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF ENHANCERS OR OPERATORS #
  write_dygraph_js(["t", "nb_functional_regions", "nb_enhancers", "nb_operators"], dico, PATH+"viewer/src/js/best_nb_enhancers_operators.js", filename, "best_nb_enhancers_operators_div", "Number of enhancers and operators", "Simulation time", "Number of sites", filled_option, "best_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SIZE OF ENHANCERS OR OPERATORS #
  write_dygraph_js(["t", "enhancer_size", "operator_size"], dico, PATH+"viewer/src/js/best_size_enhancers_operators.js", filename, "best_size_enhancers_operators_div", "Size of enhancers and operators", "Simulation time", "Size", filled_option, "best_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # TYPES OF REGIONS #
  write_dygraph_js(["t", "nb_E_regions", "nb_TF_regions", "nb_mixed_regions"], dico, PATH+"viewer/src/js/best_region_types.js", filename, "best_region_types_div", "Types of functional regions", "Simulation time", "Number by type", stacked_option, "best_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SIZE OF REGIONS BY TYPE #
  write_dygraph_js(["t", "E_region_size", "TF_region_size", "mixed_region_size"], dico, PATH+"viewer/src/js/best_region_size.js", filename, "best_region_size_div", "Size of functional regions", "Simulation time", "Mean size by type", filled_option, "best_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SIZE OF OPERONS BY TYPE #
  write_dygraph_js(["t", "E_operon_size", "TF_operon_size", "mixed_operon_size"], dico, PATH+"viewer/src/js/best_operon_size.js", filename, "best_operon_size_div", "Size of operons", "Simulation time", "Mean size by type", filled_option, "best_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # INHERITED PROTEINS FIGURES ----------------------------------------------------------------------#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 4) GENERATE MEAN INHERITED PROTEINS FIGURES           #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  filename       = PATH+"statistics/inherited_proteins_mean.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = False

  # NUMBER OF INHERITED PROTEINS #
  write_dygraph_js(["t", "nb_E", "nb_TF"], dico, PATH+"viewer/src/js/mean_number_inherited_proteins.js", filename, "mean_number_inherited_proteins_div", "Number of inherited proteins", "Simulation time", "Number of tuples by type", stacked_option, "mean_inherited_proteins_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # ENZYME TUPLE CLASSES #
  write_dygraph_js(["t", "nb_inner_enzymes", "nb_inflow_pumps", "nb_outflow_pumps"], dico, PATH+"viewer/src/js/mean_inherited_E_composition.js", filename, "mean_inherited_E_composition_div", "Inherited enzyme tuple classes", "Simulation time", "Number of tuples by class", stacked_option, "mean_inherited_proteins_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 5) GENERATE BEST INHERITED PROTEINS FIGURES           #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  filename       = PATH+"statistics/inherited_proteins_best.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = False

  counter = 1

  # NUMBER OF INHERITED PROTEINS #
  write_dygraph_js(["t", "nb_E", "nb_TF"], dico, PATH+"viewer/src/js/best_number_inherited_proteins.js", filename, "best_number_inherited_proteins_div", "Number of inherited proteins", "Simulation time", "Number of tuples by type", stacked_option, "best_inherited_proteins_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # ENZYME TUPLE CLASSES #
  write_dygraph_js(["t", "nb_inner_enzymes", "nb_inflow_pumps", "nb_outflow_pumps"], dico, PATH+"viewer/src/js/best_inherited_E_composition.js", filename, "best_inherited_E_composition_div", "Inherited enzyme tuple classes", "Simulation time", "Number of tuples by class", stacked_option, "best_inherited_proteins_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # TROPHIC NETWORK FIGURES -------------------------------------------------------------------------#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 6) GENERATE TROPHIC NETWORK PROFILE FIGURES           #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  filename       = PATH+"statistics/trophic_network_profile.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = False

  counter = 1

  # NUMBER OF GROUPS BY TROPHIC NETWORK LEVEL #
  write_dygraph_js(["t", "level0_nb_groups", "level1_nb_groups", "level2_nb_groups", "nolevel_nb_groups"], dico, PATH+"viewer/src/js/trophic_network_nb_groups.js", filename, "trophic_network_nb_groups_div", "Number of groups per trophic level", "Simulation time", "Number of groups", stacked_option, "trophic_network_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF CELLS BY TROPHIC NETWORK LEVEL #
  write_dygraph_js(["t", "level0_nb_cells", "level1_nb_cells", "level2_nb_cells", "nolevel_nb_cells"], dico, PATH+"viewer/src/js/trophic_network_nb_cells.js", filename, "trophic_network_nb_cells_div", "Number of cells per trophic level", "Simulation time", "Number of cells", stacked_option, "trophic_network_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF GROUP APPEARANCES AND EXTINCTIONS #
  write_dygraph_js(["t", "nb_group_appearances", "nb_group_extinctions"], dico, PATH+"viewer/src/js/trophic_network_group_events.js", filename, "trophic_network_group_events_div", "Number of trophic group new events", "Simulation time", "Number of groups", stacked_option, "trophic_network_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # MEAN GROUP LIFESPAN #
  write_dygraph_js(["t", "mean_group_lifespan"], dico, PATH+"viewer/src/js/trophic_network_group_lifespan.js", filename, "trophic_network_group_lifespan_div", "Mean trophic group lifespan", "Simulation time", "Lifespan", stacked_option, "trophic_network_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # MEAN GROUP TRUE DIVERSITY #
  write_dygraph_js(["t", "mean_true_diversity"], dico, PATH+"viewer/src/js/trophic_network_true_diversity.js", filename, "trophic_network_true_diversity_div", "Mean trophic group true diversity", "Simulation time", "True diversity", stacked_option, "trophic_network_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # MEAN GROUP MAXIMUM TIME DISTANCE #
  write_dygraph_js(["t", "mean_time_distance"], dico, PATH+"viewer/src/js/trophic_network_time_distance.js", filename, "trophic_network_time_distance_div", "Mean trophic group maximum time distance", "Simulation time", "Maximum time distance", stacked_option, "trophic_network_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # TREE STRUCTURE FIGURES -------------------------------------------------------------------------#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 7) GENERATE TREE STRUCTURE FIGURES                    #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  filename       = PATH+"statistics/tree_structure.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = False

  counter = 1

  # NUMBER OF NODES IN THE LINEAGE TREE #
  write_dygraph_js(["t", "lineage_nb_nodes"], dico, PATH+"viewer/src/js/lineage_nb_nodes.js", filename, "lineage_nb_nodes_div", "Number of nodes in the lineage tree", "Simulation time", "Number of nodes", stacked_option, "tree_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF NODES IN THE PHYLOGENETIC TREE #
  write_dygraph_js(["t", "phylogeny_nb_nodes"], dico, PATH+"viewer/src/js/phylogeny_nb_nodes.js", filename, "phylogeny_nb_nodes_div", "Number of nodes in the phylogenetic tree", "Simulation time", "Number of nodes", stacked_option, "tree_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1
  
  # COMMON ANCESTOR AGE #
  write_dygraph_js(["t", "phylogeny_ca_age"], dico, PATH+"viewer/src/js/common_ancestor_age.js", filename, "common_ancestor_age_div", "Common ancestor age", "Simulation time", "Age", stacked_option, "tree_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1
  
  # GLOBAL CONCENTRATIONS FIGURES -------------------------------------------------------------------------#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 8) GENERATE GLOBAL CONCENTRATIONS FIGURES             #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  filename       = PATH+"statistics/global_concentrations.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = False

  counter = 1

  # REPARTITION OF THE MATTER IN THE SYSTEM #
  write_dygraph_js(["t", "matter_pop", "matter_env"], dico, PATH+"viewer/src/js/matter_repartition.js", filename, "matter_repartition_div", "Repartition of matter in the world", "Simulation time", "Amounts", stacked_option, "global_concentrations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # FLUXES BILAN #
  write_dygraph_js(["t", "matter_inflow", "matter_outflow"], dico, PATH+"viewer/src/js/fluxes_bilan.js", filename, "fluxes_bilan_div", "Bilan of matter fluxes", "Simulation time", "Amounts", filled_option, "global_concentrations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1


