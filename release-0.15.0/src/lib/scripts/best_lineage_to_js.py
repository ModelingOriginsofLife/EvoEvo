
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
  print "=== WRITE JAVASCRIPT FILES FROM BEST LINEAGE ==="
  print "Usage: python best_lineage_to_js.py [parameters]"
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
  f = open(filename, "r")
  l = f.readline()
  l = l.strip("\n")
  l = l.split(" ")
  L = len(l)
  l = f.readline()
  while l:
    l = l.strip("\n")
    l = l.split(" ")
    if len(l) == L:
      if count < PERIOD:
        try:
          mean  += float(l[variable_dict[variable_name]])
        except:
          sys.exit()
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

  #-------------------------------------------------------#
  # 1) DEFINE OPTIONS                                     #
  #-------------------------------------------------------#
  stacked_option = ["stackedGraph: true"]
  filled_option  = ["fillGraph: true"]
  vector_length  = 500


  # PHENOTYPE FIGURES ----------------------------------------------------------------------#

  #-------------------------------------------------------#
  # 2) GENERATE BEST LINEAGE PHENOTYPE FIGURES            #
  #-------------------------------------------------------#
  filename       = PATH+"statistics/phenotype_lineage_best.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = True

  counter = 1

  # PROTEIN AMOUNTS #
  write_dygraph_js(["generation", "inherited_TF_amount", "inherited_E_amount", "TF_amount", "E_amount"], dico, PATH+"viewer/src/js/best_lineage_protein_amounts.js", filename, "best_lineage_protein_amounts_div", "Protein amounts", "Generations", "Concentration", stacked_option, "best_lineage_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # METABOLIC AMOUNTS #
  write_dygraph_js(["generation", "inherited_metabolic_amount", "metabolic_amount"], dico, PATH+"viewer/src/js/best_lineage_metabolic_amounts.js", filename, "best_lineage_metabolic_amounts_div", "Metabolic concentrations", "Generations", "Concentration", stacked_option, "best_lineage_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SCORE #
  write_dygraph_js(["generation", "min_score", "mean_score", "max_score"], dico, PATH+"viewer/src/js/best_lineage_score.js", filename, "best_lineage_score_div", "Score", "Generations", "Score", filled_option, "best_lineage_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # ENERGY #
  write_dygraph_js(["generation", "min_energy", "mean_energy", "max_energy"], dico, PATH+"viewer/src/js/best_lineage_energy.js", filename, "best_lineage_energy_div", "Energy", "Generations", "Energy", filled_option, "best_lineage_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # LIFESPAN #
  write_dygraph_js(["generation", "lifespan"], dico, PATH+"viewer/src/js/best_lineage_lifespan.js", filename, "best_lineage_lifespan_div", "Lifespan", "Generations", "Lifespan (simulation steps)", filled_option, "best_lineage_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF DIVISIONS #
  write_dygraph_js(["generation", "number_of_divisions"], dico, PATH+"viewer/src/js/best_lineage_nb_divisions.js", filename, "best_lineage_nb_divisions_div", "Number of divisions", "Generations", "Number of divisions", filled_option, "best_lineage_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # TOXICITY ACCUMULATION #
  write_dygraph_js(["generation", "toxicity"], dico, PATH+"viewer/src/js/best_lineage_toxicity.js", filename, "best_lineage_toxicity_div", "Toxicity accumulation", "Generations", "Toxicity level", filled_option, "best_lineage_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S EXCHANGES #
  write_dygraph_js(["generation", "metabolic_uptake", "metabolic_release"], dico, PATH+"viewer/src/js/best_lineage_exchanges.js", filename, "best_lineage_exchanges_div", "Exchanges with environment", "Generations", "Concentrations", stacked_option, "best_lineage_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S METABOLIC GROWTH RATE #
  write_dygraph_js(["generation", "metabolic_growth_rate"], dico, PATH+"viewer/src/js/best_lineage_metabolic_growth_rate.js", filename, "best_lineage_metabolic_growth_rate_div", "Metabolic growth rate", "Generations", "Growth rate", filled_option, "best_lineage_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S METABOLIC NODES AND EDGES #
  write_dygraph_js(["generation", "grn_nb_nodes", "grn_nb_edges"], dico, PATH+"viewer/src/js/best_lineage_grn_nodes_edges.js", filename, "best_lineage_grn_nodes_edges_div", "Genetic regulation network structure", "Generations", "Number of elements", filled_option, "best_lineage_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # CELL'S METABOLIC NODES AND EDGES #
  write_dygraph_js(["generation", "metabolic_nb_nodes", "metabolic_nb_edges"], dico, PATH+"viewer/src/js/best_lineage_metabolic_nodes_edges.js", filename, "best_lineage_metabolic_nodes_edges_div", "Metabolic network structure", "Generations", "Number of elements", filled_option, "best_lineage_phenotype_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # GENOME STRUCTURE FIGURES ----------------------------------------------------------------------#

  #-------------------------------------------------------#
  # 3) GENERATE BEST LINEAGE GENOME STRUCTURE FIGURES     #
  #-------------------------------------------------------#
  filename       = PATH+"statistics/genome_structure_lineage_best.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = True

  counter = 1

  # GENOME SIZE #
  write_dygraph_js(["generation", "genome_size", "genome_functional_size"], dico, PATH+"viewer/src/js/best_lineage_genome_size.js", filename, "best_lineage_genome_size_div", "Genome size", "Generations", "Number of tuples", filled_option, "best_lineage_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF PEARLS #
  genome_option = ["stackedGraph: true"]
  write_dygraph_js(["generation", "genome_nb_NC", "genome_nb_E", "genome_nb_TF", "genome_nb_BS", "genome_nb_P"], dico, PATH+"viewer/src/js/best_lineage_nb_pearls.js", filename, "best_lineage_nb_pearls_div", "Genome composition", "Generations", "Number of tuples by type", genome_option, "best_lineage_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # ENZYME PEARLS COMPOSITION #
  write_dygraph_js(["generation", "genome_nb_inner_enzymes", "genome_nb_inflow_pumps", "genome_nb_outflow_pumps"], dico, PATH+"viewer/src/js/best_lineage_E_composition.js", filename, "best_lineage_E_composition_div", "Enzyme tuple classes", "Generations", "Number of pearls by type", stacked_option, "best_lineage_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # REGULATION REDUNDANCY #
  write_dygraph_js(["generation", "mean_regulation_redundancy"], dico, PATH+"viewer/src/js/best_lineage_regulation_redundancy.js", filename, "best_lineage_regulation_redundancy_div", "Regulation redundancy", "Generations", "Number of redundant tuples", filled_option, "best_lineage_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # METABOLIC REDUNDANCY #
  write_dygraph_js(["generation", "mean_metabolic_redundancy"], dico, PATH+"viewer/src/js/best_lineage_metabolic_redundancy.js", filename, "best_lineage_metabolic_redundancy_div", "Metabolic redundancy", "Generations", "Number of redundant tuples", filled_option, "best_lineage_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF ENHANCERS OR OPERATORS #
  write_dygraph_js(["generation", "nb_functional_regions", "nb_enhancers", "nb_operators"], dico, PATH+"viewer/src/js/best_lineage_nb_enhancers_operators.js", filename, "best_lineage_nb_enhancers_operators_div", "Number of enhancers and operators", "Generations", "Number of sites", filled_option, "best_lineage_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SIZE OF ENHANCERS OR OPERATORS #
  write_dygraph_js(["generation", "mean_enhancer_size", "mean_operator_size"], dico, PATH+"viewer/src/js/best_lineage_size_enhancers_operators.js", filename, "best_lineage_size_enhancers_operators_div", "Size of enhancers and operators", "Generations", "Size", filled_option, "best_lineage_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # TYPES OF REGIONS #
  write_dygraph_js(["generation", "nb_E_regions", "nb_TF_regions", "nb_mixed_regions"], dico, PATH+"viewer/src/js/best_lineage_region_types.js", filename, "best_lineage_region_types_div", "Types of functional regions", "Generations", "Number by type", stacked_option, "best_lineage_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SIZE OF REGIONS BY TYPE #
  write_dygraph_js(["generation", "mean_E_region_size", "mean_TF_region_size", "mean_mixed_region_size"], dico, PATH+"viewer/src/js/best_lineage_region_size.js", filename, "best_lineage_region_size_div", "Size of functional regions", "Generations", "Mean size by type", filled_option, "best_lineage_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # SIZE OF OPERONS BY TYPE #
  write_dygraph_js(["generation", "mean_E_operon_size", "mean_TF_operon_size", "mean_mixed_operon_size"], dico, PATH+"viewer/src/js/best_lineage_operon_size.js", filename, "best_lineage_operon_size_div", "Size of operons", "Generations", "Mean size by type", filled_option, "best_lineage_genome_structure_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1
  

  # INHERITED PROTEINS FIGURES ----------------------------------------------------------------------#
  
  #-------------------------------------------------------#
  # 4) GENERATE BEST LINEAGE INHERITED PROTEINS FIGURES   #
  #-------------------------------------------------------#
  filename       = PATH+"statistics/inherited_proteins_lineage_best.txt"
  dico           = get_variable_names(filename)
  data_cumulated = False
  data_reversed  = True

  counter = 1

  # NUMBER OF INHERITED PROTEINS #
  write_dygraph_js(["generation", "inherited_nb_E", "inherited_nb_TF"], dico, PATH+"viewer/src/js/best_lineage_number_inherited_proteins.js", filename, "best_lineage_number_inherited_proteins_div", "Number of inherited proteins", "Generations", "Number of tuples by type", stacked_option, "best_lineage_inherited_proteins_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # ENZYME TUPLE CLASSES #
  write_dygraph_js(["generation", "inherited_nb_inner_enzymes", "inherited_nb_inflow_pumps", "inherited_nb_outflow_pumps"], dico, PATH+"viewer/src/js/best_lineage_inherited_E_composition.js", filename, "best_lineage_inherited_E_composition_div", "Inherited enzyme tuple classes", "Generations", "Number of tuples by class", stacked_option, "best_lineage_inherited_proteins_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1


  # FIXED MUTATIONS FIGURES ----------------------------------------------------------------------#

  #-------------------------------------------------------#
  # 5) GENERATE BEST LINEAGE FIXED MUTATIONS FIGURES      #
  #-------------------------------------------------------#
  filename       = PATH+"statistics/fixed_mutations_lineage_best.txt"
  dico           = get_variable_names(filename)
  data_cumulated = True
  data_reversed  = True

  counter = 1

  # POINT MUTATIONS #
  write_dygraph_js(["generation", "nb_NC_point_mutations", "nb_E_point_mutations", "nb_TF_point_mutations", "nb_BS_point_mutations", "nb_P_point_mutations"], dico, PATH+"viewer/src/js/best_lineage_point_mutations.js", filename, "best_lineage_point_mutations_div", "Point mutations", "Generations", "Number of events", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NC PEARL TRANSITIONS #
  write_dygraph_js(["generation", "nb_NC_to_E_transitions", "nb_NC_to_TF_transitions", "nb_NC_to_BS_transitions", "nb_NC_to_P_transitions"], dico, PATH+"viewer/src/js/best_lineage_NC_transitions.js", filename, "best_lineage_NC_transitions_div", "NC tuples transitions", "Generations", "Number of events", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # E PEARL TRANSITIONS #
  write_dygraph_js(["generation", "nb_E_to_NC_transitions", "nb_E_to_TF_transitions", "nb_E_to_BS_transitions", "nb_E_to_P_transitions"], dico, PATH+"viewer/src/js/best_lineage_E_transitions.js", filename, "best_lineage_E_transitions_div", "E tuples transitions", "Generations", "Number of events", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # TF PEARL TRANSITIONS #
  write_dygraph_js(["generation", "nb_TF_to_NC_transitions", "nb_TF_to_E_transitions", "nb_TF_to_BS_transitions", "nb_TF_to_P_transitions"], dico, PATH+"viewer/src/js/best_lineage_TF_transitions.js", filename, "best_lineage_TF_transitions_div", "TF tuples transitions", "Generations", "Number of events", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # BS PEARL TRANSITIONS #
  write_dygraph_js(["generation", "nb_BS_to_NC_transitions", "nb_BS_to_E_transitions", "nb_BS_to_TF_transitions", "nb_BS_to_P_transitions"], dico, PATH+"viewer/src/js/best_lineage_BS_transitions.js", filename, "best_lineage_BS_transitions_div", "BS tuples transitions", "Generations", "Number of events", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # P PEARL TRANSITIONS #
  write_dygraph_js(["generation", "nb_P_to_NC_transitions", "nb_P_to_E_transitions", "nb_P_to_TF_transitions", "nb_P_to_BS_transitions"], dico, PATH+"viewer/src/js/best_lineage_P_transitions.js", filename, "best_lineage_P_transitions_div", "P tuples transitions", "Generations", "Number of events", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # METABOLIC SPACE MUTATION SIZE #
  write_dygraph_js(["generation", "s_mutation_size", "p_mutation_size"], dico, PATH+"viewer/src/js/best_lineage_metabolic_space_mutation_size.js", filename, "best_lineage_metabolic_space_mutation_size_div", "Metabolic space mutation size", "Generations", "Mutation size", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # RATE CONSTANT MUTATION SIZE #
  write_dygraph_js(["generation", "km_mutation_size", "kcat_mutation_size", "BS_tag_mutation_size", "coE_tag_mutation_size", "beta_mutation_size"], dico, PATH+"viewer/src/js/best_lineage_rate_constant_mutation_size.js", filename, "best_lineage_rate_constant_mutation_size_div", "Rate constants mutation size", "Generations", "Mutation size", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # DUPLICATIONS #
  write_dygraph_js(["generation", "nb_duplicated_NC", "nb_duplicated_E", "nb_duplicated_TF", "nb_duplicated_BS", "nb_duplicated_P"], dico, PATH+"viewer/src/js/best_lineage_duplications.js", filename, "best_lineage_duplications_div", "Duplications", "Generations", "Number of events", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # DELETIONS #
  write_dygraph_js(["generation", "nb_deleted_NC", "nb_deleted_E", "nb_deleted_TF", "nb_deleted_BS", "nb_deleted_P"], dico, PATH+"viewer/src/js/best_lineage_deletions.js", filename, "best_lineage_deletions_div", "Deletions", "Generations", "Number of events", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF REARRANGEMENTS #
  write_dygraph_js(["generation", "nb_duplications", "nb_deletions", "nb_translocations", "nb_inversions"], dico, PATH+"viewer/src/js/best_lineage_nb_rearrangements.js", filename, "best_lineage_nb_rearrangements_div", "Number of rearrangements", "Generations", "Number of events", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1
  
  # REARRANGEMENTS SIZE #
  write_dygraph_js(["generation", "duplication_size", "deletion_size", "translocation_size", "inversion_size"], dico, PATH+"viewer/src/js/best_lineage_rearrangement_size.js", filename, "best_lineage_rearrangement_size_div", "Rearrangements size", "Generations", "Rearrangement size", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

  # NUMBER OF HORIZONTAL GENE TRANSFERS #
  write_dygraph_js(["generation", "nb_NC_HGT", "nb_E_HGT", "nb_TF_HGT", "nb_BS_HGT", "nb_P_HGT"], dico, PATH+"viewer/src/js/best_lineage_nb_hgt.js", filename, "best_lineage_nb_hgt_div", "Number of HGT", "Generations", "Number of HGT", stacked_option, "best_lineage_fixed_mutations_"+str(counter), vector_length, data_cumulated, data_reversed)
  counter += 1

