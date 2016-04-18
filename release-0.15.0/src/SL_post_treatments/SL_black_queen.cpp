
/**
 * \file      SL_black_queen.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      27-03-2016
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Measure the black queen effect on A/B ecotypes
 */

/****************************************************************************
 * Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * E-mail: charles.rocabert@inria.fr
 * Web: http://www.evoevo.eu/
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#include "../cmake/Config.h"

#include <iostream>
#include <sstream>
#include <cstring>
#include <sys/stat.h>
#include <assert.h>

#include "../lib/Macros.h"
#include "../lib/Parameters.h"
#include "../lib/Population.h"
#include "../lib/Tree.h"
#include "../lib/Node.h"

const std::string EXECUTABLE_NAME                = "build/bin/SL_black_queen";
const std::string DEFAULT_FILENAME               = "parameters.txt";
const std::string DEFAULT_POPULATION_PATH        = "population_to_load";
const std::string DEFAULT_LINEAGE_TREE_PATH      = "lineage_tree_to_load";
const std::string DEFAULT_PHYLOGENETIC_TREE_PATH = "phylogenetic_tree_to_load";

void readArgs( int argc, char const** argv, std::string& optional_filename, std::string& optional_population_path, std::string& optional_lineage_tree_path, std::string& optional_phylogenetic_tree_path );
void printUsage( void );
void printHeader( void );


/**
 * \brief    Main function
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \return   \e int
 */
int main( int argc, char const** argv )
{
  /*---------------------------------------*/
  /* 1) Read command line arguments        */
  /*---------------------------------------*/
  std::string optional_filename               = "";
  std::string optional_population_path        = "";
  std::string optional_lineage_tree_path      = "";
  std::string optional_phylogenetic_tree_path = "";
  readArgs(argc, argv, optional_filename, optional_population_path, optional_lineage_tree_path, optional_phylogenetic_tree_path);
  
  /*---------------------------------------*/
  /* 2) Load parameters from file          */
  /*---------------------------------------*/
  printHeader();
  Parameters* parameters = new Parameters();
  bool load_successful   = false;
  if (strcmp(optional_filename.c_str(), "") != 0)
  {
    load_successful = parameters->load_parameters_from_file(optional_filename);
  }
  else
  {
    load_successful = parameters->load_parameters_from_file(DEFAULT_FILENAME);
  }
  if (!load_successful)
  {
    std::cout << "Error during parameters loading.\n";
    exit(EXIT_FAILURE);
  }
  
  /*---------------------------------------*/
  /* 3) Load the evolved population        */
  /*---------------------------------------*/
  Population* evolved_population = NULL;
  if (strcmp(optional_population_path.c_str(), "") != 0)
  {
    gzFile pop_file    = gzopen(optional_population_path.c_str(), "r");
    evolved_population = new Population(parameters, pop_file);
    gzclose(pop_file);
  }
  else
  {
    gzFile pop_file    = gzopen(DEFAULT_POPULATION_PATH.c_str(), "r");
    evolved_population = new Population(parameters, pop_file);
    gzclose(pop_file);
  }
  
  /*---------------------------------------*/
  /* 4) Load the evolved lineage tree      */
  /*---------------------------------------*/
  Tree* evolved_lineage = NULL;
  if (strcmp(optional_lineage_tree_path.c_str(), "") != 0)
  {
    gzFile lineage_file = gzopen(optional_lineage_tree_path.c_str(), "r");
    evolved_lineage     = new Tree(parameters, evolved_population, lineage_file);
    gzclose(lineage_file);
  }
  else
  {
    gzFile lineage_file = gzopen(DEFAULT_LINEAGE_TREE_PATH.c_str(), "r");
    evolved_lineage     = new Tree(parameters, evolved_population, lineage_file);
    gzclose(lineage_file);
  }
  
  /*---------------------------------------*/
  /* 5) Load the evolved phylogenetic tree */
  /*---------------------------------------*/
  Tree* evolved_phylogeny = NULL;
  if (strcmp(optional_phylogenetic_tree_path.c_str(), "") != 0)
  {
    gzFile phylogeny_file = gzopen(optional_phylogenetic_tree_path.c_str(), "r");
    evolved_phylogeny     = new Tree(parameters, evolved_population, phylogeny_file);
    gzclose(phylogeny_file);
  }
  else
  {
    gzFile phylogeny_file = gzopen(DEFAULT_PHYLOGENETIC_TREE_PATH.c_str(), "r");
    evolved_phylogeny     = new Tree(parameters, evolved_population, phylogeny_file);
    gzclose(phylogeny_file);
  }
  
  /*---------------------------------------*/
  /* 4) Run the post-treatment             */
  /*---------------------------------------*/
  
  /*********************************/
  /* 4.1) Get the master root      */
  /*********************************/
  Node* master_root = evolved_lineage->get_node(0);
  
  /*********************************/
  /* 4.2) Get the common ancestor  */
  /*********************************/
  std::ofstream output("rootedge.txt", std::ios::out | std::ios::trunc);
  Node* CA = master_root;
  while (CA->get_number_of_children() == 1)
  {
    if (CA->get_id() != 0)
    {
      output << CA->get_replication_report()->get_birth_time() << " ";
      output << CA->get_replication_report()->get_generation() << " ";
      output << CA->get_replication_report()->get_trophic_level() << "\n";
    }
    CA = CA->get_child(0);
  }
  if (CA->get_id() != 0)
  {
    output << CA->get_replication_report()->get_birth_time() << " ";
    output << CA->get_replication_report()->get_generation() << " ";
    output << CA->get_replication_report()->get_trophic_level() << "\n";
  }
  output.close();
  
  /*********************************/
  /* 4.3) Get nascent A/B subtrees */
  /*********************************/
  Node* subA1 = CA->get_child(0);
  Node* subA2 = CA->get_child(1);
  
  /*********************************/
  /* 4.4) Climb subbranches        */
  /*********************************/
  output.open("subbranch1.txt", std::ios::out | std::ios::trunc);
  while (subA1->get_number_of_children() == 1)
  {
    output << subA1->get_replication_report()->get_birth_time() << " ";
    output << subA1->get_replication_report()->get_generation() << " ";
    output << subA1->get_replication_report()->get_trophic_level() << "\n";
    subA1 = subA1->get_child(0);
  }
  output.close();
  
  output.open("subbranch2.txt", std::ios::out | std::ios::trunc);
  while (subA2->get_number_of_children() == 1)
  {
    output << subA2->get_replication_report()->get_birth_time() << " ";
    output << subA2->get_replication_report()->get_generation() << " ";
    output << subA2->get_replication_report()->get_trophic_level() << "\n";
    subA2 = subA2->get_child(0);
  }
  output.close();
  
  /*---------------------------------------*/
  /* 5) Free the memory                    */
  /*---------------------------------------*/
  delete evolved_phylogeny;
  evolved_phylogeny = NULL;
  delete evolved_lineage;
  evolved_lineage = NULL;
  delete evolved_population;
  evolved_population = NULL;
  delete parameters;
  parameters = NULL;
  
  return EXIT_SUCCESS;
}


/**
 * \brief    Read command line arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    size_t& rep
 * \param    std::string& optional_filename
 * \param    string& optional_population_path
 * \param    std::string& optional_lineage_tree_path
 * \param    std::string& optional_phylogenetic_tree_path
 * \return   \e void
 */
void readArgs( int argc, char const** argv, std::string& optional_filename, std::string& optional_population_path, std::string& optional_lineage_tree_path, std::string& optional_phylogenetic_tree_path )
{
  for (int i = 0; i < argc; i++)
  {
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
    {
      printUsage();
      exit(EXIT_SUCCESS);
    }
    if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0)
    {
      std::cout << PACKAGE << " (" << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << ")\n";
      exit(EXIT_SUCCESS);
    }
    if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--file") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        optional_filename = argv[i+1];
      }
    }
    if (strcmp(argv[i], "-pop") == 0 || strcmp(argv[i], "--population-path") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        optional_population_path = argv[i+1];
      }
    }
    if (strcmp(argv[i], "-lineage-tree") == 0 || strcmp(argv[i], "--lineage-tree-path") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        optional_lineage_tree_path = argv[i+1];
      }
    }
    if (strcmp(argv[i], "-phylogenetic-tree") == 0 || strcmp(argv[i], "--phylogenetic-tree-path") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        optional_phylogenetic_tree_path = argv[i+1];
      }
    }
  }
}

/**
 * \brief    Print usage
 * \details  --
 * \param    void
 * \return   \e void
 */
void printUsage( void )
{
  std::cout << "\n";
  std::cout << "***************************************************************************\n";
#ifdef DEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( debug )\n";
#endif
#ifdef NDEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( release )\n";
#endif
  std::cout << " This software is dedicated to EvoEvo WP2 models development               \n";
  std::cout << "                                                                           \n";
  std::cout << " Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " E-mail: charles.rocabert@inria.fr                                         \n";
  std::cout << " Web: http://www.evoevo.eu/                                                \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "Usage: SL_black_queen -h or --help\n";
  std::cout << "   or: SL_black_queen [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -rep, --rep\n";
  std::cout << "        specify the number of repetitions (mandatory)\n";
  std::cout << "  -f, --file\n";
  std::cout << "        specify parameters file (default: parameters.txt)\n";
  std::cout << "  -pop, --population-path\n";
  std::cout << "        specify the path of the population backup to load (default: population_to_load)\n";
  std::cout << "  -lineage-tree, --lineage-tree-path\n";
  std::cout << "        specify the path of the lineage tree backup to load (default: lineage_tree_to_load)\n";
  std::cout << "  -phylogenetic-tree, --phylogenetic-tree-path\n";
  std::cout << "        specify the path of the phylogenetic tree backup to load (default: phylogenetic_tree_to_load)\n";
  std::cout << "\n";
}

/**
 * \brief    Print header
 * \details  --
 * \param    void
 * \return   \e void
 */
void printHeader( void )
{
  std::cout << "\n";
  std::cout << "***************************************************************************\n";
#ifdef DEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( debug )\n";
#endif
#ifdef NDEBUG
  std::cout << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( release )\n";
#endif
  std::cout << " This software is dedicated to EvoEvo WP2 models development               \n";
  std::cout << "                                                                           \n";
  std::cout << " Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon \n";
  std::cout << " E-mail: charles.rocabert@inria.fr                                         \n";
  std::cout << " Web: http://www.evoevo.eu/                                                \n";
  std::cout << "                                                                           \n";
  std::cout << " This program comes with ABSOLUTELY NO WARRANTY.                           \n";
  std::cout << " This is free software, and you are welcome to redistribute it under       \n";
  std::cout << " certain conditions; See the GNU General Public License for details        \n";
  std::cout << "***************************************************************************\n";
  std::cout << "\n";
}
