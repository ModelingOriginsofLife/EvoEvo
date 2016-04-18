
/**
 * \file      SL_create_competition_experiment.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      25-03-2016
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Create a long-term competition experiments with A/B ecotypes
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
#include "../lib/Environment.h"
#include "../lib/Tree.h"
#include "../lib/Simulation.h"

const std::string EXECUTABLE_NAME         = "build/bin/SL_create_competition_experiment";
const std::string DEFAULT_FILENAME        = "parameters.txt";
const std::string DEFAULT_POPULATION_PATH = "population_to_load";

void readArgs( int argc, char const** argv, std::string& optional_filename, std::string& optional_population_path );
void printUsage( void );
void printHeader( void );
void create_folders( char const** argv );
void replace_population( Simulation* simulation, Parameters* parameters, Population* evolved_population, double L_prop, bool identical, unsigned long int seed );
Simulation* create_experiment( Parameters* parameters, Population* evolved_population, double L_prop, bool identical, unsigned long int seed );

/**
 * \brief    Main function
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \return   \e int
 */
int main( int argc, char const** argv )
{
  /*--------------------------------*/
  /* 1) Read command line arguments */
  /*--------------------------------*/
  std::string optional_filename        = "";
  std::string optional_population_path = "";
  readArgs(argc, argv, optional_filename, optional_population_path);
  
  /*--------------------------------*/
  /* 2) Load parameters from file   */
  /*--------------------------------*/
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
  
  /*--------------------------------*/
  /* 3) Load the evolved population */
  /*--------------------------------*/
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
  
  /*--------------------------------*/
  /* 4) Create folders              */
  /*--------------------------------*/
  create_folders(argv);
  
  /*--------------------------------*/
  /* 5) Create experiment           */
  /*--------------------------------*/
  Prng* prng = new Prng();
  prng->set_seed(parameters->get_seed());
  
  Simulation* simulation = create_experiment(parameters, evolved_population, 0.5, false, (unsigned long int)prng->uniform(1, 10000000));
  
  /*--------------------------------*/
  /* 6) save experiment             */
  /*--------------------------------*/
  printf("> Save simulation ...\n");
  simulation->save_experiment();
  
  /*--------------------------------*/
  /* 7) save trees                  */
  /*--------------------------------*/
  printf("> Save trees ...\n");
  simulation->save_trees();
  
  /*--------------------------------*/
  /* 8) Free the memory             */
  /*--------------------------------*/
  delete simulation;
  simulation = NULL;
  delete prng;
  prng = NULL;
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
 * \param    std::string& optional_filename
 * \param    string& optional_population_path
 * \return   \e void
 */
void readArgs( int argc, char const** argv, std::string& optional_filename, std::string& optional_population_path )
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
  std::cout << "Usage: SL_create_competition_experiment -h or --help\n";
  std::cout << "   or: SL_create_competition_experiment [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -f, --file\n";
  std::cout << "        specify parameters file (default: parameters.txt)\n";
  std::cout << "  -pop, --population-path\n";
  std::cout << "        specify the path of the population backup to load (default: population_to_load)\n\n";
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

/**
 * \brief    Create folders
 * \details  --
 * \param    char const** argv
 * \return   \e void
 */
void create_folders( char const** argv )
{
  /*---------------------------------*/
  /* 1) create simulation folders    */
  /*---------------------------------*/
  system("rm -rf statistics");
  mkdir("statistics", 0777);
  system("rm -rf environment");
  mkdir("environment", 0777);
  system("rm -rf population");
  mkdir("population", 0777);
  system("rm -rf trophic_network");
  mkdir("trophic_network", 0777);
  system("rm -rf tree");
  mkdir("tree", 0777);
  system("rm -rf parameters");
  mkdir("parameters", 0777);
  system("rm -rf prng");
  mkdir("prng", 0777);
  system("rm -rf figures");
  mkdir("figures", 0777);
  
  /*---------------------------------*/
  /* 2) copy simulation viewer and   */
  /*    tracker                      */
  /*---------------------------------*/
  system("rm -rf viewer");
  std::string command = "cp -r " + std::string(argv[0]).substr(0, std::string(argv[0]).size()-EXECUTABLE_NAME.size()) + "src/lib/scripts/viewer ./";
  system(command.c_str());
  command = "cp " + std::string(argv[0]).substr(0, std::string(argv[0]).size()-EXECUTABLE_NAME.size()) + "src/lib/scripts/track_cell.py ./";
  system(command.c_str());
  
  /*---------------------------------*/
  /* 3) Write code version in a file */
  /*---------------------------------*/
  std::ofstream f("version.txt", std::ios::out | std::ios::trunc);
#ifdef DEBUG
  f << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( debug )\n";
#endif
#ifdef NDEBUG
  f << " " << PACKAGE << " " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << " ( release )\n";
#endif
  f.close();
}

/**
 * \brief    Replace the individuals of the simulation by evolved individuals in the right proportions
 * \details  --
 * \param    Simulation* simulation
 * \param    Parameters* parameters
 * \param    Population* evolved_population
 * \param    double L_prop
 * \param    bool identical
 * \param    unsigned long int seed
 * \return   \e void
 */
void replace_population( Simulation* simulation, Parameters* parameters, Population* evolved_population, double L_prop, bool identical, unsigned long int seed )
{
  assert(L_prop >= 0.0);
  assert(L_prop <= 1.0);
  if (parameters->get_width() != evolved_population->get_width() || parameters->get_height() != evolved_population->get_height())
  {
    printf("Error: evolved population grid size is different from parameters file grid size. Exit.\n");
    exit(EXIT_FAILURE);
  }
  
  /*--------------------------------------------------------------*/
  /* 1) If identical boolean is true, just replace every cells    */
  /*--------------------------------------------------------------*/
  if (identical)
  {
    for (size_t i = 0; i < parameters->get_width()*parameters->get_height(); i++)
    {
      simulation->get_population()->get_cell(i)->replace_data(evolved_population->get_cell(i));
    }
    return;
  }
  
  /*--------------------------------------------------------------*/
  /* 2) Isolate evolved identifiers thanks to their trophic level */
  /*--------------------------------------------------------------*/
  std::vector<size_t> L_indexes;
  std::vector<size_t> S_indexes;
  for (size_t i = 0; i < parameters->get_width()*parameters->get_height(); i++)
  {
    if (evolved_population->get_cell(i)->isAlive() && evolved_population->get_cell(i)->get_trophic_level() == LEVEL_1)
    {
      L_indexes.push_back(i);
    }
    else if (evolved_population->get_cell(i)->isAlive() && evolved_population->get_cell(i)->get_trophic_level() == LEVEL_2)
    {
      S_indexes.push_back(i);
    }
  }
  
  /*--------------------------------------------------------------*/
  /* 3) Replace cells in the right proportion                     */
  /*--------------------------------------------------------------*/
  Prng* prng = new Prng();
  prng->set_seed(seed);
  for (size_t i = 0; i < parameters->get_width()*parameters->get_height(); i++)
  {
    if (prng->uniform() < L_prop)
    {
      size_t draw = (size_t)prng->uniform(0, (int)L_indexes.size()-1);
      simulation->get_population()->get_cell(i)->replace_data(evolved_population->get_cell(L_indexes[draw]));
    }
    else
    {
      size_t draw = (size_t)prng->uniform(0, (int)S_indexes.size()-1);
      simulation->get_population()->get_cell(i)->replace_data(evolved_population->get_cell(S_indexes[draw]));
    }
  }
  L_indexes.clear();
  S_indexes.clear();
  delete prng;
  prng = NULL;
}

/**
 * \brief    Create an experiment
 * \details  --
 * \param    Parameters* parameters
 * \param    Population* evolved_population
 * \param    double L_prop
 * \param    bool identical
 * \param    unsigned long int seed
 * \return   \e Simulation*
 */
Simulation* create_experiment( Parameters* parameters, Population* evolved_population, double L_prop, bool identical, unsigned long int seed )
{
  /*----------------------------------*/
  /* 1) initialize PRNG seed          */
  /*----------------------------------*/
  parameters->set_seed(seed);
  parameters->get_prng()->set_seed(seed);
  
  /*----------------------------------*/
  /* 2) create a new simulation       */
  /*----------------------------------*/
  Simulation* simulation = new Simulation(parameters);
  
  /*----------------------------------*/
  /* 3) load the evolved population   */
  /*----------------------------------*/
  replace_population(simulation, parameters, evolved_population, L_prop, identical, seed);
  
  /*----------------------------------*/
  /* 4) initialize the population     */
  /*----------------------------------*/
  simulation->initialize_from_evolved_population();
  
  return simulation;
}


