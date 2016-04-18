
/**
 * \file      SL_heatmap.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      18-03-2016
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Compute the heatmap of SL state after 10*333 timesteps
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
#include <cmath>
#include <sys/stat.h>
#include <assert.h>

#include "../lib/Macros.h"
#include "../lib/Parameters.h"
#include "../lib/Population.h"
#include "../lib/Environment.h"
#include "../lib/Tree.h"
#include "../lib/Simulation.h"

const std::string EXECUTABLE_NAME         = "build/bin/SL_heatmap";
const std::string DEFAULT_FILENAME        = "parameters.txt";
const std::string DEFAULT_POPULATION_PATH = "population_to_load";

enum final_state
{
  EXTINCTION  = 0, /*!< Whole population extinct  */
  COEXISTENCE = 1, /*!< S and L coexist           */
  EXCLUSION   = 2  /*!< L fixed in the population */
};

void readArgs( int argc, char const** argv, size_t& rep, std::string& optional_filename, std::string& optional_population_path );
void printUsage( void );
void printHeader( void );
void create_folders( void );
void replace_population( Simulation* simulation, Parameters* parameters, Population* evolved_population );
Simulation* create_experiment( Parameters* parameters, Population* evolved_population, unsigned long int seed );
bool run_experiment( Simulation* simulation, size_t simulation_time );
void measure_SL_state( Parameters* parameters, Population* evolved_population, double exoConc, double envFreq, size_t REPS, final_state* states );


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
  size_t      rep                      = 0;
  std::string optional_filename        = "";
  std::string optional_population_path = "";
  readArgs(argc, argv, rep, optional_filename, optional_population_path);
  
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
  /* 4) Run the post-treatment      */
  /*--------------------------------*/
  std::ofstream results("heatmap.txt", std::ios::out | std::ios::trunc);
  results << "exo freq exclusion coexistence extinction\n";
  
  for (double exoConc = 0.0; exoConc <= 50.0; exoConc += 1.0)
  {
    for (double envFreq = -4.0; envFreq <= 0.0; envFreq += 0.5)
    {
      std::cout << "Running sim for exoConc=" << exoConc << " and envFreq=" << envFreq << " ...\n";
      final_state* states = new final_state[rep];
      measure_SL_state(parameters, evolved_population, exoConc, pow(10.0, envFreq), rep, states);
      double exclusion   = 0.0;
      double coexistence = 0.0;
      double extinction  = 0.0;
      for (size_t i = 0; i < rep; i++)
      {
        if (states[i] == EXCLUSION)
        {
          exclusion += 1.0;
        }
        else if (states[i] == COEXISTENCE)
        {
          coexistence += 1.0;
        }
        else if (states[i] == EXTINCTION)
        {
          extinction += 1.0;
        }
      }
      exclusion   /= rep;
      coexistence /= rep;
      extinction  /= rep;
      results << exoConc << " " << envFreq << " " << exclusion << " " << coexistence << " " << extinction << "\n";
      results.flush();
      delete[] states;
      states = NULL;
    }
  }
  results.close();
  
  /*--------------------------------*/
  /* 5) Free the memory             */
  /*--------------------------------*/
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
 * \return   \e void
 */
void readArgs( int argc, char const** argv, size_t& rep, std::string& optional_filename, std::string& optional_population_path )
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
    if (strcmp(argv[i], "-rep") == 0 || strcmp(argv[i], "--rep") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        rep = (size_t)atoi(argv[i+1]);
      }
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
  std::cout << "Usage: SL_heatmap -h or --help\n";
  std::cout << "   or: SL_heatmap [options]\n";
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
  std::cout << "        specify the path of the population backup to load (default: population_to_load)\n\n";
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
 * \param    void
 * \return   \e void
 */
void create_folders( void )
{
  /*------------------------------*/
  /* 1) create simulation folders */
  /*------------------------------*/
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
}

/**
 * \brief    Replace the individuals of the simulation by evolved individuals in the right proportions
 * \details  --
 * \param    Simulation* simulation
 * \param    Parameters* parameters
 * \param    Population* evolved_population
 * \return   \e void
 */
void replace_population( Simulation* simulation, Parameters* parameters, Population* evolved_population )
{
  if (parameters->get_width() != evolved_population->get_width() || parameters->get_height() != evolved_population->get_height())
  {
    printf("Error: evolved population grid size is different from parameters file grid size. Exit.\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < parameters->get_width()*parameters->get_height(); i++)
  {
    simulation->get_population()->get_cell(i)->replace_data(evolved_population->get_cell(i));
  }
  return;
}

/**
 * \brief    Create an experiment
 * \details  --
 * \param    Parameters* parameters
 * \param    Population* evolved_population
 * \param    unsigned long int seed
 * \return   \e Simulation*
 */
Simulation* create_experiment( Parameters* parameters, Population* evolved_population, unsigned long int seed )
{
  /*----------------------------------*/
  /* 1) create simulation folders     */
  /*----------------------------------*/
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
  system("rm -rf statistics");
  mkdir("statistics", 0777);
  
  /*----------------------------------*/
  /* 2) initialize PRNG seed          */
  /*----------------------------------*/
  parameters->set_seed(seed);
  parameters->get_prng()->set_seed(seed);
  
  /*----------------------------------*/
  /* 3) create a new simulation       */
  /*----------------------------------*/
  Simulation* simulation = new Simulation(parameters);
  
  /*----------------------------------*/
  /* 4) replace the population        */
  /*----------------------------------*/
  replace_population(simulation, parameters, evolved_population);
  
  /*----------------------------------*/
  /* 5) initialize the population     */
  /*----------------------------------*/
  simulation->initialize_from_evolved_population();
  
  return simulation;
}

/**
 * \brief    Run an experiment
 * \details  --
 * \param    Simulation* simulation
 * \param    size_t simulation_time
 * \return   \e bool
 */
bool run_experiment( Simulation* simulation, size_t simulation_time )
{
  /*----------------------------------*/
  /* 1) run simulation                */
  /*----------------------------------*/
  double TIME    = simulation_time;
  bool   run     = true;
  size_t step    = 1;
  bool   extinct = false;
  while (simulation->get_population()->get_time() < TIME && run)
  {
    /*--------------------------------------------*/
    /* 1.1) update simulation                     */
    /*--------------------------------------------*/
    simulation->update();
    
    /*--------------------------------------------*/
    /* 1.2) stop if popsize = 0, else write stats */
    /*--------------------------------------------*/
    if (simulation->get_population()->get_population_size() == 0)
    {
      extinct = true;
      run     = false;
    }
    else
    {
      simulation->get_statistics()->write_stats();
      simulation->get_statistics()->flush_files();
    }
    
    /*--------------------------------------------*/
    /* 1.3) update timestep                       */
    /*--------------------------------------------*/
    step++;
  }
  
  /*----------------------------------*/
  /* 2) close statistics files and    */
  /*    free memory                   */
  /*----------------------------------*/
  simulation->get_statistics()->close_files();
  
  return extinct;
}

/**
 * \brief    Measure the final state of the population with exogenous concentrartion "exoConc" and environment frequency "envFreq"
 * \details  --
 * \param    Parameters* parameter
 * \param    Population* evolved_population
 * \param    double exoConc
 * \param    double envFreq
 * \param    size_t REPS
 * \param    final_state* states
 * \return   \e void
 */
void measure_SL_state( Parameters* parameters, Population* evolved_population, double exoConc, double envFreq, size_t REPS, final_state* states )
{
  Prng* prng = new Prng();
  prng->set_seed((unsigned long int)time(NULL));
  
  std::cout << envFreq << "\n";
  
  parameters->get_environment_properties()->concentration_range.min = exoConc;
  parameters->get_environment_properties()->concentration_range.max = exoConc;
  parameters->get_environment_properties()->introduction_rate       = envFreq;
  
  size_t rep = 0;
  while (rep < REPS)
  {
    std::cout << "Rep " << rep << " ...\n";
    Simulation* simulation = create_experiment(parameters, evolved_population, (unsigned long int)prng->uniform(1, 10000000));
    replace_population(simulation, parameters, evolved_population);
    bool extinct = run_experiment(simulation, 10*333);
    
    if (extinct)
    {
      states[rep] = EXTINCTION;
    }
    else
    {
      size_t level0   = simulation->get_trophic_network()->get_nb_level_0_cells();
      size_t level1   = simulation->get_trophic_network()->get_nb_level_1_cells();
      size_t level2   = simulation->get_trophic_network()->get_nb_level_2_cells();
      size_t no_level = simulation->get_trophic_network()->get_nb_no_level_cells();
      if ((level0 > 0 || level1 > 0) && (level2 == 0 && no_level == 0))
      {
        states[rep] = EXCLUSION;
      }
      else if ((level0 == 0 && level1 == 0) && (level2 == 0 && no_level == 0))
      {
        printf("aïïïïïïï !\n");
      }
      else if ((level0 > 0 || level1 > 0) && (level2 > 0 || no_level > 0))
      {
        states[rep] = COEXISTENCE;
      }
    }
    
    delete simulation;
    simulation = NULL;
    
    rep++;
  }
  delete prng;
  prng = NULL;
}
