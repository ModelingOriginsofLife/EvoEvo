
/**
 * \file      SL_frequency_dependence.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      18-03-2016
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Measure the frequency dependent fitness by competition experiments
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

const std::string EXECUTABLE_NAME         = "build/bin/SL_frequency_dependence";
const std::string DEFAULT_FILENAME        = "parameters.txt";
const std::string DEFAULT_POPULATION_PATH = "population_to_load";

void readArgs( int argc, char const** argv, size_t& rep, std::string& optional_filename, std::string& optional_population_path );
void printUsage( void );
void printHeader( void );
void create_folders( void );
void replace_population( Simulation* simulation, Parameters* parameters, Population* evolved_population, bool identical, double L_prop, unsigned long int seed );
Simulation* create_experiment( Parameters* parameters, Population* evolved_population, unsigned long int seed, bool identical, double L_prop );
bool run_experiment( Simulation* simulation, size_t simulation_time, size_t* nb_1, size_t* nb_2 );
void run_competition_experiment( Prng* prng, Parameters* parameters, Population* evolved_population, size_t REPS, double L_prop );
void measure_frequency_dependent_fitness( Parameters* parameters, Population* evolved_population, size_t REPS );


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
  measure_frequency_dependent_fitness(parameters, evolved_population, rep);
  
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
  std::cout << "Usage: SL_frequency_dependence -h or --help\n";
  std::cout << "   or: SL_frequency_dependence [options]\n";
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
 * \param    bool identical
 * \param    double L_prop
 * \param    unsigned long int seed
 * \param    double& nb_level_1
 * \param    double& nb_level_2
 * \return   \e void
 */
void replace_population( Simulation* simulation, Parameters* parameters, Population* evolved_population, bool identical, double L_prop, unsigned long int seed, double& nb_level_1, double& nb_level_2 )
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
    nb_level_1 = 0.0;
    nb_level_2 = 0.0;
    for (size_t i = 0; i < parameters->get_width()*parameters->get_height(); i++)
    {
      simulation->get_population()->get_cell(i)->replace_data(evolved_population->get_cell(i));
      if (evolved_population->get_cell(i)->get_trophic_level() == LEVEL_0 || evolved_population->get_cell(i)->get_trophic_level() == LEVEL_1)
      {
        nb_level_1 += 1.0;
      }
      else if (evolved_population->get_cell(i)->get_trophic_level() == LEVEL_2 || evolved_population->get_cell(i)->get_trophic_level() == NO_LEVEL)
      {
        nb_level_2 += 1.0;
      }
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
  nb_level_1 = 0.0;
  nb_level_2 = 0.0;
  for (size_t i = 0; i < parameters->get_width()*parameters->get_height(); i++)
  {
    if (prng->uniform() < L_prop)
    {
      size_t draw = (size_t)prng->uniform(0, (int)L_indexes.size()-1);
      simulation->get_population()->get_cell(i)->replace_data(evolved_population->get_cell(L_indexes[draw]));
      nb_level_1 += 1.0;
    }
    else
    {
      size_t draw = (size_t)prng->uniform(0, (int)S_indexes.size()-1);
      simulation->get_population()->get_cell(i)->replace_data(evolved_population->get_cell(S_indexes[draw]));
      nb_level_2 += 1.0;
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
 * \param    unsigned long int seed
 * \param    bool identical
 * \param    double L_prop
 * \param    double& nb_level_1
 * \param    double& nb_level_2
 * \return   \e Simulation*
 */
Simulation* create_experiment( Parameters* parameters, Population* evolved_population, unsigned long int seed, bool identical, double L_prop, double& nb_level_1, double& nb_level_2 )
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
  /* 4) load the evolved population   */
  /*----------------------------------*/
  replace_population(simulation, parameters, evolved_population, identical, L_prop, seed, nb_level_1, nb_level_2);
  
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
 * \param    size_t* nb_1
 * \param    size_t* nb_2
 * \return   \e bool
 */
bool run_experiment( Simulation* simulation, size_t simulation_time, size_t* nb_1, size_t* nb_2 )
{
  /*----------------------------------*/
  /* 1) run simulation                */
  /*----------------------------------*/
  double TIME    = simulation->get_population()->get_time()+simulation_time;
  bool   run     = true;
  bool   extinct = false;
  while (simulation->get_population()->get_time() < TIME && run)
  {
    /*--------------------------------------------*/
    /* 1.1) update simulation                     */
    /*--------------------------------------------*/
    simulation->update();
    
    nb_1[simulation->get_population()->get_time()-1] = simulation->get_trophic_network()->get_nb_level_0_cells()+simulation->get_trophic_network()->get_nb_level_1_cells();
    nb_2[simulation->get_population()->get_time()-1] = simulation->get_trophic_network()->get_nb_level_2_cells()+simulation->get_trophic_network()->get_nb_no_level_cells();
    
    /*--------------------------------------------*/
    /* 1.2) stop if popsize = 0, else write stats */
    /*--------------------------------------------*/
    if (simulation->get_population()->get_population_size() == 0)
    {
      extinct = true;
      run     = false;
    }
  }
  
  /*----------------------------------*/
  /* 2) close statistics files and    */
  /*    free memory                   */
  /*----------------------------------*/
  simulation->get_statistics()->close_files();
  
  return extinct;
}

/**
 * \brief    Run a competition experiment with repetitions
 * \details  --
 * \param    Prng* prng
 * \param    Parameters* parameters
 * \param    Population* evolved_population
 * \param    size_t REPS
 * \param    double Lprop
 * \return   \e void
 */
void run_competition_experiment( Prng* prng, Parameters* parameters, Population* evolved_population, size_t REPS, double L_prop )
{
  std::cout << "Running competition experiment with initial S frequency of " << 1-L_prop << " ...\n";
  
  double* A_0          = new double[REPS];
  double* A_05         = new double[REPS];
  double* A_1          = new double[REPS];
  double* A_2          = new double[REPS];
  double* A_4          = new double[REPS];
  double* B_0          = new double[REPS];
  double* B_05         = new double[REPS];
  double* B_1          = new double[REPS];
  double* B_2          = new double[REPS];
  double* B_4          = new double[REPS];
  size_t* A_extinction = new size_t[REPS];
  size_t* B_extinction = new size_t[REPS];
  size_t** nb_1        = new size_t*[REPS];
  size_t** nb_2        = new size_t*[REPS];
  
  size_t rep = 0;
  while (rep < REPS)
  {
    std::cout << "> Rep " << rep << " ...\n";
    Simulation* simulation = NULL;
    double initial_n_1 = 0.0;
    double initial_n_2 = 0.0;
    nb_1[rep] = new size_t[4*333];
    nb_2[rep] = new size_t[4*333];
    for (size_t i = 0; i < 4*333; ++i)
    {
      nb_1[rep][i] = 0;
      nb_2[rep][i] = 0;
    }
    simulation = create_experiment(parameters, evolved_population, (unsigned long int)prng->uniform(1, 10000000), false, L_prop, initial_n_1, initial_n_2);
    
    A_0[rep] = initial_n_1;
    B_0[rep] = initial_n_2;
    
    A_05[rep] = A_1[rep] = A_2[rep] = A_4[rep] = B_05[rep] = B_1[rep] = B_2[rep] = B_4[rep] = 0.0;
    A_extinction[rep] = B_extinction[rep] = 0;
    
    bool extinct = false;
    
    /*********************/
    /* First half saison */
    /*********************/
    std::cout << "> to time 167 ...\n";
    
    extinct = run_experiment(simulation, 167, nb_1[rep], nb_2[rep]);
    assert(simulation->get_population()->get_time() == 167);
    if (!extinct)
    {
      A_05[rep] = (double)simulation->get_trophic_network()->get_nb_level_0_cells()+(double)simulation->get_trophic_network()->get_nb_level_1_cells();
      B_05[rep] = (double)simulation->get_trophic_network()->get_nb_level_2_cells()+(double)simulation->get_trophic_network()->get_nb_no_level_cells();
    }
    
    /**********************/
    /* Second half saison */
    /**********************/
    std::cout << "> to time 333 ...\n";
    
    if (!extinct)
    {
      extinct = run_experiment(simulation, 166, nb_1[rep], nb_2[rep]);
      assert(simulation->get_population()->get_time() == 333);
      if (!extinct)
      {
        A_1[rep] = (double)simulation->get_trophic_network()->get_nb_level_0_cells()+(double)simulation->get_trophic_network()->get_nb_level_1_cells();
        B_1[rep] = (double)simulation->get_trophic_network()->get_nb_level_2_cells()+(double)simulation->get_trophic_network()->get_nb_no_level_cells();
      }
    }
    
    /*****************/
    /* Second saison */
    /*****************/
    std::cout << "> to time 666 ...\n";
    
    if (!extinct)
    {
      extinct = run_experiment(simulation, 333, nb_1[rep], nb_2[rep]);
      assert(simulation->get_population()->get_time() == 666);
      if (!extinct)
      {
        A_2[rep] = (double)simulation->get_trophic_network()->get_nb_level_0_cells()+(double)simulation->get_trophic_network()->get_nb_level_1_cells();
        B_2[rep] = (double)simulation->get_trophic_network()->get_nb_level_2_cells()+(double)simulation->get_trophic_network()->get_nb_no_level_cells();
      }
    }
    
    /*****************/
    /* Fourth saison */
    /*****************/
    std::cout << "> to time 1332 ...\n";
    
    if (!extinct)
    {
      extinct = run_experiment(simulation, 666, nb_1[rep], nb_2[rep]);
      assert(simulation->get_population()->get_time() == 1332);
      if (!extinct)
      {
        A_4[rep] = (double)simulation->get_trophic_network()->get_nb_level_0_cells()+(double)simulation->get_trophic_network()->get_nb_level_1_cells();
        B_4[rep] = (double)simulation->get_trophic_network()->get_nb_level_2_cells()+(double)simulation->get_trophic_network()->get_nb_no_level_cells();
      }
    }
    
    /********************/
    /* Extinction dates */
    /********************/
    
    if (extinct)
    {
      bool A_extinct = false;
      bool B_extinct = false;
      for (size_t i = 0; i < 4*333; ++i)
      {
        if (!A_extinct && nb_1[rep][i] == 0)
        {
          A_extinction[rep] = i+1;
          A_extinct = true;
        }
        if (!B_extinct && nb_2[rep][i] == 0)
        {
          B_extinction[rep] = i+1;
          B_extinct = true;
        }
      }
    }
    
    delete simulation;
    simulation = NULL;
    
    rep++;
  }
  
  std::stringstream filename;
  filename << "./frequency_dependent_fitness_" << (1-L_prop) << ".txt";
  std::ofstream results(filename.str().c_str(), std::ios::out | std::ios::trunc);
  results << "rep A0 A05 A1 A2 A4 B0 B05 B1 B2 b4 Aext Bext\n";
  for (size_t rep = 0; rep < REPS; rep++)
  {
    results << rep+1 << " " << A_0[rep] << " " << A_05[rep] << " " << A_1[rep] << " " << A_2[rep] << " " << A_4[rep] << " " << B_0[rep] << " " << B_05[rep] << " " << B_1[rep] << " " << B_2[rep] << " " << B_4[rep] << " " << A_extinction[rep] << " " << B_extinction[rep] << "\n";
  }
  results.close();
  
  double* mean_nb_1 = new double[4*333];
  double* mean_nb_2 = new double[4*333];
  double* var_nb_1  = new double[4*333];
  double* var_nb_2  = new double[4*333];
  for (size_t t = 0; t < 4*333; ++t)
  {
    mean_nb_1[t] = 0.0;
    mean_nb_2[t] = 0.0;
    var_nb_1[t]  = 0.0;
    var_nb_2[t]  = 0.0;
    for (size_t rep = 0; rep < REPS; ++rep)
    {
      mean_nb_1[t] += nb_1[rep][t];
      mean_nb_2[t] += nb_2[rep][t];
      var_nb_1[t]  += nb_1[rep][t]*nb_1[rep][t];
      var_nb_2[t]  += nb_2[rep][t]*nb_2[rep][t];
    }
  }
  for (size_t t = 0; t < 4*333; ++t)
  {
    mean_nb_1[t] /= REPS;
    mean_nb_2[t] /= REPS;
    var_nb_1[t]  /= REPS;
    var_nb_2[t]  /= REPS;
    var_nb_1[t]  -= mean_nb_1[t]*mean_nb_1[t];
    var_nb_2[t]  -= mean_nb_2[t]*mean_nb_2[t];
  }
  
  filename.str("");
  filename << "./ecotypes_evolution_" << (1-L_prop) << ".txt";
  results.open(filename.str().c_str(), std::ios::out | std::ios::trunc);
  results << "t mean1 mean2 var1 var2\n";
  for (size_t t = 0; t < 4*333; ++t)
  {
    results << t+1 << " " << mean_nb_1[t] << " " << mean_nb_2[t] << " " << var_nb_1[t] << " " << var_nb_2[t] << "\n";
  }
  results.close();
  
  delete[] A_0;
  A_0 = NULL;
  delete[] A_05;
  A_05 = NULL;
  delete[] A_1;
  A_1 = NULL;
  delete[] A_2;
  A_2 = NULL;
  delete[] A_4;
  A_4 = NULL;
  delete[] B_0;
  B_0 = NULL;
  delete[] B_05;
  B_05 = NULL;
  delete[] B_1;
  B_1 = NULL;
  delete[] B_2;
  B_2 = NULL;
  delete[] B_4;
  B_4 = NULL;
  delete[] A_extinction;
  A_extinction = NULL;
  delete[] B_extinction;
  B_extinction = NULL;
  for (size_t rep = 0; rep < REPS; ++rep)
  {
    delete[] nb_1[rep];
    nb_1[rep] = NULL;
    delete[] nb_2[rep];
    nb_2[rep] = NULL;
  }
  delete nb_1;
  nb_1 = NULL;
  delete nb_2;
  nb_2 = NULL;
}

/**
 * \brief    Measure the frequency-dependent fitness of S/L strains
 * \details  --
 * \param    Parameters* parameter
 * \param    Population* evolved_population
 * \param    size_t REPS
 * \return   \e void
 */
void measure_frequency_dependent_fitness( Parameters* parameters, Population* evolved_population, size_t REPS )
{
  Prng* prng = new Prng();
  prng->set_seed(parameters->get_seed());
  
  /*-------------------------------------------------*/
  /* 1) evaluate S relative fitness at frequency 0.1 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, REPS, 0.9);
  
  /*-------------------------------------------------*/
  /* 2) evaluate S relative fitness at frequency 0.2 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, REPS, 0.8);
  
  /*-------------------------------------------------*/
  /* 3) evaluate S relative fitness at frequency 0.3 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, REPS, 0.7);
  
  /*-------------------------------------------------*/
  /* 4) evaluate S relative fitness at frequency 0.4 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, REPS, 0.6);
  
  /*-------------------------------------------------*/
  /* 5) evaluate S relative fitness at frequency 0.5 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, REPS, 0.5);
  
  /*-------------------------------------------------*/
  /* 6) evaluate S relative fitness at frequency 0.6 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, REPS, 0.4);
  
  /*-------------------------------------------------*/
  /* 7) evaluate S relative fitness at frequency 0.7 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, REPS, 0.3);
  
  /*-------------------------------------------------*/
  /* 8) evaluate S relative fitness at frequency 0.8 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, REPS, 0.2);
  
  /*-------------------------------------------------*/
  /* 9) evaluate S relative fitness at frequency 0.9 */
  /*-------------------------------------------------*/
  run_competition_experiment(prng, parameters, evolved_population, REPS, 0.1);
  
  delete prng;
  prng = NULL;
}
