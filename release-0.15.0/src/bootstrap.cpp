
/**
 * \file      bootstrap.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Find a good seed by bootstrap
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
#include <cstring>
#include <sys/stat.h>
#include <assert.h>

#include "./lib/Macros.h"
#include "./lib/Parameters.h"
#include "./lib/Simulation.h"
#include "./lib/Prng.h"

#if WITH_GRAPHICS_CONTEXT
  #include <SFML/Graphics.hpp>
  #include "./lib/GraphicDisplay.h"
#endif

const std::string EXECUTABLE_NAME    = "build/bin/bootstrap";
const std::string DEFAULT_FILENAME   = "parameters.txt";
const size_t DEFAULT_MINIMUM_TIME    = 500;
const size_t DEFAULT_MINIMUM_POPSIZE = 500;
const size_t DEFAULT_TRIALS          = 1000;

void readArgs( int argc, char const** argv, std::string& optional_filename, size_t& minimum_time, size_t& minimum_popsize, size_t& trials, bool& activate_graphic_display );
void printUsage( void );
void printHeader( void );
void create_folders( char const** argv );
void create_experiment( Parameters* parameters, unsigned long int seed );
#if WITH_GRAPHICS_CONTEXT
  void create_graphic_context( Parameters* parameters, sf::RenderWindow* pop_window, sf::RenderWindow* env_window );
  bool run_experiment( Parameters* parameters, size_t minimum_time, size_t minimum_popsize, bool activate_graphic_display, sf::RenderWindow* pop_window, sf::RenderWindow* env_window, std::string font_path );
#else
  bool run_experiment( Parameters* parameters, size_t minimum_time, size_t minimum_popsize );
#endif

/**
 * \brief    Main function
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \return   \e int
 */
int main( int argc, char const** argv )
{
  /*----------------------------------*/
  /* 1) read command line arguments   */
  /*----------------------------------*/
  std::string optional_filename = "";
  size_t minimum_time           = 0;
  size_t minimum_popsize        = 0;
  size_t trials                 = 0;
  bool activate_graphic_display = false;
  readArgs(argc, argv, optional_filename, minimum_time, minimum_popsize, trials, activate_graphic_display);
  if (minimum_time == 0)
  {
    minimum_time = DEFAULT_MINIMUM_TIME;
  }
  if (minimum_popsize == 0)
  {
    minimum_popsize = DEFAULT_MINIMUM_POPSIZE;
  }
  if (trials == 0)
  {
    trials = DEFAULT_TRIALS;
  }
  
  /*----------------------------------*/
  /* 2) load parameters from file     */
  /*----------------------------------*/
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
  
  /*----------------------------------*/
  /* 3) open render windows           */
  /*----------------------------------*/
#if WITH_GRAPHICS_CONTEXT
  sf::RenderWindow pop_window;
  sf::RenderWindow env_window;
  if (activate_graphic_display)
  {
    create_graphic_context(parameters, &pop_window, &env_window);
  }
  std::string font_path = std::string(argv[0]).substr(0, std::string(argv[0]).size()-EXECUTABLE_NAME.size()) + "src/lib/fonts/andale_mono.ttf";
#endif
  
  /*----------------------------------*/
  /* 4) create special prng for seed  */
  /*    drawing                       */
  /*----------------------------------*/
  Prng* seed_prng = new Prng();
  seed_prng->set_seed(parameters->get_seed());
  
  /*----------------------------------*/
  /* 5) explore seeds                 */
  /*----------------------------------*/
  unsigned long int seed  = 0;
  size_t number_of_trials = 1;
  bool success            = false;
  std::cout << "Launch exploration for " << trials << " trials, " << minimum_time << " timesteps to reach with a minimum population size of " << minimum_popsize << " cells ...\n";
  while (!success && number_of_trials <= trials)
  {
    seed = (unsigned long int)seed_prng->uniform(1, MAXIMUM_SEED);
    std::cout << "  > Trial " << number_of_trials << "/" << trials << " with seed " << seed << "\n";
    create_experiment(parameters, seed);
#if WITH_GRAPHICS_CONTEXT
    success = run_experiment(parameters, minimum_time, minimum_popsize, activate_graphic_display, &pop_window, &env_window, font_path);
#else
    success = run_experiment(parameters, minimum_time, minimum_popsize);
#endif
    number_of_trials++;
  }
  
  /*----------------------------------*/
  /* 6) if a seed is found, save      */
  /*    experiment                    */
  /*----------------------------------*/
  if (success)
  {
    std::cout << "Simulation saved with seed " << seed << "\n\n";
    create_folders(argv);
    parameters->set_seed(seed);
    parameters->get_prng()->set_seed(seed);
    Simulation* simulation = new Simulation(parameters);
    simulation->initialize(NEW_EXPERIMENT);
    simulation->save_experiment();
    simulation->save_trees();
    delete simulation;
    simulation = NULL;
    parameters->write("parameters.txt");
  }
  else
  {
    printf("No suitable seed has been found.\n\n");
    system("rm -rf environment");
    system("rm -rf population");
    system("rm -rf tree");
    system("rm -rf parameters");
    system("rm -rf prng");
  }
  
  /*----------------------------------*/
  /* 7) free memory                   */
  /*----------------------------------*/
  delete parameters;
  parameters = NULL;
  delete seed_prng;
  seed_prng = NULL;
  return EXIT_SUCCESS;
}


/**
 * \brief    Read command line arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    string& optional_file
 * \param    int& minimum_time
 * \param    size_t& minimum_popsize
 * \param    size_t& trials
 * \param    bool& activate_graphic_display
 * \return   \e void
 */
void readArgs( int argc, char const** argv, std::string& optional_filename, size_t& minimum_time, size_t& minimum_popsize, size_t& trials, bool& activate_graphic_display )
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
    if (strcmp(argv[i], "-min") == 0 || strcmp(argv[i], "--minimum-time") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        minimum_time = atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-pop") == 0 || strcmp(argv[i], "--minimum-pop-size") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        minimum_popsize = atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--trials") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        trials = atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-g") == 0 || strcmp(argv[i], "--graphics") == 0)
    {
      activate_graphic_display = true;
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
  std::cout << "Usage: bootstrap -h or --help\n";
  std::cout << "   or: bootstrap [-f param-file] [-min minimum-time] [-pop minimum-pop-size] [-t trials] [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -f, --file\n";
  std::cout << "        specify parameters file (default: parameters.txt)\n";
  std::cout << "  -min, --minimum-time\n";
  std::cout << "        specify the minimum time the new population must survive (default: 100)\n";
  std::cout << "  -pop, --minimum-pop-size\n";
  std::cout << "        specify the minimum size the new population must maintain (default: 500)\n";
  std::cout << "  -t, --trials\n";
  std::cout << "        specify the number of trials (default: 1000)\n";
  std::cout << "  -g, --graphics\n";
  std::cout << "        activate graphic display\n";
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
  system("rm -rf tree");
  mkdir("tree", 0777);
  system("rm -rf parameters");
  mkdir("parameters", 0777);
  system("rm -rf prng");
  mkdir("prng", 0777);
  system("rm -rf figures");
  mkdir("figures", 0777);
  
  /*---------------------------------*/
  /* 2) copy simulation viewer       */
  /*---------------------------------*/
  system("rm -rf viewer");
  std::string command = "cp -r " + std::string(argv[0]).substr(0, std::string(argv[0]).size()-EXECUTABLE_NAME.size()) + "src/lib/scripts/viewer ./";
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
 * \brief    Create the graphic context
 * \details  Basically, open graphic windows
 * \param    Parameters* parameters
 * \param    sf::RenderWindow* pop_window
 * \param    sf::RenderWindow* env_window
 * \return   \e void
 */
#if WITH_GRAPHICS_CONTEXT
void create_graphic_context( Parameters* parameters, sf::RenderWindow* pop_window, sf::RenderWindow* env_window )
{
  int x_inter = ((int)parameters->get_width()-1)*CELL_SPACE;
  int y_inter = ((int)parameters->get_height()-1)*CELL_SPACE;
  int width   = (int)parameters->get_width()*CELL_SCALE + x_inter + 4*SPAN + GRADIENT_SCALE + TEXT_SCALE;
  int height  = (int)parameters->get_height()*CELL_SCALE + y_inter + 2*SPAN;
  
  pop_window->create(sf::VideoMode(width, height), "Population");
  env_window->create(sf::VideoMode(width, height), "Environment");
  
  pop_window->setPosition(sf::Vector2i(100, 100));
  env_window->setPosition(sf::Vector2i(100+width+50, 100));
  
  if (FRAMERATE > 0)
  {
    pop_window->setFramerateLimit(FRAMERATE);
    env_window->setFramerateLimit(FRAMERATE);
  }
}
#endif

/**
 * \brief    Create an experiment
 * \details  Create an experiment with the seed given in parameter
 * \param    Parameters* parameters
 * \param    unsigned long int seed
 * \return   \e void
 */
void create_experiment( Parameters* parameters, unsigned long int seed )
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
  
  /*----------------------------------*/
  /* 2) initialize PRNG seed          */
  /*----------------------------------*/
  parameters->set_seed(seed);
  parameters->get_prng()->set_seed(seed);
  
  /*----------------------------------*/
  /* 3) load and initialize           */
  /*    simulation                    */
  /*----------------------------------*/
  Simulation* simulation = new Simulation(parameters);
  simulation->initialize(NEW_EXPERIMENT);
  
  /*----------------------------------*/
  /* 4) save experiment               */
  /*----------------------------------*/
  simulation->save_experiment();
  
  /*----------------------------------*/
  /* 5) save trees                    */
  /*----------------------------------*/
  simulation->save_trees();
  
  /*----------------------------------*/
  /* 6) free memory                   */
  /*----------------------------------*/
  delete simulation;
  simulation = NULL;
}

/**
 * \brief    Run the experiment
 * \details  If the population reached the minimum time, returns true, else false
 * \param    Parameters* parameters
 * \param    size_t minimum_time
 * \param    size_t minimum_popsize
 * \param    bool activate_graphic_display
 * \param    sf::RenderWindow* pop_window
 * \param    sf::RenderWindow* env_window
 * \param    std::string font_path
 * \return   \e bool
 */
#if WITH_GRAPHICS_CONTEXT
bool run_experiment( Parameters* parameters, size_t minimum_time, size_t minimum_popsize, bool activate_graphic_display, sf::RenderWindow* pop_window, sf::RenderWindow* env_window, std::string font_path )
{
  /*----------------------------------*/
  /* 1) load simulation from backup   */
  /*    and graphic display           */
  /*----------------------------------*/
  Simulation* simulation = new Simulation(parameters, 0, true);
  
  /*----------------------------------*/
  /* 2) initialize simulation         */
  /*----------------------------------*/
  simulation->initialize(FROM_BACKUP);
  
  /*----------------------------------*/
  /* 3) initialize graphic display    */
  /*----------------------------------*/
  GraphicDisplay* graphic_display = NULL;
  if (activate_graphic_display)
  {
    graphic_display = new GraphicDisplay(parameters, simulation, pop_window, env_window, font_path);
  }
  
  /*----------------------------------*/
  /* 4) run simulation                */
  /*----------------------------------*/
  double TIME    = minimum_time;
  bool   run     = true;
  size_t step    = 1;
  bool   success = false;
  while (simulation->get_population()->get_time() < TIME && run)
  {
    /*--------------------------------------------*/
    /* 4.1) update simulation                     */
    /*--------------------------------------------*/
    simulation->update();
    
    /*--------------------------------------------*/
    /* 4.2) update simulation                     */
    /*--------------------------------------------*/
    if (activate_graphic_display)
    {
      run = graphic_display->display();
    }
    
    /*--------------------------------------------*/
    /* 4.3) stop if popsize = 0, else write stats */
    /*--------------------------------------------*/
    if (simulation->get_population()->get_population_size() == 0)
    {
      run = false;
    }
    
    /*--------------------------------------------*/
    /* 4.4) update timestep                       */
    /*--------------------------------------------*/
    step++;
  }
  
  /*----------------------------------*/
  /* 5) evaluate seed                 */
  /*----------------------------------*/
  if (simulation->get_population()->get_population_size() >= minimum_popsize)
  {
    success = true;
  }
  
  /*----------------------------------*/
  /* 6) close statistics files and    */
  /*    free memory                   */
  /*----------------------------------*/
  simulation->get_statistics()->close_files();
  delete simulation;
  simulation = NULL;
  if (activate_graphic_display)
  {
    delete graphic_display;
    graphic_display = NULL;
  }
  return success;
}
#else
bool run_experiment( Parameters* parameters, size_t minimum_time, size_t minimum_popsize )
{
  /*----------------------------------*/
  /* 1) load simulation from backup   */
  /*    and graphic display           */
  /*----------------------------------*/
  Simulation* simulation = new Simulation(parameters, 0, true);
  
  /*----------------------------------*/
  /* 2) initialize simulation         */
  /*----------------------------------*/
  simulation->initialize(FROM_BACKUP);
  
  /*----------------------------------*/
  /* 4) run simulation                */
  /*----------------------------------*/
  double TIME    = minimum_time;
  bool   run     = true;
  size_t step    = 1;
  bool   success = false;
  while (simulation->get_population()->get_time() < TIME && run)
  {
    /*--------------------------------------------*/
    /* 4.1) update simulation                     */
    /*--------------------------------------------*/
    simulation->update();
    
    /*--------------------------------------------*/
    /* 4.2) stop if popsize = 0, else write stats */
    /*--------------------------------------------*/
    if (simulation->get_population()->get_population_size() == 0)
    {
      run = false;
    }
    
    /*--------------------------------------------*/
    /* 4.3) update timestep                       */
    /*--------------------------------------------*/
    step++;
  }
  
  /*----------------------------------*/
  /* 5) evaluate seed                 */
  /*----------------------------------*/
  if (simulation->get_population()->get_population_size() >= minimum_popsize)
  {
    success = true;
  }
  
  /*----------------------------------*/
  /* 6) close statistics files and    */
  /*    free memory                   */
  /*----------------------------------*/
  simulation->get_statistics()->close_files();
  delete simulation;
  simulation = NULL;
  return success;
}
#endif
