
/**
 * \file      run.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Run a simulation
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

#include "./lib/Macros.h"
#include "./lib/Parameters.h"
#include "./lib/Simulation.h"

#if WITH_GRAPHICS_CONTEXT
  #include <SFML/Graphics.hpp>
  #include "./lib/GraphicDisplay.h"
#endif

const std::string EXECUTABLE_NAME         = "build/bin/run";
const size_t      DEFAULT_BACKUP_TIME     = 0;
const size_t      DEFAULT_SIMULATION_TIME = 10000;

void readArgs( int argc, char const** argv, size_t& backup_time, size_t& simulation_time, bool& activate_graphic_display );
void printUsage( void );
void printHeader( void );
void write_simulation_data( Simulation* simulation );
void generate_figures( char const** argv, Simulation* simulation );
#if WITH_GRAPHICS_CONTEXT
  void create_graphic_context( Parameters* parameters, sf::RenderWindow* pop_window, sf::RenderWindow* env_window );
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
  /* 1) Read command line arguments   */
  /*----------------------------------*/
  size_t  backup_time           = 0;
  size_t  simulation_time       = 0;
  bool activate_graphic_display = false;
  readArgs(argc, argv, backup_time, simulation_time, activate_graphic_display);
  if (backup_time == 0)
  {
    backup_time = DEFAULT_BACKUP_TIME;
  }
  if (simulation_time == 0)
  {
    simulation_time = DEFAULT_SIMULATION_TIME;
  }
  
  /*----------------------------------*/
  /* 2) Load parameters from backup   */
  /*----------------------------------*/
  Parameters* parameters = new Parameters(backup_time);
  parameters->write("parameters.out");
  
  /*----------------------------------*/
  /* 3) Load simulation from backup   */
  /*    and graphic display           */
  /*----------------------------------*/
  Simulation* simulation = new Simulation(parameters, backup_time, true);
#if WITH_GRAPHICS_CONTEXT
  GraphicDisplay* graphic_display = NULL;
  sf::RenderWindow pop_window;
  sf::RenderWindow env_window;
  std::string font_path = std::string(argv[0]).substr(0, std::string(argv[0]).size()-EXECUTABLE_NAME.size()) + "src/lib/fonts/andale_mono.ttf";
#endif
  
  /*----------------------------------*/
  /* 4) Open render windows           */
  /*----------------------------------*/
#if WITH_GRAPHICS_CONTEXT
  if (activate_graphic_display)
  {
    create_graphic_context(parameters, &pop_window, &env_window);
    graphic_display = new GraphicDisplay(parameters, simulation, &pop_window, &env_window, font_path);
  }
#endif
  
  /*----------------------------------*/
  /* 5) Print header                  */
  /*----------------------------------*/
  printHeader();
  
  /*----------------------------------*/
  /* 6) Initialize simulation         */
  /*----------------------------------*/
  simulation->initialize( FROM_BACKUP );
  write_simulation_data(simulation);
  if (parameters->get_figures_generation_step() > 0)
  {
    generate_figures(argv, simulation);
  }
  
  /*----------------------------------*/
  /* 7) Run simulation                */
  /*----------------------------------*/
  double TIME             = backup_time + simulation_time;
  bool   run              = true;
  size_t last_exp_backup  = 0;
  size_t last_tree_backup = 0;
  
  while (simulation->get_population()->get_time() < TIME && run)
  {
    /*---------------------------------------------------*/
    /* 7.1) Save best cell state in track.txt            */
    /*---------------------------------------------------*/
    Cell* best_cell = simulation->get_population()->get_cell(0, 0);
    if (best_cell->get_birth_time() == simulation->get_population()->get_time()-1)
    {
      std::ofstream tracker_file("track.txt", std::ios::out | std::ios::trunc);
      best_cell->write_state_header(tracker_file);
      best_cell->write_current_state(tracker_file, simulation->get_population()->get_time()-1, simulation->get_environment());
      tracker_file.close();
    }
    else if (best_cell->get_birth_time() < simulation->get_population()->get_time()-1)
    {
      std::ofstream tracker_file("track.txt", std::ios::out | std::ios::app);
      best_cell->write_current_state(tracker_file, simulation->get_population()->get_time()-1, simulation->get_environment());
      tracker_file.close();
    }
    
    /*---------------------------------------------------*/
    /* 7.2) Update simulation                            */
    /*---------------------------------------------------*/
    simulation->update();
    
    /*---------------------------------------------------*/
    /* 7.3) Save experiment and/or trees                 */
    /*---------------------------------------------------*/
    if (parameters->get_experiment_backup_step() > 0 && simulation->get_population()->get_time() % parameters->get_experiment_backup_step() == 0)
    {
      simulation->save_experiment();
      last_exp_backup = simulation->get_population()->get_time();
    }
    if (parameters->get_tree_backup_step() > 0 && simulation->get_population()->get_time() % parameters->get_tree_backup_step() == 0)
    {
      simulation->save_trees();
      simulation->get_statistics()->write_last_backup_file(last_exp_backup, last_tree_backup);
      last_tree_backup = simulation->get_population()->get_time();
    }
    
    /*---------------------------------------------------*/
    /* 7.4) Display graphics                             */
    /*---------------------------------------------------*/
#if WITH_GRAPHICS_CONTEXT
    if (activate_graphic_display)
    {
      run = graphic_display->display();
    }
#endif
    
    /*---------------------------------------------------*/
    /* 7.5) Stop if popsize = 0, else write stats        */
    /*---------------------------------------------------*/
    if (simulation->get_population()->get_population_size() == 0)
    {
      run = false;
    }
    else
    {
      simulation->get_statistics()->write_stats();
    }
    
    /*---------------------------------------------------*/
    /* 7.6) Print simulation report and generate figures */
    /*---------------------------------------------------*/
    size_t step = 500;
    if (parameters->get_figures_generation_step() > 0 && parameters->get_figures_generation_step() < 500)
    {
      step = parameters->get_figures_generation_step();
    }
    if (simulation->get_population()->get_time()%step == 0 && step != parameters->get_figures_generation_step())
    {
      std::cout << "-------------- Simulation report ----------------\n";
      std::cout << " Elapsed time        : " << simulation->get_population()->get_time() << "\n";
      std::cout << " Elapsed generations : " << simulation->get_statistics()->_mean_generations << "\n";
      std::cout << " Population size     : " << simulation->get_population()->get_population_size() << "\n";
      std::cout << " Best score          : " << simulation->get_max_score() << "\n";
      std::cout << "-------------------------------------------------\n";
    }
    else if (simulation->get_population()->get_time()%step == 0 && step == parameters->get_figures_generation_step())
    {
      std::cout << "-------------- Simulation report ----------------\n";
      std::cout << " Elapsed time        : " << simulation->get_population()->get_time() << "\n";
      std::cout << " Elapsed generations : " << simulation->get_statistics()->_mean_generations << "\n";
      std::cout << " Population size     : " << simulation->get_population()->get_population_size() << "\n";
      std::cout << " Best score          : " << simulation->get_max_score() << "\n";
      std::cout << "-------------------------------------------------\n";
      simulation->get_statistics()->flush_files();
      write_simulation_data(simulation);
      generate_figures(argv, simulation);
    }
    
    /*---------------------------------------------------*/
    /* 7.7) Write signal file                            */
    /*---------------------------------------------------*/
    std::ofstream file("last_simulation_timestep.txt", std::ios::out | std::ios::trunc);
    file << simulation->get_population()->get_time() << "\n";
    file.close();
  }
  
  /*----------------------------------*/
  /* 8) final output                  */
  /*----------------------------------*/
  if (simulation->get_population()->get_population_size() > 0)
  {
    std::cout << "-------------- End of simulation ----------------\n";
    std::cout << "                                                 \n";
    std::cout << " End of simulation at " << simulation->get_population()->get_time() << " timesteps.\n";
    std::cout << "                                                 \n";
    std::cout << "-------------------------------------------------\n\n";
    if (parameters->get_figures_generation_step() > 0)
    {
      simulation->get_statistics()->flush_files();
      write_simulation_data(simulation);
      generate_figures(argv, simulation);
    }
  }
  else
  {
    std::cout << "-------------- End of simulation ----------------\n";
    std::cout << "                                                 \n";
    std::cout << " Population extinction at " << simulation->get_population()->get_time() << " timesteps.\n";
    std::cout << "                                                 \n";
    std::cout << "-------------------------------------------------\n\n";
    if (parameters->get_figures_generation_step() > 0)
    {
      simulation->get_statistics()->flush_files();
      write_simulation_data(simulation);
      generate_figures(argv, simulation);
    }
  }
  
  /*----------------------------------*/
  /* 9) close statistics files and    */
  /*    free memory                   */
  /*----------------------------------*/
  simulation->get_statistics()->close_files();
  delete simulation;
  simulation = NULL;
#if WITH_GRAPHICS_CONTEXT
  if (activate_graphic_display)
  {
    delete graphic_display;
    graphic_display = NULL;
  }
#endif
  delete parameters;
  parameters = NULL;
  
  return EXIT_SUCCESS;
}


/**
 * \brief    Read arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    size_t& backup_time
 * \param    size_t& simulation_time
 * \param    bool& activate_graphic_display
 * \return   \e void
 */
void readArgs( int argc, char const** argv, size_t& backup_time, size_t& simulation_time, bool& activate_graphic_display )
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
    if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--backup-time") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        backup_time = atoi(argv[i+1]);
      }
    }
    if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--simulation-time") == 0)
    {
      if (i+1 == argc)
      {
        std::cout << "Error: command line parameter value is missing.\n";
        exit(EXIT_FAILURE);
      }
      else
      {
        simulation_time = atoi(argv[i+1]);
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
  std::cout << "Usage: run -h or --help\n";
  std::cout << "   or: run [-b backup-time] [-t simulation-time] [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -b, --backup-time\n";
  std::cout << "        set the date of the backup to load (default: 0)\n";
  std::cout << "  -t, --simulation-time\n";
  std::cout << "        set the duration of the simulation (default: 10000)\n";
  std::cout << "  -g, --graphics\n";
  std::cout << "        activate graphic display\n\n";
  std::cout << "Statistic files content is automatically managed when a simulation is reloaded from backup to avoid data loss.\n";
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
 * \brief    Write simulation data
 * \details  --
 * \param    Simulation* simulation
 * \return   \e void
 */
void write_simulation_data( Simulation* simulation )
{
  bool extincted = simulation->get_population()->get_population_size() == 0;
  
  /*-----------------------------------------------------*/
  /* A) generate figures related to population evolution */
  /*-----------------------------------------------------*/
  if (!extincted)
  {
    simulation->get_statistics()->write_best_genome_file();
    simulation->get_statistics()->write_best_inherited_proteins_file();
    simulation->get_statistics()->write_best_genetic_regulation_network_file();
    simulation->get_statistics()->write_best_metabolic_network_file();
    simulation->get_statistics()->write_best_metabolic_amounts_file();
    simulation->get_statistics()->write_last_environment_metabolic_amounts_file();
  }
  
  /*-----------------------------------------------------*/
  /* B) generate figures related to best lineage         */
  /*-----------------------------------------------------*/
  if (!extincted && simulation->get_parameters()->get_tree_backup_step() > 0 && simulation->get_lineage_tree()->get_number_of_nodes() > 0)
  {
    simulation->get_statistics()->write_last_lineage_tree_statistics();
  }
  
  /*-----------------------------------------------------*/
  /* C) generate figures related to the phylogeny        */
  /*-----------------------------------------------------*/
  if (!extincted && simulation->get_parameters()->get_tree_backup_step() > 0 && simulation->get_phylogenetic_tree()->get_number_of_nodes() > 0)
  {
    simulation->get_statistics()->write_last_phylogenetic_tree_statistics();
  }
  
  /*-----------------------------------------------------*/
  /* D) generate figures related to the trophic network  */
  /*-----------------------------------------------------*/
  if (!extincted)
  {
    simulation->get_statistics()->write_last_trophic_network_file();
  }
}

/**
 * \brief    Generate figures
 * \details  --
 * \param    char const** argv
 * \param    Simulation* simulation
 * \return   \e void
 */
void generate_figures( char const** argv, Simulation* simulation )
{
  bool extincted = simulation->get_population()->get_population_size() == 0;
  
  /*-----------------------------------------------------*/
  /* A) generate figures related to population evolution */
  /*-----------------------------------------------------*/
  if (!extincted)
  {
    simulation->get_statistics()->plot_population_figures(std::string(argv[0]), EXECUTABLE_NAME);
  }
  
  /*-----------------------------------------------------*/
  /* B) generate figures related to best lineage         */
  /*-----------------------------------------------------*/
  if (!extincted && simulation->get_parameters()->get_tree_backup_step() > 0 && simulation->get_lineage_tree()->get_number_of_nodes() > 0)
  {
    simulation->get_statistics()->plot_lineage_tree_figures(std::string(argv[0]), EXECUTABLE_NAME);
  }
  
  /*-----------------------------------------------------*/
  /* C) generate figures related to the phylogeny        */
  /*-----------------------------------------------------*/
  if (!extincted && simulation->get_parameters()->get_tree_backup_step() > 0 && simulation->get_phylogenetic_tree()->get_number_of_nodes() > 0)
  {
    simulation->get_statistics()->plot_phylogenetic_tree_figures(std::string(argv[0]), EXECUTABLE_NAME);
  }
  
  /*-----------------------------------------------------*/
  /* D) generate figures related to the trophic network  */
  /*-----------------------------------------------------*/
  if (!extincted)
  {
    simulation->get_statistics()->plot_trophic_network_figures(std::string(argv[0]), EXECUTABLE_NAME);
  }
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
