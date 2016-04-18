
/**
 * \file      SL_run_competition_experiment.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      25-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Run a long-term competition experiments with A/B ecotypes
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
#include "../lib/Simulation.h"

const std::string EXECUTABLE_NAME         = "build/bin/SL_run_competition_experiment";
const size_t      DEFAULT_BACKUP_TIME     = 0;
const size_t      DEFAULT_SIMULATION_TIME = 10000;

void readArgs( int argc, char const** argv, size_t& backup_time, size_t& simulation_time, bool& activate_graphic_display );
void printUsage( void );
void printHeader( void );
void write_simulation_data( Simulation* simulation );
void write_SL_proportion( Simulation* simulation, std::ofstream& flux );


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
  
  /*----------------------------------*/
  /* 4) Print header                  */
  /*----------------------------------*/
  printHeader();
  
  /*----------------------------------*/
  /* 5) Initialize simulation         */
  /*----------------------------------*/
  simulation->initialize( FROM_BACKUP );
  write_simulation_data(simulation);
  std::ofstream AB_file("./statistics/AB_proportion.txt", std::ios::out | std::ios::trunc);
  write_SL_proportion(simulation, AB_file);
  
  /*----------------------------------*/
  /* 6) Run simulation                */
  /*----------------------------------*/
  double TIME             = backup_time + simulation_time;
  bool   run              = true;
  size_t last_exp_backup  = 0;
  size_t last_tree_backup = 0;
  
  while (simulation->get_population()->get_time() < TIME && run)
  {
    /*---------------------------------------------------*/
    /* 6.1) Update simulation                            */
    /*---------------------------------------------------*/
    simulation->update();
    write_SL_proportion(simulation, AB_file);
    
    /*---------------------------------------------------*/
    /* 6.2) Save experiment and/or trees                 */
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
    /* 6.3) Stop if popsize = 0, else write stats        */
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
    /* 6.4) Print simulation report and generate figures */
    /*---------------------------------------------------*/
    size_t step = 500;
    if (simulation->get_population()->get_time()%step == 0)
    {
      simulation->get_statistics()->flush_files();
      write_simulation_data(simulation);
    }
    
    /*---------------------------------------------------*/
    /* 6.5) Write signal file                            */
    /*---------------------------------------------------*/
    std::ofstream file("last_simulation_timestep.txt", std::ios::out | std::ios::trunc);
    file << simulation->get_population()->get_time() << "\n";
    file.close();
  }
  
  /*----------------------------------*/
  /* 7) final output                  */
  /*----------------------------------*/
  simulation->get_statistics()->flush_files();
  write_simulation_data(simulation);
  AB_file.close();
  
  /*----------------------------------*/
  /* 8) close statistics files and    */
  /*    free memory                   */
  /*----------------------------------*/
  simulation->get_statistics()->close_files();
  delete simulation;
  simulation = NULL;
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
  std::cout << "Usage: SL_run_competition_experiment -h or --help\n";
  std::cout << "   or: SL_run_competition_experiment [-b backup-time] [-t simulation-time] [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -b, --backup-time\n";
  std::cout << "        set the date of the backup to load (default: 0)\n";
  std::cout << "  -t, --simulation-time\n";
  std::cout << "        set the duration of the simulation (default: 10000)\n\n";
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
 * \brief    Write the proportion of A/B ecotypes
 * \details  --
 * \param    Simulation* simulation
 * \param    std::ofstream& flux
 * \return   \e void
 */
void write_SL_proportion( Simulation* simulation, std::ofstream& flux )
{
  size_t nbA = simulation->get_trophic_network()->get_nb_level_0_cells()+simulation->get_trophic_network()->get_nb_level_1_cells();;
  size_t nbB = simulation->get_trophic_network()->get_nb_level_2_cells()+simulation->get_trophic_network()->get_nb_no_level_cells();
  flux << nbA << " " << nbB << "\n";
}
