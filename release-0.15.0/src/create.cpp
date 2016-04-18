
/**
 * \file      create.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Create a simulation
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

const std::string EXECUTABLE_NAME  = "build/bin/create";
const std::string DEFAULT_FILENAME = "parameters.txt";

void readArgs( int argc, char const** argv, std::string& optional_filename, bool& random_seed );
void printUsage( void );
void printHeader( void );
void create_folders( char const** argv );


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
  bool        random_seed       = false;
  readArgs(argc, argv, optional_filename, random_seed);
  
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
  /* 3) create folders                */
  /*----------------------------------*/
  create_folders(argv);
  
  /*----------------------------------*/
  /* 4) initialize PRNG seed          */
  /*----------------------------------*/
  if (random_seed)
  {
    parameters->get_prng()->set_seed((unsigned long int)time(NULL));
    parameters->set_seed((unsigned long int)parameters->get_prng()->uniform(1, MAXIMUM_SEED));
    parameters->get_prng()->set_seed(parameters->get_seed());
    parameters->write("parameters.txt");
  }
  else
  {
    parameters->get_prng()->set_seed((unsigned long int)parameters->get_seed());
  }
  
  /*----------------------------------*/
  /* 5) load and initialize           */
  /*    simulation                    */
  /*----------------------------------*/
  printf("> Read parameters file and initialize simulation ...\n");
  Simulation* simulation = new Simulation(parameters);
  simulation->initialize(NEW_EXPERIMENT);
  
  /*----------------------------------*/
  /* 6) save experiment               */
  /*----------------------------------*/
  printf("> Save simulation ...\n");
  simulation->save_experiment();
  
  /*----------------------------------*/
  /* 7) save trees                    */
  /*----------------------------------*/
  printf("> Save trees ...\n");
  simulation->save_trees();
  
  /*----------------------------------*/
  /* 8) free memory                   */
  /*----------------------------------*/
  delete simulation;
  simulation = NULL;
  delete parameters;
  parameters = NULL;
  printf("\nSimulation creation successful.\n");
  printf("To run the simulation, do 'run -b 0 -t simulation-time [options]'.\n\n");
  return EXIT_SUCCESS;
}


/**
 * \brief    Read command line arguments
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \param    string& optional_file
 * \param    bool& random_seed
 * \return   \e void
 */
void readArgs( int argc, char const** argv, std::string& optional_filename, bool& random_seed )
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
    if (strcmp(argv[i], "-rs") == 0 || strcmp(argv[i], "--random-seed") == 0)
    {
      random_seed = true;
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
  std::cout << "Usage: create -h or --help\n";
  std::cout << "   or: create [-f param-file] [options]\n";
  std::cout << "Options are:\n";
  std::cout << "  -h, --help\n";
  std::cout << "        print this help, then exit\n";
  std::cout << "  -v, --version\n";
  std::cout << "        print the current version, then exit\n";
  std::cout << "  -f, --file\n";
  std::cout << "        specify parameters file (default: parameters.txt)\n\n";
  std::cout << "  -rs, --random-seed\n";
  std::cout << "        specify if the prng seed should be at random\n\n";
  std::cout << "Be aware that creating a simulation in a folder erases previous simulation.\n";
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
