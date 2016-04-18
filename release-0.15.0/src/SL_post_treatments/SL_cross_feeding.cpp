
/**
 * \file      SL_cross_feeding.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      18-03-2016
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Cross-feeding interactions analyzis
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

void add_metabolite( std::unordered_map<int, std::vector<unsigned long long int>>& map, int met, unsigned long long int cell_id );
void add_cells( std::unordered_map<unsigned long long int, char>& map, std::vector<unsigned long long int>& list );

/**
 * \brief    Main function
 * \details  --
 * \param    int argc
 * \param    char const** argv
 * \return   \e int
 */
int main( int argc, char const** argv )
{
  /*---------------------------------*/
  /* 1) Load parameters from file    */
  /*---------------------------------*/
  Parameters* parameters = new Parameters();
  bool load_successful = parameters->load_parameters_from_file("parameters.txt");
  if (!load_successful)
  {
    std::cout << "Error during parameters loading.\n";
    exit(EXIT_FAILURE);
  }
  
  /*---------------------------------*/
  /* 2) Load the evolved population  */
  /*---------------------------------*/
  Population* evolved_population = NULL;
  gzFile pop_file    = gzopen("./population/population_500000", "r");
  evolved_population = new Population(parameters, pop_file);
  gzclose(pop_file);
  
  /*---------------------------------*/
  /* 3) Load the evolved environment */
  /*---------------------------------*/
  Environment* evolved_environment = NULL;
  gzFile env_file     = gzopen("./environment/environment_500000", "r");
  evolved_environment = new Environment(parameters, env_file);
  gzclose(env_file);
  
  /*---------------------------------*/
  /* 4) Run the post-treatment       */
  /*---------------------------------*/
  
  /*****************************************************************************/
  /* 4.1) Get all the pumped or produced metabolites of ecotypes A and B or AB */
  /*****************************************************************************/
  
  std::unordered_map<int, std::vector<unsigned long long int>> A_pumped_map;
  std::unordered_map<int, std::vector<unsigned long long int>> B_pumped_map;
  
  std::unordered_map<int, std::vector<unsigned long long int>> A_produced_map;
  std::unordered_map<int, std::vector<unsigned long long int>> B_produced_map;
  
  double nbA = 0.0;
  double nbB = 0.0;
  double nbN = 0.0;
  for (size_t x = 0; x < parameters->get_width(); x++)
  {
    for (size_t y = 0; y < parameters->get_height(); y++)
    {
      Cell* cell = evolved_population->get_cell(x, y);
      
      /* If the cell is alive and active */
      if (cell->isActive() && cell->isAlive())
      {
        nbN += 1.0;
        if ((cell->get_trophic_level() == LEVEL_0 || cell->get_trophic_level() == LEVEL_1))
        {
          nbA += 1.0;
        }
        else if ((cell->get_trophic_level() == LEVEL_2 || cell->get_trophic_level() == NO_LEVEL))
        {
          nbB += 1.0;
        }
        cell->load_genome_in_ODE_system(evolved_environment, false, true);
        std::vector<int> pumped;
        std::vector<int> produced;
        cell->get_pumped_and_produced_metabolites(pumped, produced);
        for (size_t i = 0; i < pumped.size(); i++)
        {
          int met = pumped[i];
          if ((cell->get_trophic_level() == LEVEL_0 || cell->get_trophic_level() == LEVEL_1))
          {
            add_metabolite(A_pumped_map, met, cell->get_id());
          }
          else if ((cell->get_trophic_level() == LEVEL_2 || cell->get_trophic_level() == NO_LEVEL))
          {
            add_metabolite(B_pumped_map, met, cell->get_id());
          }
        }
        for (size_t i = 0; i < produced.size(); i++)
        {
          int met = produced[i];
          if ((cell->get_trophic_level() == LEVEL_0 || cell->get_trophic_level() == LEVEL_1))
          {
            add_metabolite(A_produced_map, met, cell->get_id());
          }
          else if ((cell->get_trophic_level() == LEVEL_2 || cell->get_trophic_level() == NO_LEVEL))
          {
            add_metabolite(B_produced_map, met, cell->get_id());
          }
        }
      }
    }
  }
  
  /*****************************************************************************/
  /* 4.2) Classify pumps depending on production and pumping in A or B ecotype */
  /*****************************************************************************/
  
  /* List of produced metabolites */
  std::unordered_map<int, char> A_no;
  std::unordered_map<int, char> B_no;
  std::unordered_map<int, char> AB_no;
  
  std::unordered_map<int, char> A_A;
  std::unordered_map<int, char> B_B;
  
  std::unordered_map<int, char> A_B;
  std::unordered_map<int, char> A_AB;
  std::unordered_map<int, char> AB_B;
  
  std::unordered_map<int, char> B_A;
  std::unordered_map<int, char> B_AB;
  std::unordered_map<int, char> AB_A;
  
  std::unordered_map<int, char> AB_AB;
  
  /* List of producer cells */
  std::unordered_map<unsigned long long int, char> A_no_producers;
  std::unordered_map<unsigned long long int, char> B_no_producers;
  std::unordered_map<unsigned long long int, char> AB_no_producers;
  
  std::unordered_map<unsigned long long int, char> A_A_producers;
  std::unordered_map<unsigned long long int, char> B_B_producers;
  
  std::unordered_map<unsigned long long int, char> A_B_producers;
  std::unordered_map<unsigned long long int, char> A_AB_producers;
  std::unordered_map<unsigned long long int, char> AB_B_producers;
  
  std::unordered_map<unsigned long long int, char> B_A_producers;
  std::unordered_map<unsigned long long int, char> B_AB_producers;
  std::unordered_map<unsigned long long int, char> AB_A_producers;
  
  std::unordered_map<unsigned long long int, char> AB_AB_producers;
  
  /* List of consumer cells */
  std::unordered_map<unsigned long long int, char> A_A_consumers;
  std::unordered_map<unsigned long long int, char> B_B_consumers;
  
  std::unordered_map<unsigned long long int, char> A_B_consumers;
  std::unordered_map<unsigned long long int, char> A_AB_consumers;
  std::unordered_map<unsigned long long int, char> AB_B_consumers;
  
  std::unordered_map<unsigned long long int, char> B_A_consumers;
  std::unordered_map<unsigned long long int, char> B_AB_consumers;
  std::unordered_map<unsigned long long int, char> AB_A_consumers;
  
  std::unordered_map<unsigned long long int, char> AB_AB_consumers;
  
  /*****************************************************************************/
  /* 4.3) Explore metabolites produced by A                                    */
  /*****************************************************************************/
  
  for (std::unordered_map<int, std::vector<unsigned long long int>>::iterator it = A_produced_map.begin(); it != A_produced_map.end(); ++it)
  {
    int met = it->first;
    
    /*----------------------------------------*/
    /* A) If met is exclusively produced by A */
    /*----------------------------------------*/
    if (B_produced_map.find(met) == B_produced_map.end())
    {
      if (A_pumped_map.find(met) == A_pumped_map.end() && B_pumped_map.find(met) == B_pumped_map.end())
      {
        A_no[met] = 1;
        add_cells(A_no_producers, it->second);
      }
      else if (A_pumped_map.find(met) != A_pumped_map.end() && B_pumped_map.find(met) == B_pumped_map.end())
      {
        A_A[met] = 1;
        add_cells(A_A_producers, it->second);
        add_cells(A_A_consumers, A_pumped_map[met]);
      }
      else if (A_pumped_map.find(met) == A_pumped_map.end() && B_pumped_map.find(met) != B_pumped_map.end())
      {
        A_B[met] = 1;
        add_cells(A_B_producers, it->second);
        add_cells(A_B_consumers, B_pumped_map[met]);
      }
      else if (A_pumped_map.find(met) != A_pumped_map.end() && B_pumped_map.find(met) != B_pumped_map.end())
      {
        A_AB[met] = 1;
        add_cells(A_AB_producers, it->second);
        add_cells(A_AB_consumers, A_pumped_map[met]);
        add_cells(A_AB_consumers, B_pumped_map[met]);
      }
    }
    
    /*----------------------------------------*/
    /* B) If met is produced by A and by B    */
    /*----------------------------------------*/
    if (B_produced_map.find(met) != B_produced_map.end())
    {
      if (A_pumped_map.find(met) == A_pumped_map.end() && B_pumped_map.find(met) == B_pumped_map.end())
      {
        AB_no[met] = 1;
        add_cells(AB_no_producers, it->second);
        add_cells(AB_no_producers, B_produced_map[met]);
      }
      else if (A_pumped_map.find(met) != A_pumped_map.end() && B_pumped_map.find(met) == B_pumped_map.end())
      {
        AB_A[met] = 1;
        add_cells(AB_A_producers, it->second);
        add_cells(AB_A_producers, B_produced_map[met]);
        add_cells(AB_A_consumers, A_pumped_map[met]);
      }
      else if (A_pumped_map.find(met) == A_pumped_map.end() && B_pumped_map.find(met) != B_pumped_map.end())
      {
        AB_B[met] = 1;
        add_cells(AB_B_producers, it->second);
        add_cells(AB_B_producers, B_produced_map[met]);
        add_cells(AB_B_consumers, B_pumped_map[met]);
      }
      else if (A_pumped_map.find(met) != A_pumped_map.end() && B_pumped_map.find(met) != B_pumped_map.end())
      {
        AB_AB[met] = 1;
        add_cells(AB_AB_producers, it->second);
        add_cells(AB_AB_producers, B_produced_map[met]);
        add_cells(AB_AB_consumers, A_pumped_map[met]);
        add_cells(AB_AB_consumers, B_pumped_map[met]);
      }
    }
  }
  
  /*************************************/
  /* Explore metabolites produced by B */
  /*************************************/
  
  for (std::unordered_map<int, std::vector<unsigned long long int>>::iterator it = B_produced_map.begin(); it != B_produced_map.end(); ++it)
  {
    int met = it->first;
    
    /*----------------------------------------*/
    /* A) If met is exclusively produced by B */
    /*----------------------------------------*/
    if (A_produced_map.find(met) == A_produced_map.end())
    {
      if (B_pumped_map.find(met) == B_pumped_map.end() && A_pumped_map.find(met) == A_pumped_map.end())
      {
        B_no[met] = 1;
        add_cells(B_no_producers, it->second);
      }
      else if (B_pumped_map.find(met) != B_pumped_map.end() && A_pumped_map.find(met) == A_pumped_map.end())
      {
        B_B[met] = 1;
        add_cells(B_B_producers, it->second);
        add_cells(B_B_consumers, B_pumped_map[met]);
      }
      else if (B_pumped_map.find(met) == B_pumped_map.end() && A_pumped_map.find(met) != A_pumped_map.end())
      {
        B_A[met] = 1;
        add_cells(B_A_producers, it->second);
        add_cells(B_A_consumers, A_pumped_map[met]);
      }
      else if (B_pumped_map.find(met) != B_pumped_map.end() && A_pumped_map.find(met) != A_pumped_map.end())
      {
        B_AB[met] = 1;
        add_cells(B_AB_producers, it->second);
        add_cells(B_AB_consumers, A_pumped_map[met]);
        add_cells(B_AB_consumers, B_pumped_map[met]);
      }
    }
    
    /*----------------------------------------*/
    /* B) If met is produced by B and by A    */
    /*----------------------------------------*/
    if (A_produced_map.find(met) != A_produced_map.end())
    {
      if (B_pumped_map.find(met) == B_pumped_map.end() && A_pumped_map.find(met) == A_pumped_map.end())
      {
        AB_no[met] = 1;
        add_cells(AB_no_producers, it->second);
        add_cells(AB_no_producers, A_produced_map[met]);
      }
      else if (B_pumped_map.find(met) != B_pumped_map.end() && A_pumped_map.find(met) == A_pumped_map.end())
      {
        AB_B[met] = 1;
        add_cells(AB_B_producers, it->second);
        add_cells(AB_B_producers, A_produced_map[met]);
        add_cells(AB_B_consumers, B_pumped_map[met]);
      }
      else if (B_pumped_map.find(met) == B_pumped_map.end() && A_pumped_map.find(met) != A_pumped_map.end())
      {
        AB_A[met] = 1;
        add_cells(AB_A_producers, it->second);
        add_cells(AB_A_producers, A_produced_map[met]);
        add_cells(AB_A_consumers, A_pumped_map[met]);
      }
      else if (B_pumped_map.find(met) != B_pumped_map.end() && A_pumped_map.find(met) != A_pumped_map.end())
      {
        AB_AB[met] = 1;
        add_cells(AB_AB_producers, it->second);
        add_cells(AB_AB_producers, A_produced_map[met]);
        add_cells(AB_AB_consumers, B_pumped_map[met]);
        add_cells(AB_AB_consumers, A_pumped_map[met]);
      }
    }
  }
  
  std::ofstream results("cross_feeding_produced_met.txt", std::ios::out | std::ios::trunc);
  results << "A_no B_no AB_no A_A B_B A_B A_AB AB_B B_A B_AB AB_A AB_AB\n";
  results << A_no.size() << " ";
  results << B_no.size() << " ";
  results << AB_no.size() << " ";
  results << A_A.size() << " ";
  results << B_B.size() << " ";
  results << A_B.size() << " ";
  results << A_AB.size() << " ";
  results << AB_B.size() << " ";
  results << B_A.size() << " ";
  results << B_AB.size() << " ";
  results << AB_A.size() << " ";
  results << AB_AB.size() << "\n";
  results.close();
  
  results.open("cross_feeding_nb_producers.txt", std::ios::out | std::ios::trunc);
  results << "A_no B_no AB_no A_A B_B A_B A_AB AB_B B_A B_AB AB_A AB_AB\n";
  results << (double)A_no_producers.size()/nbA << " ";
  results << (double)B_no_producers.size()/nbB << " ";
  results << (double)AB_no_producers.size()/nbN << " ";
  results << (double)A_A_producers.size()/nbA << " ";
  results << (double)B_B_producers.size()/nbB << " ";
  results << (double)A_B_producers.size()/nbA << " ";
  results << (double)A_AB_producers.size()/nbA << " ";
  results << (double)AB_B_producers.size()/nbN << " ";
  results << (double)B_A_producers.size()/nbB << " ";
  results << (double)B_AB_producers.size()/nbB << " ";
  results << (double)AB_A_producers.size()/nbN << " ";
  results << (double)AB_AB_producers.size()/nbN << "\n";
  results.close();
  
  results.open("cross_feeding_nb_consumers.txt", std::ios::out | std::ios::trunc);
  results << "A_no B_no AB_no A_A B_B A_B A_AB AB_B B_A B_AB AB_A AB_AB\n";
  results << 0.0 << " ";
  results << 0.0 << " ";
  results << 0.0 << " ";
  results << (double)A_A_consumers.size()/nbA << " ";
  results << (double)B_B_consumers.size()/nbB << " ";
  results << (double)A_B_consumers.size()/nbA << " ";
  results << (double)A_AB_consumers.size()/nbA << " ";
  results << (double)AB_B_consumers.size()/nbN << " ";
  results << (double)B_A_consumers.size()/nbB << " ";
  results << (double)B_AB_consumers.size()/nbB << " ";
  results << (double)AB_A_consumers.size()/nbN << " ";
  results << (double)AB_AB_consumers.size()/nbN << "\n";
  results.close();
  
  /*---------------------------------*/
  /* 5) Free the memory              */
  /*---------------------------------*/
  delete evolved_population;
  evolved_population = NULL;
  delete evolved_environment;
  evolved_environment = NULL;
  delete parameters;
  parameters = NULL;
  
  return EXIT_SUCCESS;
}


/**
 * \brief    Add metabolite to a map
 * \details  --
 * \param    std::unordered_map<int, std::vector<unsigned long long int>>& map
 * \param    int met
 * \param    unsigned long long int cell_id
 * \return   \e void
 */
void add_metabolite( std::unordered_map<int, std::vector<unsigned long long int>>& map, int met, unsigned long long int cell_id )
{
  if (map.find(met) == map.end())
  {
    map[met] = std::vector<unsigned long long int>();
    map[met].push_back(cell_id);
  }
  else
  {
    map[met].push_back(cell_id);
  }
}

/**
 * \brief    Transfer cell identifiers from vector to map
 * \details  --
 * \param    std::unordered_map<unsigned long long int, char>& map
 * \param    std::vector<unsigned long long int>& list
 * \return   \e void
 */
void add_cells( std::unordered_map<unsigned long long int, char>& map, std::vector<unsigned long long int>& list )
{
  for (size_t i = 0; i < list.size(); i++)
  {
    map[list[i]] = 1;
  }
}
