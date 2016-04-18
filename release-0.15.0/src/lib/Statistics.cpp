
/**
 * \file      Statistics.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Statistics class definition
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

#include "Statistics.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \param    Population* population
 * \param    Environment* environment
 * \param    TrophicNetwork* trophic_network
 * \param    Tree* lineage_tree
 * \param    Tree* phylogenetic_tree
 * \param    bool clean_statistic_files
 * \return   \e void
 */
Statistics::Statistics( Parameters* parameters, Population* population, Environment* environment, TrophicNetwork* trophic_network, Tree* lineage_tree, Tree* phylogenetic_tree, bool clean_statistic_files )
{
  /*------------------------------------------------------------------ simulation variables */
  
  _parameters        = parameters;
  _population        = population;
  _environment       = environment;
  _trophic_network   = trophic_network;
  _lineage_tree      = lineage_tree;
  _phylogenetic_tree = phylogenetic_tree;
  
  /*------------------------------------------------------------------ open files */
  
  if (clean_statistic_files)
  {
    if (_population->get_population_size() > 0)
    {
      clean_files();
    }
    open_files(clean_statistic_files);
  }
  else
  {
    open_files(clean_statistic_files);
  }
  
  /*------------------------------------------------------------------ initialize variables */
  
  init_variables();
  
}

/*----------------------------
 * DESTRUCTORS
 *----------------------------*/

/**
 * \brief    Destructor
 * \details  --
 * \param    void
 * \return   \e void
 */
Statistics::~Statistics( void )
{
  /* NOTHING TO DELETE */
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    CLean files
 * \details  Remove lines older than population age
 * \param    void
 * \return   \e void
 */
void Statistics::clean_files( void )
{
  clean_phenotype_mean_file();
  clean_genome_structure_mean_file();
  clean_inherited_proteins_mean_file();
  clean_phenotype_var_file();
  clean_genome_structure_var_file();
  clean_inherited_proteins_var_file();
  clean_phenotype_best_file();
  clean_genome_structure_best_file();
  clean_inherited_proteins_best_file();
  clean_environment_metabolic_amounts_file();
  clean_trophic_network_profile_file();
  clean_level_0_file();
  clean_level_1_file();
  clean_level_2_file();
  clean_no_level_file();
  clean_global_concentrations_file();
  clean_tree_structure_file();
  clean_best_id_file();
  chmod("./statistics/tmp.txt", 0777);
  std::remove("./statistics/tmp.txt");
}

/**
 * \brief    Open files
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::open_files( bool clean_statistic_files )
{
  if (clean_statistic_files)
  {
    if (_population->get_time() == 0)
    {
      /*------------------------------------------------------------------ MEAN statistical variables */
      
      _phenotype_mean_file.open("./statistics/phenotype_mean.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/phenotype_mean.txt", 0777);
      _genome_structure_mean_file.open("./statistics/genome_structure_mean.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/genome_structure_mean.txt", 0777);
      _inherited_proteins_mean_file.open("./statistics/inherited_proteins_mean.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/_inherited_proteins_mean.txt", 0777);
      
      /*------------------------------------------------------------------ VARIANCE statistical variables */
      
      _phenotype_var_file.open("./statistics/phenotype_var.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/phenotype_var.txt", 0777);
      _genome_structure_var_file.open("./statistics/genome_structure_var.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/genome_structure_var.txt", 0777);
      _inherited_proteins_var_file.open("./statistics/inherited_proteins_var.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/_inherited_proteins_var.txt", 0777);
      
      /*------------------------------------------------------------------ BEST statistical variables */
      
      _phenotype_best_file.open("./statistics/phenotype_best.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/phenotype_best.txt", 0777);
      _genome_structure_best_file.open("./statistics/genome_structure_best.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/genome_structure_best.txt", 0777);
      _inherited_proteins_best_file.open("./statistics/inherited_proteins_best.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/_inherited_proteins_best.txt", 0777);
      _best_id_file.open("./statistics/best_id.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/best_id.txt", 0777);
      
      /*------------------------------------------------------------------ ENVIRONMENT statistical variables */
      
      _environment_metabolic_amounts_file.open("./statistics/environment_metabolic_amounts.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/environment_metabolic_amounts.txt", 0777);
      
      /*------------------------------------------------------------------ TROPHIC NETWORK statistical variables */
      
      _trophic_network_profile_file.open("./statistics/trophic_network_profile.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/trophic_network_profile.txt", 0777);
      _level_0_statistics_file.open("./statistics/level_0.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/level_0.txt", 0777);
      _level_1_statistics_file.open("./statistics/level_1.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/level_1.txt", 0777);
      _level_2_statistics_file.open("./statistics/level_2.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/level_2.txt", 0777);
      _no_level_statistics_file.open("./statistics/no_level.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/no_level.txt", 0777);
      
      /*------------------------------------------------------------------ GLOBAL CONCENTRATIONS statistical variables */
      
      _global_concentrations_file.open("./statistics/global_concentrations.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/global_concentrations.txt", 0777);
      
      /*------------------------------------------------------------------ TREE STRUCTURE statistical variables */
      
      _tree_structure_file.open("./statistics/tree_structure.txt", std::ios::out | std::ios::trunc);
      chmod("./statistics/tree_structure.txt", 0777);
      
      write_headers();
    }
    else
    {
      /*------------------------------------------------------------------ MEAN statistical variables */
      
      _phenotype_mean_file.open("./statistics/phenotype_mean.txt", std::ios::out | std::ios::app);
      chmod("./statistics/phenotype_mean.txt", 0777);
      _genome_structure_mean_file.open("./statistics/genome_structure_mean.txt", std::ios::out | std::ios::app);
      chmod("./statistics/genome_structure_mean.txt", 0777);
      _inherited_proteins_mean_file.open("./statistics/inherited_proteins_mean.txt", std::ios::out | std::ios::app);
      chmod("./statistics/_inherited_proteins_mean.txt", 0777);
      
      /*------------------------------------------------------------------ VARIANCE statistical variables */
      
      _phenotype_var_file.open("./statistics/phenotype_var.txt", std::ios::out | std::ios::app);
      chmod("./statistics/phenotype_var.txt", 0777);
      _genome_structure_var_file.open("./statistics/genome_structure_var.txt", std::ios::out | std::ios::app);
      chmod("./statistics/genome_structure_var.txt", 0777);
      _inherited_proteins_var_file.open("./statistics/inherited_proteins_var.txt", std::ios::out | std::ios::app);
      chmod("./statistics/_inherited_proteins_var.txt", 0777);
      
      /*------------------------------------------------------------------ BEST statistical variables */
      
      _phenotype_best_file.open("./statistics/phenotype_best.txt", std::ios::out | std::ios::app);
      chmod("./statistics/phenotype_best.txt", 0777);
      _genome_structure_best_file.open("./statistics/genome_structure_best.txt", std::ios::out | std::ios::app);
      chmod("./statistics/genome_structure_best.txt", 0777);
      _inherited_proteins_best_file.open("./statistics/inherited_proteins_best.txt", std::ios::out | std::ios::app);
      chmod("./statistics/_inherited_proteins_best.txt", 0777);
      _best_id_file.open("./statistics/best_id.txt", std::ios::out | std::ios::app);
      chmod("./statistics/best_id.txt", 0777);
      
      /*------------------------------------------------------------------ ENVIRONMENT statistical variables */
      
      _environment_metabolic_amounts_file.open("./statistics/environment_metabolic_amounts.txt", std::ios::out | std::ios::app);
      chmod("./statistics/environment_metabolic_amounts.txt", 0777);
      
      /*------------------------------------------------------------------ TROPHIC NETWORK statistical variables */
      
      _trophic_network_profile_file.open("./statistics/trophic_network_profile.txt", std::ios::out | std::ios::app);
      chmod("./statistics/trophic_network_profile.txt", 0777);
      _level_0_statistics_file.open("./statistics/level_0.txt", std::ios::out | std::ios::app);
      chmod("./statistics/level_0.txt", 0777);
      _level_1_statistics_file.open("./statistics/level_1.txt", std::ios::out | std::ios::app);
      chmod("./statistics/level_1.txt", 0777);
      _level_2_statistics_file.open("./statistics/level_2.txt", std::ios::out | std::ios::app);
      chmod("./statistics/level_2.txt", 0777);
      _no_level_statistics_file.open("./statistics/no_level.txt", std::ios::out | std::ios::app);
      chmod("./statistics/no_level.txt", 0777);
      
      /*------------------------------------------------------------------ GLOBAL CONCENTRATIONS statistical variables */
      
      _global_concentrations_file.open("./statistics/global_concentrations.txt", std::ios::out | std::ios::app);
      chmod("./statistics/global_concentrations.txt", 0777);
      
      /*------------------------------------------------------------------ TREE STRUCTURE statistical variables */
      
      _tree_structure_file.open("./statistics/tree_structure.txt", std::ios::out | std::ios::app);
      chmod("./statistics/tree_structure.txt", 0777);
    }
  }
  else if (!clean_statistic_files)
  {
    /*------------------------------------------------------------------ MEAN statistical variables */
    
    _phenotype_mean_file.open("./statistics/phenotype_mean.txt", std::ios::out | std::ios::app);
    chmod("./statistics/phenotype_mean.txt", 0777);
    _genome_structure_mean_file.open("./statistics/genome_structure_mean.txt", std::ios::out | std::ios::app);
    chmod("./statistics/genome_structure_mean.txt", 0777);
    _inherited_proteins_mean_file.open("./statistics/inherited_proteins_mean.txt", std::ios::out | std::ios::app);
    chmod("./statistics/_inherited_proteins_mean.txt", 0777);
    
    /*------------------------------------------------------------------ VARIANCE statistical variables */
    
    _phenotype_var_file.open("./statistics/phenotype_var.txt", std::ios::out | std::ios::app);
    chmod("./statistics/phenotype_var.txt", 0777);
    _genome_structure_var_file.open("./statistics/genome_structure_var.txt", std::ios::out | std::ios::app);
    chmod("./statistics/genome_structure_var.txt", 0777);
    _inherited_proteins_var_file.open("./statistics/inherited_proteins_var.txt", std::ios::out | std::ios::app);
    chmod("./statistics/_inherited_proteins_var.txt", 0777);
    
    /*------------------------------------------------------------------ BEST statistical variables */
    
    _phenotype_best_file.open("./statistics/phenotype_best.txt", std::ios::out | std::ios::app);
    chmod("./statistics/phenotype_best.txt", 0777);
    _genome_structure_best_file.open("./statistics/genome_structure_best.txt", std::ios::out | std::ios::app);
    chmod("./statistics/genome_structure_best.txt", 0777);
    _inherited_proteins_best_file.open("./statistics/inherited_proteins_best.txt", std::ios::out | std::ios::app);
    chmod("./statistics/_inherited_proteins_best.txt", 0777);
    _best_id_file.open("./statistics/best_id.txt", std::ios::out | std::ios::app);
    chmod("./statistics/best_id.txt", 0777);
    
    /*------------------------------------------------------------------ ENVIRONMENT statistical variables */
    
    _environment_metabolic_amounts_file.open("./statistics/environment_metabolic_amounts.txt", std::ios::out | std::ios::app);
    chmod("./statistics/environment_metabolic_amounts.txt", 0777);
    
    /*------------------------------------------------------------------ TROPHIC NETWORK statistical variables */
    
    _trophic_network_profile_file.open("./statistics/trophic_network_profile.txt", std::ios::out | std::ios::app);
    chmod("./statistics/trophic_network_profile.txt", 0777);
    _level_0_statistics_file.open("./statistics/level_0.txt", std::ios::out | std::ios::app);
    chmod("./statistics/level_0.txt", 0777);
    _level_1_statistics_file.open("./statistics/level_1.txt", std::ios::out | std::ios::app);
    chmod("./statistics/level_1.txt", 0777);
    _level_2_statistics_file.open("./statistics/level_2.txt", std::ios::out | std::ios::app);
    chmod("./statistics/level_2.txt", 0777);
    _no_level_statistics_file.open("./statistics/no_level.txt", std::ios::out | std::ios::app);
    chmod("./statistics/no_level.txt", 0777);
    
    /*------------------------------------------------------------------ GLOBAL CONCENTRATIONS statistical variables */
    
    _global_concentrations_file.open("./statistics/global_concentrations.txt", std::ios::out | std::ios::app);
    chmod("./statistics/global_concentrations.txt", 0777);
    
    /*------------------------------------------------------------------ TREE STRUCTURE statistical variables */
    
    _tree_structure_file.open("./statistics/tree_structure.txt", std::ios::out | std::ios::app);
    chmod("./statistics/tree_structure.txt", 0777);
  }
}

/**
 * \brief    Write files headers
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_headers( void )
{
  /*------------------------------------------------------------------ MEAN statistical variables */
  
  write_phenotype_mean_file_header();
  write_genome_structure_mean_file_header();
  write_inherited_proteins_mean_file_header();
  
  /*------------------------------------------------------------------ VARIANCE statistical variables */
  
  write_phenotype_var_file_header();
  write_genome_structure_var_file_header();
  write_inherited_proteins_var_file_header();
  
  /*------------------------------------------------------------------ BEST statistical variables */
  
  write_phenotype_best_file_header();
  write_genome_structure_best_file_header();
  write_inherited_proteins_best_file_header();
  write_best_id_file_header();
  
  /*------------------------------------------------------------------ TROPHIC NETWORK statistical variables */
  
  write_trophic_network_profile_file_header();
  write_level_0_file_header();
  write_level_1_file_header();
  write_level_2_file_header();
  write_no_level_file_header();
  
  /*------------------------------------------------------------------ GLOBAL CONCENTRATIONS statistical variables */
  
  write_global_concentrations_file_header();
  
  /*------------------------------------------------------------------ TREE STRUCTURE statistical variables */
  
  write_tree_structure_file_header();
}

/**
 * \brief    Write statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_stats( void )
{
  /*------------------------------------------------------------------ MEAN statistical variables */
  
  write_phenotype_mean_file_stats();
  write_genome_structure_mean_file_stats();
  write_inherited_proteins_mean_file_stats();
  
  /*------------------------------------------------------------------ VARIANCE statistical variables */
  
  write_phenotype_var_file_stats();
  write_genome_structure_var_file_stats();
  write_inherited_proteins_var_file_stats();
  
  /*------------------------------------------------------------------ BEST statistical variables */
  
  write_phenotype_best_file_stats();
  write_genome_structure_best_file_stats();
  write_inherited_proteins_best_file_stats();
  write_best_id_file();
  
  /*------------------------------------------------------------------ ENVIRONMENT statistical variables */
  
  write_environment_metabolic_amounts_stats();
  
  /*------------------------------------------------------------------ TROPHIC NETWORK statistical variables */
  
  write_trophic_network_profile_stats();
  write_level_0_stats();
  write_level_1_stats();
  write_level_2_stats();
  write_no_level_stats();
  
  /*------------------------------------------------------------------ GLOBAL CONCENTRATIONS statistical variables */
  
  write_global_concentrations_stats();
  
  /*------------------------------------------------------------------ TREE STRUCTURE statistical variables */
  
  write_tree_structure_stats();
}

/**
 * \brief    Flush files
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::flush_files( void )
{
  /*------------------------------------------------------------------ MEAN statistical variables */
  
  _phenotype_mean_file.flush();
  _genome_structure_mean_file.flush();
  _inherited_proteins_mean_file.flush();
  
  /*------------------------------------------------------------------ VARIANCE statistical variables */
  
  _phenotype_var_file.flush();
  _genome_structure_var_file.flush();
  _inherited_proteins_var_file.flush();
  
  /*------------------------------------------------------------------ BEST statistical variables */
  
  _phenotype_best_file.flush();
  _genome_structure_best_file.flush();
  _inherited_proteins_best_file.flush();
  _best_id_file.flush();
  
  /*------------------------------------------------------------------ ENVIRONMENT statistical variables */
  
  _environment_metabolic_amounts_file.flush();
  
  /*------------------------------------------------------------------ TROPHIC NETWORK statistical variables */
  
  _trophic_network_profile_file.flush();
  _level_0_statistics_file.flush();
  _level_1_statistics_file.flush();
  _level_2_statistics_file.flush();
  _no_level_statistics_file.flush();
  
  /*------------------------------------------------------------------ GLOBAL CONCENTRATIONS statistical variables */
  
  _global_concentrations_file.flush();
  
  /*------------------------------------------------------------------ TREE STRUCTURE statistical variables */
  
  _tree_structure_file.flush();
}

/**
 *\brief     Close files
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::close_files( void )
{
  /*------------------------------------------------------------------ MEAN statistical variables */
  
  _phenotype_mean_file.close();
  _genome_structure_mean_file.close();
  _inherited_proteins_mean_file.close();
  
  /*------------------------------------------------------------------ VARIANCE statistical variables */
  
  _phenotype_var_file.close();
  _genome_structure_var_file.close();
  _inherited_proteins_var_file.close();
  
  /*------------------------------------------------------------------ BEST statistical variables */
  
  _phenotype_best_file.close();
  _genome_structure_best_file.close();
  _inherited_proteins_best_file.close();
  _best_id_file.close();
  
  /*------------------------------------------------------------------ ENVIRONMENT statistical variables */
  
  _environment_metabolic_amounts_file.close();
  
  /*------------------------------------------------------------------ TROPHIC NETWORK statistical variables */
  
  _trophic_network_profile_file.close();
  _level_0_statistics_file.close();
  _level_1_statistics_file.close();
  _level_2_statistics_file.close();
  _no_level_statistics_file.close();
  
  /*------------------------------------------------------------------ GLOBAL CONCENTRATIONS statistical variables */
  
  _global_concentrations_file.close();
  
  /*------------------------------------------------------------------ TREE STRUCTURE statistical variables */
  
  _tree_structure_file.close();
}

/**
 * \brief    Initialize statistic variables
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::init_variables( void )
{
  /*------------------------------------------------------------------ MEAN statistical variables */
  
  /* PHENOTYPE */
  _mean_generations                = 0.0;
  _mean_inherited_TF_amount        = 0.0;
  _mean_inherited_E_amount         = 0.0;
  _mean_TF_amount                  = 0.0;
  _mean_E_amount                   = 0.0;
  _mean_inherited_metabolic_amount = 0.0;
  _mean_metabolic_amount           = 0.0;
  _mean_energy                     = 0.0;
  _mean_score                      = 0.0;
  _mean_lifespan                   = 0.0;
  _mean_number_of_divisions        = 0.0;
  _mean_toxicity                   = 0.0;
  _mean_metabolic_uptake           = 0.0;
  _mean_metabolic_release          = 0.0;
  _mean_metabolic_growth_rate      = 0.0;
  _mean_Dmetabolic_growth_rate     = 0.0;
  _mean_grn_nb_nodes               = 0.0;
  _mean_grn_nb_edges               = 0.0;
  _mean_metabolic_nb_nodes         = 0.0;
  _mean_metabolic_nb_edges         = 0.0;
  _mean_regulation_redundancy      = 0.0;
  _mean_metabolic_redundancy       = 0.0;
  _mean_cumulated_error            = 0.0;
  
  /* GENOME STRUCTURE */
  _mean_genome_size                   = 0.0;
  _mean_functional_size               = 0.0;
  _mean_genome_nb_NC                  = 0.0;
  _mean_genome_nb_E                   = 0.0;
  _mean_genome_nb_TF                  = 0.0;
  _mean_genome_nb_BS                  = 0.0;
  _mean_genome_nb_P                   = 0.0;
  _mean_genome_nb_inner_enzymes       = 0.0;
  _mean_genome_nb_inflow_pumps        = 0.0;
  _mean_genome_nb_outflow_pumps       = 0.0;
  _mean_genome_nb_functional_regions  = 0.0;
  _mean_genome_nb_enhancers           = 0.0;
  _mean_genome_nb_operators           = 0.0;
  _mean_genome_nb_E_regions           = 0.0;
  _mean_genome_nb_TF_regions          = 0.0;
  _mean_genome_nb_mixed_regions       = 0.0;
  _mean_genome_functional_region_size = 0.0;
  _mean_genome_E_region_size          = 0.0;
  _mean_genome_TF_region_size         = 0.0;
  _mean_genome_mixed_region_size      = 0.0;
  _mean_genome_enhancer_size          = 0.0;
  _mean_genome_operator_size          = 0.0;
  _mean_genome_operon_size            = 0.0;
  _mean_genome_E_operon_size          = 0.0;
  _mean_genome_TF_operon_size         = 0.0;
  _mean_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _mean_inherited_size             = 0.0;
  _mean_inherited_nb_E             = 0.0;
  _mean_inherited_nb_TF            = 0.0;
  _mean_inherited_nb_inner_enzymes = 0.0;
  _mean_inherited_nb_inflow_pumps  = 0.0;
  _mean_inherited_nb_outflow_pumps = 0.0;
  
  /*------------------------------------------------------------------ VARIANCE statistical variables */
  
  /* PHENOTYPE */
  _var_generations                = 0.0;
  _var_inherited_TF_amount        = 0.0;
  _var_inherited_E_amount         = 0.0;
  _var_TF_amount                  = 0.0;
  _var_E_amount                   = 0.0;
  _var_inherited_metabolic_amount = 0.0;
  _var_metabolic_amount           = 0.0;
  _var_energy                     = 0.0;
  _var_score                      = 0.0;
  _var_lifespan                   = 0.0;
  _var_number_of_divisions        = 0.0;
  _var_toxicity                   = 0.0;
  _var_metabolic_uptake           = 0.0;
  _var_metabolic_release          = 0.0;
  _var_metabolic_growth_rate      = 0.0;
  _var_Dmetabolic_growth_rate     = 0.0;
  _var_grn_nb_nodes               = 0.0;
  _var_grn_nb_edges               = 0.0;
  _var_metabolic_nb_nodes         = 0.0;
  _var_metabolic_nb_edges         = 0.0;
  _var_regulation_redundancy      = 0.0;
  _var_metabolic_redundancy       = 0.0;
  _var_cumulated_error            = 0.0;
  
  /* GENOME STRUCTURE */
  _var_genome_size                   = 0.0;
  _var_functional_size               = 0.0;
  _var_genome_nb_NC                  = 0.0;
  _var_genome_nb_E                   = 0.0;
  _var_genome_nb_TF                  = 0.0;
  _var_genome_nb_BS                  = 0.0;
  _var_genome_nb_P                   = 0.0;
  _var_genome_nb_inner_enzymes       = 0.0;
  _var_genome_nb_inflow_pumps        = 0.0;
  _var_genome_nb_outflow_pumps       = 0.0;
  _var_genome_nb_functional_regions  = 0.0;
  _var_genome_nb_enhancers           = 0.0;
  _var_genome_nb_operators           = 0.0;
  _var_genome_nb_E_regions           = 0.0;
  _var_genome_nb_TF_regions          = 0.0;
  _var_genome_nb_mixed_regions       = 0.0;
  _var_genome_functional_region_size = 0.0;
  _var_genome_E_region_size          = 0.0;
  _var_genome_TF_region_size         = 0.0;
  _var_genome_mixed_region_size      = 0.0;
  _var_genome_enhancer_size          = 0.0;
  _var_genome_operator_size          = 0.0;
  _var_genome_operon_size            = 0.0;
  _var_genome_E_operon_size          = 0.0;
  _var_genome_TF_operon_size         = 0.0;
  _var_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _var_inherited_size             = 0.0;
  _var_inherited_nb_E             = 0.0;
  _var_inherited_nb_TF            = 0.0;
  _var_inherited_nb_inner_enzymes = 0.0;
  _var_inherited_nb_inflow_pumps  = 0.0;
  _var_inherited_nb_outflow_pumps = 0.0;
  
  /*------------------------------------------------------------------ BEST statistical variables */
  
  /* PHENOTYPE */
  _best_id                         = 0.0;
  _best_generations                = 0.0;
  _best_inherited_TF_amount        = 0.0;
  _best_inherited_E_amount         = 0.0;
  _best_TF_amount                  = 0.0;
  _best_E_amount                   = 0.0;
  _best_inherited_metabolic_amount = 0.0;
  _best_metabolic_amount           = 0.0;
  _best_energy                     = 0.0;
  _best_score                      = 0.0;
  _best_lifespan                   = 0.0;
  _best_number_of_divisions        = 0.0;
  _best_toxicity                   = 0.0;
  _best_metabolic_uptake           = 0.0;
  _best_metabolic_release          = 0.0;
  _best_metabolic_growth_rate      = 0.0;
  _best_Dmetabolic_growth_rate     = 0.0;
  _best_grn_nb_nodes               = 0.0;
  _best_grn_nb_edges               = 0.0;
  _best_metabolic_nb_nodes         = 0.0;
  _best_metabolic_nb_edges         = 0.0;
  _best_regulation_redundancy      = 0.0;
  _best_metabolic_redundancy       = 0.0;
  _best_trophic_level              = 0.0;
  _best_cumulated_error            = 0.0;
  
  /* GENOME STRUCTURE */
  _best_genome_size                   = 0.0;
  _best_functional_size               = 0.0;
  _best_genome_nb_NC                  = 0.0;
  _best_genome_nb_E                   = 0.0;
  _best_genome_nb_TF                  = 0.0;
  _best_genome_nb_BS                  = 0.0;
  _best_genome_nb_P                   = 0.0;
  _best_genome_nb_inner_enzymes       = 0.0;
  _best_genome_nb_inflow_pumps        = 0.0;
  _best_genome_nb_outflow_pumps       = 0.0;
  _best_genome_nb_functional_regions  = 0.0;
  _best_genome_nb_enhancers           = 0.0;
  _best_genome_nb_operators           = 0.0;
  _best_genome_nb_E_regions           = 0.0;
  _best_genome_nb_TF_regions          = 0.0;
  _best_genome_nb_mixed_regions       = 0.0;
  _best_genome_functional_region_size = 0.0;
  _best_genome_E_region_size          = 0.0;
  _best_genome_TF_region_size         = 0.0;
  _best_genome_mixed_region_size      = 0.0;
  _best_genome_enhancer_size          = 0.0;
  _best_genome_operator_size          = 0.0;
  _best_genome_operon_size            = 0.0;
  _best_genome_E_operon_size          = 0.0;
  _best_genome_TF_operon_size         = 0.0;
  _best_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _best_inherited_size             = 0.0;
  _best_inherited_nb_E             = 0.0;
  _best_inherited_nb_TF            = 0.0;
  _best_inherited_nb_inner_enzymes = 0.0;
  _best_inherited_nb_inflow_pumps  = 0.0;
  _best_inherited_nb_outflow_pumps = 0.0;
}

/**
 *\brief     Add individual to statistics
 * \details  --
 * \param    size_t pos
 * \return   \e void
 */
void Statistics::add_individual( size_t pos )
{
  Cell& cell                            = *_population->get_cell(pos);
  ReplicationReport& replication_report = *(cell.get_replication_report());
  
  /*------------------------------------------------------------------ MEAN statistical variables */
  
  /* PHENOTYPE */
  _mean_generations                += cell.get_generation();
  _mean_inherited_TF_amount        += cell.get_inherited_TF_amount();
  _mean_inherited_E_amount         += cell.get_inherited_E_amount();
  _mean_TF_amount                  += cell.get_TF_amount();
  _mean_E_amount                   += cell.get_E_amount();
  _mean_inherited_metabolic_amount += cell.get_inherited_metabolic_amount();
  _mean_metabolic_amount           += cell.get_species_list()->get_amount();
  _mean_energy                     += cell.get_energy();
  _mean_score                      += cell.get_score();
  _mean_lifespan                   += cell.get_lifespan();
  _mean_number_of_divisions        += cell.get_number_of_divisions();
  _mean_toxicity                   += cell.get_toxicity();
  _mean_metabolic_uptake           += cell.get_metabolic_uptake();
  _mean_metabolic_release          += cell.get_metabolic_release();
  _mean_metabolic_growth_rate      += cell.get_metabolic_growth_rate();
  _mean_Dmetabolic_growth_rate     += cell.get_Dmetabolic_growth_rate();
  _mean_grn_nb_nodes               += cell.get_grn_nb_nodes();
  _mean_grn_nb_edges               += cell.get_grn_nb_edges();
  _mean_metabolic_nb_nodes         += cell.get_metabolic_nb_nodes();
  _mean_metabolic_nb_edges         += cell.get_metabolic_nb_edges();
  _mean_regulation_redundancy      += replication_report.get_mean_regulation_redundancy();
  _mean_metabolic_redundancy       += replication_report.get_mean_metabolic_redundancy();
  _mean_cumulated_error            += cell.get_ode()->get_cumulated_error();
  
  /* GENOME STRUCTURE */
  _mean_genome_size                   += cell.get_genome()->get_size();
  _mean_functional_size               += replication_report.get_genome_functional_size();
  _mean_genome_nb_NC                  += cell.get_genome()->get_nb_NC();
  _mean_genome_nb_E                   += cell.get_genome()->get_nb_E();
  _mean_genome_nb_TF                  += cell.get_genome()->get_nb_TF();
  _mean_genome_nb_BS                  += cell.get_genome()->get_nb_BS();
  _mean_genome_nb_P                   += cell.get_genome()->get_nb_P();
  _mean_genome_nb_inner_enzymes       += cell.get_genome()->get_nb_inner_enzymes();
  _mean_genome_nb_inflow_pumps        += cell.get_genome()->get_nb_inflow_pumps();
  _mean_genome_nb_outflow_pumps       += cell.get_genome()->get_nb_outflow_pumps();
  _mean_genome_nb_functional_regions  += replication_report.get_nb_functional_regions();
  _mean_genome_nb_enhancers           += replication_report.get_nb_enhancers();
  _mean_genome_nb_operators           += replication_report.get_nb_operators();
  _mean_genome_nb_E_regions           += replication_report.get_nb_E_regions();
  _mean_genome_nb_TF_regions          += replication_report.get_nb_TF_regions();
  _mean_genome_nb_mixed_regions       += replication_report.get_nb_mixed_regions();
  _mean_genome_functional_region_size += replication_report.get_mean_functional_region_size();
  _mean_genome_E_region_size          += replication_report.get_mean_E_region_size();
  _mean_genome_TF_region_size         += replication_report.get_mean_TF_region_size();
  _mean_genome_mixed_region_size      += replication_report.get_mean_mixed_region_size();
  _mean_genome_enhancer_size          += replication_report.get_mean_enhancer_size();
  _mean_genome_operator_size          += replication_report.get_mean_operator_size();
  _mean_genome_operon_size            += replication_report.get_mean_operon_size();
  _mean_genome_E_operon_size          += replication_report.get_mean_E_operon_size();
  _mean_genome_TF_operon_size         += replication_report.get_mean_TF_operon_size();
  _mean_genome_mixed_operon_size      += replication_report.get_mean_mixed_operon_size();
  
  /* INHERITED STRUCTURE */
  if (_parameters->get_enzymatic_inheritance())
  {
    _mean_inherited_size             += cell.get_inherited_proteins()->get_size();
    _mean_inherited_nb_E             += cell.get_inherited_proteins()->get_nb_E();
    _mean_inherited_nb_TF            += cell.get_inherited_proteins()->get_nb_TF();
    _mean_inherited_nb_inner_enzymes += cell.get_inherited_proteins()->get_nb_inner_enzymes();
    _mean_inherited_nb_inflow_pumps  += cell.get_inherited_proteins()->get_nb_inflow_pumps();
    _mean_inherited_nb_outflow_pumps += cell.get_inherited_proteins()->get_nb_outflow_pumps();
  }
  
  /*------------------------------------------------------------------ VARIANCE statistical variables */
  
  /* PHENOTYPE */
  _var_generations                += cell.get_generation()*cell.get_generation();
  _var_inherited_TF_amount        += cell.get_inherited_TF_amount()*cell.get_inherited_TF_amount();
  _var_inherited_E_amount         += cell.get_inherited_E_amount()*cell.get_inherited_E_amount();
  _var_TF_amount                  += cell.get_TF_amount()*cell.get_TF_amount();
  _var_E_amount                   += cell.get_E_amount()*cell.get_E_amount();
  _var_inherited_metabolic_amount += cell.get_inherited_metabolic_amount()*cell.get_inherited_metabolic_amount();
  _var_metabolic_amount           += cell.get_species_list()->get_amount()*cell.get_species_list()->get_amount();
  _var_energy                     += cell.get_energy()*cell.get_energy();
  _var_score                      += cell.get_score()*cell.get_score();
  _var_lifespan                   += cell.get_lifespan()*cell.get_lifespan();
  _var_number_of_divisions        += cell.get_number_of_divisions()*cell.get_number_of_divisions();
  _var_toxicity                   += cell.get_toxicity()*cell.get_toxicity();
  _var_metabolic_uptake           += cell.get_metabolic_uptake()*cell.get_metabolic_uptake();
  _var_metabolic_release          += cell.get_metabolic_release()*cell.get_metabolic_release();
  _var_metabolic_growth_rate      += cell.get_metabolic_growth_rate()*cell.get_metabolic_growth_rate();
  _var_Dmetabolic_growth_rate     += cell.get_Dmetabolic_growth_rate()*cell.get_Dmetabolic_growth_rate();
  _var_grn_nb_nodes               += cell.get_grn_nb_nodes()*cell.get_grn_nb_nodes();
  _var_grn_nb_edges               += cell.get_grn_nb_edges()*cell.get_grn_nb_edges();
  _var_metabolic_nb_nodes         += cell.get_metabolic_nb_nodes()*cell.get_metabolic_nb_nodes();
  _var_metabolic_nb_edges         += cell.get_metabolic_nb_edges()*cell.get_metabolic_nb_edges();
  _var_regulation_redundancy      += replication_report.get_mean_regulation_redundancy()*replication_report.get_mean_regulation_redundancy();
  _var_metabolic_redundancy       += replication_report.get_mean_metabolic_redundancy()*replication_report.get_mean_metabolic_redundancy();
  _var_cumulated_error            += cell.get_ode()->get_cumulated_error()*cell.get_ode()->get_cumulated_error();
  
  /* GENOME STRUCTURE */
  _var_genome_size                   += cell.get_genome()->get_size()*cell.get_genome()->get_size();
  _var_functional_size               += replication_report.get_genome_functional_size()*replication_report.get_genome_functional_size();
  _var_genome_nb_NC                  += cell.get_genome()->get_nb_NC()*cell.get_genome()->get_nb_NC();
  _var_genome_nb_E                   += cell.get_genome()->get_nb_E()*cell.get_genome()->get_nb_E();
  _var_genome_nb_TF                  += cell.get_genome()->get_nb_TF()*cell.get_genome()->get_nb_TF();
  _var_genome_nb_BS                  += cell.get_genome()->get_nb_BS()*cell.get_genome()->get_nb_BS();
  _var_genome_nb_P                   += cell.get_genome()->get_nb_P()*cell.get_genome()->get_nb_P();
  _var_genome_nb_inner_enzymes       += cell.get_genome()->get_nb_inner_enzymes()*cell.get_genome()->get_nb_inner_enzymes();
  _var_genome_nb_inflow_pumps        += cell.get_genome()->get_nb_inflow_pumps()*cell.get_genome()->get_nb_inflow_pumps();
  _var_genome_nb_outflow_pumps       += cell.get_genome()->get_nb_outflow_pumps()*cell.get_genome()->get_nb_outflow_pumps();
  _var_genome_nb_functional_regions  += replication_report.get_nb_functional_regions()*replication_report.get_nb_functional_regions();
  _var_genome_nb_enhancers           += replication_report.get_nb_enhancers()*replication_report.get_nb_enhancers();
  _var_genome_nb_operators           += replication_report.get_nb_operators()*replication_report.get_nb_operators();
  _var_genome_nb_E_regions           += replication_report.get_nb_E_regions()*replication_report.get_nb_E_regions();
  _var_genome_nb_TF_regions          += replication_report.get_nb_TF_regions()*replication_report.get_nb_TF_regions();
  _var_genome_nb_mixed_regions       += replication_report.get_nb_mixed_regions()*replication_report.get_nb_mixed_regions();
  _var_genome_functional_region_size += replication_report.get_mean_functional_region_size()*replication_report.get_mean_functional_region_size();
  _var_genome_E_region_size          += replication_report.get_mean_E_region_size()*replication_report.get_mean_E_region_size();
  _var_genome_TF_region_size         += replication_report.get_mean_TF_region_size()*replication_report.get_mean_TF_region_size();
  _var_genome_mixed_region_size      += replication_report.get_mean_mixed_region_size()*replication_report.get_mean_mixed_region_size();
  _var_genome_enhancer_size          += replication_report.get_mean_enhancer_size()*replication_report.get_mean_enhancer_size();
  _var_genome_operator_size          += replication_report.get_mean_operator_size()*replication_report.get_mean_operator_size();
  _var_genome_operon_size            += replication_report.get_mean_operon_size()*replication_report.get_mean_operon_size();
  _var_genome_E_operon_size          += replication_report.get_mean_E_operon_size()*replication_report.get_mean_E_operon_size();
  _var_genome_TF_operon_size         += replication_report.get_mean_TF_operon_size()*replication_report.get_mean_TF_operon_size();
  _var_genome_mixed_operon_size      += replication_report.get_mean_mixed_operon_size()*replication_report.get_mean_mixed_operon_size();
  
  /* INHERITED STRUCTURE */
  if (_parameters->get_enzymatic_inheritance())
  {
    _var_inherited_size             += cell.get_inherited_proteins()->get_size()*cell.get_inherited_proteins()->get_size();
    _var_inherited_nb_E             += cell.get_inherited_proteins()->get_nb_E()*cell.get_inherited_proteins()->get_nb_E();
    _var_inherited_nb_TF            += cell.get_inherited_proteins()->get_nb_TF()*cell.get_inherited_proteins()->get_nb_TF();
    _var_inherited_nb_inner_enzymes += cell.get_inherited_proteins()->get_nb_inner_enzymes()*cell.get_inherited_proteins()->get_nb_inner_enzymes();
    _var_inherited_nb_inflow_pumps  += cell.get_inherited_proteins()->get_nb_inflow_pumps()*cell.get_inherited_proteins()->get_nb_inflow_pumps();
    _var_inherited_nb_outflow_pumps += cell.get_inherited_proteins()->get_nb_outflow_pumps()*cell.get_inherited_proteins()->get_nb_outflow_pumps();
  }
}

/**
 *\brief     Add best individual to statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::add_best_individual( void )
{
  Cell& cell                            = *_population->get_best_cell();
  ReplicationReport& replication_report = *(cell.get_replication_report());
  
  /* PHENOTYPE */
  _best_id                         = cell.get_id();
  _best_generations                = cell.get_generation();
  _best_inherited_TF_amount        = cell.get_inherited_TF_amount();
  _best_inherited_E_amount         = cell.get_inherited_E_amount();
  _best_TF_amount                  = cell.get_TF_amount();
  _best_E_amount                   = cell.get_E_amount();
  _best_inherited_metabolic_amount = cell.get_inherited_metabolic_amount();
  _best_metabolic_amount           = cell.get_species_list()->get_amount();
  _best_energy                     = cell.get_energy();
  _best_score                      = cell.get_score();
  _best_lifespan                   = cell.get_lifespan();
  _best_number_of_divisions        = cell.get_number_of_divisions();
  _best_toxicity                   = cell.get_toxicity();
  _best_metabolic_uptake           = cell.get_metabolic_uptake();
  _best_metabolic_release          = cell.get_metabolic_release();
  _best_metabolic_growth_rate      = cell.get_metabolic_growth_rate();
  _best_Dmetabolic_growth_rate     = cell.get_Dmetabolic_growth_rate();
  _best_grn_nb_nodes               = cell.get_grn_nb_nodes();
  _best_grn_nb_edges               = cell.get_grn_nb_edges();
  _best_metabolic_nb_nodes         = cell.get_metabolic_nb_nodes();
  _best_metabolic_nb_edges         = cell.get_metabolic_nb_edges();
  _best_regulation_redundancy      = replication_report.get_mean_regulation_redundancy();
  _best_metabolic_redundancy       = replication_report.get_mean_metabolic_redundancy();
  _best_trophic_level              = cell.get_trophic_level();
  _best_cumulated_error            = cell.get_ode()->get_cumulated_error();
  
  /* GENOME STRUCTURE */
  _best_genome_size                   = cell.get_genome()->get_size();
  _best_functional_size               = replication_report.get_genome_functional_size();
  _best_genome_nb_NC                  = cell.get_genome()->get_nb_NC();
  _best_genome_nb_E                   = cell.get_genome()->get_nb_E();
  _best_genome_nb_TF                  = cell.get_genome()->get_nb_TF();
  _best_genome_nb_BS                  = cell.get_genome()->get_nb_BS();
  _best_genome_nb_P                   = cell.get_genome()->get_nb_P();
  _best_genome_nb_inner_enzymes       = cell.get_genome()->get_nb_inner_enzymes();
  _best_genome_nb_inflow_pumps        = cell.get_genome()->get_nb_inflow_pumps();
  _best_genome_nb_outflow_pumps       = cell.get_genome()->get_nb_outflow_pumps();
  _best_genome_nb_functional_regions  = replication_report.get_nb_functional_regions();
  _best_genome_nb_enhancers           = replication_report.get_nb_enhancers();
  _best_genome_nb_operators           = replication_report.get_nb_operators();
  _best_genome_nb_E_regions           = replication_report.get_nb_E_regions();
  _best_genome_nb_TF_regions          = replication_report.get_nb_TF_regions();
  _best_genome_nb_mixed_regions       = replication_report.get_nb_mixed_regions();
  _best_genome_functional_region_size = replication_report.get_mean_functional_region_size();
  _best_genome_E_region_size          = replication_report.get_mean_E_region_size();
  _best_genome_TF_region_size         = replication_report.get_mean_TF_region_size();
  _best_genome_mixed_region_size      = replication_report.get_mean_mixed_region_size();
  _best_genome_enhancer_size          = replication_report.get_mean_enhancer_size();
  _best_genome_operator_size          = replication_report.get_mean_operator_size();
  _best_genome_operon_size            = replication_report.get_mean_operon_size();
  _best_genome_E_operon_size          = replication_report.get_mean_E_operon_size();
  _best_genome_TF_operon_size         = replication_report.get_mean_TF_operon_size();
  _best_genome_mixed_operon_size      = replication_report.get_mean_mixed_operon_size();
  
  /* INHERITED STRUCTURE */
  if (_parameters->get_enzymatic_inheritance())
  {
    _best_inherited_size             = cell.get_inherited_proteins()->get_size();
    _best_inherited_nb_E             = cell.get_inherited_proteins()->get_nb_E();
    _best_inherited_nb_TF            = cell.get_inherited_proteins()->get_nb_TF();
    _best_inherited_nb_inner_enzymes = cell.get_inherited_proteins()->get_nb_inner_enzymes();
    _best_inherited_nb_inflow_pumps  = cell.get_inherited_proteins()->get_nb_inflow_pumps();
    _best_inherited_nb_outflow_pumps = cell.get_inherited_proteins()->get_nb_outflow_pumps();
  }
}

/**
 * \brief    Compute means and variances
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::compute_mean_and_var( void )
{
  if (_population->get_population_size() == 0)
  {
    init_variables();
  }
  else
  {
    double pop_size = _population->get_population_size();
    
    /*------------------------------------------------------------------ MEAN statistical variables */
    
    /* PHENOTYPE */
    _mean_generations                /= pop_size;
    _mean_inherited_TF_amount        /= pop_size;
    _mean_inherited_E_amount         /= pop_size;
    _mean_TF_amount                  /= pop_size;
    _mean_E_amount                   /= pop_size;
    _mean_inherited_metabolic_amount /= pop_size;
    _mean_metabolic_amount           /= pop_size;
    _mean_energy                     /= pop_size;
    _mean_score                      /= pop_size;
    _mean_lifespan                   /= pop_size;
    _mean_number_of_divisions        /= pop_size;
    _mean_toxicity                   /= pop_size;
    _mean_metabolic_uptake           /= pop_size;
    _mean_metabolic_release          /= pop_size;
    _mean_metabolic_growth_rate      /= pop_size;
    _mean_Dmetabolic_growth_rate     /= pop_size;
    _mean_grn_nb_nodes               /= pop_size;
    _mean_grn_nb_edges               /= pop_size;
    _mean_metabolic_nb_nodes         /= pop_size;
    _mean_metabolic_nb_edges         /= pop_size;
    _mean_regulation_redundancy      /= pop_size;
    _mean_metabolic_redundancy       /= pop_size;
    _mean_cumulated_error            /= pop_size;
    
    /* GENOME STRUCTURE */
    _mean_genome_size                   /= pop_size;
    _mean_functional_size               /= pop_size;
    _mean_genome_nb_NC                  /= pop_size;
    _mean_genome_nb_E                   /= pop_size;
    _mean_genome_nb_TF                  /= pop_size;
    _mean_genome_nb_BS                  /= pop_size;
    _mean_genome_nb_P                   /= pop_size;
    _mean_genome_nb_inner_enzymes       /= pop_size;
    _mean_genome_nb_inflow_pumps        /= pop_size;
    _mean_genome_nb_outflow_pumps       /= pop_size;
    _mean_genome_nb_functional_regions  /= pop_size;
    _mean_genome_nb_enhancers           /= pop_size;
    _mean_genome_nb_operators           /= pop_size;
    _mean_genome_nb_E_regions           /= pop_size;
    _mean_genome_nb_TF_regions          /= pop_size;
    _mean_genome_nb_mixed_regions       /= pop_size;
    _mean_genome_functional_region_size /= pop_size;
    _mean_genome_E_region_size          /= pop_size;
    _mean_genome_TF_region_size         /= pop_size;
    _mean_genome_mixed_region_size      /= pop_size;
    _mean_genome_enhancer_size          /= pop_size;
    _mean_genome_operator_size          /= pop_size;
    _mean_genome_operon_size            /= pop_size;
    _mean_genome_E_operon_size          /= pop_size;
    _mean_genome_TF_operon_size         /= pop_size;
    _mean_genome_mixed_operon_size      /= pop_size;
    
    /* INHERITED STRUCTURE */
    if (_parameters->get_enzymatic_inheritance())
    {
      _mean_inherited_size             /= pop_size;
      _mean_inherited_nb_E             /= pop_size;
      _mean_inherited_nb_TF            /= pop_size;
      _mean_inherited_nb_inner_enzymes /= pop_size;
      _mean_inherited_nb_inflow_pumps  /= pop_size;
      _mean_inherited_nb_outflow_pumps /= pop_size;
    }
    
    /*------------------------------------------------------------------ VARIANCE statistical variables */
    
    /* PHENOTYPE */
    _var_generations                /= pop_size;
    _var_inherited_TF_amount        /= pop_size;
    _var_inherited_E_amount         /= pop_size;
    _var_TF_amount                  /= pop_size;
    _var_E_amount                   /= pop_size;
    _var_inherited_metabolic_amount /= pop_size;
    _var_metabolic_amount           /= pop_size;
    _var_energy                     /= pop_size;
    _var_score                      /= pop_size;
    _var_lifespan                   /= pop_size;
    _var_number_of_divisions        /= pop_size;
    _var_toxicity                   /= pop_size;
    _var_metabolic_uptake           /= pop_size;
    _var_metabolic_release          /= pop_size;
    _var_metabolic_growth_rate      /= pop_size;
    _var_Dmetabolic_growth_rate     /= pop_size;
    _var_grn_nb_nodes               /= pop_size;
    _var_grn_nb_edges               /= pop_size;
    _var_metabolic_nb_nodes         /= pop_size;
    _var_metabolic_nb_edges         /= pop_size;
    _var_regulation_redundancy      /= pop_size;
    _var_metabolic_redundancy       /= pop_size;
    _var_cumulated_error            /= pop_size;
    
    /* GENOME STRUCTURE */
    _var_genome_size                   /= pop_size;
    _var_functional_size               /= pop_size;
    _var_genome_nb_NC                  /= pop_size;
    _var_genome_nb_E                   /= pop_size;
    _var_genome_nb_TF                  /= pop_size;
    _var_genome_nb_BS                  /= pop_size;
    _var_genome_nb_P                   /= pop_size;
    _var_genome_nb_inner_enzymes       /= pop_size;
    _var_genome_nb_inflow_pumps        /= pop_size;
    _var_genome_nb_outflow_pumps       /= pop_size;
    _var_genome_nb_functional_regions  /= pop_size;
    _var_genome_nb_enhancers           /= pop_size;
    _var_genome_nb_operators           /= pop_size;
    _var_genome_nb_E_regions           /= pop_size;
    _var_genome_nb_TF_regions          /= pop_size;
    _var_genome_nb_mixed_regions       /= pop_size;
    _var_genome_functional_region_size /= pop_size;
    _var_genome_E_region_size          /= pop_size;
    _var_genome_TF_region_size         /= pop_size;
    _var_genome_mixed_region_size      /= pop_size;
    _var_genome_enhancer_size          /= pop_size;
    _var_genome_operator_size          /= pop_size;
    _var_genome_operon_size            /= pop_size;
    _var_genome_E_operon_size          /= pop_size;
    _var_genome_TF_operon_size         /= pop_size;
    _var_genome_mixed_operon_size      /= pop_size;
    
    /* INHERITED STRUCTURE */
    if (_parameters->get_enzymatic_inheritance())
    {
      _var_inherited_size             /= pop_size;
      _var_inherited_nb_E             /= pop_size;
      _var_inherited_nb_TF            /= pop_size;
      _var_inherited_nb_inner_enzymes /= pop_size;
      _var_inherited_nb_inflow_pumps  /= pop_size;
      _var_inherited_nb_outflow_pumps /= pop_size;
    }
    
    /* PHENOTYPE */
    _var_generations                -= _mean_generations*_mean_generations;
    _var_inherited_TF_amount        -= _mean_inherited_TF_amount*_mean_inherited_TF_amount;
    _var_inherited_E_amount         -= _mean_inherited_E_amount*_mean_inherited_E_amount;
    _var_TF_amount                  -= _mean_TF_amount*_mean_TF_amount;
    _var_E_amount                   -= _mean_E_amount*_mean_E_amount;
    _var_inherited_metabolic_amount -= _mean_inherited_metabolic_amount*_mean_inherited_metabolic_amount;
    _var_metabolic_amount           -= _mean_metabolic_amount*_mean_metabolic_amount;
    _var_energy                     -= _mean_energy*_mean_energy;
    _var_score                      -= _mean_score*_mean_score;
    _var_lifespan                   -= _mean_lifespan*_mean_lifespan;
    _var_number_of_divisions        -= _mean_number_of_divisions*_mean_number_of_divisions;
    _var_toxicity                   -= _mean_toxicity*_mean_toxicity;
    _var_metabolic_uptake           -= _mean_metabolic_uptake*_mean_metabolic_uptake;
    _var_metabolic_release          -= _mean_metabolic_release*_mean_metabolic_release;
    _var_metabolic_growth_rate      -= _mean_metabolic_growth_rate*_mean_metabolic_growth_rate;
    _var_Dmetabolic_growth_rate     -= _mean_Dmetabolic_growth_rate*_mean_Dmetabolic_growth_rate;
    _var_grn_nb_nodes               -= _mean_grn_nb_nodes*_mean_grn_nb_nodes;
    _var_grn_nb_edges               -= _mean_grn_nb_edges*_mean_grn_nb_edges;
    _var_metabolic_nb_nodes         -= _mean_metabolic_nb_nodes*_mean_metabolic_nb_nodes;
    _var_metabolic_nb_edges         -= _mean_metabolic_nb_edges*_mean_metabolic_nb_edges;
    _var_regulation_redundancy      -= _mean_regulation_redundancy*_mean_regulation_redundancy;
    _var_metabolic_redundancy       -= _mean_metabolic_redundancy*_mean_metabolic_redundancy;
    _var_cumulated_error            -= _mean_cumulated_error*_mean_cumulated_error;
    
    /* GENOME STRUCTURE */
    _var_genome_size                   -= _mean_genome_size*_mean_genome_size;
    _var_functional_size               -= _mean_functional_size*_mean_functional_size;
    _var_genome_nb_NC                  -= _mean_genome_nb_NC*_mean_genome_nb_NC;
    _var_genome_nb_E                   -= _mean_genome_nb_E*_mean_genome_nb_E;
    _var_genome_nb_TF                  -= _mean_genome_nb_TF*_mean_genome_nb_TF;
    _var_genome_nb_BS                  -= _mean_genome_nb_BS*_mean_genome_nb_BS;
    _var_genome_nb_P                   -= _mean_genome_nb_P*_mean_genome_nb_P;
    _var_genome_nb_inner_enzymes       -= _mean_genome_nb_inner_enzymes*_mean_genome_nb_inner_enzymes;
    _var_genome_nb_inflow_pumps        -= _mean_genome_nb_inflow_pumps*_mean_genome_nb_inflow_pumps;
    _var_genome_nb_outflow_pumps       -= _mean_genome_nb_outflow_pumps*_mean_genome_nb_outflow_pumps;
    _var_genome_nb_functional_regions  -= _mean_genome_nb_functional_regions*_mean_genome_nb_functional_regions;
    _var_genome_nb_enhancers           -= _mean_genome_nb_enhancers*_mean_genome_nb_enhancers;
    _var_genome_nb_operators           -= _mean_genome_nb_operators*_mean_genome_nb_operators;
    _var_genome_nb_E_regions           -= _mean_genome_nb_E_regions*_mean_genome_nb_E_regions;
    _var_genome_nb_TF_regions          -= _mean_genome_nb_TF_regions*_mean_genome_nb_TF_regions;
    _var_genome_nb_mixed_regions       -= _mean_genome_nb_mixed_regions*_mean_genome_nb_mixed_regions;
    _var_genome_functional_region_size -= _mean_genome_functional_region_size*_mean_genome_functional_region_size;
    _var_genome_E_region_size          -= _mean_genome_E_region_size*_mean_genome_E_region_size;
    _var_genome_TF_region_size         -= _mean_genome_TF_region_size*_mean_genome_TF_region_size;
    _var_genome_mixed_region_size      -= _mean_genome_mixed_region_size*_mean_genome_mixed_region_size;
    _var_genome_enhancer_size          -= _mean_genome_enhancer_size*_mean_genome_enhancer_size;
    _var_genome_operator_size          -= _mean_genome_operator_size*_mean_genome_operator_size;
    _var_genome_operon_size            -= _mean_genome_operon_size*_mean_genome_operon_size;
    _var_genome_E_operon_size          -= _mean_genome_E_operon_size*_mean_genome_E_operon_size;
    _var_genome_TF_operon_size         -= _mean_genome_TF_operon_size*_mean_genome_TF_operon_size;
    _var_genome_mixed_operon_size      -= _mean_genome_mixed_operon_size*_mean_genome_mixed_operon_size;
    
    /* INHERITED STRUCTURE */
    _var_inherited_size             -= _mean_inherited_size*_mean_inherited_size;
    _var_inherited_nb_E             -= _mean_inherited_nb_E*_mean_inherited_nb_E;
    _var_inherited_nb_TF            -= _mean_inherited_nb_TF*_mean_inherited_nb_TF;
    _var_inherited_nb_inner_enzymes -= _mean_inherited_nb_inner_enzymes*_mean_inherited_nb_inner_enzymes;
    _var_inherited_nb_inflow_pumps  -= _mean_inherited_nb_inflow_pumps*_mean_inherited_nb_inflow_pumps;
    _var_inherited_nb_outflow_pumps -= _mean_inherited_nb_outflow_pumps*_mean_inherited_nb_outflow_pumps;
  }
}

/**
 * \brief    Write best genome
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_best_genome_file( void )
{
  _best_genome_file.open("./statistics/best_genome.txt", std::ios::out | std::ios::trunc);
  _population->get_best_cell()->write_genome(_best_genome_file);
  _best_genome_file.close();
}

/**
 * \brief    Write best inherited proteins
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_best_inherited_proteins_file( void )
{
  _best_inherited_proteins_file.open("./statistics/best_inherited_proteins.txt", std::ios::out | std::ios::trunc);
  _population->get_best_cell()->write_inherited_proteins(_best_inherited_proteins_file);
  _best_inherited_proteins_file.close();
}

/**
 * \brief    Write best genetic regulation network
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_best_genetic_regulation_network_file( void )
{
  _best_grn_nodes_file.open("./statistics/best_grn_nodes.txt", std::ios::out | std::ios::trunc);
  _best_grn_edges_file.open("./statistics/best_grn_edges.txt", std::ios::out | std::ios::trunc);
  _population->get_best_cell()->write_genetic_regulation_network(_best_grn_nodes_file, _best_grn_edges_file);
  _best_grn_nodes_file.close();
  _best_grn_edges_file.close();
}

/**
 * \brief    Write best metabolic network
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_best_metabolic_network_file( void )
{
  _best_metabolic_nodes_file.open("./statistics/best_metabolic_nodes.txt", std::ios::out | std::ios::trunc);
  _best_metabolic_edges_file.open("./statistics/best_metabolic_edges.txt", std::ios::out | std::ios::trunc);
  _population->get_best_cell()->write_metabolic_network(_best_metabolic_nodes_file, _best_metabolic_edges_file);
  _best_metabolic_nodes_file.close();
  _best_metabolic_edges_file.close();
}

/**
 * \brief    Write best metabolic amounts
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_best_metabolic_amounts_file( void )
{
  _best_metabolic_amounts_file.open("./statistics/best_metabolic_amounts.txt", std::ios::out | std::ios::trunc);
  _population->get_best_cell()->write_metabolic_amounts(_best_metabolic_amounts_file);
  _best_metabolic_amounts_file.close();
}

/**
 * \brief    Write last environment metabolic state vector
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_last_environment_metabolic_amounts_file( void )
{
  _environment->write_metabolic_state_vector("./statistics/last_environment_metabolic_amounts.txt");
  _environment->write_local_metabolic_state_vector("./statistics/last_local_environment_metabolic_amounts.txt", _population->get_best_position());
}

/**
 * \brief    Write last trophic network
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_last_trophic_network_file( void )
{
  _trophic_network->write_trophic_network("./statistics/trophic_network_nodes.txt", "./statistics/trophic_network_edges.txt");
}

/**
 * \brief    Write last lineage tree statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_last_lineage_tree_statistics( void )
{
  _lineage_tree->write_lineage_statistics("lineage_best", _lineage_tree->get_best_alive_node()->get_id());
}

/**
 * \brief    Write last phylogenetic tree statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_last_phylogenetic_tree_statistics( void )
{
  _phylogenetic_tree->write_newick_tree("./statistics/phylogenetic_tree.phb");
  _phylogenetic_tree->write_trophic_data("./statistics/phylogenetic_tree_trophic_data.txt");
}

/**
 * \brief    Write last backup
 * \details  --
 * \param    size_t last_exp_backup
 * \param    size_t last_tree_backup
 * \return   \e void
 */
void Statistics::write_last_backup_file( size_t last_exp_backup, size_t last_tree_backup )
{
  _last_backup_file.open("./statistics/last_backup.txt", std::ios::out | std::ios::trunc);
  _last_backup_file << " " << last_exp_backup << " " << last_tree_backup << "\n";
  _last_backup_file.close();
}

/**
 * \brief    Plot population figures
 * \details  Run Rscripts generating population figures in ./figures/ folder
 * \param    std::string app_path
 * \param    std::string executable_name
 * \return   \e void
 */
void Statistics::plot_population_figures( std::string app_path, std::string executable_name )
{
  /*--------------------------------------------*/
  /* A) write javascript scripts                */
  /*--------------------------------------------*/
  std::string command = "python " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/statistics_to_js.py -p ./ > /dev/null &";
  system(command.c_str());
  command = "python " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/genome_to_js.py -p ./ > /dev/null &";
  system(command.c_str());
  command = "python " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/inherited_proteins_to_js.py -p ./ > /dev/null &";
  system(command.c_str());
  command = "python " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/networks_to_js.py -p ./ > /dev/null &";
  system(command.c_str());
  command = "python " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/environment_to_js.py -p ./ > /dev/null &";
  system(command.c_str());
  
  /*--------------------------------------------*/
  /* B) plot best genome structure              */
  /*--------------------------------------------*/
#ifdef DEBUG
  command = "Rscript " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/plot_best_genome_structure.R ./ > /dev/null &";
  system(command.c_str());
#endif
#ifdef NDEBUG
  command = "Rscript " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/plot_best_genome_structure.R ./ -no-warnings > /dev/null &";
  system(command.c_str());
#endif
  
  /*--------------------------------------------*/
  /* C) plot best metabolic state vector        */
  /*--------------------------------------------*/
#ifdef DEBUG
  command = "Rscript " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/plot_best_phenotype.R ./ > /dev/null &";
  system(command.c_str());
#endif
#ifdef NDEBUG
  command = "Rscript " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/plot_best_phenotype.R ./ -no-warnings > /dev/null &";
  system(command.c_str());
#endif
}

/**
 * \brief    Plot lineage tree figures
 * \details  Run Rscripts generating lineage tree figures in ./figures/ folder
 * \param    std::string app_path
 * \param    std::string executable_name
 * \return   \e void
 */
void Statistics::plot_lineage_tree_figures( std::string app_path, std::string executable_name )
{
  std::string command = "python " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/best_lineage_to_js.py -p ./ > /dev/null &";
  system(command.c_str());
}

/**
 * \brief    Plot phylogenetic tree figures
 * \details  Run Rscripts generating phylogenetic tree figures in ./figures/ folder
 * \param    std::string app_path
 * \param    std::string executable_name
 * \return   \e void
 */
void Statistics::plot_phylogenetic_tree_figures( std::string app_path, std::string executable_name )
{
#ifdef DEBUG
  std::string command = "Rscript " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/plot_phylogenetic_tree.R ./ > /dev/null &";
  system(command.c_str());
  command = "Rscript " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/plot_phylogenetic_tree_with_trophic_level.R ./ > /dev/null &";
  system(command.c_str());
  command = "Rscript " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/plot_phylogenetic_tree_with_trophic_group.R ./ > /dev/null &";
  system(command.c_str());
#endif
#ifdef NDEBUG
  std::string command = "Rscript " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/plot_phylogenetic_tree.R ./ -no-warnings > /dev/null &";
  system(command.c_str());
  command = "Rscript " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/plot_phylogenetic_tree_with_trophic_level.R ./ -no-warnings > /dev/null &";
  system(command.c_str());
  command = "Rscript " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/plot_phylogenetic_tree_with_trophic_group.R ./ -no-warnings > /dev/null &";
  system(command.c_str());
#endif
}

/**
 * \brief    Plot trophic network figures
 * \details  Run Rscripts generating trophic network figures in ./figures/ folder
 * \param    std::string app_path
 * \param    std::string executable_name
 * \return   \e void
 */
void Statistics::plot_trophic_network_figures( std::string app_path, std::string executable_name )
{
  std::string command = "python " + app_path.substr(0, app_path.size()-executable_name.size()) + "src/lib/scripts/trophic_network_to_js.py -p ./ > /dev/null &";
  system(command.c_str());
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Clean phenotype mean file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_phenotype_mean_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/phenotype_mean.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/phenotype_mean.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _phenotype_mean_file.open("./statistics/phenotype_mean.txt", std::ios::out | std::ios::trunc);
    write_phenotype_mean_file_header();
    std::string line;
    while(getline(file, line))
    {
      _phenotype_mean_file << line << "\n";
    }
    _phenotype_mean_file.close();
    file.close();
  }
}

/**
 * \brief    Clean genome structure mean file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_genome_structure_mean_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/genome_structure_mean.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/genome_structure_mean.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _genome_structure_mean_file.open("./statistics/genome_structure_mean.txt", std::ios::out | std::ios::trunc);
    write_genome_structure_mean_file_header();
    std::string line;
    while(getline(file, line))
    {
      _genome_structure_mean_file << line << "\n";
    }
    _genome_structure_mean_file.close();
    file.close();
  }
}

/**
 * \brief    Clean inherited proteins mean file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_inherited_proteins_mean_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/inherited_proteins_mean.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/inherited_proteins_mean.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _inherited_proteins_mean_file.open("./statistics/inherited_proteins_mean.txt", std::ios::out | std::ios::trunc);
    write_inherited_proteins_mean_file_header();
    std::string line;
    while(getline(file, line))
    {
      _inherited_proteins_mean_file << line << "\n";
    }
    _inherited_proteins_mean_file.close();
    file.close();
  }
}

/**
 * \brief    Clean phenotype variance file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_phenotype_var_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/phenotype_var.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/phenotype_var.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _phenotype_var_file.open("./statistics/phenotype_var.txt", std::ios::out | std::ios::trunc);
    write_phenotype_var_file_header();
    std::string line;
    while(getline(file, line))
    {
      _phenotype_var_file << line << "\n";
    }
    _phenotype_var_file.close();
    file.close();
  }
}

/**
 * \brief    Clean genome structure variance file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_genome_structure_var_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/genome_structure_var.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/genome_structure_var.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _genome_structure_var_file.open("./statistics/genome_structure_var.txt", std::ios::out | std::ios::trunc);
    write_genome_structure_var_file_header();
    std::string line;
    while(getline(file, line))
    {
      _genome_structure_var_file << line << "\n";
    }
    _genome_structure_var_file.close();
    file.close();
  }
}

/**
 * \brief    Clean inherited proteins variance file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_inherited_proteins_var_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/inherited_proteins_var.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/inherited_proteins_var.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _inherited_proteins_var_file.open("./statistics/inherited_proteins_var.txt", std::ios::out | std::ios::trunc);
    write_inherited_proteins_var_file_header();
    std::string line;
    while(getline(file, line))
    {
      _inherited_proteins_var_file << line << "\n";
    }
    _inherited_proteins_var_file.close();
    file.close();
  }
}

/**
 * \brief    Clean phenotype best file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_phenotype_best_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/phenotype_best.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/phenotype_best.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _phenotype_best_file.open("./statistics/phenotype_best.txt", std::ios::out | std::ios::trunc);
    write_phenotype_best_file_header();
    std::string line;
    while(getline(file, line))
    {
      _phenotype_best_file << line << "\n";
    }
    _phenotype_best_file.close();
    file.close();
  }
}

/**
 * \brief    Clean genome structure best file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_genome_structure_best_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/genome_structure_best.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/genome_structure_best.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _genome_structure_best_file.open("./statistics/genome_structure_best.txt", std::ios::out | std::ios::trunc);
    write_genome_structure_best_file_header();
    std::string line;
    while(getline(file, line))
    {
      _genome_structure_best_file << line << "\n";
    }
    _genome_structure_best_file.close();
    file.close();
  }
}

/**
 * \brief    Clean inherited proteins best file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_inherited_proteins_best_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/inherited_proteins_best.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/inherited_proteins_best.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _inherited_proteins_best_file.open("./statistics/inherited_proteins_best.txt", std::ios::out | std::ios::trunc);
    write_inherited_proteins_best_file_header();
    std::string line;
    while(getline(file, line))
    {
      _inherited_proteins_best_file << line << "\n";
    }
    _inherited_proteins_best_file.close();
    file.close();
  }
}

/**
 * \brief    Clean environment metabolic amounts file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_environment_metabolic_amounts_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/environment_metabolic_amounts.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/environment_metabolic_amounts.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _environment_metabolic_amounts_file.open("./statistics/environment_metabolic_amounts.txt", std::ios::out | std::ios::trunc);
    std::string line;
    while(getline(file, line))
    {
      _environment_metabolic_amounts_file << line << "\n";
    }
    _environment_metabolic_amounts_file.close();
    file.close();
  }
}

/**
 * \brief    Clean trophic network profile file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_trophic_network_profile_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/trophic_network_profile.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/trophic_network_profile.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _trophic_network_profile_file.open("./statistics/trophic_network_profile.txt", std::ios::out | std::ios::trunc);
    write_trophic_network_profile_file_header();
    std::string line;
    while(getline(file, line))
    {
      _trophic_network_profile_file << line << "\n";
    }
    _trophic_network_profile_file.close();
    file.close();
  }
}

/**
 * \brief    Clean level 0 statistics file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_level_0_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/level_0.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/level_0.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _level_0_statistics_file.open("./statistics/level_0.txt", std::ios::out | std::ios::trunc);
    write_level_0_file_header();
    std::string line;
    while(getline(file, line))
    {
      _level_0_statistics_file << line << "\n";
    }
    _level_0_statistics_file.close();
    file.close();
  }
}

/**
 * \brief    Clean level 1 statistics file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_level_1_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/level_1.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/level_1.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _level_1_statistics_file.open("./statistics/level_1.txt", std::ios::out | std::ios::trunc);
    write_level_1_file_header();
    std::string line;
    while(getline(file, line))
    {
      _level_1_statistics_file << line << "\n";
    }
    _level_1_statistics_file.close();
    file.close();
  }
}

/**
 * \brief    Clean level 2 statistics file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_level_2_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/level_2.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/level_2.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _level_2_statistics_file.open("./statistics/level_2.txt", std::ios::out | std::ios::trunc);
    write_level_2_file_header();
    std::string line;
    while(getline(file, line))
    {
      _level_2_statistics_file << line << "\n";
    }
    _level_2_statistics_file.close();
    file.close();
  }
}

/**
 * \brief    Clean no level statistics file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_no_level_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/no_level.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/no_level.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _no_level_statistics_file.open("./statistics/no_level.txt", std::ios::out | std::ios::trunc);
    write_no_level_file_header();
    std::string line;
    while(getline(file, line))
    {
      _no_level_statistics_file << line << "\n";
    }
    _no_level_statistics_file.close();
    file.close();
  }
}

/**
 * \brief    Clean global concentrations file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_global_concentrations_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/global_concentrations.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/global_concentrations.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _global_concentrations_file.open("./statistics/global_concentrations.txt", std::ios::out | std::ios::trunc);
    write_global_concentrations_file_header();
    std::string line;
    while(getline(file, line))
    {
      _global_concentrations_file << line << "\n";
    }
    _global_concentrations_file.close();
    file.close();
  }
}

/**
 * \brief    Clean tree structure file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_tree_structure_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/tree_structure.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tree_structure.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _tree_structure_file.open("./statistics/tree_structure.txt", std::ios::out | std::ios::trunc);
    write_tree_structure_file_header();
    std::string line;
    while(getline(file, line))
    {
      _tree_structure_file << line << "\n";
    }
    _tree_structure_file.close();
    file.close();
  }
}

/**
 * \brief    Clean best identifier file
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::clean_best_id_file( void )
{
  /* A) copy statistics in temporary file until current time -----------------*/
  
  std::ifstream file("./statistics/best_id.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/best_id.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    std::ofstream tmp("./statistics/tmp.txt", std::ios::out | std::ios::trunc);
    std::string   line;
    size_t        time;
    size_t        counter = 0;
    while(getline(file, line))
    {
      if (counter > 0)
      {
        std::stringstream flux;
        flux.str(line.c_str());
        flux >> time;
        if (time <= _population->get_time())
        {
          tmp << line << "\n";
        }
        else
        {
          break;
        }
      }
      counter++;
    }
    tmp.close();
    file.close();
  }
  
  /* B) recover statistics file ------------------------------------------------*/
  
  file.open("./statistics/tmp.txt", std::ios::in);
  if(!file)
  {
    std::cout << "Error: ./statistics/tmp.txt file not found.\n";
    exit(EXIT_FAILURE);
  }
  else
  {
    _best_id_file.open("./statistics/best_id.txt", std::ios::out | std::ios::trunc);
    write_best_id_file_header();
    std::string line;
    while(getline(file, line))
    {
      _best_id_file << line << "\n";
    }
    _best_id_file.close();
    file.close();
  }
}

/**
 * \brief    Write phenotype mean header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_phenotype_mean_file_header( void )
{
  _phenotype_mean_file << "t" << " " <<
  "generations" << " " <<
  "population_size" << " " <<
  "growth_rate" << " " <<
  "inherited_TF_amount" << " " <<
  "inherited_E_amount" << " " <<
  "TF_amount" << " " <<
  "E_amount" << " " <<
  "inherited_metabolic_amount" << " " <<
  "metabolic_amount" << " " <<
  "energy" << " " <<
  "score" << " " <<
  "lifespan" << " " <<
  "number_of_divisions" << " " <<
  "toxicity" << " " <<
  "metabolic_uptake" << " " <<
  "metabolic_release" << " " <<
  "metabolic_growth_rate" << " " <<
  "Dmetabolic_growth_rate" << " " <<
  "grn_nb_nodes" << " " <<
  "grn_nb_edges" << " " <<
  "metabolic_nb_nodes" << " " <<
  "metabolic_nb_edges" << " " <<
  "regulation_redundancy" << " " <<
  "metabolic_redundancy" << " " <<
  "cumulated_error" << "\n";
}

/**
 * \brief    Write genome structure mean header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_genome_structure_mean_file_header( void )
{
  _genome_structure_mean_file << "t" << " " <<
  "genome_size" << " " <<
  "functional_size" << " " <<
  "nb_NC" << " " <<
  "nb_E" << " " <<
  "nb_TF" << " " <<
  "nb_BS" << " " <<
  "nb_P" << " " <<
  "nb_inner_enzymes" << " " <<
  "nb_inflow_pumps" << " " <<
  "nb_outflow_pumps" << " " <<
  "nb_functional_regions" << " " <<
  "nb_enhancers" << " " <<
  "nb_operators" << " " <<
  "nb_E_regions" << " " <<
  "nb_TF_regions" << " " <<
  "nb_mixed_regions" << " " <<
  "functional_region_size" << " " <<
  "E_region_size" << " " <<
  "TF_region_size" << " " <<
  "mixed_region_size" << " " <<
  "enhancer_size" << " " <<
  "operator_size" << " " <<
  "operon_size" << " " <<
  "E_operon_size" << " " <<
  "TF_operon_size" << " " <<
  "mixed_operon_size" << "\n";
}

/**
 * \brief    Write inherited proteins mean header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_inherited_proteins_mean_file_header( void )
{
  _inherited_proteins_mean_file << "t" << " " <<
  "inherited_size" << " " <<
  "nb_E" << " " <<
  "nb_TF" << " " <<
  "nb_inner_enzymes" << " " <<
  "nb_inflow_pumps" << " " <<
  "nb_outflow_pumps" << "\n";
}

/**
 * \brief    Write phenotype variance header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_phenotype_var_file_header( void )
{
  _phenotype_var_file << "t" << " " <<
  "generations" << " " <<
  "inherited_TF_amount" << " " <<
  "inherited_E_amount" << " " <<
  "TF_amount" << " " <<
  "E_amount" << " " <<
  "inherited_metabolic_amount" << " " <<
  "metabolic_amount" << " " <<
  "energy" << " " <<
  "score" << " " <<
  "lifespan" << " " <<
  "number_of_divisions" << " " <<
  "toxicity" << " " <<
  "metabolic_uptake" << " " <<
  "metabolic_release" << " " <<
  "metabolic_growth_rate" << " " <<
  "Dmetabolic_growth_rate" << " " <<
  "grn_nb_nodes" << " " <<
  "grn_nb_edges" << " " <<
  "metabolic_nb_nodes" << " " <<
  "metabolic_nb_edges" << " " <<
  "regulation_redundancy" << " " <<
  "metabolic_redundancy" << " " <<
  "cumulated_error" << "\n";
}

/**
 * \brief    Write genome structure variance header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_genome_structure_var_file_header( void )
{
  _genome_structure_var_file << "t" << " " <<
  "genome_size" << " " <<
  "functional_size" << " " <<
  "nb_NC" << " " <<
  "nb_E" << " " <<
  "nb_TF" << " " <<
  "nb_BS" << " " <<
  "nb_P" << " " <<
  "nb_inner_enzymes" << " " <<
  "nb_inflow_pumps" << " " <<
  "nb_outflow_pumps" << " " <<
  "nb_functional_regions" << " " <<
  "nb_enhancers" << " " <<
  "nb_operators" << " " <<
  "nb_E_regions" << " " <<
  "nb_TF_regions" << " " <<
  "nb_mixed_regions" << " " <<
  "functional_region_size" << " " <<
  "E_region_size" << " " <<
  "TF_region_size" << " " <<
  "mixed_region_size" << " " <<
  "enhancer_size" << " " <<
  "operator_size" << " " <<
  "operon_size" << " " <<
  "E_operon_size" << " " <<
  "TF_operon_size" << " " <<
  "mixed_operon_size" << "\n";
}

/**
 * \brief    Write inherited proteins variance header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_inherited_proteins_var_file_header( void )
{
  _inherited_proteins_var_file << "t" << " " <<
  "inherited_size" << " " <<
  "nb_E" << " " <<
  "nb_TF" << " " <<
  "nb_inner_enzymes" << " " <<
  "nb_inflow_pumps" << " " <<
  "nb_outflow_pumps" << "\n";
}

/**
 * \brief    Write phenotype best header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_phenotype_best_file_header( void )
{
  _phenotype_best_file << "t" << " " <<
  "generations" << " " <<
  "inherited_TF_amount" << " " <<
  "inherited_E_amount" << " " <<
  "TF_amount" << " " <<
  "E_amount" << " " <<
  "inherited_metabolic_amount" << " " <<
  "metabolic_amount" << " " <<
  "energy" << " " <<
  "score" << " " <<
  "lifespan" << " " <<
  "number_of_divisions" << " " <<
  "toxicity" << " " <<
  "metabolic_uptake" << " " <<
  "metabolic_release" << " " <<
  "metabolic_growth_rate" << " " <<
  "Dmetabolic_growth_rate" << " " <<
  "grn_nb_nodes" << " " <<
  "grn_nb_edges" << " " <<
  "metabolic_nb_nodes" << " " <<
  "metabolic_nb_edges" << " " <<
  "regulation_redundancy" << " " <<
  "metabolic_redundancy" << " " <<
  "trophic_level" << " " <<
  "cumulated_error" << "\n";
}

/**
 * \brief    Write genome structure best header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_genome_structure_best_file_header( void )
{
  _genome_structure_best_file << "t" << " " <<
  "genome_size" << " " <<
  "functional_size" << " " <<
  "nb_NC" << " " <<
  "nb_E" << " " <<
  "nb_TF" << " " <<
  "nb_BS" << " " <<
  "nb_P" << " " <<
  "nb_inner_enzymes" << " " <<
  "nb_inflow_pumps" << " " <<
  "nb_outflow_pumps" << " " <<
  "nb_functional_regions" << " " <<
  "nb_enhancers" << " " <<
  "nb_operators" << " " <<
  "nb_E_regions" << " " <<
  "nb_TF_regions" << " " <<
  "nb_mixed_regions" << " " <<
  "functional_region_size" << " " <<
  "E_region_size" << " " <<
  "TF_region_size" << " " <<
  "mixed_region_size" << " " <<
  "enhancer_size" << " " <<
  "operator_size" << " " <<
  "operon_size" << " " <<
  "E_operon_size" << " " <<
  "TF_operon_size" << " " <<
  "mixed_operon_size" << "\n";
}

/**
 * \brief    Write inherited proteins best header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_inherited_proteins_best_file_header( void )
{
  _inherited_proteins_best_file << "t" << " " <<
  "inherited_size" << " " <<
  "nb_E" << " " <<
  "nb_TF" << " " <<
  "nb_inner_enzymes" << " " <<
  "nb_inflow_pumps" << " " <<
  "nb_outflow_pumps" << "\n";
}

/**
 * \brief    Write trophic network profile header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_trophic_network_profile_file_header( void )
{
  _trophic_network_profile_file << "t" << " " <<
  "level0_nb_cells" << " " <<
  "level1_nb_cells" << " " <<
  "level2_nb_cells" << " " <<
  "nolevel_nb_cells" << " " <<
  "level0_nb_groups" << " " <<
  "level1_nb_groups" << " " <<
  "level2_nb_groups" << " " <<
  "nolevel_nb_groups" << " " <<
  "nb_group_appearances" << " " <<
  "nb_group_extinctions" << " " <<
  "mean_group_lifespan" << " " <<
  "min_true_diversity" << " " <<
  "min_time_distance" << " " <<
  "min_generation_distance" << " " <<
  "min_point_mutation_distance" << " " <<
  "min_hgt_distance" << " " <<
  "min_duplication_distance" << " " <<
  "min_deletion_distance" << " " <<
  "min_inversion_distance" << " " <<
  "min_translocation_distance" << " " <<
  "mean_true_diversity" << " " <<
  "mean_time_distance" << " " <<
  "mean_generation_distance" << " " <<
  "mean_point_mutation_distance" << " " <<
  "mean_hgt_distance" << " " <<
  "mean_duplication_distance" << " " <<
  "mean_deletion_distance" << " " <<
  "mean_inversion_distance" << " " <<
  "mean_translocation_distance" << " " <<
  "max_true_diversity" << " " <<
  "max_time_distance" << " " <<
  "max_generation_distance" << " " <<
  "max_point_mutation_distance" << " " <<
  "max_hgt_distance" << " " <<
  "max_duplication_distance" << " " <<
  "max_deletion_distance" << " " <<
  "max_inversion_distance" << " " <<
  "max_translocation_distance" << "\n";
}

/**
 * \brief    Write level 0 header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_level_0_file_header( void )
{
  _level_0_statistics_file << "t" << " " <<
  
  /* GENETIC DIVERSITY */
  "true_diversity" << " " <<
  "max_time_distance" << " " <<
  "max_generation_distance" << " " <<
  "max_point_mutation_distance" << " " <<
  "max_hgt_distance" << " " <<
  "max_duplication_distance" << " " <<
  "max_deletion_distance" << " " <<
  "max_inversion_distance" << " " <<
  "max_translocation_distance" << " " <<
  
  /* PHENOTYPE */
  "generations" << " " <<
  "inherited_TF_amount" << " " <<
  "inherited_E_amount" << " " <<
  "TF_amount" << " " <<
  "amount" << " " <<
  "inherited_metabolic_amount" << " " <<
  "metabolic_amount" << " " <<
  "energy" << " " <<
  "score" << " " <<
  "lifespan" << " " <<
  "number_of_divisions" << " " <<
  "toxicity" << " " <<
  "metabolic_uptake" << " " <<
  "metabolic_release" << " " <<
  "metabolic_growth_rate" << " " <<
  "Dmetabolic_growth_rate" << " " <<
  "grn_nb_nodes" << " " <<
  "grn_nb_edges" << " " <<
  "metabolic_nb_nodes" << " " <<
  "metabolic_nb_edges" << " " <<
  "regulation_redundancy" << " " <<
  "metabolic_redundancy" << " " <<
  
  /* GENOME STRUCTURE */
  "genome_size" << " " <<
  "functional_size" << " " <<
  "genome_nb_NC" << " " <<
  "genome_nb_E" << " " <<
  "genome_nb_TF" << " " <<
  "genome_nb_BS" << " " <<
  "genome_nb_P" << " " <<
  "genome_nb_inner_enzymes" << " " <<
  "genome_nb_inflow_pumps" << " " <<
  "genome_nb_outflow_pumps" << " " <<
  "genome_nb_functional_regions" << " " <<
  "genome_nb_enhancers" << " " <<
  "genome_nb_operators" << " " <<
  "genome_nb_E_regions" << " " <<
  "genome_nb_TF_regions" << " " <<
  "genome_nb_mixed_regions" << " " <<
  "genome_functional_region_size" << " " <<
  "genome_E_region_size" << " " <<
  "genome_TF_region_size" << " " <<
  "genome_mixed_region_size" << " " <<
  "genome_enhancer_size" << " " <<
  "genome_operator_size" << " " <<
  "genome_operon_size" << " " <<
  "genome_E_operon_size" << " " <<
  "genome_TF_operon_size" << " " <<
  "genome_mixed_operon_size" << " " <<
  
  /* INHERITED STRUCTURE */
  "inherited_size" << " " <<
  "inherited_nb_E" << " " <<
  "inherited_nb_TF" << " " <<
  "inherited_nb_inner_enzymes" << " " <<
  "inherited_nb_inflow_pumps"  << " " <<
  "inherited_nb_outflow_pumps" << "\n";
}

/**
 * \brief    Write level 1 header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_level_1_file_header( void )
{
  _level_1_statistics_file << "t" << " " <<
  
  /* GENETIC DIVERSITY */
  "true_diversity" << " " <<
  "max_time_distance" << " " <<
  "max_generation_distance" << " " <<
  "max_point_mutation_distance" << " " <<
  "max_hgt_distance" << " " <<
  "max_duplication_distance" << " " <<
  "max_deletion_distance" << " " <<
  "max_inversion_distance" << " " <<
  "max_translocation_distance" << " " <<
  
  /* PHENOTYPE */
  "generations" << " " <<
  "inherited_TF_amount" << " " <<
  "inherited_E_amount" << " " <<
  "TF_amount" << " " <<
  "amount" << " " <<
  "inherited_metabolic_amount" << " " <<
  "metabolic_amount" << " " <<
  "energy" << " " <<
  "score" << " " <<
  "lifespan" << " " <<
  "number_of_divisions" << " " <<
  "toxicity" << " " <<
  "metabolic_uptake" << " " <<
  "metabolic_release" << " " <<
  "metabolic_growth_rate" << " " <<
  "Dmetabolic_growth_rate" << " " <<
  "grn_nb_nodes" << " " <<
  "grn_nb_edges" << " " <<
  "metabolic_nb_nodes" << " " <<
  "metabolic_nb_edges" << " " <<
  "regulation_redundancy" << " " <<
  "metabolic_redundancy" << " " <<
  
  /* GENOME STRUCTURE */
  "genome_size" << " " <<
  "functional_size" << " " <<
  "genome_nb_NC" << " " <<
  "genome_nb_E" << " " <<
  "genome_nb_TF" << " " <<
  "genome_nb_BS" << " " <<
  "genome_nb_P" << " " <<
  "genome_nb_inner_enzymes" << " " <<
  "genome_nb_inflow_pumps" << " " <<
  "genome_nb_outflow_pumps" << " " <<
  "genome_nb_functional_regions" << " " <<
  "genome_nb_enhancers" << " " <<
  "genome_nb_operators" << " " <<
  "genome_nb_E_regions" << " " <<
  "genome_nb_TF_regions" << " " <<
  "genome_nb_mixed_regions" << " " <<
  "genome_functional_region_size" << " " <<
  "genome_E_region_size" << " " <<
  "genome_TF_region_size" << " " <<
  "genome_mixed_region_size" << " " <<
  "genome_enhancer_size" << " " <<
  "genome_operator_size" << " " <<
  "genome_operon_size" << " " <<
  "genome_E_operon_size" << " " <<
  "genome_TF_operon_size" << " " <<
  "genome_mixed_operon_size" << " " <<
  
  /* INHERITED STRUCTURE */
  "inherited_size" << " " <<
  "inherited_nb_E" << " " <<
  "inherited_nb_TF" << " " <<
  "inherited_nb_inner_enzymes" << " " <<
  "inherited_nb_inflow_pumps"  << " " <<
  "inherited_nb_outflow_pumps" << "\n";
}

/**
 * \brief    Write level 2 header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_level_2_file_header( void )
{
  _level_2_statistics_file << "t" << " " <<
  
  /* GENETIC DIVERSITY */
  "true_diversity" << " " <<
  "max_time_distance" << " " <<
  "max_generation_distance" << " " <<
  "max_point_mutation_distance" << " " <<
  "max_hgt_distance" << " " <<
  "max_duplication_distance" << " " <<
  "max_deletion_distance" << " " <<
  "max_inversion_distance" << " " <<
  "max_translocation_distance" << " " <<
  
  /* PHENOTYPE */
  "generations" << " " <<
  "inherited_TF_amount" << " " <<
  "inherited_E_amount" << " " <<
  "TF_amount" << " " <<
  "amount" << " " <<
  "inherited_metabolic_amount" << " " <<
  "metabolic_amount" << " " <<
  "energy" << " " <<
  "score" << " " <<
  "lifespan" << " " <<
  "number_of_divisions" << " " <<
  "toxicity" << " " <<
  "metabolic_uptake" << " " <<
  "metabolic_release" << " " <<
  "metabolic_growth_rate" << " " <<
  "Dmetabolic_growth_rate" << " " <<
  "grn_nb_nodes" << " " <<
  "grn_nb_edges" << " " <<
  "metabolic_nb_nodes" << " " <<
  "metabolic_nb_edges" << " " <<
  "regulation_redundancy" << " " <<
  "metabolic_redundancy" << " " <<
  
  /* GENOME STRUCTURE */
  "genome_size" << " " <<
  "functional_size" << " " <<
  "genome_nb_NC" << " " <<
  "genome_nb_E" << " " <<
  "genome_nb_TF" << " " <<
  "genome_nb_BS" << " " <<
  "genome_nb_P" << " " <<
  "genome_nb_inner_enzymes" << " " <<
  "genome_nb_inflow_pumps" << " " <<
  "genome_nb_outflow_pumps" << " " <<
  "genome_nb_functional_regions" << " " <<
  "genome_nb_enhancers" << " " <<
  "genome_nb_operators" << " " <<
  "genome_nb_E_regions" << " " <<
  "genome_nb_TF_regions" << " " <<
  "genome_nb_mixed_regions" << " " <<
  "genome_functional_region_size" << " " <<
  "genome_E_region_size" << " " <<
  "genome_TF_region_size" << " " <<
  "genome_mixed_region_size" << " " <<
  "genome_enhancer_size" << " " <<
  "genome_operator_size" << " " <<
  "genome_operon_size" << " " <<
  "genome_E_operon_size" << " " <<
  "genome_TF_operon_size" << " " <<
  "genome_mixed_operon_size" << " " <<
  
  /* INHERITED STRUCTURE */
  "inherited_size" << " " <<
  "inherited_nb_E" << " " <<
  "inherited_nb_TF" << " " <<
  "inherited_nb_inner_enzymes" << " " <<
  "inherited_nb_inflow_pumps"  << " " <<
  "inherited_nb_outflow_pumps" << "\n";
}

/**
 * \brief    Write no level header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_no_level_file_header( void )
{
  _no_level_statistics_file << "t" << " " <<
  
  /* GENETIC DIVERSITY */
  "true_diversity" << " " <<
  "max_time_distance" << " " <<
  "max_generation_distance" << " " <<
  "max_point_mutation_distance" << " " <<
  "max_hgt_distance" << " " <<
  "max_duplication_distance" << " " <<
  "max_deletion_distance" << " " <<
  "max_inversion_distance" << " " <<
  "max_translocation_distance" << " " <<
  
  /* PHENOTYPE */
  "generations" << " " <<
  "inherited_TF_amount" << " " <<
  "inherited_E_amount" << " " <<
  "TF_amount" << " " <<
  "amount" << " " <<
  "inherited_metabolic_amount" << " " <<
  "metabolic_amount" << " " <<
  "energy" << " " <<
  "score" << " " <<
  "lifespan" << " " <<
  "number_of_divisions" << " " <<
  "toxicity" << " " <<
  "metabolic_uptake" << " " <<
  "metabolic_release" << " " <<
  "metabolic_growth_rate" << " " <<
  "Dmetabolic_growth_rate" << " " <<
  "grn_nb_nodes" << " " <<
  "grn_nb_edges" << " " <<
  "metabolic_nb_nodes" << " " <<
  "metabolic_nb_edges" << " " <<
  "regulation_redundancy" << " " <<
  "metabolic_redundancy" << " " <<
  
  /* GENOME STRUCTURE */
  "genome_size" << " " <<
  "functional_size" << " " <<
  "genome_nb_NC" << " " <<
  "genome_nb_E" << " " <<
  "genome_nb_TF" << " " <<
  "genome_nb_BS" << " " <<
  "genome_nb_P" << " " <<
  "genome_nb_inner_enzymes" << " " <<
  "genome_nb_inflow_pumps" << " " <<
  "genome_nb_outflow_pumps" << " " <<
  "genome_nb_functional_regions" << " " <<
  "genome_nb_enhancers" << " " <<
  "genome_nb_operators" << " " <<
  "genome_nb_E_regions" << " " <<
  "genome_nb_TF_regions" << " " <<
  "genome_nb_mixed_regions" << " " <<
  "genome_functional_region_size" << " " <<
  "genome_E_region_size" << " " <<
  "genome_TF_region_size" << " " <<
  "genome_mixed_region_size" << " " <<
  "genome_enhancer_size" << " " <<
  "genome_operator_size" << " " <<
  "genome_operon_size" << " " <<
  "genome_E_operon_size" << " " <<
  "genome_TF_operon_size" << " " <<
  "genome_mixed_operon_size" << " " <<
  
  /* INHERITED STRUCTURE */
  "inherited_size" << " " <<
  "inherited_nb_E" << " " <<
  "inherited_nb_TF" << " " <<
  "inherited_nb_inner_enzymes" << " " <<
  "inherited_nb_inflow_pumps"  << " " <<
  "inherited_nb_outflow_pumps" << "\n";
}

/**
 * \brief    Write global concentrations header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_global_concentrations_file_header( void )
{
  _global_concentrations_file << "t" << " " <<
  "matter_bilan" << " " <<
  "matter_inflow" << " " <<
  "matter_outflow" << " " <<
  "matter_pop" << " " <<
  "matter_env" << "\n";
}

/**
 * \brief    Write tree structure header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_tree_structure_file_header( void )
{
  _tree_structure_file << "t" << " " <<
  "lineage_nb_nodes" << " " <<
  "lineage_nb_roots" << " " <<
  "lineage_nb_normal" << " " <<
  "lineage_nb_dead" << " " <<
  "lineage_nb_alive" << " " <<
  "phylogeny_nb_nodes" << " " <<
  "phylogeny_nb_roots" << " " <<
  "phylogeny_nb_normal" << " " <<
  "phylogeny_nb_dead" << " " <<
  "phylogeny_nb_alive" << " " <<
  "phylogeny_ca_age" << "\n";
}

/**
 * \brief    Write best identifier header
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_best_id_file_header( void )
{
  _best_id_file << "t" << " " << "generation" << " " << "best_id" << "\n";
}

/**
 * \brief    Write phenotype mean statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_phenotype_mean_file_stats( void )
{
  _phenotype_mean_file << _population->get_time() << " " <<
  _mean_generations << " " <<
  _population->get_population_size() << " " <<
  _population->get_growth_rate() << " " <<
  _mean_inherited_TF_amount << " " <<
  _mean_inherited_E_amount << " " <<
  _mean_TF_amount << " " <<
  _mean_E_amount << " " <<
  _mean_inherited_metabolic_amount << " " <<
  _mean_metabolic_amount << " " <<
  _mean_energy << " " <<
  _mean_score << " " <<
  _mean_lifespan << " " <<
  _mean_number_of_divisions << " " <<
  _mean_toxicity << " " <<
  _mean_metabolic_uptake << " " <<
  _mean_metabolic_release << " " <<
  _mean_metabolic_growth_rate << " " <<
  _mean_Dmetabolic_growth_rate << " " <<
  _mean_grn_nb_nodes << " " <<
  _mean_grn_nb_edges << " " <<
  _mean_metabolic_nb_nodes << " " <<
  _mean_metabolic_nb_edges << " " <<
  _mean_regulation_redundancy << " " <<
  _mean_metabolic_redundancy << " " <<
  _mean_cumulated_error << "\n";
}

/**
 * \brief    Write genome structure mean statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_genome_structure_mean_file_stats( void )
{
  _genome_structure_mean_file << _population->get_time() << " " <<
  _mean_genome_size << " " <<
  _mean_functional_size << " " <<
  _mean_genome_nb_NC << " " <<
  _mean_genome_nb_E << " " <<
  _mean_genome_nb_TF << " " <<
  _mean_genome_nb_BS << " " <<
  _mean_genome_nb_P << " " <<
  _mean_genome_nb_inner_enzymes << " " <<
  _mean_genome_nb_inflow_pumps << " " <<
  _mean_genome_nb_outflow_pumps << " " <<
  _mean_genome_nb_functional_regions << " " <<
  _mean_genome_nb_enhancers << " " <<
  _mean_genome_nb_operators << " " <<
  _mean_genome_nb_E_regions << " " <<
  _mean_genome_nb_TF_regions << " " <<
  _mean_genome_nb_mixed_regions << " " <<
  _mean_genome_functional_region_size << " " <<
  _mean_genome_E_region_size << " " <<
  _mean_genome_TF_region_size << " " <<
  _mean_genome_mixed_region_size << " " <<
  _mean_genome_enhancer_size << " " <<
  _mean_genome_operator_size << " " <<
  _mean_genome_operon_size << " " <<
  _mean_genome_E_operon_size << " " <<
  _mean_genome_TF_operon_size << " " <<
  _mean_genome_mixed_operon_size << "\n";
}

/**
 * \brief    Write inherited proteins mean statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_inherited_proteins_mean_file_stats( void )
{
  _inherited_proteins_mean_file << _population->get_time() << " " <<
  _mean_inherited_size << " " <<
  _mean_inherited_nb_E << " " <<
  _mean_inherited_nb_TF << " " <<
  _mean_inherited_nb_inner_enzymes << " " <<
  _mean_inherited_nb_inflow_pumps << " " <<
  _mean_inherited_nb_outflow_pumps << "\n";
}

/**
 * \brief    Write phenotype variance statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_phenotype_var_file_stats( void )
{
  _phenotype_var_file << _population->get_time() << " " <<
  _var_generations << " " <<
  _var_inherited_TF_amount << " " <<
  _var_inherited_E_amount << " " <<
  _var_TF_amount << " " <<
  _var_E_amount << " " <<
  _var_inherited_metabolic_amount << " " <<
  _var_metabolic_amount << " " <<
  _var_energy << " " <<
  _var_score << " " <<
  _var_lifespan << " " <<
  _var_number_of_divisions << " " <<
  _var_toxicity << " " <<
  _var_metabolic_uptake << " " <<
  _var_metabolic_release << " " <<
  _var_metabolic_growth_rate << " " <<
  _var_Dmetabolic_growth_rate << " " <<
  _var_grn_nb_nodes << " " <<
  _var_grn_nb_edges << " " <<
  _var_metabolic_nb_nodes << " " <<
  _var_metabolic_nb_edges << " " <<
  _var_regulation_redundancy << " " <<
  _var_metabolic_redundancy << " " <<
  _var_cumulated_error << "\n";
}

/**
 * \brief    Write genome structure variance statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_genome_structure_var_file_stats( void )
{
  _genome_structure_var_file << _population->get_time() << " " <<
  _var_genome_size << " " <<
  _var_functional_size << " " <<
  _var_genome_nb_NC << " " <<
  _var_genome_nb_E << " " <<
  _var_genome_nb_TF << " " <<
  _var_genome_nb_BS << " " <<
  _var_genome_nb_P << " " <<
  _var_genome_nb_inner_enzymes << " " <<
  _var_genome_nb_inflow_pumps << " " <<
  _var_genome_nb_outflow_pumps << " " <<
  _var_genome_nb_functional_regions << " " <<
  _var_genome_nb_enhancers << " " <<
  _var_genome_nb_operators << " " <<
  _var_genome_nb_E_regions << " " <<
  _var_genome_nb_TF_regions << " " <<
  _var_genome_nb_mixed_regions << " " <<
  _var_genome_functional_region_size << " " <<
  _var_genome_E_region_size << " " <<
  _var_genome_TF_region_size << " " <<
  _var_genome_mixed_region_size << " " <<
  _var_genome_enhancer_size << " " <<
  _var_genome_operator_size << " " <<
  _var_genome_operon_size << " " <<
  _var_genome_E_operon_size << " " <<
  _var_genome_TF_operon_size << " " <<
  _var_genome_mixed_operon_size << "\n";
}

/**
 * \brief    Write inherited proteins variance statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_inherited_proteins_var_file_stats( void )
{
  _inherited_proteins_var_file << _population->get_time() << " " <<
  _var_inherited_size << " " <<
  _var_inherited_nb_E << " " <<
  _var_inherited_nb_TF << " " <<
  _var_inherited_nb_inner_enzymes << " " <<
  _var_inherited_nb_inflow_pumps << " " <<
  _var_inherited_nb_outflow_pumps << "\n";
}

/**
 * \brief    Write best phenotype statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_phenotype_best_file_stats( void )
{
  _phenotype_best_file << _population->get_time() << " " <<
  _best_generations << " " <<
  _best_inherited_TF_amount << " " <<
  _best_inherited_E_amount << " " <<
  _best_TF_amount << " " <<
  _best_E_amount << " " <<
  _best_inherited_metabolic_amount << " " <<
  _best_metabolic_amount << " " <<
  _best_energy << " " <<
  _best_score << " " <<
  _best_lifespan << " " <<
  _best_number_of_divisions << " " <<
  _best_toxicity << " " <<
  _best_metabolic_uptake << " " <<
  _best_metabolic_release << " " <<
  _best_metabolic_growth_rate << " " <<
  _best_Dmetabolic_growth_rate << " " <<
  _best_grn_nb_nodes << " " <<
  _best_grn_nb_edges << " " <<
  _best_metabolic_nb_nodes << " " <<
  _best_metabolic_nb_edges << " " <<
  _best_regulation_redundancy << " " <<
  _best_metabolic_redundancy << " " <<
  _best_trophic_level << " " <<
  _best_cumulated_error << "\n";
}

/**
 * \brief    Write best genome structure statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_genome_structure_best_file_stats( void )
{
  _genome_structure_best_file << _population->get_time() << " " <<
  _best_genome_size << " " <<
  _best_functional_size << " " <<
  _best_genome_nb_NC << " " <<
  _best_genome_nb_E << " " <<
  _best_genome_nb_TF << " " <<
  _best_genome_nb_BS << " " <<
  _best_genome_nb_P << " " <<
  _best_genome_nb_inner_enzymes << " " <<
  _best_genome_nb_inflow_pumps << " " <<
  _best_genome_nb_outflow_pumps << " " <<
  _best_genome_nb_functional_regions << " " <<
  _best_genome_nb_enhancers << " " <<
  _best_genome_nb_operators << " " <<
  _best_genome_nb_E_regions << " " <<
  _best_genome_nb_TF_regions << " " <<
  _best_genome_nb_mixed_regions << " " <<
  _best_genome_functional_region_size << " " <<
  _best_genome_E_region_size << " " <<
  _best_genome_TF_region_size << " " <<
  _best_genome_mixed_region_size << " " <<
  _best_genome_enhancer_size << " " <<
  _best_genome_operator_size << " " <<
  _best_genome_operon_size << " " <<
  _best_genome_E_operon_size << " " <<
  _best_genome_TF_operon_size << " " <<
  _best_genome_mixed_operon_size << "\n";
}

/**
 * \brief    Write best inherited proteins statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_inherited_proteins_best_file_stats( void )
{
  _inherited_proteins_best_file << _population->get_time() << " " <<
  _best_inherited_size << " " <<
  _best_inherited_nb_E << " " <<
  _best_inherited_nb_TF << " " <<
  _best_inherited_nb_inner_enzymes << " " <<
  _best_inherited_nb_inflow_pumps << " " <<
  _best_inherited_nb_outflow_pumps << "\n";
}

/**
 * \brief    Write environment metabolic amounts statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_environment_metabolic_amounts_stats( void )
{
  _environment_metabolic_amounts_file << _population->get_time();
  for (int i = 0; i < (int)_environment->get_X()->get_size(); i++)
  {
    _environment_metabolic_amounts_file << " " << _environment->get_X()->get(i+1);
  }
  _environment_metabolic_amounts_file << "\n";
}

/**
 * \brief    Write trophic network profile statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_trophic_network_profile_stats( void )
{
  _trophic_network_profile_file << _population->get_time() << " " <<
  _trophic_network->get_nb_level_0_cells() << " " <<
  _trophic_network->get_nb_level_1_cells() << " " <<
  _trophic_network->get_nb_level_2_cells() << " " <<
  _trophic_network->get_nb_no_level_cells() << " " <<
  _trophic_network->get_nb_level_0_groups() << " " <<
  _trophic_network->get_nb_level_1_groups() << " " <<
  _trophic_network->get_nb_level_2_groups() << " " <<
  _trophic_network->get_nb_no_level_groups() << " " <<
  _trophic_network->get_nb_group_appearances() << " " <<
  _trophic_network->get_nb_group_extinctions() << " " <<
  _trophic_network->get_mean_group_lifespan() << " " <<
  _trophic_network->get_min_true_diversity() << " " <<
  _trophic_network->get_min_max_time_distance() << " " <<
  _trophic_network->get_min_max_generation_distance() << " " <<
  _trophic_network->get_min_max_point_mutation_distance() << " " <<
  _trophic_network->get_min_max_hgt_distance() << " " <<
  _trophic_network->get_min_max_duplication_distance() << " " <<
  _trophic_network->get_min_max_deletion_distance() << " " <<
  _trophic_network->get_min_max_inversion_distance() << " " <<
  _trophic_network->get_min_max_translocation_distance() << " " <<
  _trophic_network->get_mean_true_diversity() << " " <<
  _trophic_network->get_mean_max_time_distance() << " " <<
  _trophic_network->get_mean_max_generation_distance() << " " <<
  _trophic_network->get_mean_max_point_mutation_distance() << " " <<
  _trophic_network->get_mean_max_hgt_distance() << " " <<
  _trophic_network->get_mean_max_duplication_distance() << " " <<
  _trophic_network->get_mean_max_deletion_distance() << " " <<
  _trophic_network->get_mean_max_inversion_distance() << " " <<
  _trophic_network->get_mean_max_translocation_distance() << " " <<
  _trophic_network->get_max_true_diversity() << " " <<
  _trophic_network->get_max_max_time_distance() << " " <<
  _trophic_network->get_max_max_generation_distance() << " " <<
  _trophic_network->get_max_max_point_mutation_distance() << " " <<
  _trophic_network->get_max_max_hgt_distance() << " " <<
  _trophic_network->get_max_max_duplication_distance() << " " <<
  _trophic_network->get_max_max_deletion_distance() << " " <<
  _trophic_network->get_max_max_inversion_distance() << " " <<
  _trophic_network->get_max_max_translocation_distance() << "\n";
}

/**
 * \brief    Write level 0 statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_level_0_stats( void )
{
  _level_0_statistics_file << _population->get_time() << " " <<
  
  /* GENETIC DIVERSITY */
  _trophic_network->get_level_0_true_diversity() << " " <<
  _trophic_network->get_level_0_max_time_distance() << " " <<
  _trophic_network->get_level_0_max_generation_distance() << " " <<
  _trophic_network->get_level_0_max_point_mutation_distance() << " " <<
  _trophic_network->get_level_0_max_hgt_distance() << " " <<
  _trophic_network->get_level_0_max_duplication_distance() << " " <<
  _trophic_network->get_level_0_max_deletion_distance() << " " <<
  _trophic_network->get_level_0_max_inversion_distance() << " " <<
  _trophic_network->get_level_0_max_translocation_distance() << " " <<
  
  /* PHENOTYPE */
  _trophic_network->get_level_0_generations() << " " <<
  _trophic_network->get_level_0_inherited_TF_amount() << " " <<
  _trophic_network->get_level_0_inherited_E_amount() << " " <<
  _trophic_network->get_level_0_TF_amount() << " " <<
  _trophic_network->get_level_0_E_amount() << " " <<
  _trophic_network->get_level_0_inherited_metabolic_amount() << " " <<
  _trophic_network->get_level_0_metabolic_amount() << " " <<
  _trophic_network->get_level_0_energy() << " " <<
  _trophic_network->get_level_0_score() << " " <<
  _trophic_network->get_level_0_lifespan() << " " <<
  _trophic_network->get_level_0_number_of_divisions() << " " <<
  _trophic_network->get_level_0_toxicity() << " " <<
  _trophic_network->get_level_0_metabolic_uptake() << " " <<
  _trophic_network->get_level_0_metabolic_release() << " " <<
  _trophic_network->get_level_0_metabolic_growth_rate() << " " <<
  _trophic_network->get_level_0_Dmetabolic_growth_rate() << " " <<
  _trophic_network->get_level_0_grn_nb_nodes() << " " <<
  _trophic_network->get_level_0_grn_nb_edges() << " " <<
  _trophic_network->get_level_0_metabolic_nb_nodes() << " " <<
  _trophic_network->get_level_0_metabolic_nb_edges() << " " <<
  _trophic_network->get_level_0_regulation_redundancy() << " " <<
  _trophic_network->get_level_0_metabolic_redundancy() << " " <<
  
  /* GENOME STRUCTURE */
  _trophic_network->get_level_0_genome_size() << " " <<
  _trophic_network->get_level_0_functional_size() << " " <<
  _trophic_network->get_level_0_genome_nb_NC() << " " <<
  _trophic_network->get_level_0_genome_nb_E() << " " <<
  _trophic_network->get_level_0_genome_nb_TF() << " " <<
  _trophic_network->get_level_0_genome_nb_BS() << " " <<
  _trophic_network->get_level_0_genome_nb_P() << " " <<
  _trophic_network->get_level_0_genome_nb_inner_enzymes() << " " <<
  _trophic_network->get_level_0_genome_nb_inflow_pumps() << " " <<
  _trophic_network->get_level_0_genome_nb_outflow_pumps() << " " <<
  _trophic_network->get_level_0_genome_nb_functional_regions() << " " <<
  _trophic_network->get_level_0_genome_nb_enhancers() << " " <<
  _trophic_network->get_level_0_genome_nb_operators() << " " <<
  _trophic_network->get_level_0_genome_nb_E_regions() << " " <<
  _trophic_network->get_level_0_genome_nb_TF_regions() << " " <<
  _trophic_network->get_level_0_genome_nb_mixed_regions() << " " <<
  _trophic_network->get_level_0_genome_functional_region_size() << " " <<
  _trophic_network->get_level_0_genome_E_region_size() << " " <<
  _trophic_network->get_level_0_genome_TF_region_size() << " " <<
  _trophic_network->get_level_0_genome_mixed_region_size() << " " <<
  _trophic_network->get_level_0_genome_enhancer_size() << " " <<
  _trophic_network->get_level_0_genome_operator_size() << " " <<
  _trophic_network->get_level_0_genome_operon_size() << " " <<
  _trophic_network->get_level_0_genome_E_operon_size() << " " <<
  _trophic_network->get_level_0_genome_TF_operon_size() << " " <<
  _trophic_network->get_level_0_genome_mixed_operon_size() << " " <<
  
  /* INHERITED STRUCTURE */
  _trophic_network->get_level_0_inherited_size() << " " <<
  _trophic_network->get_level_0_inherited_nb_E() << " " <<
  _trophic_network->get_level_0_inherited_nb_TF() << " " <<
  _trophic_network->get_level_0_inherited_nb_inner_enzymes() << " " <<
  _trophic_network->get_level_0_inherited_nb_inflow_pumps() << " " <<
  _trophic_network->get_level_0_inherited_nb_outflow_pumps() << "\n";
}

/**
 * \brief    Write level 1 statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_level_1_stats( void )
{
  _level_1_statistics_file << _population->get_time() << " " <<
  
  /* GENETIC DIVERSITY */
  _trophic_network->get_level_1_true_diversity() << " " <<
  _trophic_network->get_level_1_max_time_distance() << " " <<
  _trophic_network->get_level_1_max_generation_distance() << " " <<
  _trophic_network->get_level_1_max_point_mutation_distance() << " " <<
  _trophic_network->get_level_1_max_hgt_distance() << " " <<
  _trophic_network->get_level_1_max_duplication_distance() << " " <<
  _trophic_network->get_level_1_max_deletion_distance() << " " <<
  _trophic_network->get_level_1_max_inversion_distance() << " " <<
  _trophic_network->get_level_1_max_translocation_distance() << " " <<
  
  /* PHENOTYPE */
  _trophic_network->get_level_1_generations() << " " <<
  _trophic_network->get_level_1_inherited_TF_amount() << " " <<
  _trophic_network->get_level_1_inherited_E_amount() << " " <<
  _trophic_network->get_level_1_TF_amount() << " " <<
  _trophic_network->get_level_1_E_amount() << " " <<
  _trophic_network->get_level_1_inherited_metabolic_amount() << " " <<
  _trophic_network->get_level_1_metabolic_amount() << " " <<
  _trophic_network->get_level_1_energy() << " " <<
  _trophic_network->get_level_1_score() << " " <<
  _trophic_network->get_level_1_lifespan() << " " <<
  _trophic_network->get_level_1_number_of_divisions() << " " <<
  _trophic_network->get_level_1_toxicity() << " " <<
  _trophic_network->get_level_1_metabolic_uptake() << " " <<
  _trophic_network->get_level_1_metabolic_release() << " " <<
  _trophic_network->get_level_1_metabolic_growth_rate() << " " <<
  _trophic_network->get_level_1_Dmetabolic_growth_rate() << " " <<
  _trophic_network->get_level_1_grn_nb_nodes() << " " <<
  _trophic_network->get_level_1_grn_nb_edges() << " " <<
  _trophic_network->get_level_1_metabolic_nb_nodes() << " " <<
  _trophic_network->get_level_1_metabolic_nb_edges() << " " <<
  _trophic_network->get_level_1_regulation_redundancy() << " " <<
  _trophic_network->get_level_1_metabolic_redundancy() << " " <<
  
  /* GENOME STRUCTURE */
  _trophic_network->get_level_1_genome_size() << " " <<
  _trophic_network->get_level_1_functional_size() << " " <<
  _trophic_network->get_level_1_genome_nb_NC() << " " <<
  _trophic_network->get_level_1_genome_nb_E() << " " <<
  _trophic_network->get_level_1_genome_nb_TF() << " " <<
  _trophic_network->get_level_1_genome_nb_BS() << " " <<
  _trophic_network->get_level_1_genome_nb_P() << " " <<
  _trophic_network->get_level_1_genome_nb_inner_enzymes() << " " <<
  _trophic_network->get_level_1_genome_nb_inflow_pumps() << " " <<
  _trophic_network->get_level_1_genome_nb_outflow_pumps() << " " <<
  _trophic_network->get_level_1_genome_nb_functional_regions() << " " <<
  _trophic_network->get_level_1_genome_nb_enhancers() << " " <<
  _trophic_network->get_level_1_genome_nb_operators() << " " <<
  _trophic_network->get_level_1_genome_nb_E_regions() << " " <<
  _trophic_network->get_level_1_genome_nb_TF_regions() << " " <<
  _trophic_network->get_level_1_genome_nb_mixed_regions() << " " <<
  _trophic_network->get_level_1_genome_functional_region_size() << " " <<
  _trophic_network->get_level_1_genome_E_region_size() << " " <<
  _trophic_network->get_level_1_genome_TF_region_size() << " " <<
  _trophic_network->get_level_1_genome_mixed_region_size() << " " <<
  _trophic_network->get_level_1_genome_enhancer_size() << " " <<
  _trophic_network->get_level_1_genome_operator_size() << " " <<
  _trophic_network->get_level_1_genome_operon_size() << " " <<
  _trophic_network->get_level_1_genome_E_operon_size() << " " <<
  _trophic_network->get_level_1_genome_TF_operon_size() << " " <<
  _trophic_network->get_level_1_genome_mixed_operon_size() << " " <<
  
  /* INHERITED STRUCTURE */
  _trophic_network->get_level_1_inherited_size() << " " <<
  _trophic_network->get_level_1_inherited_nb_E() << " " <<
  _trophic_network->get_level_1_inherited_nb_TF() << " " <<
  _trophic_network->get_level_1_inherited_nb_inner_enzymes() << " " <<
  _trophic_network->get_level_1_inherited_nb_inflow_pumps() << " " <<
  _trophic_network->get_level_1_inherited_nb_outflow_pumps() << "\n";
}

/**
 * \brief    Write level 2 statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_level_2_stats( void )
{
  _level_2_statistics_file << _population->get_time() << " " <<
  
  /* GENETIC DIVERSITY */
  _trophic_network->get_level_2_true_diversity() << " " <<
  _trophic_network->get_level_2_max_time_distance() << " " <<
  _trophic_network->get_level_2_max_generation_distance() << " " <<
  _trophic_network->get_level_2_max_point_mutation_distance() << " " <<
  _trophic_network->get_level_2_max_hgt_distance() << " " <<
  _trophic_network->get_level_2_max_duplication_distance() << " " <<
  _trophic_network->get_level_2_max_deletion_distance() << " " <<
  _trophic_network->get_level_2_max_inversion_distance() << " " <<
  _trophic_network->get_level_2_max_translocation_distance() << " " <<
  
  /* PHENOTYPE */
  _trophic_network->get_level_2_generations() << " " <<
  _trophic_network->get_level_2_inherited_TF_amount() << " " <<
  _trophic_network->get_level_2_inherited_E_amount() << " " <<
  _trophic_network->get_level_2_TF_amount() << " " <<
  _trophic_network->get_level_2_E_amount() << " " <<
  _trophic_network->get_level_2_inherited_metabolic_amount() << " " <<
  _trophic_network->get_level_2_metabolic_amount() << " " <<
  _trophic_network->get_level_2_energy() << " " <<
  _trophic_network->get_level_2_score() << " " <<
  _trophic_network->get_level_2_lifespan() << " " <<
  _trophic_network->get_level_2_number_of_divisions() << " " <<
  _trophic_network->get_level_2_toxicity() << " " <<
  _trophic_network->get_level_2_metabolic_uptake() << " " <<
  _trophic_network->get_level_2_metabolic_release() << " " <<
  _trophic_network->get_level_2_metabolic_growth_rate() << " " <<
  _trophic_network->get_level_2_Dmetabolic_growth_rate() << " " <<
  _trophic_network->get_level_2_grn_nb_nodes() << " " <<
  _trophic_network->get_level_2_grn_nb_edges() << " " <<
  _trophic_network->get_level_2_metabolic_nb_nodes() << " " <<
  _trophic_network->get_level_2_metabolic_nb_edges() << " " <<
  _trophic_network->get_level_2_regulation_redundancy() << " " <<
  _trophic_network->get_level_2_metabolic_redundancy() << " " <<
  
  /* GENOME STRUCTURE */
  _trophic_network->get_level_2_genome_size() << " " <<
  _trophic_network->get_level_2_functional_size() << " " <<
  _trophic_network->get_level_2_genome_nb_NC() << " " <<
  _trophic_network->get_level_2_genome_nb_E() << " " <<
  _trophic_network->get_level_2_genome_nb_TF() << " " <<
  _trophic_network->get_level_2_genome_nb_BS() << " " <<
  _trophic_network->get_level_2_genome_nb_P() << " " <<
  _trophic_network->get_level_2_genome_nb_inner_enzymes() << " " <<
  _trophic_network->get_level_2_genome_nb_inflow_pumps() << " " <<
  _trophic_network->get_level_2_genome_nb_outflow_pumps() << " " <<
  _trophic_network->get_level_2_genome_nb_functional_regions() << " " <<
  _trophic_network->get_level_2_genome_nb_enhancers() << " " <<
  _trophic_network->get_level_2_genome_nb_operators() << " " <<
  _trophic_network->get_level_2_genome_nb_E_regions() << " " <<
  _trophic_network->get_level_2_genome_nb_TF_regions() << " " <<
  _trophic_network->get_level_2_genome_nb_mixed_regions() << " " <<
  _trophic_network->get_level_2_genome_functional_region_size() << " " <<
  _trophic_network->get_level_2_genome_E_region_size() << " " <<
  _trophic_network->get_level_2_genome_TF_region_size() << " " <<
  _trophic_network->get_level_2_genome_mixed_region_size() << " " <<
  _trophic_network->get_level_2_genome_enhancer_size() << " " <<
  _trophic_network->get_level_2_genome_operator_size() << " " <<
  _trophic_network->get_level_2_genome_operon_size() << " " <<
  _trophic_network->get_level_2_genome_E_operon_size() << " " <<
  _trophic_network->get_level_2_genome_TF_operon_size() << " " <<
  _trophic_network->get_level_2_genome_mixed_operon_size() << " " <<
  
  /* INHERITED STRUCTURE */
  _trophic_network->get_level_2_inherited_size() << " " <<
  _trophic_network->get_level_2_inherited_nb_E() << " " <<
  _trophic_network->get_level_2_inherited_nb_TF() << " " <<
  _trophic_network->get_level_2_inherited_nb_inner_enzymes() << " " <<
  _trophic_network->get_level_2_inherited_nb_inflow_pumps() << " " <<
  _trophic_network->get_level_2_inherited_nb_outflow_pumps() << "\n";
}

/**
 * \brief    Write no level statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_no_level_stats( void )
{
  _no_level_statistics_file << _population->get_time() << " " <<
  
  /* GENETIC DIVERSITY */
  _trophic_network->get_no_level_true_diversity() << " " <<
  _trophic_network->get_no_level_max_time_distance() << " " <<
  _trophic_network->get_no_level_max_generation_distance() << " " <<
  _trophic_network->get_no_level_max_point_mutation_distance() << " " <<
  _trophic_network->get_no_level_max_hgt_distance() << " " <<
  _trophic_network->get_no_level_max_duplication_distance() << " " <<
  _trophic_network->get_no_level_max_deletion_distance() << " " <<
  _trophic_network->get_no_level_max_inversion_distance() << " " <<
  _trophic_network->get_no_level_max_translocation_distance() << " " <<
  
  /* PHENOTYPE */
  _trophic_network->get_no_level_generations() << " " <<
  _trophic_network->get_no_level_inherited_TF_amount() << " " <<
  _trophic_network->get_no_level_inherited_E_amount() << " " <<
  _trophic_network->get_no_level_TF_amount() << " " <<
  _trophic_network->get_no_level_E_amount() << " " <<
  _trophic_network->get_no_level_inherited_metabolic_amount() << " " <<
  _trophic_network->get_no_level_metabolic_amount() << " " <<
  _trophic_network->get_no_level_energy() << " " <<
  _trophic_network->get_no_level_score() << " " <<
  _trophic_network->get_no_level_lifespan() << " " <<
  _trophic_network->get_no_level_number_of_divisions() << " " <<
  _trophic_network->get_no_level_toxicity() << " " <<
  _trophic_network->get_no_level_metabolic_uptake() << " " <<
  _trophic_network->get_no_level_metabolic_release() << " " <<
  _trophic_network->get_no_level_metabolic_growth_rate() << " " <<
  _trophic_network->get_no_level_Dmetabolic_growth_rate() << " " <<
  _trophic_network->get_no_level_grn_nb_nodes() << " " <<
  _trophic_network->get_no_level_grn_nb_edges() << " " <<
  _trophic_network->get_no_level_metabolic_nb_nodes() << " " <<
  _trophic_network->get_no_level_metabolic_nb_edges() << " " <<
  _trophic_network->get_no_level_regulation_redundancy() << " " <<
  _trophic_network->get_no_level_metabolic_redundancy() << " " <<
  
  /* GENOME STRUCTURE */
  _trophic_network->get_no_level_genome_size() << " " <<
  _trophic_network->get_no_level_functional_size() << " " <<
  _trophic_network->get_no_level_genome_nb_NC() << " " <<
  _trophic_network->get_no_level_genome_nb_E() << " " <<
  _trophic_network->get_no_level_genome_nb_TF() << " " <<
  _trophic_network->get_no_level_genome_nb_BS() << " " <<
  _trophic_network->get_no_level_genome_nb_P() << " " <<
  _trophic_network->get_no_level_genome_nb_inner_enzymes() << " " <<
  _trophic_network->get_no_level_genome_nb_inflow_pumps() << " " <<
  _trophic_network->get_no_level_genome_nb_outflow_pumps() << " " <<
  _trophic_network->get_no_level_genome_nb_functional_regions() << " " <<
  _trophic_network->get_no_level_genome_nb_enhancers() << " " <<
  _trophic_network->get_no_level_genome_nb_operators() << " " <<
  _trophic_network->get_no_level_genome_nb_E_regions() << " " <<
  _trophic_network->get_no_level_genome_nb_TF_regions() << " " <<
  _trophic_network->get_no_level_genome_nb_mixed_regions() << " " <<
  _trophic_network->get_no_level_genome_functional_region_size() << " " <<
  _trophic_network->get_no_level_genome_E_region_size() << " " <<
  _trophic_network->get_no_level_genome_TF_region_size() << " " <<
  _trophic_network->get_no_level_genome_mixed_region_size() << " " <<
  _trophic_network->get_no_level_genome_enhancer_size() << " " <<
  _trophic_network->get_no_level_genome_operator_size() << " " <<
  _trophic_network->get_no_level_genome_operon_size() << " " <<
  _trophic_network->get_no_level_genome_E_operon_size() << " " <<
  _trophic_network->get_no_level_genome_TF_operon_size() << " " <<
  _trophic_network->get_no_level_genome_mixed_operon_size() << " " <<
  
  /* INHERITED STRUCTURE */
  _trophic_network->get_no_level_inherited_size() << " " <<
  _trophic_network->get_no_level_inherited_nb_E() << " " <<
  _trophic_network->get_no_level_inherited_nb_TF() << " " <<
  _trophic_network->get_no_level_inherited_nb_inner_enzymes() << " " <<
  _trophic_network->get_no_level_inherited_nb_inflow_pumps() << " " <<
  _trophic_network->get_no_level_inherited_nb_outflow_pumps() << "\n";
}

/**
 * \brief    Write global concentrations statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_global_concentrations_stats( void )
{
  double pop_amount = _population->get_total_amount();
  double env_amount = _environment->get_total_amount();
  _global_concentrations_file << _population->get_time() << " " <<
  pop_amount+env_amount << " " <<
  _environment->get_inflowing_amount() << " " <<
  _environment->get_outflowing_amount() << " " <<
  pop_amount << " " <<
  env_amount << "\n";
}

/**
 * \brief    Write tree structure data
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_tree_structure_stats( void )
{
  /*----------------------------------*/
  /* 1) Explore the lineage tree      */
  /*----------------------------------*/
  size_t lineage_nb_nodes  = 0;
  size_t lineage_nb_roots  = 0;
  size_t lineage_nb_normal = 0;
  size_t lineage_nb_dead   = 0;
  size_t lineage_nb_alive  = 0;
  Node* current_node = _lineage_tree->get_first_node();
  while (current_node != NULL)
  {
    lineage_nb_nodes++;
    if (current_node->isRoot())
    {
      lineage_nb_roots++;
    }
    else if (current_node->isNormal())
    {
      lineage_nb_normal++;
    }
    if (current_node->isDead())
    {
      lineage_nb_dead++;
    }
    else if (current_node->isAlive())
    {
      lineage_nb_alive++;
    }
    current_node = _lineage_tree->get_next_node();
  }
  
  /*----------------------------------*/
  /* 2) Explore the phylogenetic tree */
  /*----------------------------------*/
  size_t phylogeny_nb_nodes  = 0;
  size_t phylogeny_nb_roots  = 0;
  size_t phylogeny_nb_normal = 0;
  size_t phylogeny_nb_dead   = 0;
  size_t phylogeny_nb_alive  = 0;
  current_node = _phylogenetic_tree->get_first_node();
  while (current_node != NULL)
  {
    phylogeny_nb_nodes++;
    if (current_node->isRoot())
    {
      phylogeny_nb_roots++;
    }
    else if (current_node->isNormal())
    {
      phylogeny_nb_normal++;
    }
    if (current_node->isDead())
    {
      phylogeny_nb_dead++;
    }
    else if (current_node->isAlive())
    {
      phylogeny_nb_alive++;
    }
    current_node = _phylogenetic_tree->get_next_node();
  }
  
  /*----------------------------------*/
  /* 3) Write data                    */
  /*----------------------------------*/
  _tree_structure_file << _population->get_time() << " " <<
  lineage_nb_nodes << " " <<
  lineage_nb_roots << " " <<
  lineage_nb_normal << " " <<
  lineage_nb_dead << " " <<
  lineage_nb_alive << " " <<
  phylogeny_nb_nodes << " " <<
  phylogeny_nb_roots << " " <<
  phylogeny_nb_normal << " " <<
  phylogeny_nb_dead << " " <<
  phylogeny_nb_alive << " " <<
  _phylogenetic_tree->get_common_ancestor_age() << "\n";
}

/**
 * \brief    Write best identifier
 * \details  --
 * \param    void
 * \return   \e void
 */
void Statistics::write_best_id_file( void )
{
  _best_id_file << _population->get_time() << " " << _best_generations << " " << _best_id << "\n";
}
