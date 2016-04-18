
/**
 * \file      TrophicNetwork.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      10-09-2015
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     TrophicNetwork class definition
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

#include "TrophicNetwork.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \param    Population* population
 * \param    Environment* environment
 * \return   \e void
 */
TrophicNetwork::TrophicNetwork( Parameters* parameters, Population* population, Environment* environment )
{
  /*------------------------------------------------------------------ Trophic network attributes */
  
  _parameters  = parameters;
  _population  = population;
  _environment = environment;
  _current_id  = 0;
  _group_map.clear();
  _iterator = _group_map.begin();
  
  /*------------------------------------------------------------------ Trophic network statistics */
  
  _nb_level_0_groups  = 0;
  _nb_level_1_groups  = 0;
  _nb_level_2_groups  = 0;
  _nb_no_level_groups = 0;
  
  _nb_level_0_cells  = 0;
  _nb_level_1_cells  = 0;
  _nb_level_2_cells  = 0;
  _nb_no_level_cells = 0;
  
  _nb_group_appearances = 0;
  _nb_group_extinctions = 0;
  
  _mean_group_lifespan = 0.0;
  
  _min_true_diversity              = 0.0;
  _min_max_time_distance           = 0.0;
  _min_max_generation_distance     = 0.0;
  _min_max_point_mutation_distance = 0.0;
  _min_max_hgt_distance            = 0.0;
  _min_max_duplication_distance    = 0.0;
  _min_max_deletion_distance       = 0.0;
  _min_max_inversion_distance      = 0.0;
  _min_max_translocation_distance  = 0.0;
  
  _mean_true_diversity              = 0.0;
  _mean_max_time_distance           = 0.0;
  _mean_max_generation_distance     = 0.0;
  _mean_max_point_mutation_distance = 0.0;
  _mean_max_hgt_distance            = 0.0;
  _mean_max_duplication_distance    = 0.0;
  _mean_max_deletion_distance       = 0.0;
  _mean_max_inversion_distance      = 0.0;
  _mean_max_translocation_distance  = 0.0;
  
  _max_true_diversity              = 0.0;
  _max_max_time_distance           = 0.0;
  _max_max_generation_distance     = 0.0;
  _max_max_point_mutation_distance = 0.0;
  _max_max_hgt_distance            = 0.0;
  _max_max_duplication_distance    = 0.0;
  _max_max_deletion_distance       = 0.0;
  _max_max_inversion_distance      = 0.0;
  _max_max_translocation_distance  = 0.0;
  
  /*------------------------------------------------------------------ LEVEL 0 statistics */
  
  /* GENETIC DIVERSITY */
  _level_0_true_diversity              = 0.0;
  _level_0_max_time_distance           = 0.0;
  _level_0_max_generation_distance     = 0.0;
  _level_0_max_point_mutation_distance = 0.0;
  _level_0_max_hgt_distance            = 0.0;
  _level_0_max_duplication_distance    = 0.0;
  _level_0_max_deletion_distance       = 0.0;
  _level_0_max_inversion_distance      = 0.0;
  _level_0_max_translocation_distance  = 0.0;
  
  /* PHENOTYPE */
  _level_0_generations                = 0.0;
  _level_0_inherited_TF_amount        = 0.0;
  _level_0_inherited_E_amount         = 0.0;
  _level_0_TF_amount                  = 0.0;
  _level_0_E_amount                   = 0.0;
  _level_0_inherited_metabolic_amount = 0.0;
  _level_0_metabolic_amount           = 0.0;
  _level_0_energy                     = 0.0;
  _level_0_score                      = 0.0;
  _level_0_lifespan                   = 0.0;
  _level_0_number_of_divisions        = 0.0;
  _level_0_toxicity                   = 0.0;
  _level_0_metabolic_uptake           = 0.0;
  _level_0_metabolic_release          = 0.0;
  _level_0_metabolic_growth_rate      = 0.0;
  _level_0_Dmetabolic_growth_rate     = 0.0;
  _level_0_grn_nb_nodes               = 0.0;
  _level_0_grn_nb_edges               = 0.0;
  _level_0_metabolic_nb_nodes         = 0.0;
  _level_0_metabolic_nb_edges         = 0.0;
  _level_0_regulation_redundancy      = 0.0;
  _level_0_metabolic_redundancy       = 0.0;
  
  /* GENOME STRUCTURE */
  _level_0_genome_size                   = 0.0;
  _level_0_functional_size               = 0.0;
  _level_0_genome_nb_NC                  = 0.0;
  _level_0_genome_nb_E                   = 0.0;
  _level_0_genome_nb_TF                  = 0.0;
  _level_0_genome_nb_BS                  = 0.0;
  _level_0_genome_nb_P                   = 0.0;
  _level_0_genome_nb_inner_enzymes       = 0.0;
  _level_0_genome_nb_inflow_pumps        = 0.0;
  _level_0_genome_nb_outflow_pumps       = 0.0;
  _level_0_genome_nb_functional_regions  = 0.0;
  _level_0_genome_nb_enhancers           = 0.0;
  _level_0_genome_nb_operators           = 0.0;
  _level_0_genome_nb_E_regions           = 0.0;
  _level_0_genome_nb_TF_regions          = 0.0;
  _level_0_genome_nb_mixed_regions       = 0.0;
  _level_0_genome_functional_region_size = 0.0;
  _level_0_genome_E_region_size          = 0.0;
  _level_0_genome_TF_region_size         = 0.0;
  _level_0_genome_mixed_region_size      = 0.0;
  _level_0_genome_enhancer_size          = 0.0;
  _level_0_genome_operator_size          = 0.0;
  _level_0_genome_operon_size            = 0.0;
  _level_0_genome_E_operon_size          = 0.0;
  _level_0_genome_TF_operon_size         = 0.0;
  _level_0_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _level_0_inherited_size             = 0.0;
  _level_0_inherited_nb_E             = 0.0;
  _level_0_inherited_nb_TF            = 0.0;
  _level_0_inherited_nb_inner_enzymes = 0.0;
  _level_0_inherited_nb_inflow_pumps  = 0.0;
  _level_0_inherited_nb_outflow_pumps = 0.0;
  
  /*------------------------------------------------------------------ LEVEL 1 statistics */
  
  /* GENETIC DIVERSITY */
  _level_1_true_diversity              = 0.0;
  _level_1_max_time_distance           = 0.0;
  _level_1_max_generation_distance     = 0.0;
  _level_1_max_point_mutation_distance = 0.0;
  _level_1_max_hgt_distance            = 0.0;
  _level_1_max_duplication_distance    = 0.0;
  _level_1_max_deletion_distance       = 0.0;
  _level_1_max_inversion_distance      = 0.0;
  _level_1_max_translocation_distance  = 0.0;
  
  /* PHENOTYPE */
  _level_1_generations                = 0.0;
  _level_1_inherited_TF_amount        = 0.0;
  _level_1_inherited_E_amount         = 0.0;
  _level_1_TF_amount                  = 0.0;
  _level_1_E_amount                   = 0.0;
  _level_1_inherited_metabolic_amount = 0.0;
  _level_1_metabolic_amount           = 0.0;
  _level_1_energy                     = 0.0;
  _level_1_score                      = 0.0;
  _level_1_lifespan                   = 0.0;
  _level_1_number_of_divisions        = 0.0;
  _level_1_toxicity                   = 0.0;
  _level_1_metabolic_uptake           = 0.0;
  _level_1_metabolic_release          = 0.0;
  _level_1_metabolic_growth_rate      = 0.0;
  _level_1_Dmetabolic_growth_rate     = 0.0;
  _level_1_grn_nb_nodes               = 0.0;
  _level_1_grn_nb_edges               = 0.0;
  _level_1_metabolic_nb_nodes         = 0.0;
  _level_1_metabolic_nb_edges         = 0.0;
  _level_1_regulation_redundancy      = 0.0;
  _level_1_metabolic_redundancy       = 0.0;
  
  /* GENOME STRUCTURE */
  _level_1_genome_size                   = 0.0;
  _level_1_functional_size               = 0.0;
  _level_1_genome_nb_NC                  = 0.0;
  _level_1_genome_nb_E                   = 0.0;
  _level_1_genome_nb_TF                  = 0.0;
  _level_1_genome_nb_BS                  = 0.0;
  _level_1_genome_nb_P                   = 0.0;
  _level_1_genome_nb_inner_enzymes       = 0.0;
  _level_1_genome_nb_inflow_pumps        = 0.0;
  _level_1_genome_nb_outflow_pumps       = 0.0;
  _level_1_genome_nb_functional_regions  = 0.0;
  _level_1_genome_nb_enhancers           = 0.0;
  _level_1_genome_nb_operators           = 0.0;
  _level_1_genome_nb_E_regions           = 0.0;
  _level_1_genome_nb_TF_regions          = 0.0;
  _level_1_genome_nb_mixed_regions       = 0.0;
  _level_1_genome_functional_region_size = 0.0;
  _level_1_genome_E_region_size          = 0.0;
  _level_1_genome_TF_region_size         = 0.0;
  _level_1_genome_mixed_region_size      = 0.0;
  _level_1_genome_enhancer_size          = 0.0;
  _level_1_genome_operator_size          = 0.0;
  _level_1_genome_operon_size            = 0.0;
  _level_1_genome_E_operon_size          = 0.0;
  _level_1_genome_TF_operon_size         = 0.0;
  _level_1_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _level_1_inherited_size             = 0.0;
  _level_1_inherited_nb_E             = 0.0;
  _level_1_inherited_nb_TF            = 0.0;
  _level_1_inherited_nb_inner_enzymes = 0.0;
  _level_1_inherited_nb_inflow_pumps  = 0.0;
  _level_1_inherited_nb_outflow_pumps = 0.0;
  
  /*------------------------------------------------------------------ LEVEL 2 statistics */
  
  /* GENETIC DIVERSITY */
  _level_2_true_diversity              = 0.0;
  _level_2_max_time_distance           = 0.0;
  _level_2_max_generation_distance     = 0.0;
  _level_2_max_point_mutation_distance = 0.0;
  _level_2_max_hgt_distance            = 0.0;
  _level_2_max_duplication_distance    = 0.0;
  _level_2_max_deletion_distance       = 0.0;
  _level_2_max_inversion_distance      = 0.0;
  _level_2_max_translocation_distance  = 0.0;
  
  /* PHENOTYPE */
  _level_2_generations                = 0.0;
  _level_2_inherited_TF_amount        = 0.0;
  _level_2_inherited_E_amount         = 0.0;
  _level_2_TF_amount                  = 0.0;
  _level_2_E_amount                   = 0.0;
  _level_2_inherited_metabolic_amount = 0.0;
  _level_2_metabolic_amount           = 0.0;
  _level_2_energy                     = 0.0;
  _level_2_score                      = 0.0;
  _level_2_lifespan                   = 0.0;
  _level_2_number_of_divisions        = 0.0;
  _level_2_toxicity                   = 0.0;
  _level_2_metabolic_uptake           = 0.0;
  _level_2_metabolic_release          = 0.0;
  _level_2_metabolic_growth_rate      = 0.0;
  _level_2_Dmetabolic_growth_rate     = 0.0;
  _level_2_grn_nb_nodes               = 0.0;
  _level_2_grn_nb_edges               = 0.0;
  _level_2_metabolic_nb_nodes         = 0.0;
  _level_2_metabolic_nb_edges         = 0.0;
  _level_2_regulation_redundancy      = 0.0;
  _level_2_metabolic_redundancy       = 0.0;
  
  /* GENOME STRUCTURE */
  _level_2_genome_size                   = 0.0;
  _level_2_functional_size               = 0.0;
  _level_2_genome_nb_NC                  = 0.0;
  _level_2_genome_nb_E                   = 0.0;
  _level_2_genome_nb_TF                  = 0.0;
  _level_2_genome_nb_BS                  = 0.0;
  _level_2_genome_nb_P                   = 0.0;
  _level_2_genome_nb_inner_enzymes       = 0.0;
  _level_2_genome_nb_inflow_pumps        = 0.0;
  _level_2_genome_nb_outflow_pumps       = 0.0;
  _level_2_genome_nb_functional_regions  = 0.0;
  _level_2_genome_nb_enhancers           = 0.0;
  _level_2_genome_nb_operators           = 0.0;
  _level_2_genome_nb_E_regions           = 0.0;
  _level_2_genome_nb_TF_regions          = 0.0;
  _level_2_genome_nb_mixed_regions       = 0.0;
  _level_2_genome_functional_region_size = 0.0;
  _level_2_genome_E_region_size          = 0.0;
  _level_2_genome_TF_region_size         = 0.0;
  _level_2_genome_mixed_region_size      = 0.0;
  _level_2_genome_enhancer_size          = 0.0;
  _level_2_genome_operator_size          = 0.0;
  _level_2_genome_operon_size            = 0.0;
  _level_2_genome_E_operon_size          = 0.0;
  _level_2_genome_TF_operon_size         = 0.0;
  _level_2_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _level_2_inherited_size             = 0.0;
  _level_2_inherited_nb_E             = 0.0;
  _level_2_inherited_nb_TF            = 0.0;
  _level_2_inherited_nb_inner_enzymes = 0.0;
  _level_2_inherited_nb_inflow_pumps  = 0.0;
  _level_2_inherited_nb_outflow_pumps = 0.0;

  /*------------------------------------------------------------------ NO LEVEL statistics */
  
  /* GENETIC DIVERSITY */
  _no_level_true_diversity              = 0.0;
  _no_level_max_time_distance           = 0.0;
  _no_level_max_generation_distance     = 0.0;
  _no_level_max_point_mutation_distance = 0.0;
  _no_level_max_hgt_distance            = 0.0;
  _no_level_max_duplication_distance    = 0.0;
  _no_level_max_deletion_distance       = 0.0;
  _no_level_max_inversion_distance      = 0.0;
  _no_level_max_translocation_distance  = 0.0;
  
  /* PHENOTYPE */
  _no_level_generations                = 0.0;
  _no_level_inherited_TF_amount        = 0.0;
  _no_level_inherited_E_amount         = 0.0;
  _no_level_TF_amount                  = 0.0;
  _no_level_E_amount                   = 0.0;
  _no_level_inherited_metabolic_amount = 0.0;
  _no_level_metabolic_amount           = 0.0;
  _no_level_energy                     = 0.0;
  _no_level_score                      = 0.0;
  _no_level_lifespan                   = 0.0;
  _no_level_number_of_divisions        = 0.0;
  _no_level_toxicity                   = 0.0;
  _no_level_metabolic_uptake           = 0.0;
  _no_level_metabolic_release          = 0.0;
  _no_level_metabolic_growth_rate      = 0.0;
  _no_level_Dmetabolic_growth_rate     = 0.0;
  _no_level_grn_nb_nodes               = 0.0;
  _no_level_grn_nb_edges               = 0.0;
  _no_level_metabolic_nb_nodes         = 0.0;
  _no_level_metabolic_nb_edges         = 0.0;
  _no_level_regulation_redundancy      = 0.0;
  _no_level_metabolic_redundancy       = 0.0;
  
  /* GENOME STRUCTURE */
  _no_level_genome_size                   = 0.0;
  _no_level_functional_size               = 0.0;
  _no_level_genome_nb_NC                  = 0.0;
  _no_level_genome_nb_E                   = 0.0;
  _no_level_genome_nb_TF                  = 0.0;
  _no_level_genome_nb_BS                  = 0.0;
  _no_level_genome_nb_P                   = 0.0;
  _no_level_genome_nb_inner_enzymes       = 0.0;
  _no_level_genome_nb_inflow_pumps        = 0.0;
  _no_level_genome_nb_outflow_pumps       = 0.0;
  _no_level_genome_nb_functional_regions  = 0.0;
  _no_level_genome_nb_enhancers           = 0.0;
  _no_level_genome_nb_operators           = 0.0;
  _no_level_genome_nb_E_regions           = 0.0;
  _no_level_genome_nb_TF_regions          = 0.0;
  _no_level_genome_nb_mixed_regions       = 0.0;
  _no_level_genome_functional_region_size = 0.0;
  _no_level_genome_E_region_size          = 0.0;
  _no_level_genome_TF_region_size         = 0.0;
  _no_level_genome_mixed_region_size      = 0.0;
  _no_level_genome_enhancer_size          = 0.0;
  _no_level_genome_operator_size          = 0.0;
  _no_level_genome_operon_size            = 0.0;
  _no_level_genome_E_operon_size          = 0.0;
  _no_level_genome_TF_operon_size         = 0.0;
  _no_level_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _no_level_inherited_size             = 0.0;
  _no_level_inherited_nb_E             = 0.0;
  _no_level_inherited_nb_TF            = 0.0;
  _no_level_inherited_nb_inner_enzymes = 0.0;
  _no_level_inherited_nb_inflow_pumps  = 0.0;
  _no_level_inherited_nb_outflow_pumps = 0.0;
  
}

/**
 * \brief    Constructor from backup file
 * \details  --
 * \param    Parameters* parameters
 * \param    Population* population
 * \param    Environment* environment
 * \param    gzFile backup_file
 * \return   \e void
 */
TrophicNetwork::TrophicNetwork( Parameters* parameters, Population* population, Environment* environment, gzFile backup_file )
{
  /*------------------------------------------------------------------ Trophic network attributes */
  
  _parameters  = parameters;
  _population  = population;
  _environment = environment;
  gzread( backup_file, &_current_id, sizeof(_current_id) );
  size_t n = 0;
  gzread( backup_file, &n, sizeof(n) );
  _group_map.clear();
  for (size_t i = 0; i < n; i++)
  {
    TrophicGroup* current_group = new TrophicGroup(backup_file);
    _group_map[current_group->get_identifier()] = current_group;
  }
  _iterator = _group_map.begin();
  
  /*------------------------------------------------------------------ Trophic network statistics */
  
  gzread( backup_file, &_nb_level_0_groups,  sizeof(_nb_level_0_groups) );
  gzread( backup_file, &_nb_level_1_groups,  sizeof(_nb_level_1_groups) );
  gzread( backup_file, &_nb_level_2_groups,  sizeof(_nb_level_2_groups) );
  gzread( backup_file, &_nb_no_level_groups, sizeof(_nb_no_level_groups) );
  
  gzread( backup_file, &_nb_level_0_cells,  sizeof(_nb_level_0_cells) );
  gzread( backup_file, &_nb_level_1_cells,  sizeof(_nb_level_1_cells) );
  gzread( backup_file, &_nb_level_2_cells,  sizeof(_nb_level_2_cells) );
  gzread( backup_file, &_nb_no_level_cells, sizeof(_nb_no_level_cells) );
  
  gzread( backup_file, &_nb_group_appearances, sizeof(_nb_group_appearances) );
  gzread( backup_file, &_nb_group_extinctions, sizeof(_nb_group_extinctions) );
  
  gzread( backup_file, &_mean_group_lifespan, sizeof(_mean_group_lifespan) );
  
  gzread( backup_file, &_min_true_diversity,              sizeof(_min_true_diversity) );
  gzread( backup_file, &_min_max_time_distance,           sizeof(_min_max_time_distance) );
  gzread( backup_file, &_min_max_generation_distance,     sizeof(_min_max_generation_distance) );
  gzread( backup_file, &_min_max_point_mutation_distance, sizeof(_min_max_point_mutation_distance) );
  gzread( backup_file, &_min_max_hgt_distance,            sizeof(_min_max_hgt_distance) );
  gzread( backup_file, &_min_max_duplication_distance,    sizeof(_min_max_duplication_distance) );
  gzread( backup_file, &_min_max_deletion_distance,       sizeof(_min_max_deletion_distance) );
  gzread( backup_file, &_min_max_inversion_distance,      sizeof(_min_max_inversion_distance) );
  gzread( backup_file, &_min_max_translocation_distance,  sizeof(_min_max_translocation_distance) );
  
  gzread( backup_file, &_mean_true_diversity,              sizeof(_mean_true_diversity) );
  gzread( backup_file, &_mean_max_time_distance,           sizeof(_mean_max_time_distance) );
  gzread( backup_file, &_mean_max_generation_distance,     sizeof(_mean_max_generation_distance) );
  gzread( backup_file, &_mean_max_point_mutation_distance, sizeof(_mean_max_point_mutation_distance) );
  gzread( backup_file, &_mean_max_hgt_distance,            sizeof(_mean_max_hgt_distance) );
  gzread( backup_file, &_mean_max_duplication_distance,    sizeof(_mean_max_duplication_distance) );
  gzread( backup_file, &_mean_max_deletion_distance,       sizeof(_mean_max_deletion_distance) );
  gzread( backup_file, &_mean_max_inversion_distance,      sizeof(_mean_max_inversion_distance) );
  gzread( backup_file, &_mean_max_translocation_distance,  sizeof(_mean_max_translocation_distance) );
  
  gzread( backup_file, &_max_true_diversity,              sizeof(_max_true_diversity) );
  gzread( backup_file, &_max_max_time_distance,           sizeof(_max_max_time_distance) );
  gzread( backup_file, &_max_max_generation_distance,     sizeof(_max_max_generation_distance) );
  gzread( backup_file, &_max_max_point_mutation_distance, sizeof(_max_max_point_mutation_distance) );
  gzread( backup_file, &_max_max_hgt_distance,            sizeof(_max_max_hgt_distance) );
  gzread( backup_file, &_max_max_duplication_distance,    sizeof(_max_max_duplication_distance) );
  gzread( backup_file, &_max_max_deletion_distance,       sizeof(_max_max_deletion_distance) );
  gzread( backup_file, &_max_max_inversion_distance,      sizeof(_max_max_inversion_distance) );
  gzread( backup_file, &_max_max_translocation_distance,  sizeof(_max_max_translocation_distance) );
  
  /*------------------------------------------------------------------ LEVEL 0 statistics */
  
  /* GENETIC DIVERSITY */
  gzread( backup_file, &_level_0_true_diversity,              sizeof(_level_0_true_diversity) );
  gzread( backup_file, &_level_0_max_time_distance,           sizeof(_level_0_max_time_distance) );
  gzread( backup_file, &_level_0_max_generation_distance,     sizeof(_level_0_max_generation_distance) );
  gzread( backup_file, &_level_0_max_point_mutation_distance, sizeof(_level_0_max_point_mutation_distance) );
  gzread( backup_file, &_level_0_max_hgt_distance,            sizeof(_level_0_max_hgt_distance) );
  gzread( backup_file, &_level_0_max_duplication_distance,    sizeof(_level_0_max_duplication_distance) );
  gzread( backup_file, &_level_0_max_deletion_distance,       sizeof(_level_0_max_deletion_distance) );
  gzread( backup_file, &_level_0_max_inversion_distance,      sizeof(_level_0_max_inversion_distance) );
  gzread( backup_file, &_level_0_max_translocation_distance,  sizeof(_level_0_max_translocation_distance) );
  
  /* PHENOTYPE */
  gzread( backup_file, &_level_0_generations,                sizeof(_level_0_generations) );
  gzread( backup_file, &_level_0_inherited_TF_amount,        sizeof(_level_0_inherited_TF_amount) );
  gzread( backup_file, &_level_0_inherited_E_amount,         sizeof(_level_0_inherited_E_amount) );
  gzread( backup_file, &_level_0_TF_amount,                  sizeof(_level_0_TF_amount) );
  gzread( backup_file, &_level_0_E_amount,                   sizeof(_level_0_E_amount) );
  gzread( backup_file, &_level_0_inherited_metabolic_amount, sizeof(_level_0_inherited_metabolic_amount) );
  gzread( backup_file, &_level_0_metabolic_amount,           sizeof(_level_0_metabolic_amount) );
  gzread( backup_file, &_level_0_energy,                     sizeof(_level_0_energy) );
  gzread( backup_file, &_level_0_score,                      sizeof(_level_0_score) );
  gzread( backup_file, &_level_0_lifespan,                   sizeof(_level_0_lifespan) );
  gzread( backup_file, &_level_0_number_of_divisions,        sizeof(_level_0_number_of_divisions) );
  gzread( backup_file, &_level_0_toxicity,                   sizeof(_level_0_toxicity) );
  gzread( backup_file, &_level_0_metabolic_uptake,           sizeof(_level_0_metabolic_uptake) );
  gzread( backup_file, &_level_0_metabolic_release,          sizeof(_level_0_metabolic_release) );
  gzread( backup_file, &_level_0_metabolic_growth_rate,      sizeof(_level_0_metabolic_growth_rate) );
  gzread( backup_file, &_level_0_Dmetabolic_growth_rate,     sizeof(_level_0_Dmetabolic_growth_rate) );
  gzread( backup_file, &_level_0_grn_nb_nodes,               sizeof(_level_0_grn_nb_nodes) );
  gzread( backup_file, &_level_0_grn_nb_edges,               sizeof(_level_0_grn_nb_edges) );
  gzread( backup_file, &_level_0_metabolic_nb_nodes,         sizeof(_level_0_metabolic_nb_nodes) );
  gzread( backup_file, &_level_0_metabolic_nb_edges,         sizeof(_level_0_metabolic_nb_edges) );
  gzread( backup_file, &_level_0_regulation_redundancy,      sizeof(_level_0_regulation_redundancy) );
  gzread( backup_file, &_level_0_metabolic_redundancy,       sizeof(_level_0_metabolic_redundancy) );
  
  /* GENOME STRUCTURE */
  gzread( backup_file, &_level_0_genome_size,                   sizeof(_level_0_genome_size) );
  gzread( backup_file, &_level_0_functional_size,               sizeof(_level_0_functional_size) );
  gzread( backup_file, &_level_0_genome_nb_NC,                  sizeof(_level_0_genome_nb_NC) );
  gzread( backup_file, &_level_0_genome_nb_E,                   sizeof(_level_0_genome_nb_E) );
  gzread( backup_file, &_level_0_genome_nb_TF,                  sizeof(_level_0_genome_nb_TF) );
  gzread( backup_file, &_level_0_genome_nb_BS,                  sizeof(_level_0_genome_nb_BS) );
  gzread( backup_file, &_level_0_genome_nb_P,                   sizeof(_level_0_genome_nb_P) );
  gzread( backup_file, &_level_0_genome_nb_inner_enzymes,       sizeof(_level_0_genome_nb_inner_enzymes) );
  gzread( backup_file, &_level_0_genome_nb_inflow_pumps,        sizeof(_level_0_genome_nb_inflow_pumps) );
  gzread( backup_file, &_level_0_genome_nb_outflow_pumps,       sizeof(_level_0_genome_nb_outflow_pumps) );
  gzread( backup_file, &_level_0_genome_nb_functional_regions,  sizeof(_level_0_genome_nb_functional_regions) );
  gzread( backup_file, &_level_0_genome_nb_enhancers,           sizeof(_level_0_genome_nb_enhancers) );
  gzread( backup_file, &_level_0_genome_nb_operators,           sizeof(_level_0_genome_nb_operators) );
  gzread( backup_file, &_level_0_genome_nb_E_regions,           sizeof(_level_0_genome_nb_E_regions) );
  gzread( backup_file, &_level_0_genome_nb_TF_regions,          sizeof(_level_0_genome_nb_TF_regions) );
  gzread( backup_file, &_level_0_genome_nb_mixed_regions,       sizeof(_level_0_genome_nb_mixed_regions) );
  gzread( backup_file, &_level_0_genome_functional_region_size, sizeof(_level_0_genome_functional_region_size) );
  gzread( backup_file, &_level_0_genome_E_region_size,          sizeof(_level_0_genome_E_region_size) );
  gzread( backup_file, &_level_0_genome_TF_region_size,         sizeof(_level_0_genome_TF_region_size) );
  gzread( backup_file, &_level_0_genome_mixed_region_size,      sizeof(_level_0_genome_mixed_region_size) );
  gzread( backup_file, &_level_0_genome_enhancer_size,          sizeof(_level_0_genome_enhancer_size) );
  gzread( backup_file, &_level_0_genome_operator_size,          sizeof(_level_0_genome_operator_size) );
  gzread( backup_file, &_level_0_genome_operon_size,            sizeof(_level_0_genome_operon_size) );
  gzread( backup_file, &_level_0_genome_E_operon_size,          sizeof(_level_0_genome_E_operon_size) );
  gzread( backup_file, &_level_0_genome_TF_operon_size,         sizeof(_level_0_genome_TF_operon_size) );
  gzread( backup_file, &_level_0_genome_mixed_operon_size,      sizeof(_level_0_genome_mixed_operon_size) );
  
  /* INHERITED STRUCTURE */
  gzread( backup_file, &_level_0_inherited_size,             sizeof(_level_0_inherited_size) );
  gzread( backup_file, &_level_0_inherited_nb_E,             sizeof(_level_0_inherited_nb_E) );
  gzread( backup_file, &_level_0_inherited_nb_TF,            sizeof(_level_0_inherited_nb_TF) );
  gzread( backup_file, &_level_0_inherited_nb_inner_enzymes, sizeof(_level_0_inherited_nb_inner_enzymes) );
  gzread( backup_file, &_level_0_inherited_nb_inflow_pumps,  sizeof(_level_0_inherited_nb_inflow_pumps) );
  gzread( backup_file, &_level_0_inherited_nb_outflow_pumps, sizeof(_level_0_inherited_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ LEVEL 1 statistics */
  
  /* GENETIC DIVERSITY */
  gzread( backup_file, &_level_1_true_diversity,              sizeof(_level_1_true_diversity) );
  gzread( backup_file, &_level_1_max_time_distance,           sizeof(_level_1_max_time_distance) );
  gzread( backup_file, &_level_1_max_generation_distance,     sizeof(_level_1_max_generation_distance) );
  gzread( backup_file, &_level_1_max_point_mutation_distance, sizeof(_level_1_max_point_mutation_distance) );
  gzread( backup_file, &_level_1_max_hgt_distance,            sizeof(_level_1_max_hgt_distance) );
  gzread( backup_file, &_level_1_max_duplication_distance,    sizeof(_level_1_max_duplication_distance) );
  gzread( backup_file, &_level_1_max_deletion_distance,       sizeof(_level_1_max_deletion_distance) );
  gzread( backup_file, &_level_1_max_inversion_distance,      sizeof(_level_1_max_inversion_distance) );
  gzread( backup_file, &_level_1_max_translocation_distance,  sizeof(_level_1_max_translocation_distance) );
  
  /* PHENOTYPE */
  gzread( backup_file, &_level_1_generations,                sizeof(_level_1_generations) );
  gzread( backup_file, &_level_1_inherited_TF_amount,        sizeof(_level_1_inherited_TF_amount) );
  gzread( backup_file, &_level_1_inherited_E_amount,         sizeof(_level_1_inherited_E_amount) );
  gzread( backup_file, &_level_1_TF_amount,                  sizeof(_level_1_TF_amount) );
  gzread( backup_file, &_level_1_E_amount,                   sizeof(_level_1_E_amount) );
  gzread( backup_file, &_level_1_inherited_metabolic_amount, sizeof(_level_1_inherited_metabolic_amount) );
  gzread( backup_file, &_level_1_metabolic_amount,           sizeof(_level_1_metabolic_amount) );
  gzread( backup_file, &_level_1_energy,                     sizeof(_level_1_energy) );
  gzread( backup_file, &_level_1_score,                      sizeof(_level_1_score) );
  gzread( backup_file, &_level_1_lifespan,                   sizeof(_level_1_lifespan) );
  gzread( backup_file, &_level_1_number_of_divisions,        sizeof(_level_1_number_of_divisions) );
  gzread( backup_file, &_level_1_toxicity,                   sizeof(_level_1_toxicity) );
  gzread( backup_file, &_level_1_metabolic_uptake,           sizeof(_level_1_metabolic_uptake) );
  gzread( backup_file, &_level_1_metabolic_release,          sizeof(_level_1_metabolic_release) );
  gzread( backup_file, &_level_1_metabolic_growth_rate,      sizeof(_level_1_metabolic_growth_rate) );
  gzread( backup_file, &_level_1_Dmetabolic_growth_rate,     sizeof(_level_1_Dmetabolic_growth_rate) );
  gzread( backup_file, &_level_1_grn_nb_nodes,               sizeof(_level_1_grn_nb_nodes) );
  gzread( backup_file, &_level_1_grn_nb_edges,               sizeof(_level_1_grn_nb_edges) );
  gzread( backup_file, &_level_1_metabolic_nb_nodes,         sizeof(_level_1_metabolic_nb_nodes) );
  gzread( backup_file, &_level_1_metabolic_nb_edges,         sizeof(_level_1_metabolic_nb_edges) );
  gzread( backup_file, &_level_1_regulation_redundancy,      sizeof(_level_1_regulation_redundancy) );
  gzread( backup_file, &_level_1_metabolic_redundancy,       sizeof(_level_1_metabolic_redundancy) );
  
  /* GENOME STRUCTURE */
  gzread( backup_file, &_level_1_genome_size,                   sizeof(_level_1_genome_size) );
  gzread( backup_file, &_level_1_functional_size,               sizeof(_level_1_functional_size) );
  gzread( backup_file, &_level_1_genome_nb_NC,                  sizeof(_level_1_genome_nb_NC) );
  gzread( backup_file, &_level_1_genome_nb_E,                   sizeof(_level_1_genome_nb_E) );
  gzread( backup_file, &_level_1_genome_nb_TF,                  sizeof(_level_1_genome_nb_TF) );
  gzread( backup_file, &_level_1_genome_nb_BS,                  sizeof(_level_1_genome_nb_BS) );
  gzread( backup_file, &_level_1_genome_nb_P,                   sizeof(_level_1_genome_nb_P) );
  gzread( backup_file, &_level_1_genome_nb_inner_enzymes,       sizeof(_level_1_genome_nb_inner_enzymes) );
  gzread( backup_file, &_level_1_genome_nb_inflow_pumps,        sizeof(_level_1_genome_nb_inflow_pumps) );
  gzread( backup_file, &_level_1_genome_nb_outflow_pumps,       sizeof(_level_1_genome_nb_outflow_pumps) );
  gzread( backup_file, &_level_1_genome_nb_functional_regions,  sizeof(_level_1_genome_nb_functional_regions) );
  gzread( backup_file, &_level_1_genome_nb_enhancers,           sizeof(_level_1_genome_nb_enhancers) );
  gzread( backup_file, &_level_1_genome_nb_operators,           sizeof(_level_1_genome_nb_operators) );
  gzread( backup_file, &_level_1_genome_nb_E_regions,           sizeof(_level_1_genome_nb_E_regions) );
  gzread( backup_file, &_level_1_genome_nb_TF_regions,          sizeof(_level_1_genome_nb_TF_regions) );
  gzread( backup_file, &_level_1_genome_nb_mixed_regions,       sizeof(_level_1_genome_nb_mixed_regions) );
  gzread( backup_file, &_level_1_genome_functional_region_size, sizeof(_level_1_genome_functional_region_size) );
  gzread( backup_file, &_level_1_genome_E_region_size,          sizeof(_level_1_genome_E_region_size) );
  gzread( backup_file, &_level_1_genome_TF_region_size,         sizeof(_level_1_genome_TF_region_size) );
  gzread( backup_file, &_level_1_genome_mixed_region_size,      sizeof(_level_1_genome_mixed_region_size) );
  gzread( backup_file, &_level_1_genome_enhancer_size,          sizeof(_level_1_genome_enhancer_size) );
  gzread( backup_file, &_level_1_genome_operator_size,          sizeof(_level_1_genome_operator_size) );
  gzread( backup_file, &_level_1_genome_operon_size,            sizeof(_level_1_genome_operon_size) );
  gzread( backup_file, &_level_1_genome_E_operon_size,          sizeof(_level_1_genome_E_operon_size) );
  gzread( backup_file, &_level_1_genome_TF_operon_size,         sizeof(_level_1_genome_TF_operon_size) );
  gzread( backup_file, &_level_1_genome_mixed_operon_size,      sizeof(_level_1_genome_mixed_operon_size) );
  
  /* INHERITED STRUCTURE */
  gzread( backup_file, &_level_1_inherited_size,             sizeof(_level_1_inherited_size) );
  gzread( backup_file, &_level_1_inherited_nb_E,             sizeof(_level_1_inherited_nb_E) );
  gzread( backup_file, &_level_1_inherited_nb_TF,            sizeof(_level_1_inherited_nb_TF) );
  gzread( backup_file, &_level_1_inherited_nb_inner_enzymes, sizeof(_level_1_inherited_nb_inner_enzymes) );
  gzread( backup_file, &_level_1_inherited_nb_inflow_pumps,  sizeof(_level_1_inherited_nb_inflow_pumps) );
  gzread( backup_file, &_level_1_inherited_nb_outflow_pumps, sizeof(_level_1_inherited_nb_outflow_pumps) );

  /*------------------------------------------------------------------ LEVEL 2 statistics */
  
  /* GENETIC DIVERSITY */
  gzread( backup_file, &_level_2_true_diversity,              sizeof(_level_2_true_diversity) );
  gzread( backup_file, &_level_2_max_time_distance,           sizeof(_level_2_max_time_distance) );
  gzread( backup_file, &_level_2_max_generation_distance,     sizeof(_level_2_max_generation_distance) );
  gzread( backup_file, &_level_2_max_point_mutation_distance, sizeof(_level_2_max_point_mutation_distance) );
  gzread( backup_file, &_level_2_max_hgt_distance,            sizeof(_level_2_max_hgt_distance) );
  gzread( backup_file, &_level_2_max_duplication_distance,    sizeof(_level_2_max_duplication_distance) );
  gzread( backup_file, &_level_2_max_deletion_distance,       sizeof(_level_2_max_deletion_distance) );
  gzread( backup_file, &_level_2_max_inversion_distance,      sizeof(_level_2_max_inversion_distance) );
  gzread( backup_file, &_level_2_max_translocation_distance,  sizeof(_level_2_max_translocation_distance) );
  
  /* PHENOTYPE */
  gzread( backup_file, &_level_2_generations,                sizeof(_level_2_generations) );
  gzread( backup_file, &_level_2_inherited_TF_amount,        sizeof(_level_2_inherited_TF_amount) );
  gzread( backup_file, &_level_2_inherited_E_amount,         sizeof(_level_2_inherited_E_amount) );
  gzread( backup_file, &_level_2_TF_amount,                  sizeof(_level_2_TF_amount) );
  gzread( backup_file, &_level_2_E_amount,                   sizeof(_level_2_E_amount) );
  gzread( backup_file, &_level_2_inherited_metabolic_amount, sizeof(_level_2_inherited_metabolic_amount) );
  gzread( backup_file, &_level_2_metabolic_amount,           sizeof(_level_2_metabolic_amount) );
  gzread( backup_file, &_level_2_energy,                     sizeof(_level_2_energy) );
  gzread( backup_file, &_level_2_score,                      sizeof(_level_2_score) );
  gzread( backup_file, &_level_2_lifespan,                   sizeof(_level_2_lifespan) );
  gzread( backup_file, &_level_2_number_of_divisions,        sizeof(_level_2_number_of_divisions) );
  gzread( backup_file, &_level_2_toxicity,                   sizeof(_level_2_toxicity) );
  gzread( backup_file, &_level_2_metabolic_uptake,           sizeof(_level_2_metabolic_uptake) );
  gzread( backup_file, &_level_2_metabolic_release,          sizeof(_level_2_metabolic_release) );
  gzread( backup_file, &_level_2_metabolic_growth_rate,      sizeof(_level_2_metabolic_growth_rate) );
  gzread( backup_file, &_level_2_Dmetabolic_growth_rate,     sizeof(_level_2_Dmetabolic_growth_rate) );
  gzread( backup_file, &_level_2_grn_nb_nodes,               sizeof(_level_2_grn_nb_nodes) );
  gzread( backup_file, &_level_2_grn_nb_edges,               sizeof(_level_2_grn_nb_edges) );
  gzread( backup_file, &_level_2_metabolic_nb_nodes,         sizeof(_level_2_metabolic_nb_nodes) );
  gzread( backup_file, &_level_2_metabolic_nb_edges,         sizeof(_level_2_metabolic_nb_edges) );
  gzread( backup_file, &_level_2_regulation_redundancy,      sizeof(_level_2_regulation_redundancy) );
  gzread( backup_file, &_level_2_metabolic_redundancy,       sizeof(_level_2_metabolic_redundancy) );
  
  /* GENOME STRUCTURE */
  gzread( backup_file, &_level_2_genome_size,                   sizeof(_level_2_genome_size) );
  gzread( backup_file, &_level_2_functional_size,               sizeof(_level_2_functional_size) );
  gzread( backup_file, &_level_2_genome_nb_NC,                  sizeof(_level_2_genome_nb_NC) );
  gzread( backup_file, &_level_2_genome_nb_E,                   sizeof(_level_2_genome_nb_E) );
  gzread( backup_file, &_level_2_genome_nb_TF,                  sizeof(_level_2_genome_nb_TF) );
  gzread( backup_file, &_level_2_genome_nb_BS,                  sizeof(_level_2_genome_nb_BS) );
  gzread( backup_file, &_level_2_genome_nb_P,                   sizeof(_level_2_genome_nb_P) );
  gzread( backup_file, &_level_2_genome_nb_inner_enzymes,       sizeof(_level_2_genome_nb_inner_enzymes) );
  gzread( backup_file, &_level_2_genome_nb_inflow_pumps,        sizeof(_level_2_genome_nb_inflow_pumps) );
  gzread( backup_file, &_level_2_genome_nb_outflow_pumps,       sizeof(_level_2_genome_nb_outflow_pumps) );
  gzread( backup_file, &_level_2_genome_nb_functional_regions,  sizeof(_level_2_genome_nb_functional_regions) );
  gzread( backup_file, &_level_2_genome_nb_enhancers,           sizeof(_level_2_genome_nb_enhancers) );
  gzread( backup_file, &_level_2_genome_nb_operators,           sizeof(_level_2_genome_nb_operators) );
  gzread( backup_file, &_level_2_genome_nb_E_regions,           sizeof(_level_2_genome_nb_E_regions) );
  gzread( backup_file, &_level_2_genome_nb_TF_regions,          sizeof(_level_2_genome_nb_TF_regions) );
  gzread( backup_file, &_level_2_genome_nb_mixed_regions,       sizeof(_level_2_genome_nb_mixed_regions) );
  gzread( backup_file, &_level_2_genome_functional_region_size, sizeof(_level_2_genome_functional_region_size) );
  gzread( backup_file, &_level_2_genome_E_region_size,          sizeof(_level_2_genome_E_region_size) );
  gzread( backup_file, &_level_2_genome_TF_region_size,         sizeof(_level_2_genome_TF_region_size) );
  gzread( backup_file, &_level_2_genome_mixed_region_size,      sizeof(_level_2_genome_mixed_region_size) );
  gzread( backup_file, &_level_2_genome_enhancer_size,          sizeof(_level_2_genome_enhancer_size) );
  gzread( backup_file, &_level_2_genome_operator_size,          sizeof(_level_2_genome_operator_size) );
  gzread( backup_file, &_level_2_genome_operon_size,            sizeof(_level_2_genome_operon_size) );
  gzread( backup_file, &_level_2_genome_E_operon_size,          sizeof(_level_2_genome_E_operon_size) );
  gzread( backup_file, &_level_2_genome_TF_operon_size,         sizeof(_level_2_genome_TF_operon_size) );
  gzread( backup_file, &_level_2_genome_mixed_operon_size,      sizeof(_level_2_genome_mixed_operon_size) );
  
  /* INHERITED STRUCTURE */
  gzread( backup_file, &_level_2_inherited_size,             sizeof(_level_2_inherited_size) );
  gzread( backup_file, &_level_2_inherited_nb_E,             sizeof(_level_2_inherited_nb_E) );
  gzread( backup_file, &_level_2_inherited_nb_TF,            sizeof(_level_2_inherited_nb_TF) );
  gzread( backup_file, &_level_2_inherited_nb_inner_enzymes, sizeof(_level_2_inherited_nb_inner_enzymes) );
  gzread( backup_file, &_level_2_inherited_nb_inflow_pumps,  sizeof(_level_2_inherited_nb_inflow_pumps) );
  gzread( backup_file, &_level_2_inherited_nb_outflow_pumps, sizeof(_level_2_inherited_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ NO LEVEL statistics */
  
  /* GENETIC DIVERSITY */
  gzread( backup_file, &_no_level_true_diversity,              sizeof(_no_level_true_diversity) );
  gzread( backup_file, &_no_level_max_time_distance,           sizeof(_no_level_max_time_distance) );
  gzread( backup_file, &_no_level_max_generation_distance,     sizeof(_no_level_max_generation_distance) );
  gzread( backup_file, &_no_level_max_point_mutation_distance, sizeof(_no_level_max_point_mutation_distance) );
  gzread( backup_file, &_no_level_max_hgt_distance,            sizeof(_no_level_max_hgt_distance) );
  gzread( backup_file, &_no_level_max_duplication_distance,    sizeof(_no_level_max_duplication_distance) );
  gzread( backup_file, &_no_level_max_deletion_distance,       sizeof(_no_level_max_deletion_distance) );
  gzread( backup_file, &_no_level_max_inversion_distance,      sizeof(_no_level_max_inversion_distance) );
  gzread( backup_file, &_no_level_max_translocation_distance,  sizeof(_no_level_max_translocation_distance) );
  
  /* PHENOTYPE */
  gzread( backup_file, &_no_level_generations,                sizeof(_no_level_generations) );
  gzread( backup_file, &_no_level_inherited_TF_amount,        sizeof(_no_level_inherited_TF_amount) );
  gzread( backup_file, &_no_level_inherited_E_amount,         sizeof(_no_level_inherited_E_amount) );
  gzread( backup_file, &_no_level_TF_amount,                  sizeof(_no_level_TF_amount) );
  gzread( backup_file, &_no_level_E_amount,                   sizeof(_no_level_E_amount) );
  gzread( backup_file, &_no_level_inherited_metabolic_amount, sizeof(_no_level_inherited_metabolic_amount) );
  gzread( backup_file, &_no_level_metabolic_amount,           sizeof(_no_level_metabolic_amount) );
  gzread( backup_file, &_no_level_energy,                     sizeof(_no_level_energy) );
  gzread( backup_file, &_no_level_score,                      sizeof(_no_level_score) );
  gzread( backup_file, &_no_level_lifespan,                   sizeof(_no_level_lifespan) );
  gzread( backup_file, &_no_level_number_of_divisions,        sizeof(_no_level_number_of_divisions) );
  gzread( backup_file, &_no_level_toxicity,                   sizeof(_no_level_toxicity) );
  gzread( backup_file, &_no_level_metabolic_uptake,           sizeof(_no_level_metabolic_uptake) );
  gzread( backup_file, &_no_level_metabolic_release,          sizeof(_no_level_metabolic_release) );
  gzread( backup_file, &_no_level_metabolic_growth_rate,      sizeof(_no_level_metabolic_growth_rate) );
  gzread( backup_file, &_no_level_Dmetabolic_growth_rate,     sizeof(_no_level_Dmetabolic_growth_rate) );
  gzread( backup_file, &_no_level_grn_nb_nodes,               sizeof(_no_level_grn_nb_nodes) );
  gzread( backup_file, &_no_level_grn_nb_edges,               sizeof(_no_level_grn_nb_edges) );
  gzread( backup_file, &_no_level_metabolic_nb_nodes,         sizeof(_no_level_metabolic_nb_nodes) );
  gzread( backup_file, &_no_level_metabolic_nb_edges,         sizeof(_no_level_metabolic_nb_edges) );
  gzread( backup_file, &_no_level_regulation_redundancy,      sizeof(_no_level_regulation_redundancy) );
  gzread( backup_file, &_no_level_metabolic_redundancy,       sizeof(_no_level_metabolic_redundancy) );
  
  /* GENOME STRUCTURE */
  gzread( backup_file, &_no_level_genome_size,                   sizeof(_no_level_genome_size) );
  gzread( backup_file, &_no_level_functional_size,               sizeof(_no_level_functional_size) );
  gzread( backup_file, &_no_level_genome_nb_NC,                  sizeof(_no_level_genome_nb_NC) );
  gzread( backup_file, &_no_level_genome_nb_E,                   sizeof(_no_level_genome_nb_E) );
  gzread( backup_file, &_no_level_genome_nb_TF,                  sizeof(_no_level_genome_nb_TF) );
  gzread( backup_file, &_no_level_genome_nb_BS,                  sizeof(_no_level_genome_nb_BS) );
  gzread( backup_file, &_no_level_genome_nb_P,                   sizeof(_no_level_genome_nb_P) );
  gzread( backup_file, &_no_level_genome_nb_inner_enzymes,       sizeof(_no_level_genome_nb_inner_enzymes) );
  gzread( backup_file, &_no_level_genome_nb_inflow_pumps,        sizeof(_no_level_genome_nb_inflow_pumps) );
  gzread( backup_file, &_no_level_genome_nb_outflow_pumps,       sizeof(_no_level_genome_nb_outflow_pumps) );
  gzread( backup_file, &_no_level_genome_nb_functional_regions,  sizeof(_no_level_genome_nb_functional_regions) );
  gzread( backup_file, &_no_level_genome_nb_enhancers,           sizeof(_no_level_genome_nb_enhancers) );
  gzread( backup_file, &_no_level_genome_nb_operators,           sizeof(_no_level_genome_nb_operators) );
  gzread( backup_file, &_no_level_genome_nb_E_regions,           sizeof(_no_level_genome_nb_E_regions) );
  gzread( backup_file, &_no_level_genome_nb_TF_regions,          sizeof(_no_level_genome_nb_TF_regions) );
  gzread( backup_file, &_no_level_genome_nb_mixed_regions,       sizeof(_no_level_genome_nb_mixed_regions) );
  gzread( backup_file, &_no_level_genome_functional_region_size, sizeof(_no_level_genome_functional_region_size) );
  gzread( backup_file, &_no_level_genome_E_region_size,          sizeof(_no_level_genome_E_region_size) );
  gzread( backup_file, &_no_level_genome_TF_region_size,         sizeof(_no_level_genome_TF_region_size) );
  gzread( backup_file, &_no_level_genome_mixed_region_size,      sizeof(_no_level_genome_mixed_region_size) );
  gzread( backup_file, &_no_level_genome_enhancer_size,          sizeof(_no_level_genome_enhancer_size) );
  gzread( backup_file, &_no_level_genome_operator_size,          sizeof(_no_level_genome_operator_size) );
  gzread( backup_file, &_no_level_genome_operon_size,            sizeof(_no_level_genome_operon_size) );
  gzread( backup_file, &_no_level_genome_E_operon_size,          sizeof(_no_level_genome_E_operon_size) );
  gzread( backup_file, &_no_level_genome_TF_operon_size,         sizeof(_no_level_genome_TF_operon_size) );
  gzread( backup_file, &_no_level_genome_mixed_operon_size,      sizeof(_no_level_genome_mixed_operon_size) );
  
  /* INHERITED STRUCTURE */
  gzread( backup_file, &_no_level_inherited_size,             sizeof(_no_level_inherited_size) );
  gzread( backup_file, &_no_level_inherited_nb_E,             sizeof(_no_level_inherited_nb_E) );
  gzread( backup_file, &_no_level_inherited_nb_TF,            sizeof(_no_level_inherited_nb_TF) );
  gzread( backup_file, &_no_level_inherited_nb_inner_enzymes, sizeof(_no_level_inherited_nb_inner_enzymes) );
  gzread( backup_file, &_no_level_inherited_nb_inflow_pumps,  sizeof(_no_level_inherited_nb_inflow_pumps) );
  gzread( backup_file, &_no_level_inherited_nb_outflow_pumps, sizeof(_no_level_inherited_nb_outflow_pumps) );

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
TrophicNetwork::~TrophicNetwork( void )
{
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    delete _iterator->second;
    _iterator->second = NULL;
  }
  _group_map.clear();
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void TrophicNetwork::save( gzFile backup_file )
{
  /*------------------------------------------------------------------ Trophic network attributes */
  
  gzwrite( backup_file, &_current_id, sizeof(_current_id) );
  size_t n = _group_map.size();
  gzwrite( backup_file, &n, sizeof(n) );
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    _iterator->second->save(backup_file);
  }
  
  /*------------------------------------------------------------------ Trophic network statistics */
  
  gzwrite( backup_file, &_nb_level_0_groups,  sizeof(_nb_level_0_groups) );
  gzwrite( backup_file, &_nb_level_1_groups,  sizeof(_nb_level_1_groups) );
  gzwrite( backup_file, &_nb_level_2_groups,  sizeof(_nb_level_2_groups) );
  gzwrite( backup_file, &_nb_no_level_groups, sizeof(_nb_no_level_groups) );
  
  gzwrite( backup_file, &_nb_level_0_cells,  sizeof(_nb_level_0_cells) );
  gzwrite( backup_file, &_nb_level_1_cells,  sizeof(_nb_level_1_cells) );
  gzwrite( backup_file, &_nb_level_2_cells,  sizeof(_nb_level_2_cells) );
  gzwrite( backup_file, &_nb_no_level_cells, sizeof(_nb_no_level_cells) );
  
  gzwrite( backup_file, &_nb_group_appearances, sizeof(_nb_group_appearances) );
  gzwrite( backup_file, &_nb_group_extinctions, sizeof(_nb_group_extinctions) );
  
  gzwrite( backup_file, &_mean_group_lifespan, sizeof(_mean_group_lifespan) );
  
  gzwrite( backup_file, &_min_true_diversity,              sizeof(_min_true_diversity) );
  gzwrite( backup_file, &_min_max_time_distance,           sizeof(_min_max_time_distance) );
  gzwrite( backup_file, &_min_max_generation_distance,     sizeof(_min_max_generation_distance) );
  gzwrite( backup_file, &_min_max_point_mutation_distance, sizeof(_min_max_point_mutation_distance) );
  gzwrite( backup_file, &_min_max_hgt_distance,            sizeof(_min_max_hgt_distance) );
  gzwrite( backup_file, &_min_max_duplication_distance,    sizeof(_min_max_duplication_distance) );
  gzwrite( backup_file, &_min_max_deletion_distance,       sizeof(_min_max_deletion_distance) );
  gzwrite( backup_file, &_min_max_inversion_distance,      sizeof(_min_max_inversion_distance) );
  gzwrite( backup_file, &_min_max_translocation_distance,  sizeof(_min_max_translocation_distance) );
  
  gzwrite( backup_file, &_mean_true_diversity,              sizeof(_mean_true_diversity) );
  gzwrite( backup_file, &_mean_max_time_distance,           sizeof(_mean_max_time_distance) );
  gzwrite( backup_file, &_mean_max_generation_distance,     sizeof(_mean_max_generation_distance) );
  gzwrite( backup_file, &_mean_max_point_mutation_distance, sizeof(_mean_max_point_mutation_distance) );
  gzwrite( backup_file, &_mean_max_hgt_distance,            sizeof(_mean_max_hgt_distance) );
  gzwrite( backup_file, &_mean_max_duplication_distance,    sizeof(_mean_max_duplication_distance) );
  gzwrite( backup_file, &_mean_max_deletion_distance,       sizeof(_mean_max_deletion_distance) );
  gzwrite( backup_file, &_mean_max_inversion_distance,      sizeof(_mean_max_inversion_distance) );
  gzwrite( backup_file, &_mean_max_translocation_distance,  sizeof(_mean_max_translocation_distance) );
  
  gzwrite( backup_file, &_max_true_diversity,              sizeof(_max_true_diversity) );
  gzwrite( backup_file, &_max_max_time_distance,           sizeof(_max_max_time_distance) );
  gzwrite( backup_file, &_max_max_generation_distance,     sizeof(_max_max_generation_distance) );
  gzwrite( backup_file, &_max_max_point_mutation_distance, sizeof(_max_max_point_mutation_distance) );
  gzwrite( backup_file, &_max_max_hgt_distance,            sizeof(_max_max_hgt_distance) );
  gzwrite( backup_file, &_max_max_duplication_distance,    sizeof(_max_max_duplication_distance) );
  gzwrite( backup_file, &_max_max_deletion_distance,       sizeof(_max_max_deletion_distance) );
  gzwrite( backup_file, &_max_max_inversion_distance,      sizeof(_max_max_inversion_distance) );
  gzwrite( backup_file, &_max_max_translocation_distance,  sizeof(_max_max_translocation_distance) );
  
  /*------------------------------------------------------------------ LEVEL 0 statistics */
  
  /* GENETIC DIVERSITY */
  gzwrite( backup_file, &_level_0_true_diversity,              sizeof(_level_0_true_diversity) );
  gzwrite( backup_file, &_level_0_max_time_distance,           sizeof(_level_0_max_time_distance) );
  gzwrite( backup_file, &_level_0_max_generation_distance,     sizeof(_level_0_max_generation_distance) );
  gzwrite( backup_file, &_level_0_max_point_mutation_distance, sizeof(_level_0_max_point_mutation_distance) );
  gzwrite( backup_file, &_level_0_max_hgt_distance,            sizeof(_level_0_max_hgt_distance) );
  gzwrite( backup_file, &_level_0_max_duplication_distance,    sizeof(_level_0_max_duplication_distance) );
  gzwrite( backup_file, &_level_0_max_deletion_distance,       sizeof(_level_0_max_deletion_distance) );
  gzwrite( backup_file, &_level_0_max_inversion_distance,      sizeof(_level_0_max_inversion_distance) );
  gzwrite( backup_file, &_level_0_max_translocation_distance,  sizeof(_level_0_max_translocation_distance) );
  
  /* PHENOTYPE */
  gzwrite( backup_file, &_level_0_generations,                sizeof(_level_0_generations) );
  gzwrite( backup_file, &_level_0_inherited_TF_amount,        sizeof(_level_0_inherited_TF_amount) );
  gzwrite( backup_file, &_level_0_inherited_E_amount,         sizeof(_level_0_inherited_E_amount) );
  gzwrite( backup_file, &_level_0_TF_amount,                  sizeof(_level_0_TF_amount) );
  gzwrite( backup_file, &_level_0_E_amount,                   sizeof(_level_0_E_amount) );
  gzwrite( backup_file, &_level_0_inherited_metabolic_amount, sizeof(_level_0_inherited_metabolic_amount) );
  gzwrite( backup_file, &_level_0_metabolic_amount,           sizeof(_level_0_metabolic_amount) );
  gzwrite( backup_file, &_level_0_energy,                     sizeof(_level_0_energy) );
  gzwrite( backup_file, &_level_0_score,                      sizeof(_level_0_score) );
  gzwrite( backup_file, &_level_0_lifespan,                   sizeof(_level_0_lifespan) );
  gzwrite( backup_file, &_level_0_number_of_divisions,        sizeof(_level_0_number_of_divisions) );
  gzwrite( backup_file, &_level_0_toxicity,                   sizeof(_level_0_toxicity) );
  gzwrite( backup_file, &_level_0_metabolic_uptake,           sizeof(_level_0_metabolic_uptake) );
  gzwrite( backup_file, &_level_0_metabolic_release,          sizeof(_level_0_metabolic_release) );
  gzwrite( backup_file, &_level_0_metabolic_growth_rate,      sizeof(_level_0_metabolic_growth_rate) );
  gzwrite( backup_file, &_level_0_Dmetabolic_growth_rate,     sizeof(_level_0_Dmetabolic_growth_rate) );
  gzwrite( backup_file, &_level_0_grn_nb_nodes,               sizeof(_level_0_grn_nb_nodes) );
  gzwrite( backup_file, &_level_0_grn_nb_edges,               sizeof(_level_0_grn_nb_edges) );
  gzwrite( backup_file, &_level_0_metabolic_nb_nodes,         sizeof(_level_0_metabolic_nb_nodes) );
  gzwrite( backup_file, &_level_0_metabolic_nb_edges,         sizeof(_level_0_metabolic_nb_edges) );
  gzwrite( backup_file, &_level_0_regulation_redundancy,      sizeof(_level_0_regulation_redundancy) );
  gzwrite( backup_file, &_level_0_metabolic_redundancy,       sizeof(_level_0_metabolic_redundancy) );
  
  /* GENOME STRUCTURE */
  gzwrite( backup_file, &_level_0_genome_size,                   sizeof(_level_0_genome_size) );
  gzwrite( backup_file, &_level_0_functional_size,               sizeof(_level_0_functional_size) );
  gzwrite( backup_file, &_level_0_genome_nb_NC,                  sizeof(_level_0_genome_nb_NC) );
  gzwrite( backup_file, &_level_0_genome_nb_E,                   sizeof(_level_0_genome_nb_E) );
  gzwrite( backup_file, &_level_0_genome_nb_TF,                  sizeof(_level_0_genome_nb_TF) );
  gzwrite( backup_file, &_level_0_genome_nb_BS,                  sizeof(_level_0_genome_nb_BS) );
  gzwrite( backup_file, &_level_0_genome_nb_P,                   sizeof(_level_0_genome_nb_P) );
  gzwrite( backup_file, &_level_0_genome_nb_inner_enzymes,       sizeof(_level_0_genome_nb_inner_enzymes) );
  gzwrite( backup_file, &_level_0_genome_nb_inflow_pumps,        sizeof(_level_0_genome_nb_inflow_pumps) );
  gzwrite( backup_file, &_level_0_genome_nb_outflow_pumps,       sizeof(_level_0_genome_nb_outflow_pumps) );
  gzwrite( backup_file, &_level_0_genome_nb_functional_regions,  sizeof(_level_0_genome_nb_functional_regions) );
  gzwrite( backup_file, &_level_0_genome_nb_enhancers,           sizeof(_level_0_genome_nb_enhancers) );
  gzwrite( backup_file, &_level_0_genome_nb_operators,           sizeof(_level_0_genome_nb_operators) );
  gzwrite( backup_file, &_level_0_genome_nb_E_regions,           sizeof(_level_0_genome_nb_E_regions) );
  gzwrite( backup_file, &_level_0_genome_nb_TF_regions,          sizeof(_level_0_genome_nb_TF_regions) );
  gzwrite( backup_file, &_level_0_genome_nb_mixed_regions,       sizeof(_level_0_genome_nb_mixed_regions) );
  gzwrite( backup_file, &_level_0_genome_functional_region_size, sizeof(_level_0_genome_functional_region_size) );
  gzwrite( backup_file, &_level_0_genome_E_region_size,          sizeof(_level_0_genome_E_region_size) );
  gzwrite( backup_file, &_level_0_genome_TF_region_size,         sizeof(_level_0_genome_TF_region_size) );
  gzwrite( backup_file, &_level_0_genome_mixed_region_size,      sizeof(_level_0_genome_mixed_region_size) );
  gzwrite( backup_file, &_level_0_genome_enhancer_size,          sizeof(_level_0_genome_enhancer_size) );
  gzwrite( backup_file, &_level_0_genome_operator_size,          sizeof(_level_0_genome_operator_size) );
  gzwrite( backup_file, &_level_0_genome_operon_size,            sizeof(_level_0_genome_operon_size) );
  gzwrite( backup_file, &_level_0_genome_E_operon_size,          sizeof(_level_0_genome_E_operon_size) );
  gzwrite( backup_file, &_level_0_genome_TF_operon_size,         sizeof(_level_0_genome_TF_operon_size) );
  gzwrite( backup_file, &_level_0_genome_mixed_operon_size,      sizeof(_level_0_genome_mixed_operon_size) );
  
  /* INHERITED STRUCTURE */
  gzwrite( backup_file, &_level_0_inherited_size,             sizeof(_level_0_inherited_size) );
  gzwrite( backup_file, &_level_0_inherited_nb_E,             sizeof(_level_0_inherited_nb_E) );
  gzwrite( backup_file, &_level_0_inherited_nb_TF,            sizeof(_level_0_inherited_nb_TF) );
  gzwrite( backup_file, &_level_0_inherited_nb_inner_enzymes, sizeof(_level_0_inherited_nb_inner_enzymes) );
  gzwrite( backup_file, &_level_0_inherited_nb_inflow_pumps,  sizeof(_level_0_inherited_nb_inflow_pumps) );
  gzwrite( backup_file, &_level_0_inherited_nb_outflow_pumps, sizeof(_level_0_inherited_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ LEVEL 1 statistics */
  
  /* GENETIC DIVERSITY */
  gzwrite( backup_file, &_level_1_true_diversity,              sizeof(_level_1_true_diversity) );
  gzwrite( backup_file, &_level_1_max_time_distance,           sizeof(_level_1_max_time_distance) );
  gzwrite( backup_file, &_level_1_max_generation_distance,     sizeof(_level_1_max_generation_distance) );
  gzwrite( backup_file, &_level_1_max_point_mutation_distance, sizeof(_level_1_max_point_mutation_distance) );
  gzwrite( backup_file, &_level_1_max_hgt_distance,            sizeof(_level_1_max_hgt_distance) );
  gzwrite( backup_file, &_level_1_max_duplication_distance,    sizeof(_level_1_max_duplication_distance) );
  gzwrite( backup_file, &_level_1_max_deletion_distance,       sizeof(_level_1_max_deletion_distance) );
  gzwrite( backup_file, &_level_1_max_inversion_distance,      sizeof(_level_1_max_inversion_distance) );
  gzwrite( backup_file, &_level_1_max_translocation_distance,  sizeof(_level_1_max_translocation_distance) );
  
  /* PHENOTYPE */
  gzwrite( backup_file, &_level_1_generations,                sizeof(_level_1_generations) );
  gzwrite( backup_file, &_level_1_inherited_TF_amount,        sizeof(_level_1_inherited_TF_amount) );
  gzwrite( backup_file, &_level_1_inherited_E_amount,         sizeof(_level_1_inherited_E_amount) );
  gzwrite( backup_file, &_level_1_TF_amount,                  sizeof(_level_1_TF_amount) );
  gzwrite( backup_file, &_level_1_E_amount,                   sizeof(_level_1_E_amount) );
  gzwrite( backup_file, &_level_1_inherited_metabolic_amount, sizeof(_level_1_inherited_metabolic_amount) );
  gzwrite( backup_file, &_level_1_metabolic_amount,           sizeof(_level_1_metabolic_amount) );
  gzwrite( backup_file, &_level_1_energy,                     sizeof(_level_1_energy) );
  gzwrite( backup_file, &_level_1_score,                      sizeof(_level_1_score) );
  gzwrite( backup_file, &_level_1_lifespan,                   sizeof(_level_1_lifespan) );
  gzwrite( backup_file, &_level_1_number_of_divisions,        sizeof(_level_1_number_of_divisions) );
  gzwrite( backup_file, &_level_1_toxicity,                   sizeof(_level_1_toxicity) );
  gzwrite( backup_file, &_level_1_metabolic_uptake,           sizeof(_level_1_metabolic_uptake) );
  gzwrite( backup_file, &_level_1_metabolic_release,          sizeof(_level_1_metabolic_release) );
  gzwrite( backup_file, &_level_1_metabolic_growth_rate,      sizeof(_level_1_metabolic_growth_rate) );
  gzwrite( backup_file, &_level_1_Dmetabolic_growth_rate,     sizeof(_level_1_Dmetabolic_growth_rate) );
  gzwrite( backup_file, &_level_1_grn_nb_nodes,               sizeof(_level_1_grn_nb_nodes) );
  gzwrite( backup_file, &_level_1_grn_nb_edges,               sizeof(_level_1_grn_nb_edges) );
  gzwrite( backup_file, &_level_1_metabolic_nb_nodes,         sizeof(_level_1_metabolic_nb_nodes) );
  gzwrite( backup_file, &_level_1_metabolic_nb_edges,         sizeof(_level_1_metabolic_nb_edges) );
  gzwrite( backup_file, &_level_1_regulation_redundancy,      sizeof(_level_1_regulation_redundancy) );
  gzwrite( backup_file, &_level_1_metabolic_redundancy,       sizeof(_level_1_metabolic_redundancy) );
  
  /* GENOME STRUCTURE */
  gzwrite( backup_file, &_level_1_genome_size,                   sizeof(_level_1_genome_size) );
  gzwrite( backup_file, &_level_1_functional_size,               sizeof(_level_1_functional_size) );
  gzwrite( backup_file, &_level_1_genome_nb_NC,                  sizeof(_level_1_genome_nb_NC) );
  gzwrite( backup_file, &_level_1_genome_nb_E,                   sizeof(_level_1_genome_nb_E) );
  gzwrite( backup_file, &_level_1_genome_nb_TF,                  sizeof(_level_1_genome_nb_TF) );
  gzwrite( backup_file, &_level_1_genome_nb_BS,                  sizeof(_level_1_genome_nb_BS) );
  gzwrite( backup_file, &_level_1_genome_nb_P,                   sizeof(_level_1_genome_nb_P) );
  gzwrite( backup_file, &_level_1_genome_nb_inner_enzymes,       sizeof(_level_1_genome_nb_inner_enzymes) );
  gzwrite( backup_file, &_level_1_genome_nb_inflow_pumps,        sizeof(_level_1_genome_nb_inflow_pumps) );
  gzwrite( backup_file, &_level_1_genome_nb_outflow_pumps,       sizeof(_level_1_genome_nb_outflow_pumps) );
  gzwrite( backup_file, &_level_1_genome_nb_functional_regions,  sizeof(_level_1_genome_nb_functional_regions) );
  gzwrite( backup_file, &_level_1_genome_nb_enhancers,           sizeof(_level_1_genome_nb_enhancers) );
  gzwrite( backup_file, &_level_1_genome_nb_operators,           sizeof(_level_1_genome_nb_operators) );
  gzwrite( backup_file, &_level_1_genome_nb_E_regions,           sizeof(_level_1_genome_nb_E_regions) );
  gzwrite( backup_file, &_level_1_genome_nb_TF_regions,          sizeof(_level_1_genome_nb_TF_regions) );
  gzwrite( backup_file, &_level_1_genome_nb_mixed_regions,       sizeof(_level_1_genome_nb_mixed_regions) );
  gzwrite( backup_file, &_level_1_genome_functional_region_size, sizeof(_level_1_genome_functional_region_size) );
  gzwrite( backup_file, &_level_1_genome_E_region_size,          sizeof(_level_1_genome_E_region_size) );
  gzwrite( backup_file, &_level_1_genome_TF_region_size,         sizeof(_level_1_genome_TF_region_size) );
  gzwrite( backup_file, &_level_1_genome_mixed_region_size,      sizeof(_level_1_genome_mixed_region_size) );
  gzwrite( backup_file, &_level_1_genome_enhancer_size,          sizeof(_level_1_genome_enhancer_size) );
  gzwrite( backup_file, &_level_1_genome_operator_size,          sizeof(_level_1_genome_operator_size) );
  gzwrite( backup_file, &_level_1_genome_operon_size,            sizeof(_level_1_genome_operon_size) );
  gzwrite( backup_file, &_level_1_genome_E_operon_size,          sizeof(_level_1_genome_E_operon_size) );
  gzwrite( backup_file, &_level_1_genome_TF_operon_size,         sizeof(_level_1_genome_TF_operon_size) );
  gzwrite( backup_file, &_level_1_genome_mixed_operon_size,      sizeof(_level_1_genome_mixed_operon_size) );
  
  /* INHERITED STRUCTURE */
  gzwrite( backup_file, &_level_1_inherited_size,             sizeof(_level_1_inherited_size) );
  gzwrite( backup_file, &_level_1_inherited_nb_E,             sizeof(_level_1_inherited_nb_E) );
  gzwrite( backup_file, &_level_1_inherited_nb_TF,            sizeof(_level_1_inherited_nb_TF) );
  gzwrite( backup_file, &_level_1_inherited_nb_inner_enzymes, sizeof(_level_1_inherited_nb_inner_enzymes) );
  gzwrite( backup_file, &_level_1_inherited_nb_inflow_pumps,  sizeof(_level_1_inherited_nb_inflow_pumps) );
  gzwrite( backup_file, &_level_1_inherited_nb_outflow_pumps, sizeof(_level_1_inherited_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ LEVEL 2 statistics */
  
  /* GENETIC DIVERSITY */
  gzwrite( backup_file, &_level_2_true_diversity,              sizeof(_level_2_true_diversity) );
  gzwrite( backup_file, &_level_2_max_time_distance,           sizeof(_level_2_max_time_distance) );
  gzwrite( backup_file, &_level_2_max_generation_distance,     sizeof(_level_2_max_generation_distance) );
  gzwrite( backup_file, &_level_2_max_point_mutation_distance, sizeof(_level_2_max_point_mutation_distance) );
  gzwrite( backup_file, &_level_2_max_hgt_distance,            sizeof(_level_2_max_hgt_distance) );
  gzwrite( backup_file, &_level_2_max_duplication_distance,    sizeof(_level_2_max_duplication_distance) );
  gzwrite( backup_file, &_level_2_max_deletion_distance,       sizeof(_level_2_max_deletion_distance) );
  gzwrite( backup_file, &_level_2_max_inversion_distance,      sizeof(_level_2_max_inversion_distance) );
  gzwrite( backup_file, &_level_2_max_translocation_distance,  sizeof(_level_2_max_translocation_distance) );
  
  /* PHENOTYPE */
  gzwrite( backup_file, &_level_2_generations,                sizeof(_level_2_generations) );
  gzwrite( backup_file, &_level_2_inherited_TF_amount,        sizeof(_level_2_inherited_TF_amount) );
  gzwrite( backup_file, &_level_2_inherited_E_amount,         sizeof(_level_2_inherited_E_amount) );
  gzwrite( backup_file, &_level_2_TF_amount,                  sizeof(_level_2_TF_amount) );
  gzwrite( backup_file, &_level_2_E_amount,                   sizeof(_level_2_E_amount) );
  gzwrite( backup_file, &_level_2_inherited_metabolic_amount, sizeof(_level_2_inherited_metabolic_amount) );
  gzwrite( backup_file, &_level_2_metabolic_amount,           sizeof(_level_2_metabolic_amount) );
  gzwrite( backup_file, &_level_2_energy,                     sizeof(_level_2_energy) );
  gzwrite( backup_file, &_level_2_score,                      sizeof(_level_2_score) );
  gzwrite( backup_file, &_level_2_lifespan,                   sizeof(_level_2_lifespan) );
  gzwrite( backup_file, &_level_2_number_of_divisions,        sizeof(_level_2_number_of_divisions) );
  gzwrite( backup_file, &_level_2_toxicity,                   sizeof(_level_2_toxicity) );
  gzwrite( backup_file, &_level_2_metabolic_uptake,           sizeof(_level_2_metabolic_uptake) );
  gzwrite( backup_file, &_level_2_metabolic_release,          sizeof(_level_2_metabolic_release) );
  gzwrite( backup_file, &_level_2_metabolic_growth_rate,      sizeof(_level_2_metabolic_growth_rate) );
  gzwrite( backup_file, &_level_2_Dmetabolic_growth_rate,     sizeof(_level_2_Dmetabolic_growth_rate) );
  gzwrite( backup_file, &_level_2_grn_nb_nodes,               sizeof(_level_2_grn_nb_nodes) );
  gzwrite( backup_file, &_level_2_grn_nb_edges,               sizeof(_level_2_grn_nb_edges) );
  gzwrite( backup_file, &_level_2_metabolic_nb_nodes,         sizeof(_level_2_metabolic_nb_nodes) );
  gzwrite( backup_file, &_level_2_metabolic_nb_edges,         sizeof(_level_2_metabolic_nb_edges) );
  gzwrite( backup_file, &_level_2_regulation_redundancy,      sizeof(_level_2_regulation_redundancy) );
  gzwrite( backup_file, &_level_2_metabolic_redundancy,       sizeof(_level_2_metabolic_redundancy) );
  
  /* GENOME STRUCTURE */
  gzwrite( backup_file, &_level_2_genome_size,                   sizeof(_level_2_genome_size) );
  gzwrite( backup_file, &_level_2_functional_size,               sizeof(_level_2_functional_size) );
  gzwrite( backup_file, &_level_2_genome_nb_NC,                  sizeof(_level_2_genome_nb_NC) );
  gzwrite( backup_file, &_level_2_genome_nb_E,                   sizeof(_level_2_genome_nb_E) );
  gzwrite( backup_file, &_level_2_genome_nb_TF,                  sizeof(_level_2_genome_nb_TF) );
  gzwrite( backup_file, &_level_2_genome_nb_BS,                  sizeof(_level_2_genome_nb_BS) );
  gzwrite( backup_file, &_level_2_genome_nb_P,                   sizeof(_level_2_genome_nb_P) );
  gzwrite( backup_file, &_level_2_genome_nb_inner_enzymes,       sizeof(_level_2_genome_nb_inner_enzymes) );
  gzwrite( backup_file, &_level_2_genome_nb_inflow_pumps,        sizeof(_level_2_genome_nb_inflow_pumps) );
  gzwrite( backup_file, &_level_2_genome_nb_outflow_pumps,       sizeof(_level_2_genome_nb_outflow_pumps) );
  gzwrite( backup_file, &_level_2_genome_nb_functional_regions,  sizeof(_level_2_genome_nb_functional_regions) );
  gzwrite( backup_file, &_level_2_genome_nb_enhancers,           sizeof(_level_2_genome_nb_enhancers) );
  gzwrite( backup_file, &_level_2_genome_nb_operators,           sizeof(_level_2_genome_nb_operators) );
  gzwrite( backup_file, &_level_2_genome_nb_E_regions,           sizeof(_level_2_genome_nb_E_regions) );
  gzwrite( backup_file, &_level_2_genome_nb_TF_regions,          sizeof(_level_2_genome_nb_TF_regions) );
  gzwrite( backup_file, &_level_2_genome_nb_mixed_regions,       sizeof(_level_2_genome_nb_mixed_regions) );
  gzwrite( backup_file, &_level_2_genome_functional_region_size, sizeof(_level_2_genome_functional_region_size) );
  gzwrite( backup_file, &_level_2_genome_E_region_size,          sizeof(_level_2_genome_E_region_size) );
  gzwrite( backup_file, &_level_2_genome_TF_region_size,         sizeof(_level_2_genome_TF_region_size) );
  gzwrite( backup_file, &_level_2_genome_mixed_region_size,      sizeof(_level_2_genome_mixed_region_size) );
  gzwrite( backup_file, &_level_2_genome_enhancer_size,          sizeof(_level_2_genome_enhancer_size) );
  gzwrite( backup_file, &_level_2_genome_operator_size,          sizeof(_level_2_genome_operator_size) );
  gzwrite( backup_file, &_level_2_genome_operon_size,            sizeof(_level_2_genome_operon_size) );
  gzwrite( backup_file, &_level_2_genome_E_operon_size,          sizeof(_level_2_genome_E_operon_size) );
  gzwrite( backup_file, &_level_2_genome_TF_operon_size,         sizeof(_level_2_genome_TF_operon_size) );
  gzwrite( backup_file, &_level_2_genome_mixed_operon_size,      sizeof(_level_2_genome_mixed_operon_size) );
  
  /* INHERITED STRUCTURE */
  gzwrite( backup_file, &_level_2_inherited_size,             sizeof(_level_2_inherited_size) );
  gzwrite( backup_file, &_level_2_inherited_nb_E,             sizeof(_level_2_inherited_nb_E) );
  gzwrite( backup_file, &_level_2_inherited_nb_TF,            sizeof(_level_2_inherited_nb_TF) );
  gzwrite( backup_file, &_level_2_inherited_nb_inner_enzymes, sizeof(_level_2_inherited_nb_inner_enzymes) );
  gzwrite( backup_file, &_level_2_inherited_nb_inflow_pumps,  sizeof(_level_2_inherited_nb_inflow_pumps) );
  gzwrite( backup_file, &_level_2_inherited_nb_outflow_pumps, sizeof(_level_2_inherited_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ NO LEVEL statistics */
  
  /* GENETIC DIVERSITY */
  gzwrite( backup_file, &_no_level_true_diversity,              sizeof(_no_level_true_diversity) );
  gzwrite( backup_file, &_no_level_max_time_distance,           sizeof(_no_level_max_time_distance) );
  gzwrite( backup_file, &_no_level_max_generation_distance,     sizeof(_no_level_max_generation_distance) );
  gzwrite( backup_file, &_no_level_max_point_mutation_distance, sizeof(_no_level_max_point_mutation_distance) );
  gzwrite( backup_file, &_no_level_max_hgt_distance,            sizeof(_no_level_max_hgt_distance) );
  gzwrite( backup_file, &_no_level_max_duplication_distance,    sizeof(_no_level_max_duplication_distance) );
  gzwrite( backup_file, &_no_level_max_deletion_distance,       sizeof(_no_level_max_deletion_distance) );
  gzwrite( backup_file, &_no_level_max_inversion_distance,      sizeof(_no_level_max_inversion_distance) );
  gzwrite( backup_file, &_no_level_max_translocation_distance,  sizeof(_no_level_max_translocation_distance) );
  
  /* PHENOTYPE */
  gzwrite( backup_file, &_no_level_generations,                sizeof(_no_level_generations) );
  gzwrite( backup_file, &_no_level_inherited_TF_amount,        sizeof(_no_level_inherited_TF_amount) );
  gzwrite( backup_file, &_no_level_inherited_E_amount,         sizeof(_no_level_inherited_E_amount) );
  gzwrite( backup_file, &_no_level_TF_amount,                  sizeof(_no_level_TF_amount) );
  gzwrite( backup_file, &_no_level_E_amount,                   sizeof(_no_level_E_amount) );
  gzwrite( backup_file, &_no_level_inherited_metabolic_amount, sizeof(_no_level_inherited_metabolic_amount) );
  gzwrite( backup_file, &_no_level_metabolic_amount,           sizeof(_no_level_metabolic_amount) );
  gzwrite( backup_file, &_no_level_energy,                     sizeof(_no_level_energy) );
  gzwrite( backup_file, &_no_level_score,                      sizeof(_no_level_score) );
  gzwrite( backup_file, &_no_level_lifespan,                   sizeof(_no_level_lifespan) );
  gzwrite( backup_file, &_no_level_number_of_divisions,        sizeof(_no_level_number_of_divisions) );
  gzwrite( backup_file, &_no_level_toxicity,                   sizeof(_no_level_toxicity) );
  gzwrite( backup_file, &_no_level_metabolic_uptake,           sizeof(_no_level_metabolic_uptake) );
  gzwrite( backup_file, &_no_level_metabolic_release,          sizeof(_no_level_metabolic_release) );
  gzwrite( backup_file, &_no_level_metabolic_growth_rate,      sizeof(_no_level_metabolic_growth_rate) );
  gzwrite( backup_file, &_no_level_Dmetabolic_growth_rate,     sizeof(_no_level_Dmetabolic_growth_rate) );
  gzwrite( backup_file, &_no_level_grn_nb_nodes,               sizeof(_no_level_grn_nb_nodes) );
  gzwrite( backup_file, &_no_level_grn_nb_edges,               sizeof(_no_level_grn_nb_edges) );
  gzwrite( backup_file, &_no_level_metabolic_nb_nodes,         sizeof(_no_level_metabolic_nb_nodes) );
  gzwrite( backup_file, &_no_level_metabolic_nb_edges,         sizeof(_no_level_metabolic_nb_edges) );
  gzwrite( backup_file, &_no_level_regulation_redundancy,      sizeof(_no_level_regulation_redundancy) );
  gzwrite( backup_file, &_no_level_metabolic_redundancy,       sizeof(_no_level_metabolic_redundancy) );
  
  /* GENOME STRUCTURE */
  gzwrite( backup_file, &_no_level_genome_size,                   sizeof(_no_level_genome_size) );
  gzwrite( backup_file, &_no_level_functional_size,               sizeof(_no_level_functional_size) );
  gzwrite( backup_file, &_no_level_genome_nb_NC,                  sizeof(_no_level_genome_nb_NC) );
  gzwrite( backup_file, &_no_level_genome_nb_E,                   sizeof(_no_level_genome_nb_E) );
  gzwrite( backup_file, &_no_level_genome_nb_TF,                  sizeof(_no_level_genome_nb_TF) );
  gzwrite( backup_file, &_no_level_genome_nb_BS,                  sizeof(_no_level_genome_nb_BS) );
  gzwrite( backup_file, &_no_level_genome_nb_P,                   sizeof(_no_level_genome_nb_P) );
  gzwrite( backup_file, &_no_level_genome_nb_inner_enzymes,       sizeof(_no_level_genome_nb_inner_enzymes) );
  gzwrite( backup_file, &_no_level_genome_nb_inflow_pumps,        sizeof(_no_level_genome_nb_inflow_pumps) );
  gzwrite( backup_file, &_no_level_genome_nb_outflow_pumps,       sizeof(_no_level_genome_nb_outflow_pumps) );
  gzwrite( backup_file, &_no_level_genome_nb_functional_regions,  sizeof(_no_level_genome_nb_functional_regions) );
  gzwrite( backup_file, &_no_level_genome_nb_enhancers,           sizeof(_no_level_genome_nb_enhancers) );
  gzwrite( backup_file, &_no_level_genome_nb_operators,           sizeof(_no_level_genome_nb_operators) );
  gzwrite( backup_file, &_no_level_genome_nb_E_regions,           sizeof(_no_level_genome_nb_E_regions) );
  gzwrite( backup_file, &_no_level_genome_nb_TF_regions,          sizeof(_no_level_genome_nb_TF_regions) );
  gzwrite( backup_file, &_no_level_genome_nb_mixed_regions,       sizeof(_no_level_genome_nb_mixed_regions) );
  gzwrite( backup_file, &_no_level_genome_functional_region_size, sizeof(_no_level_genome_functional_region_size) );
  gzwrite( backup_file, &_no_level_genome_E_region_size,          sizeof(_no_level_genome_E_region_size) );
  gzwrite( backup_file, &_no_level_genome_TF_region_size,         sizeof(_no_level_genome_TF_region_size) );
  gzwrite( backup_file, &_no_level_genome_mixed_region_size,      sizeof(_no_level_genome_mixed_region_size) );
  gzwrite( backup_file, &_no_level_genome_enhancer_size,          sizeof(_no_level_genome_enhancer_size) );
  gzwrite( backup_file, &_no_level_genome_operator_size,          sizeof(_no_level_genome_operator_size) );
  gzwrite( backup_file, &_no_level_genome_operon_size,            sizeof(_no_level_genome_operon_size) );
  gzwrite( backup_file, &_no_level_genome_E_operon_size,          sizeof(_no_level_genome_E_operon_size) );
  gzwrite( backup_file, &_no_level_genome_TF_operon_size,         sizeof(_no_level_genome_TF_operon_size) );
  gzwrite( backup_file, &_no_level_genome_mixed_operon_size,      sizeof(_no_level_genome_mixed_operon_size) );
  
  /* INHERITED STRUCTURE */
  gzwrite( backup_file, &_no_level_inherited_size,             sizeof(_no_level_inherited_size) );
  gzwrite( backup_file, &_no_level_inherited_nb_E,             sizeof(_no_level_inherited_nb_E) );
  gzwrite( backup_file, &_no_level_inherited_nb_TF,            sizeof(_no_level_inherited_nb_TF) );
  gzwrite( backup_file, &_no_level_inherited_nb_inner_enzymes, sizeof(_no_level_inherited_nb_inner_enzymes) );
  gzwrite( backup_file, &_no_level_inherited_nb_inflow_pumps,  sizeof(_no_level_inherited_nb_inflow_pumps) );
  gzwrite( backup_file, &_no_level_inherited_nb_outflow_pumps, sizeof(_no_level_inherited_nb_outflow_pumps) );
  
}

/**
 * \brief    Initialize the trophic network
 * \details  This method basically creates the environment group
 * \param    void
 * \return   \e void
 */
void TrophicNetwork::initialize_trophic_network( void )
{
  size_t N    = _environment->get_species_lists_size();
  _current_id = 0;
  
  /*-----------------------------------------------------------------*/
  /* 1) Initialize the environment trophic group and profile strings */
  /*-----------------------------------------------------------------*/
  TrophicGroup* env_group = new TrophicGroup(get_new_id(), _population->get_time(), 0.0, 0.0, 0.0);
  std::string trophic_profile    = "";
  std::string production_profile = "";
  std::string uptake_profile     = "";
  std::string release_profile    = "";
  
  /*-----------------------------------------------------------------*/
  /* 2) Build the environment production profile                     */
  /*-----------------------------------------------------------------*/
  if (_parameters->get_environment_properties()->species_tag_range.law != UNIFORM)
  {
    printf("TrophicNetwork::initialize_trophic_network() method is not valid if environment species tag is different that 'UNIFORM'. Exit simulation.\n");
    exit(EXIT_FAILURE);
  }
  for (size_t i = 0; i < N; i++)
  {
    if ((i+1) >= _parameters->get_environment_properties()->species_tag_range.min && (i+1) <= _parameters->get_environment_properties()->species_tag_range.max)
    {
      production_profile += "1";
    }
    else
    {
      production_profile += "0";
    }
  }
  
  /*-----------------------------------------------------------------*/
  /* 3) Build the environment uptake and release profiles            */
  /*-----------------------------------------------------------------*/
  for (size_t i = 0; i < N; i++)
  {
    uptake_profile  += "0";
    release_profile += "0";
  }
  
  /*-----------------------------------------------------------------*/
  /* 4) Build environment trophic profile                            */
  /*-----------------------------------------------------------------*/
  trophic_profile = production_profile+uptake_profile+release_profile;
  
  /*-----------------------------------------------------------------*/
  /* 5) Update environment group                                     */
  /*-----------------------------------------------------------------*/
  env_group->set_trophic_profile(trophic_profile);
  env_group->set_production_profile(production_profile);
  env_group->set_uptake_profile(uptake_profile);
  env_group->set_release_profile(release_profile);
  
  /*-----------------------------------------------------------------*/
  /* 6) Add the environment node to the trophic network              */
  /*-----------------------------------------------------------------*/
  _group_map.clear();
  _group_map[env_group->get_identifier()] = env_group;
  _nb_no_level_groups++;
  _nb_group_appearances++;
}

/**
 * \brief    Load the population in the trophic network
 * \details  --
 * \param    void
 * \return   \e void
 */
void TrophicNetwork::load_population( void )
{
  size_t N = _environment->get_species_lists_size();
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Clear statistics and relationships */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    _iterator->second->update_profile(N);
    _iterator->second->clear();
  }
  _nb_level_0_groups    = 0;
  _nb_level_1_groups    = 0;
  _nb_level_2_groups    = 0;
  _nb_no_level_groups   = 0;
  _nb_level_0_cells     = 0;
  _nb_level_1_cells     = 0;
  _nb_level_2_cells     = 0;
  _nb_no_level_cells    = 0;
  _nb_group_appearances = 0;
  _nb_group_extinctions = 0;
  _mean_group_lifespan  = 0.0;
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Explore cell profiles              */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  std::string env_profile = _group_map[0]->get_production_profile();
  for (size_t pos = 0; pos < _population->get_width()*_population->get_height(); pos++)
  {
    Cell* cell = _population->get_cell(pos);
    if (cell->isAlive())
    {
      double*        Xcell = cell->get_species_list()->get_X();
      double*        Xenv  = _environment->get_X()->get_X();
      reaction_list* rlist = cell->get_ode()->get_reaction_list();
      
      std::string trophic_profile    = "";
      std::string production_profile = "";
      std::string uptake_profile     = "";
      std::string release_profile    = "";
      
      /*-----------------------------------------------------*/
      /* 2.1) Build the uptake profile                       */
      /*-----------------------------------------------------*/
      bool* profile = new bool[N];
      for (size_t i = 0; i < N; i++)
      {
        profile[i] = false;
      }
      for (size_t i = 0; i < rlist->metabolic_N; i++)
      {
        int s = rlist->metabolic_s[i];
        /* If there is by-products in the environment (Xenv) */
        /* or exogeneous nutrient (env_profile), the uptake  */
        /* profile is 1. Else it is 0.                       */
        if (rlist->metabolic_type[i] == INFLOWING_PUMP_ACTIVITY && (Xenv[s-1] > MINIMUM_CONCENTRATION || env_profile[s-1] == '1'))
        {
          profile[s-1] = true;
        }
      }
      /* Build the uptake profile */
      for (size_t i = 0; i < N; i++)
      {
        if (profile[i])
        {
          uptake_profile += "1";
        }
        else
        {
          uptake_profile += "0";
        }
      }
      
      /*-----------------------------------------------------*/
      /* 2.2) Build the release profile                      */
      /*-----------------------------------------------------*/
      for (size_t i = 0; i < N; i++)
      {
        profile[i] = false;
      }
      for (size_t i = 0; i < rlist->metabolic_N; i++)
      {
        int s = rlist->metabolic_s[i];
        /* If the release pump is active, set the profile at 1 */
        if (rlist->metabolic_type[i] == OUTFLOWING_PUMP_ACTIVITY && Xcell[s-1] > MINIMUM_CONCENTRATION)
        {
          profile[s-1] = true;
        }
      }
      /* Build the release profile */
      for (size_t i = 0; i < N; i++)
      {
        if (profile[i])
        {
          release_profile += "1";
        }
        else
        {
          release_profile += "0";
        }
      }
      
      /*-----------------------------------------------------*/
      /* 2.3) Build the production profile                   */
      /*-----------------------------------------------------*/
      for (size_t i = 0; i < N; i++)
      {
        profile[i] = false;
      }
      for (size_t i = 0; i < rlist->metabolic_N; i++)
      {
        int p = rlist->metabolic_p[i];
        /* If the reaction produces a not uptaken metabolite, */
        /* set the production profile to 1                    */
        if (rlist->metabolic_type[i] == CATALYTIC_ACTIVITY && uptake_profile[p-1] == '0' && Xcell[p-1] > MINIMUM_CONCENTRATION )
        {
          profile[p-1] = true;
        }
      }
      /* Build the production profile */
      for (size_t i = 0; i < N; i++)
      {
        if (profile[i])
        {
          production_profile += "1";
        }
        else
        {
          production_profile += "0";
        }
      }
      delete[] profile;
      profile = NULL;
      
      /*-----------------------------------------------------*/
      /* 2.4) Build node profile string                      */
      /*-----------------------------------------------------*/
      trophic_profile = production_profile+uptake_profile+release_profile;
      
      /*-----------------------------------------------------*/
      /* 2.5) Evaluate if the node pre-exists in the network */
      /*-----------------------------------------------------*/
      unsigned long long int group_id;
      
      /* If the group already exists, update its statistics */
      
      if (groupExists(trophic_profile, group_id))
      {
        TrophicGroup* group = _group_map[group_id];
        group->add_alive_cell(_parameters, cell);
        cell->set_trophic_group(group_id);
        cell->set_red_color(group->get_red_color());
        cell->set_green_color(group->get_green_color());
        cell->set_blue_color(group->get_blue_color());
      }
      
      /* If the group does not exist, create it */
      
      else
      {
        double red   = _parameters->get_prng()->uniform()*(255.0-50.0)+50.0;
        double green = _parameters->get_prng()->uniform()*(255.0-50.0)+50.0;
        double blue  = _parameters->get_prng()->uniform()*(255.0-50.0)+50.0;
        TrophicGroup* group = new TrophicGroup(get_new_id(), _population->get_time(), red, green, blue);
        group->set_trophic_profile(trophic_profile);
        group->set_production_profile(production_profile);
        group->set_uptake_profile(uptake_profile);
        group->set_release_profile(release_profile);
        group->add_alive_cell(_parameters, cell);
        _group_map[group->get_identifier()] = group;
        _nb_group_appearances++;
        cell->set_trophic_group(group->get_identifier());
        cell->set_red_color(group->get_red_color());
        cell->set_green_color(group->get_green_color());
        cell->set_blue_color(group->get_blue_color());
      }
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) Build trophic network links        */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  std::unordered_map<unsigned long long int, TrophicGroup*>::iterator it1;
  std::unordered_map<unsigned long long int, TrophicGroup*>::iterator it2;
  for (it1 = _group_map.begin(); it1 != _group_map.end(); ++it1)
  {
    TrophicGroup* group_i = it1->second;
    bool level_0 = false;
    bool level_2 = false;
    for (it2 = _group_map.begin(); it2 != _group_map.end(); ++it2)
    {
      if (it1->first != it2->first)
      {
        TrophicGroup* group_j = it2->second;
        /*----------------------------------------------------------------------*/
        /* For each pump in the uptake profile of the group i, if the group j   */
        /* produces or releases the same metabolite, save it in necrophagy link */
        /* or active release link of i                                          */
        /*----------------------------------------------------------------------*/
        bool necrophagy     = false;
        bool active_release = false;
        for (size_t met = 0; met < N; met++)
        {
          /*------------------------------------------------*/
          /* evaluate the production profile of the group j */
          /*------------------------------------------------*/
          if (group_i->get_uptake_profile()[met] == '1' && group_j->get_production_profile()[met] == '1')
          {
            necrophagy = true;
            if (group_j->get_identifier() == 0)
            {
              level_0 = true;
            }
            else
            {
              level_2 = true;
            }
          }
          /*------------------------------------------------*/
          /* evaluate the release profile of the group j    */
          /*------------------------------------------------*/
          if (group_i->get_uptake_profile()[met] == '1' && group_j->get_release_profile()[met] == '1')
          {
            active_release = true;
            level_2 = true;
          }
        }
        if (necrophagy && !active_release)
        {
          group_i->add_necrophagy_link(it2->first);
        }
        else if (necrophagy && active_release)
        {
          group_i->add_active_release_link(it2->first);
        }
      }
    }
    /*----------------------------------------*/
    /* Evaluate the trophic level the group i */
    /*----------------------------------------*/
    if (!level_0 && !level_2)
    {
      group_i->set_trophic_level(NO_LEVEL);
    }
    else if (level_0 && !level_2)
    {
      group_i->set_trophic_level(LEVEL_0);
    }
    else if (level_0 && level_2)
    {
      group_i->set_trophic_level(LEVEL_1);
    }
    else if (!level_0 && level_2)
    {
      group_i->set_trophic_level(LEVEL_2);
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4) Delete empty groups                */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  std::vector<unsigned long long int> to_remove;
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    if (_iterator->second->get_number_of_cells() == 0 && _iterator->first != 0)
    {
      to_remove.push_back(_iterator->first);
    }
  }
  for (size_t i = 0; i < to_remove.size(); i++)
  {
    delete _group_map[to_remove[i]];
    _group_map[to_remove[i]] = NULL;
    _group_map.erase(to_remove[i]);
  }
  _nb_group_extinctions += to_remove.size();
  to_remove.clear();
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 5) Update trophic network statistics  */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    if (_iterator->first != 0)
    {
      _iterator->second->compute_mean();
      _iterator->second->update_lifespan(_population->get_time());
      _mean_group_lifespan += _iterator->second->get_lifespan();
      if (_iterator->second->get_trophic_level() == LEVEL_0)
      {
        _nb_level_0_groups++;
        _nb_level_0_cells += _iterator->second->get_number_of_cells();
      }
      else if (_iterator->second->get_trophic_level() == LEVEL_1)
      {
        _nb_level_1_groups++;
        _nb_level_1_cells += _iterator->second->get_number_of_cells();
      }
      else if (_iterator->second->get_trophic_level() == LEVEL_2)
      {
        _nb_level_2_groups++;
        _nb_level_2_cells += _iterator->second->get_number_of_cells();
      }
      else if (_iterator->second->get_trophic_level() == NO_LEVEL)
      {
        _nb_no_level_groups++;
        _nb_no_level_cells += _iterator->second->get_number_of_cells();
      }
    }
  }
  if (_group_map.size() > 1)
  {
    _mean_group_lifespan /= (double)(_group_map.size()-1);
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 6) Update cells trophic level         */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  for (size_t pos = 0; pos < _population->get_width()*_population->get_height(); pos++)
  {
    Cell* cell = _population->get_cell(pos);
    if (cell->isAlive())
    {
      cell->set_trophic_level(_group_map[cell->get_trophic_group()]->get_trophic_level());
    }
  }
}

/**
 * \brief    Compute diversity statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void TrophicNetwork::compute_diversity_statistics( void )
{
  _min_true_diversity              = 0.0;
  _min_max_time_distance           = 0.0;
  _min_max_generation_distance     = 0.0;
  _min_max_point_mutation_distance = 0.0;
  _min_max_hgt_distance            = 0.0;
  _min_max_duplication_distance    = 0.0;
  _min_max_deletion_distance       = 0.0;
  _min_max_inversion_distance      = 0.0;
  _min_max_translocation_distance  = 0.0;
  
  _mean_true_diversity              = 0.0;
  _mean_max_time_distance           = 0.0;
  _mean_max_generation_distance     = 0.0;
  _mean_max_point_mutation_distance = 0.0;
  _mean_max_hgt_distance            = 0.0;
  _mean_max_duplication_distance    = 0.0;
  _mean_max_deletion_distance       = 0.0;
  _mean_max_inversion_distance      = 0.0;
  _mean_max_translocation_distance  = 0.0;
  
  _max_true_diversity              = 0.0;
  _max_max_time_distance           = 0.0;
  _max_max_generation_distance     = 0.0;
  _max_max_point_mutation_distance = 0.0;
  _max_max_hgt_distance            = 0.0;
  _max_max_duplication_distance    = 0.0;
  _max_max_deletion_distance       = 0.0;
  _max_max_inversion_distance      = 0.0;
  _max_max_translocation_distance  = 0.0;
  
  _iterator = _group_map.begin();
  while (_iterator != _group_map.end())
  {
    if (_iterator->first != 0)
    {
      /* Compute mean values */
      _mean_true_diversity              += _iterator->second->get_true_diversity();
      _mean_max_time_distance           += _iterator->second->get_max_time_distance();
      _mean_max_generation_distance     += _iterator->second->get_max_generation_distance();
      _mean_max_point_mutation_distance += _iterator->second->get_max_point_mutation_distance();
      _mean_max_hgt_distance            += _iterator->second->get_max_hgt_distance();
      _mean_max_duplication_distance    += _iterator->second->get_max_duplication_distance();
      _mean_max_deletion_distance       += _iterator->second->get_max_deletion_distance();
      _mean_max_inversion_distance      += _iterator->second->get_max_inversion_distance();
      _mean_max_translocation_distance  += _iterator->second->get_max_translocation_distance();
      
      /* Compute minimum values */
      if (_min_true_diversity > _iterator->second->get_true_diversity())
      {
        _min_true_diversity = _iterator->second->get_true_diversity();
      }
      if (_min_max_time_distance > _iterator->second->get_max_time_distance())
      {
        _min_max_time_distance = _iterator->second->get_max_time_distance();
      }
      if (_min_max_generation_distance > _iterator->second->get_max_generation_distance())
      {
        _min_max_generation_distance = _iterator->second->get_max_generation_distance();
      }
      if (_min_max_point_mutation_distance > _iterator->second->get_max_point_mutation_distance())
      {
        _min_max_point_mutation_distance = _iterator->second->get_max_point_mutation_distance();
      }
      if (_min_max_hgt_distance > _iterator->second->get_max_hgt_distance())
      {
        _min_max_hgt_distance = _iterator->second->get_max_hgt_distance();
      }
      if (_min_max_duplication_distance > _iterator->second->get_max_duplication_distance())
      {
        _min_max_duplication_distance = _iterator->second->get_max_duplication_distance();
      }
      if (_min_max_deletion_distance > _iterator->second->get_max_deletion_distance())
      {
        _min_max_deletion_distance = _iterator->second->get_max_deletion_distance();
      }
      if (_min_max_inversion_distance > _iterator->second->get_max_inversion_distance())
      {
        _min_max_inversion_distance = _iterator->second->get_max_inversion_distance();
      }
      if (_min_max_translocation_distance > _iterator->second->get_max_translocation_distance())
      {
        _min_max_translocation_distance = _iterator->second->get_max_translocation_distance();
      }
      
      /* Compute maximum values */
      if (_max_true_diversity < _iterator->second->get_true_diversity())
      {
        _max_true_diversity = _iterator->second->get_true_diversity();
      }
      if (_max_max_time_distance < _iterator->second->get_max_time_distance())
      {
        _max_max_time_distance = _iterator->second->get_max_time_distance();
      }
      if (_max_max_generation_distance < _iterator->second->get_max_generation_distance())
      {
        _max_max_generation_distance = _iterator->second->get_max_generation_distance();
      }
      if (_max_max_point_mutation_distance < _iterator->second->get_max_point_mutation_distance())
      {
        _max_max_point_mutation_distance = _iterator->second->get_max_point_mutation_distance();
      }
      if (_max_max_hgt_distance < _iterator->second->get_max_hgt_distance())
      {
        _max_max_hgt_distance = _iterator->second->get_max_hgt_distance();
      }
      if (_max_max_duplication_distance < _iterator->second->get_max_duplication_distance())
      {
        _max_max_duplication_distance = _iterator->second->get_max_duplication_distance();
      }
      if (_max_max_deletion_distance < _iterator->second->get_max_deletion_distance())
      {
        _max_max_deletion_distance = _iterator->second->get_max_deletion_distance();
      }
      if (_max_max_inversion_distance < _iterator->second->get_max_inversion_distance())
      {
        _max_max_inversion_distance = _iterator->second->get_max_inversion_distance();
      }
      if (_max_max_translocation_distance < _iterator->second->get_max_translocation_distance())
      {
        _max_max_translocation_distance = _iterator->second->get_max_translocation_distance();
      }
    }
    _iterator++;
  }
  
  /* Compute mean values */
  if (_group_map.size() > 1)
  {
    _mean_true_diversity              /= (double)_group_map.size()-1.0;
    _mean_max_time_distance           /= (double)_group_map.size()-1.0;
    _mean_max_generation_distance     /= (double)_group_map.size()-1.0;
    _mean_max_point_mutation_distance /= (double)_group_map.size()-1.0;
    _mean_max_hgt_distance            /= (double)_group_map.size()-1.0;
    _mean_max_duplication_distance    /= (double)_group_map.size()-1.0;
    _mean_max_deletion_distance       /= (double)_group_map.size()-1.0;
    _mean_max_inversion_distance      /= (double)_group_map.size()-1.0;
    _mean_max_translocation_distance  /= (double)_group_map.size()-1.0;
  }
  else
  {
    _mean_true_diversity              = 0.0;
    _mean_max_time_distance           = 0.0;
    _mean_max_generation_distance     = 0.0;
    _mean_max_point_mutation_distance = 0.0;
    _mean_max_hgt_distance            = 0.0;
    _mean_max_duplication_distance    = 0.0;
    _mean_max_deletion_distance       = 0.0;
    _mean_max_inversion_distance      = 0.0;
    _mean_max_translocation_distance  = 0.0;
  }
}

/**
 * \brief    Compute level statistics
 * \details  --
 * \param    void
 * \return   \e void
 */
void TrophicNetwork::compute_level_statistics( void )
{
  /*------------------------------------------------------------------ LEVEL 0 statistics */
  
  /* GENETIC DIVERSITY */
  _level_0_true_diversity              = 0.0;
  _level_0_max_time_distance           = 0.0;
  _level_0_max_generation_distance     = 0.0;
  _level_0_max_point_mutation_distance = 0.0;
  _level_0_max_hgt_distance            = 0.0;
  _level_0_max_duplication_distance    = 0.0;
  _level_0_max_deletion_distance       = 0.0;
  _level_0_max_inversion_distance      = 0.0;
  _level_0_max_translocation_distance  = 0.0;
  
  /* PHENOTYPE */
  _level_0_generations                = 0.0;
  _level_0_inherited_TF_amount        = 0.0;
  _level_0_inherited_E_amount         = 0.0;
  _level_0_TF_amount                  = 0.0;
  _level_0_E_amount                   = 0.0;
  _level_0_inherited_metabolic_amount = 0.0;
  _level_0_metabolic_amount           = 0.0;
  _level_0_energy                     = 0.0;
  _level_0_score                      = 0.0;
  _level_0_lifespan                   = 0.0;
  _level_0_number_of_divisions        = 0.0;
  _level_0_toxicity                   = 0.0;
  _level_0_metabolic_uptake           = 0.0;
  _level_0_metabolic_release          = 0.0;
  _level_0_metabolic_growth_rate      = 0.0;
  _level_0_Dmetabolic_growth_rate     = 0.0;
  _level_0_grn_nb_nodes               = 0.0;
  _level_0_grn_nb_edges               = 0.0;
  _level_0_metabolic_nb_nodes         = 0.0;
  _level_0_metabolic_nb_edges         = 0.0;
  _level_0_regulation_redundancy      = 0.0;
  _level_0_metabolic_redundancy       = 0.0;
  
  /* GENOME STRUCTURE */
  _level_0_genome_size                   = 0.0;
  _level_0_functional_size               = 0.0;
  _level_0_genome_nb_NC                  = 0.0;
  _level_0_genome_nb_E                   = 0.0;
  _level_0_genome_nb_TF                  = 0.0;
  _level_0_genome_nb_BS                  = 0.0;
  _level_0_genome_nb_P                   = 0.0;
  _level_0_genome_nb_inner_enzymes       = 0.0;
  _level_0_genome_nb_inflow_pumps        = 0.0;
  _level_0_genome_nb_outflow_pumps       = 0.0;
  _level_0_genome_nb_functional_regions  = 0.0;
  _level_0_genome_nb_enhancers           = 0.0;
  _level_0_genome_nb_operators           = 0.0;
  _level_0_genome_nb_E_regions           = 0.0;
  _level_0_genome_nb_TF_regions          = 0.0;
  _level_0_genome_nb_mixed_regions       = 0.0;
  _level_0_genome_functional_region_size = 0.0;
  _level_0_genome_E_region_size          = 0.0;
  _level_0_genome_TF_region_size         = 0.0;
  _level_0_genome_mixed_region_size      = 0.0;
  _level_0_genome_enhancer_size          = 0.0;
  _level_0_genome_operator_size          = 0.0;
  _level_0_genome_operon_size            = 0.0;
  _level_0_genome_E_operon_size          = 0.0;
  _level_0_genome_TF_operon_size         = 0.0;
  _level_0_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _level_0_inherited_size             = 0.0;
  _level_0_inherited_nb_E             = 0.0;
  _level_0_inherited_nb_TF            = 0.0;
  _level_0_inherited_nb_inner_enzymes = 0.0;
  _level_0_inherited_nb_inflow_pumps  = 0.0;
  _level_0_inherited_nb_outflow_pumps = 0.0;
  
  /*------------------------------------------------------------------ LEVEL 1 statistics */
  
  /* GENETIC DIVERSITY */
  _level_1_true_diversity              = 0.0;
  _level_1_max_time_distance           = 0.0;
  _level_1_max_generation_distance     = 0.0;
  _level_1_max_point_mutation_distance = 0.0;
  _level_1_max_hgt_distance            = 0.0;
  _level_1_max_duplication_distance    = 0.0;
  _level_1_max_deletion_distance       = 0.0;
  _level_1_max_inversion_distance      = 0.0;
  _level_1_max_translocation_distance  = 0.0;
  
  /* PHENOTYPE */
  _level_1_generations                = 0.0;
  _level_1_inherited_TF_amount        = 0.0;
  _level_1_inherited_E_amount         = 0.0;
  _level_1_TF_amount                  = 0.0;
  _level_1_E_amount                   = 0.0;
  _level_1_inherited_metabolic_amount = 0.0;
  _level_1_metabolic_amount           = 0.0;
  _level_1_energy                     = 0.0;
  _level_1_score                      = 0.0;
  _level_1_lifespan                   = 0.0;
  _level_1_number_of_divisions        = 0.0;
  _level_1_toxicity                   = 0.0;
  _level_1_metabolic_uptake           = 0.0;
  _level_1_metabolic_release          = 0.0;
  _level_1_metabolic_growth_rate      = 0.0;
  _level_1_Dmetabolic_growth_rate     = 0.0;
  _level_1_grn_nb_nodes               = 0.0;
  _level_1_grn_nb_edges               = 0.0;
  _level_1_metabolic_nb_nodes         = 0.0;
  _level_1_metabolic_nb_edges         = 0.0;
  _level_1_regulation_redundancy      = 0.0;
  _level_1_metabolic_redundancy       = 0.0;
  
  /* GENOME STRUCTURE */
  _level_1_genome_size                   = 0.0;
  _level_1_functional_size               = 0.0;
  _level_1_genome_nb_NC                  = 0.0;
  _level_1_genome_nb_E                   = 0.0;
  _level_1_genome_nb_TF                  = 0.0;
  _level_1_genome_nb_BS                  = 0.0;
  _level_1_genome_nb_P                   = 0.0;
  _level_1_genome_nb_inner_enzymes       = 0.0;
  _level_1_genome_nb_inflow_pumps        = 0.0;
  _level_1_genome_nb_outflow_pumps       = 0.0;
  _level_1_genome_nb_functional_regions  = 0.0;
  _level_1_genome_nb_enhancers           = 0.0;
  _level_1_genome_nb_operators           = 0.0;
  _level_1_genome_nb_E_regions           = 0.0;
  _level_1_genome_nb_TF_regions          = 0.0;
  _level_1_genome_nb_mixed_regions       = 0.0;
  _level_1_genome_functional_region_size = 0.0;
  _level_1_genome_E_region_size          = 0.0;
  _level_1_genome_TF_region_size         = 0.0;
  _level_1_genome_mixed_region_size      = 0.0;
  _level_1_genome_enhancer_size          = 0.0;
  _level_1_genome_operator_size          = 0.0;
  _level_1_genome_operon_size            = 0.0;
  _level_1_genome_E_operon_size          = 0.0;
  _level_1_genome_TF_operon_size         = 0.0;
  _level_1_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _level_1_inherited_size             = 0.0;
  _level_1_inherited_nb_E             = 0.0;
  _level_1_inherited_nb_TF            = 0.0;
  _level_1_inherited_nb_inner_enzymes = 0.0;
  _level_1_inherited_nb_inflow_pumps  = 0.0;
  _level_1_inherited_nb_outflow_pumps = 0.0;
  
  /*------------------------------------------------------------------ LEVEL 2 statistics */
  
  /* GENETIC DIVERSITY */
  _level_2_true_diversity              = 0.0;
  _level_2_max_time_distance           = 0.0;
  _level_2_max_generation_distance     = 0.0;
  _level_2_max_point_mutation_distance = 0.0;
  _level_2_max_hgt_distance            = 0.0;
  _level_2_max_duplication_distance    = 0.0;
  _level_2_max_deletion_distance       = 0.0;
  _level_2_max_inversion_distance      = 0.0;
  _level_2_max_translocation_distance  = 0.0;
  
  /* PHENOTYPE */
  _level_2_generations                = 0.0;
  _level_2_inherited_TF_amount        = 0.0;
  _level_2_inherited_E_amount         = 0.0;
  _level_2_TF_amount                  = 0.0;
  _level_2_E_amount                   = 0.0;
  _level_2_inherited_metabolic_amount = 0.0;
  _level_2_metabolic_amount           = 0.0;
  _level_2_energy                     = 0.0;
  _level_2_score                      = 0.0;
  _level_2_lifespan                   = 0.0;
  _level_2_number_of_divisions        = 0.0;
  _level_2_toxicity                   = 0.0;
  _level_2_metabolic_uptake           = 0.0;
  _level_2_metabolic_release          = 0.0;
  _level_2_metabolic_growth_rate      = 0.0;
  _level_2_Dmetabolic_growth_rate     = 0.0;
  _level_2_grn_nb_nodes               = 0.0;
  _level_2_grn_nb_edges               = 0.0;
  _level_2_metabolic_nb_nodes         = 0.0;
  _level_2_metabolic_nb_edges         = 0.0;
  _level_2_regulation_redundancy      = 0.0;
  _level_2_metabolic_redundancy       = 0.0;
  
  /* GENOME STRUCTURE */
  _level_2_genome_size                   = 0.0;
  _level_2_functional_size               = 0.0;
  _level_2_genome_nb_NC                  = 0.0;
  _level_2_genome_nb_E                   = 0.0;
  _level_2_genome_nb_TF                  = 0.0;
  _level_2_genome_nb_BS                  = 0.0;
  _level_2_genome_nb_P                   = 0.0;
  _level_2_genome_nb_inner_enzymes       = 0.0;
  _level_2_genome_nb_inflow_pumps        = 0.0;
  _level_2_genome_nb_outflow_pumps       = 0.0;
  _level_2_genome_nb_functional_regions  = 0.0;
  _level_2_genome_nb_enhancers           = 0.0;
  _level_2_genome_nb_operators           = 0.0;
  _level_2_genome_nb_E_regions           = 0.0;
  _level_2_genome_nb_TF_regions          = 0.0;
  _level_2_genome_nb_mixed_regions       = 0.0;
  _level_2_genome_functional_region_size = 0.0;
  _level_2_genome_E_region_size          = 0.0;
  _level_2_genome_TF_region_size         = 0.0;
  _level_2_genome_mixed_region_size      = 0.0;
  _level_2_genome_enhancer_size          = 0.0;
  _level_2_genome_operator_size          = 0.0;
  _level_2_genome_operon_size            = 0.0;
  _level_2_genome_E_operon_size          = 0.0;
  _level_2_genome_TF_operon_size         = 0.0;
  _level_2_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _level_2_inherited_size             = 0.0;
  _level_2_inherited_nb_E             = 0.0;
  _level_2_inherited_nb_TF            = 0.0;
  _level_2_inherited_nb_inner_enzymes = 0.0;
  _level_2_inherited_nb_inflow_pumps  = 0.0;
  _level_2_inherited_nb_outflow_pumps = 0.0;
  
  /*------------------------------------------------------------------ NO LEVEL statistics */
  
  /* GENETIC DIVERSITY */
  _no_level_true_diversity              = 0.0;
  _no_level_max_time_distance           = 0.0;
  _no_level_max_generation_distance     = 0.0;
  _no_level_max_point_mutation_distance = 0.0;
  _no_level_max_hgt_distance            = 0.0;
  _no_level_max_duplication_distance    = 0.0;
  _no_level_max_deletion_distance       = 0.0;
  _no_level_max_inversion_distance      = 0.0;
  _no_level_max_translocation_distance  = 0.0;
  
  /* PHENOTYPE */
  _no_level_generations                = 0.0;
  _no_level_inherited_TF_amount        = 0.0;
  _no_level_inherited_E_amount         = 0.0;
  _no_level_TF_amount                  = 0.0;
  _no_level_E_amount                   = 0.0;
  _no_level_inherited_metabolic_amount = 0.0;
  _no_level_metabolic_amount           = 0.0;
  _no_level_energy                     = 0.0;
  _no_level_score                      = 0.0;
  _no_level_lifespan                   = 0.0;
  _no_level_number_of_divisions        = 0.0;
  _no_level_toxicity                   = 0.0;
  _no_level_metabolic_uptake           = 0.0;
  _no_level_metabolic_release          = 0.0;
  _no_level_metabolic_growth_rate      = 0.0;
  _no_level_Dmetabolic_growth_rate     = 0.0;
  _no_level_grn_nb_nodes               = 0.0;
  _no_level_grn_nb_edges               = 0.0;
  _no_level_metabolic_nb_nodes         = 0.0;
  _no_level_metabolic_nb_edges         = 0.0;
  _no_level_regulation_redundancy      = 0.0;
  _no_level_metabolic_redundancy       = 0.0;
  
  /* GENOME STRUCTURE */
  _no_level_genome_size                   = 0.0;
  _no_level_functional_size               = 0.0;
  _no_level_genome_nb_NC                  = 0.0;
  _no_level_genome_nb_E                   = 0.0;
  _no_level_genome_nb_TF                  = 0.0;
  _no_level_genome_nb_BS                  = 0.0;
  _no_level_genome_nb_P                   = 0.0;
  _no_level_genome_nb_inner_enzymes       = 0.0;
  _no_level_genome_nb_inflow_pumps        = 0.0;
  _no_level_genome_nb_outflow_pumps       = 0.0;
  _no_level_genome_nb_functional_regions  = 0.0;
  _no_level_genome_nb_enhancers           = 0.0;
  _no_level_genome_nb_operators           = 0.0;
  _no_level_genome_nb_E_regions           = 0.0;
  _no_level_genome_nb_TF_regions          = 0.0;
  _no_level_genome_nb_mixed_regions       = 0.0;
  _no_level_genome_functional_region_size = 0.0;
  _no_level_genome_E_region_size          = 0.0;
  _no_level_genome_TF_region_size         = 0.0;
  _no_level_genome_mixed_region_size      = 0.0;
  _no_level_genome_enhancer_size          = 0.0;
  _no_level_genome_operator_size          = 0.0;
  _no_level_genome_operon_size            = 0.0;
  _no_level_genome_E_operon_size          = 0.0;
  _no_level_genome_TF_operon_size         = 0.0;
  _no_level_genome_mixed_operon_size      = 0.0;
  
  /* INHERITED STRUCTURE */
  _no_level_inherited_size             = 0.0;
  _no_level_inherited_nb_E             = 0.0;
  _no_level_inherited_nb_TF            = 0.0;
  _no_level_inherited_nb_inner_enzymes = 0.0;
  _no_level_inherited_nb_inflow_pumps  = 0.0;
  _no_level_inherited_nb_outflow_pumps = 0.0;
  
  /*------------------------------------------------------------------ explore the network */
  
  _iterator = _group_map.begin();
  while (_iterator != _group_map.end())
  {
    if (_iterator->first != 0)
    {
      /*-------------------------------------*/
      /* A) If the group is a level 0 group  */
      /*-------------------------------------*/
      if (_iterator->second->get_trophic_level() == LEVEL_0)
      {
        /* GENETIC DIVERSITY */
        _level_0_true_diversity += _iterator->second->get_true_diversity();
        if (_level_0_max_time_distance < _iterator->second->get_max_time_distance())
        {
          _level_0_max_time_distance = _iterator->second->get_max_time_distance();
        }
        if (_level_0_max_generation_distance < _iterator->second->get_max_generation_distance())
        {
          _level_0_max_generation_distance = _iterator->second->get_max_generation_distance();
        }
        if (_level_0_max_point_mutation_distance < _iterator->second->get_max_point_mutation_distance())
        {
          _level_0_max_point_mutation_distance = _iterator->second->get_max_point_mutation_distance();
        }
        if (_level_0_max_hgt_distance < _iterator->second->get_max_hgt_distance())
        {
          _level_0_max_hgt_distance = _iterator->second->get_max_hgt_distance();
        }
        if (_level_0_max_duplication_distance < _iterator->second->get_max_duplication_distance())
        {
          _level_0_max_duplication_distance = _iterator->second->get_max_duplication_distance();
        }
        if (_level_0_max_deletion_distance < _iterator->second->get_max_deletion_distance())
        {
          _level_0_max_deletion_distance = _iterator->second->get_max_deletion_distance();
        }
        if (_level_0_max_inversion_distance < _iterator->second->get_max_inversion_distance())
        {
          _level_0_max_inversion_distance = _iterator->second->get_max_inversion_distance();
        }
        if (_level_0_max_translocation_distance < _iterator->second->get_max_translocation_distance())
        {
          _level_0_max_translocation_distance = _iterator->second->get_max_translocation_distance();
        }
        
        /* PHENOTYPE */
        _level_0_generations                += _iterator->second->get_mean_generations();
        _level_0_inherited_TF_amount        += _iterator->second->get_mean_inherited_TF_amount();
        _level_0_inherited_E_amount         += _iterator->second->get_mean_inherited_E_amount();
        _level_0_TF_amount                  += _iterator->second->get_mean_TF_amount();
        _level_0_E_amount                   += _iterator->second->get_mean_E_amount();
        _level_0_inherited_metabolic_amount += _iterator->second->get_mean_inherited_metabolic_amount();
        _level_0_metabolic_amount           += _iterator->second->get_mean_metabolic_amount();
        _level_0_energy                     += _iterator->second->get_mean_energy();
        _level_0_score                      += _iterator->second->get_mean_score();
        _level_0_lifespan                   += _iterator->second->get_mean_lifespan();
        _level_0_number_of_divisions        += _iterator->second->get_mean_number_of_divisions();
        _level_0_toxicity                   += _iterator->second->get_mean_toxicity();
        _level_0_metabolic_uptake           += _iterator->second->get_mean_metabolic_uptake();
        _level_0_metabolic_release          += _iterator->second->get_mean_metabolic_release();
        _level_0_metabolic_growth_rate      += _iterator->second->get_mean_metabolic_growth_rate();
        _level_0_Dmetabolic_growth_rate     += _iterator->second->get_mean_Dmetabolic_growth_rate();
        _level_0_grn_nb_nodes               += _iterator->second->get_mean_grn_nb_nodes();
        _level_0_grn_nb_edges               += _iterator->second->get_mean_grn_nb_edges();
        _level_0_metabolic_nb_nodes         += _iterator->second->get_mean_metabolic_nb_nodes();
        _level_0_metabolic_nb_edges         += _iterator->second->get_mean_metabolic_nb_edges();
        _level_0_regulation_redundancy      += _iterator->second->get_mean_regulation_redundancy();
        _level_0_metabolic_redundancy       += _iterator->second->get_mean_metabolic_redundancy();
        
        /* GENOME STRUCTURE */
        _level_0_genome_size                   += _iterator->second->get_mean_genome_size();
        _level_0_functional_size               += _iterator->second->get_mean_functional_size();
        _level_0_genome_nb_NC                  += _iterator->second->get_mean_genome_nb_NC();
        _level_0_genome_nb_E                   += _iterator->second->get_mean_genome_nb_E();
        _level_0_genome_nb_TF                  += _iterator->second->get_mean_genome_nb_TF();
        _level_0_genome_nb_BS                  += _iterator->second->get_mean_genome_nb_BS();
        _level_0_genome_nb_P                   += _iterator->second->get_mean_genome_nb_P();
        _level_0_genome_nb_inner_enzymes       += _iterator->second->get_mean_genome_nb_inner_enzymes();
        _level_0_genome_nb_inflow_pumps        += _iterator->second->get_mean_genome_nb_inflow_pumps();
        _level_0_genome_nb_outflow_pumps       += _iterator->second->get_mean_genome_nb_outflow_pumps();
        _level_0_genome_nb_functional_regions  += _iterator->second->get_mean_genome_nb_functional_regions();
        _level_0_genome_nb_enhancers           += _iterator->second->get_mean_genome_nb_enhancers();
        _level_0_genome_nb_operators           += _iterator->second->get_mean_genome_nb_operators();
        _level_0_genome_nb_E_regions           += _iterator->second->get_mean_genome_nb_E_regions();
        _level_0_genome_nb_TF_regions          += _iterator->second->get_mean_genome_nb_TF_regions();
        _level_0_genome_nb_mixed_regions       += _iterator->second->get_mean_genome_nb_mixed_regions();
        _level_0_genome_functional_region_size += _iterator->second->get_mean_genome_functional_region_size();
        _level_0_genome_E_region_size          += _iterator->second->get_mean_genome_E_region_size();
        _level_0_genome_TF_region_size         += _iterator->second->get_mean_genome_TF_region_size();
        _level_0_genome_mixed_region_size      += _iterator->second->get_mean_genome_mixed_region_size();
        _level_0_genome_enhancer_size          += _iterator->second->get_mean_genome_enhancer_size();
        _level_0_genome_operator_size          += _iterator->second->get_mean_genome_operator_size();
        _level_0_genome_operon_size            += _iterator->second->get_mean_genome_operon_size();
        _level_0_genome_E_operon_size          += _iterator->second->get_mean_genome_E_operon_size();
        _level_0_genome_TF_operon_size         += _iterator->second->get_mean_genome_TF_operon_size();
        _level_0_genome_mixed_operon_size      += _iterator->second->get_mean_genome_mixed_operon_size();
        
        /* INHERITED STRUCTURE */
        _level_0_inherited_size             += _iterator->second->get_mean_inherited_size();
        _level_0_inherited_nb_E             += _iterator->second->get_mean_inherited_nb_E();
        _level_0_inherited_nb_TF            += _iterator->second->get_mean_inherited_nb_TF();
        _level_0_inherited_nb_inner_enzymes += _iterator->second->get_mean_inherited_nb_inner_enzymes();
        _level_0_inherited_nb_inflow_pumps  += _iterator->second->get_mean_inherited_nb_inflow_pumps();
        _level_0_inherited_nb_outflow_pumps += _iterator->second->get_mean_inherited_nb_outflow_pumps();
      }
      /*-------------------------------------*/
      /* B) If the group is a level 1 group  */
      /*-------------------------------------*/
      if (_iterator->second->get_trophic_level() == LEVEL_1)
      {
        /* GENETIC DIVERSITY */
        _level_1_true_diversity += _iterator->second->get_true_diversity();
        if (_level_1_max_time_distance < _iterator->second->get_max_time_distance())
        {
          _level_1_max_time_distance = _iterator->second->get_max_time_distance();
        }
        if (_level_1_max_generation_distance < _iterator->second->get_max_generation_distance())
        {
          _level_1_max_generation_distance = _iterator->second->get_max_generation_distance();
        }
        if (_level_1_max_point_mutation_distance < _iterator->second->get_max_point_mutation_distance())
        {
          _level_1_max_point_mutation_distance = _iterator->second->get_max_point_mutation_distance();
        }
        if (_level_1_max_hgt_distance < _iterator->second->get_max_hgt_distance())
        {
          _level_1_max_hgt_distance = _iterator->second->get_max_hgt_distance();
        }
        if (_level_1_max_duplication_distance < _iterator->second->get_max_duplication_distance())
        {
          _level_1_max_duplication_distance = _iterator->second->get_max_duplication_distance();
        }
        if (_level_1_max_deletion_distance < _iterator->second->get_max_deletion_distance())
        {
          _level_1_max_deletion_distance = _iterator->second->get_max_deletion_distance();
        }
        if (_level_1_max_inversion_distance < _iterator->second->get_max_inversion_distance())
        {
          _level_1_max_inversion_distance = _iterator->second->get_max_inversion_distance();
        }
        if (_level_1_max_translocation_distance < _iterator->second->get_max_translocation_distance())
        {
          _level_1_max_translocation_distance = _iterator->second->get_max_translocation_distance();
        }
        
        /* PHENOTYPE */
        _level_1_generations                += _iterator->second->get_mean_generations();
        _level_1_inherited_TF_amount        += _iterator->second->get_mean_inherited_TF_amount();
        _level_1_inherited_E_amount         += _iterator->second->get_mean_inherited_E_amount();
        _level_1_TF_amount                  += _iterator->second->get_mean_TF_amount();
        _level_1_E_amount                   += _iterator->second->get_mean_E_amount();
        _level_1_inherited_metabolic_amount += _iterator->second->get_mean_inherited_metabolic_amount();
        _level_1_metabolic_amount           += _iterator->second->get_mean_metabolic_amount();
        _level_1_energy                     += _iterator->second->get_mean_energy();
        _level_1_score                      += _iterator->second->get_mean_score();
        _level_1_lifespan                   += _iterator->second->get_mean_lifespan();
        _level_1_number_of_divisions        += _iterator->second->get_mean_number_of_divisions();
        _level_1_toxicity                   += _iterator->second->get_mean_toxicity();
        _level_1_metabolic_uptake           += _iterator->second->get_mean_metabolic_uptake();
        _level_1_metabolic_release          += _iterator->second->get_mean_metabolic_release();
        _level_1_metabolic_growth_rate      += _iterator->second->get_mean_metabolic_growth_rate();
        _level_1_Dmetabolic_growth_rate     += _iterator->second->get_mean_Dmetabolic_growth_rate();
        _level_1_grn_nb_nodes               += _iterator->second->get_mean_grn_nb_nodes();
        _level_1_grn_nb_edges               += _iterator->second->get_mean_grn_nb_edges();
        _level_1_metabolic_nb_nodes         += _iterator->second->get_mean_metabolic_nb_nodes();
        _level_1_metabolic_nb_edges         += _iterator->second->get_mean_metabolic_nb_edges();
        _level_1_regulation_redundancy      += _iterator->second->get_mean_regulation_redundancy();
        _level_1_metabolic_redundancy       += _iterator->second->get_mean_metabolic_redundancy();
        
        /* GENOME STRUCTURE */
        _level_1_genome_size                   += _iterator->second->get_mean_genome_size();
        _level_1_functional_size               += _iterator->second->get_mean_functional_size();
        _level_1_genome_nb_NC                  += _iterator->second->get_mean_genome_nb_NC();
        _level_1_genome_nb_E                   += _iterator->second->get_mean_genome_nb_E();
        _level_1_genome_nb_TF                  += _iterator->second->get_mean_genome_nb_TF();
        _level_1_genome_nb_BS                  += _iterator->second->get_mean_genome_nb_BS();
        _level_1_genome_nb_P                   += _iterator->second->get_mean_genome_nb_P();
        _level_1_genome_nb_inner_enzymes       += _iterator->second->get_mean_genome_nb_inner_enzymes();
        _level_1_genome_nb_inflow_pumps        += _iterator->second->get_mean_genome_nb_inflow_pumps();
        _level_1_genome_nb_outflow_pumps       += _iterator->second->get_mean_genome_nb_outflow_pumps();
        _level_1_genome_nb_functional_regions  += _iterator->second->get_mean_genome_nb_functional_regions();
        _level_1_genome_nb_enhancers           += _iterator->second->get_mean_genome_nb_enhancers();
        _level_1_genome_nb_operators           += _iterator->second->get_mean_genome_nb_operators();
        _level_1_genome_nb_E_regions           += _iterator->second->get_mean_genome_nb_E_regions();
        _level_1_genome_nb_TF_regions          += _iterator->second->get_mean_genome_nb_TF_regions();
        _level_1_genome_nb_mixed_regions       += _iterator->second->get_mean_genome_nb_mixed_regions();
        _level_1_genome_functional_region_size += _iterator->second->get_mean_genome_functional_region_size();
        _level_1_genome_E_region_size          += _iterator->second->get_mean_genome_E_region_size();
        _level_1_genome_TF_region_size         += _iterator->second->get_mean_genome_TF_region_size();
        _level_1_genome_mixed_region_size      += _iterator->second->get_mean_genome_mixed_region_size();
        _level_1_genome_enhancer_size          += _iterator->second->get_mean_genome_enhancer_size();
        _level_1_genome_operator_size          += _iterator->second->get_mean_genome_operator_size();
        _level_1_genome_operon_size            += _iterator->second->get_mean_genome_operon_size();
        _level_1_genome_E_operon_size          += _iterator->second->get_mean_genome_E_operon_size();
        _level_1_genome_TF_operon_size         += _iterator->second->get_mean_genome_TF_operon_size();
        _level_1_genome_mixed_operon_size      += _iterator->second->get_mean_genome_mixed_operon_size();
        
        /* INHERITED STRUCTURE */
        _level_1_inherited_size             += _iterator->second->get_mean_inherited_size();
        _level_1_inherited_nb_E             += _iterator->second->get_mean_inherited_nb_E();
        _level_1_inherited_nb_TF            += _iterator->second->get_mean_inherited_nb_TF();
        _level_1_inherited_nb_inner_enzymes += _iterator->second->get_mean_inherited_nb_inner_enzymes();
        _level_1_inherited_nb_inflow_pumps  += _iterator->second->get_mean_inherited_nb_inflow_pumps();
        _level_1_inherited_nb_outflow_pumps += _iterator->second->get_mean_inherited_nb_outflow_pumps();
      }
      /*-------------------------------------*/
      /* C) If the group is a level 2 group  */
      /*-------------------------------------*/
      if (_iterator->second->get_trophic_level() == LEVEL_2)
      {
        /* GENETIC DIVERSITY */
        _level_2_true_diversity += _iterator->second->get_true_diversity();
        if (_level_2_max_time_distance < _iterator->second->get_max_time_distance())
        {
          _level_2_max_time_distance = _iterator->second->get_max_time_distance();
        }
        if (_level_2_max_generation_distance < _iterator->second->get_max_generation_distance())
        {
          _level_2_max_generation_distance = _iterator->second->get_max_generation_distance();
        }
        if (_level_2_max_point_mutation_distance < _iterator->second->get_max_point_mutation_distance())
        {
          _level_2_max_point_mutation_distance = _iterator->second->get_max_point_mutation_distance();
        }
        if (_level_2_max_hgt_distance < _iterator->second->get_max_hgt_distance())
        {
          _level_2_max_hgt_distance = _iterator->second->get_max_hgt_distance();
        }
        if (_level_2_max_duplication_distance < _iterator->second->get_max_duplication_distance())
        {
          _level_2_max_duplication_distance = _iterator->second->get_max_duplication_distance();
        }
        if (_level_2_max_deletion_distance < _iterator->second->get_max_deletion_distance())
        {
          _level_2_max_deletion_distance = _iterator->second->get_max_deletion_distance();
        }
        if (_level_2_max_inversion_distance < _iterator->second->get_max_inversion_distance())
        {
          _level_2_max_inversion_distance = _iterator->second->get_max_inversion_distance();
        }
        if (_level_2_max_translocation_distance < _iterator->second->get_max_translocation_distance())
        {
          _level_2_max_translocation_distance = _iterator->second->get_max_translocation_distance();
        }
        
        /* PHENOTYPE */
        _level_2_generations                += _iterator->second->get_mean_generations();
        _level_2_inherited_TF_amount        += _iterator->second->get_mean_inherited_TF_amount();
        _level_2_inherited_E_amount         += _iterator->second->get_mean_inherited_E_amount();
        _level_2_TF_amount                  += _iterator->second->get_mean_TF_amount();
        _level_2_E_amount                   += _iterator->second->get_mean_E_amount();
        _level_2_inherited_metabolic_amount += _iterator->second->get_mean_inherited_metabolic_amount();
        _level_2_metabolic_amount           += _iterator->second->get_mean_metabolic_amount();
        _level_2_energy                     += _iterator->second->get_mean_energy();
        _level_2_score                      += _iterator->second->get_mean_score();
        _level_2_lifespan                   += _iterator->second->get_mean_lifespan();
        _level_2_number_of_divisions        += _iterator->second->get_mean_number_of_divisions();
        _level_2_toxicity                   += _iterator->second->get_mean_toxicity();
        _level_2_metabolic_uptake           += _iterator->second->get_mean_metabolic_uptake();
        _level_2_metabolic_release          += _iterator->second->get_mean_metabolic_release();
        _level_2_metabolic_growth_rate      += _iterator->second->get_mean_metabolic_growth_rate();
        _level_2_Dmetabolic_growth_rate     += _iterator->second->get_mean_Dmetabolic_growth_rate();
        _level_2_grn_nb_nodes               += _iterator->second->get_mean_grn_nb_nodes();
        _level_2_grn_nb_edges               += _iterator->second->get_mean_grn_nb_edges();
        _level_2_metabolic_nb_nodes         += _iterator->second->get_mean_metabolic_nb_nodes();
        _level_2_metabolic_nb_edges         += _iterator->second->get_mean_metabolic_nb_edges();
        _level_2_regulation_redundancy      += _iterator->second->get_mean_regulation_redundancy();
        _level_2_metabolic_redundancy       += _iterator->second->get_mean_metabolic_redundancy();
        
        /* GENOME STRUCTURE */
        _level_2_genome_size                   += _iterator->second->get_mean_genome_size();
        _level_2_functional_size               += _iterator->second->get_mean_functional_size();
        _level_2_genome_nb_NC                  += _iterator->second->get_mean_genome_nb_NC();
        _level_2_genome_nb_E                   += _iterator->second->get_mean_genome_nb_E();
        _level_2_genome_nb_TF                  += _iterator->second->get_mean_genome_nb_TF();
        _level_2_genome_nb_BS                  += _iterator->second->get_mean_genome_nb_BS();
        _level_2_genome_nb_P                   += _iterator->second->get_mean_genome_nb_P();
        _level_2_genome_nb_inner_enzymes       += _iterator->second->get_mean_genome_nb_inner_enzymes();
        _level_2_genome_nb_inflow_pumps        += _iterator->second->get_mean_genome_nb_inflow_pumps();
        _level_2_genome_nb_outflow_pumps       += _iterator->second->get_mean_genome_nb_outflow_pumps();
        _level_2_genome_nb_functional_regions  += _iterator->second->get_mean_genome_nb_functional_regions();
        _level_2_genome_nb_enhancers           += _iterator->second->get_mean_genome_nb_enhancers();
        _level_2_genome_nb_operators           += _iterator->second->get_mean_genome_nb_operators();
        _level_2_genome_nb_E_regions           += _iterator->second->get_mean_genome_nb_E_regions();
        _level_2_genome_nb_TF_regions          += _iterator->second->get_mean_genome_nb_TF_regions();
        _level_2_genome_nb_mixed_regions       += _iterator->second->get_mean_genome_nb_mixed_regions();
        _level_2_genome_functional_region_size += _iterator->second->get_mean_genome_functional_region_size();
        _level_2_genome_E_region_size          += _iterator->second->get_mean_genome_E_region_size();
        _level_2_genome_TF_region_size         += _iterator->second->get_mean_genome_TF_region_size();
        _level_2_genome_mixed_region_size      += _iterator->second->get_mean_genome_mixed_region_size();
        _level_2_genome_enhancer_size          += _iterator->second->get_mean_genome_enhancer_size();
        _level_2_genome_operator_size          += _iterator->second->get_mean_genome_operator_size();
        _level_2_genome_operon_size            += _iterator->second->get_mean_genome_operon_size();
        _level_2_genome_E_operon_size          += _iterator->second->get_mean_genome_E_operon_size();
        _level_2_genome_TF_operon_size         += _iterator->second->get_mean_genome_TF_operon_size();
        _level_2_genome_mixed_operon_size      += _iterator->second->get_mean_genome_mixed_operon_size();
        
        /* INHERITED STRUCTURE */
        _level_2_inherited_size             += _iterator->second->get_mean_inherited_size();
        _level_2_inherited_nb_E             += _iterator->second->get_mean_inherited_nb_E();
        _level_2_inherited_nb_TF            += _iterator->second->get_mean_inherited_nb_TF();
        _level_2_inherited_nb_inner_enzymes += _iterator->second->get_mean_inherited_nb_inner_enzymes();
        _level_2_inherited_nb_inflow_pumps  += _iterator->second->get_mean_inherited_nb_inflow_pumps();
        _level_2_inherited_nb_outflow_pumps += _iterator->second->get_mean_inherited_nb_outflow_pumps();
      }
      /*-------------------------------------*/
      /* D) If the group is a no level group */
      /*-------------------------------------*/
      if (_iterator->second->get_trophic_level() == NO_LEVEL)
      {
        /* GENETIC DIVERSITY */
        _no_level_true_diversity += _iterator->second->get_true_diversity();
        if (_no_level_max_time_distance < _iterator->second->get_max_time_distance())
        {
          _no_level_max_time_distance = _iterator->second->get_max_time_distance();
        }
        if (_no_level_max_generation_distance < _iterator->second->get_max_generation_distance())
        {
          _no_level_max_generation_distance = _iterator->second->get_max_generation_distance();
        }
        if (_no_level_max_point_mutation_distance < _iterator->second->get_max_point_mutation_distance())
        {
          _no_level_max_point_mutation_distance = _iterator->second->get_max_point_mutation_distance();
        }
        if (_no_level_max_hgt_distance < _iterator->second->get_max_hgt_distance())
        {
          _no_level_max_hgt_distance = _iterator->second->get_max_hgt_distance();
        }
        if (_no_level_max_duplication_distance < _iterator->second->get_max_duplication_distance())
        {
          _no_level_max_duplication_distance = _iterator->second->get_max_duplication_distance();
        }
        if (_no_level_max_deletion_distance < _iterator->second->get_max_deletion_distance())
        {
          _no_level_max_deletion_distance = _iterator->second->get_max_deletion_distance();
        }
        if (_no_level_max_inversion_distance < _iterator->second->get_max_inversion_distance())
        {
          _no_level_max_inversion_distance = _iterator->second->get_max_inversion_distance();
        }
        if (_no_level_max_translocation_distance < _iterator->second->get_max_translocation_distance())
        {
          _no_level_max_translocation_distance = _iterator->second->get_max_translocation_distance();
        }
        
        /* PHENOTYPE */
        _no_level_generations                += _iterator->second->get_mean_generations();
        _no_level_inherited_TF_amount        += _iterator->second->get_mean_inherited_TF_amount();
        _no_level_inherited_E_amount         += _iterator->second->get_mean_inherited_E_amount();
        _no_level_TF_amount                  += _iterator->second->get_mean_TF_amount();
        _no_level_E_amount                   += _iterator->second->get_mean_E_amount();
        _no_level_inherited_metabolic_amount += _iterator->second->get_mean_inherited_metabolic_amount();
        _no_level_metabolic_amount           += _iterator->second->get_mean_metabolic_amount();
        _no_level_energy                     += _iterator->second->get_mean_energy();
        _no_level_score                      += _iterator->second->get_mean_score();
        _no_level_lifespan                   += _iterator->second->get_mean_lifespan();
        _no_level_number_of_divisions        += _iterator->second->get_mean_number_of_divisions();
        _no_level_toxicity                   += _iterator->second->get_mean_toxicity();
        _no_level_metabolic_uptake           += _iterator->second->get_mean_metabolic_uptake();
        _no_level_metabolic_release          += _iterator->second->get_mean_metabolic_release();
        _no_level_metabolic_growth_rate      += _iterator->second->get_mean_metabolic_growth_rate();
        _no_level_Dmetabolic_growth_rate     += _iterator->second->get_mean_Dmetabolic_growth_rate();
        _no_level_grn_nb_nodes               += _iterator->second->get_mean_grn_nb_nodes();
        _no_level_grn_nb_edges               += _iterator->second->get_mean_grn_nb_edges();
        _no_level_metabolic_nb_nodes         += _iterator->second->get_mean_metabolic_nb_nodes();
        _no_level_metabolic_nb_edges         += _iterator->second->get_mean_metabolic_nb_edges();
        _no_level_regulation_redundancy      += _iterator->second->get_mean_regulation_redundancy();
        _no_level_metabolic_redundancy       += _iterator->second->get_mean_metabolic_redundancy();
        
        /* GENOME STRUCTURE */
        _no_level_genome_size                   += _iterator->second->get_mean_genome_size();
        _no_level_functional_size               += _iterator->second->get_mean_functional_size();
        _no_level_genome_nb_NC                  += _iterator->second->get_mean_genome_nb_NC();
        _no_level_genome_nb_E                   += _iterator->second->get_mean_genome_nb_E();
        _no_level_genome_nb_TF                  += _iterator->second->get_mean_genome_nb_TF();
        _no_level_genome_nb_BS                  += _iterator->second->get_mean_genome_nb_BS();
        _no_level_genome_nb_P                   += _iterator->second->get_mean_genome_nb_P();
        _no_level_genome_nb_inner_enzymes       += _iterator->second->get_mean_genome_nb_inner_enzymes();
        _no_level_genome_nb_inflow_pumps        += _iterator->second->get_mean_genome_nb_inflow_pumps();
        _no_level_genome_nb_outflow_pumps       += _iterator->second->get_mean_genome_nb_outflow_pumps();
        _no_level_genome_nb_functional_regions  += _iterator->second->get_mean_genome_nb_functional_regions();
        _no_level_genome_nb_enhancers           += _iterator->second->get_mean_genome_nb_enhancers();
        _no_level_genome_nb_operators           += _iterator->second->get_mean_genome_nb_operators();
        _no_level_genome_nb_E_regions           += _iterator->second->get_mean_genome_nb_E_regions();
        _no_level_genome_nb_TF_regions          += _iterator->second->get_mean_genome_nb_TF_regions();
        _no_level_genome_nb_mixed_regions       += _iterator->second->get_mean_genome_nb_mixed_regions();
        _no_level_genome_functional_region_size += _iterator->second->get_mean_genome_functional_region_size();
        _no_level_genome_E_region_size          += _iterator->second->get_mean_genome_E_region_size();
        _no_level_genome_TF_region_size         += _iterator->second->get_mean_genome_TF_region_size();
        _no_level_genome_mixed_region_size      += _iterator->second->get_mean_genome_mixed_region_size();
        _no_level_genome_enhancer_size          += _iterator->second->get_mean_genome_enhancer_size();
        _no_level_genome_operator_size          += _iterator->second->get_mean_genome_operator_size();
        _no_level_genome_operon_size            += _iterator->second->get_mean_genome_operon_size();
        _no_level_genome_E_operon_size          += _iterator->second->get_mean_genome_E_operon_size();
        _no_level_genome_TF_operon_size         += _iterator->second->get_mean_genome_TF_operon_size();
        _no_level_genome_mixed_operon_size      += _iterator->second->get_mean_genome_mixed_operon_size();
        
        /* INHERITED STRUCTURE */
        _no_level_inherited_size             += _iterator->second->get_mean_inherited_size();
        _no_level_inherited_nb_E             += _iterator->second->get_mean_inherited_nb_E();
        _no_level_inherited_nb_TF            += _iterator->second->get_mean_inherited_nb_TF();
        _no_level_inherited_nb_inner_enzymes += _iterator->second->get_mean_inherited_nb_inner_enzymes();
        _no_level_inherited_nb_inflow_pumps  += _iterator->second->get_mean_inherited_nb_inflow_pumps();
        _no_level_inherited_nb_outflow_pumps += _iterator->second->get_mean_inherited_nb_outflow_pumps();
      }
    }
    _iterator++;
  }
  
  /*------------------------------------------------------------------ Compute mean values */
  
  if (_nb_level_0_cells > 0)
  {
    /* GENETIC DIVERSITY */
    _level_0_true_diversity              /= (double)_nb_level_0_cells;
    _level_0_max_time_distance           /= (double)_nb_level_0_cells;
    _level_0_max_generation_distance     /= (double)_nb_level_0_cells;
    _level_0_max_point_mutation_distance /= (double)_nb_level_0_cells;
    _level_0_max_hgt_distance            /= (double)_nb_level_0_cells;
    _level_0_max_duplication_distance    /= (double)_nb_level_0_cells;
    _level_0_max_deletion_distance       /= (double)_nb_level_0_cells;
    _level_0_max_inversion_distance      /= (double)_nb_level_0_cells;
    _level_0_max_translocation_distance  /= (double)_nb_level_0_cells;
    
    /* PHENOTYPE */
    _level_0_generations                /= (double)_nb_level_0_cells;
    _level_0_inherited_TF_amount        /= (double)_nb_level_0_cells;
    _level_0_inherited_E_amount         /= (double)_nb_level_0_cells;
    _level_0_TF_amount                  /= (double)_nb_level_0_cells;
    _level_0_E_amount                   /= (double)_nb_level_0_cells;
    _level_0_inherited_metabolic_amount /= (double)_nb_level_0_cells;
    _level_0_metabolic_amount           /= (double)_nb_level_0_cells;
    _level_0_energy                     /= (double)_nb_level_0_cells;
    _level_0_score                      /= (double)_nb_level_0_cells;
    _level_0_lifespan                   /= (double)_nb_level_0_cells;
    _level_0_number_of_divisions        /= (double)_nb_level_0_cells;
    _level_0_toxicity                   /= (double)_nb_level_0_cells;
    _level_0_metabolic_uptake           /= (double)_nb_level_0_cells;
    _level_0_metabolic_release          /= (double)_nb_level_0_cells;
    _level_0_metabolic_growth_rate      /= (double)_nb_level_0_cells;
    _level_0_Dmetabolic_growth_rate     /= (double)_nb_level_0_cells;
    _level_0_grn_nb_nodes               /= (double)_nb_level_0_cells;
    _level_0_grn_nb_edges               /= (double)_nb_level_0_cells;
    _level_0_metabolic_nb_nodes         /= (double)_nb_level_0_cells;
    _level_0_metabolic_nb_edges         /= (double)_nb_level_0_cells;
    _level_0_regulation_redundancy      /= (double)_nb_level_0_cells;
    _level_0_metabolic_redundancy       /= (double)_nb_level_0_cells;
    
    /* GENOME STRUCTURE */
    _level_0_genome_size                   /= (double)_nb_level_0_cells;
    _level_0_functional_size               /= (double)_nb_level_0_cells;
    _level_0_genome_nb_NC                  /= (double)_nb_level_0_cells;
    _level_0_genome_nb_E                   /= (double)_nb_level_0_cells;
    _level_0_genome_nb_TF                  /= (double)_nb_level_0_cells;
    _level_0_genome_nb_BS                  /= (double)_nb_level_0_cells;
    _level_0_genome_nb_P                   /= (double)_nb_level_0_cells;
    _level_0_genome_nb_inner_enzymes       /= (double)_nb_level_0_cells;
    _level_0_genome_nb_inflow_pumps        /= (double)_nb_level_0_cells;
    _level_0_genome_nb_outflow_pumps       /= (double)_nb_level_0_cells;
    _level_0_genome_nb_functional_regions  /= (double)_nb_level_0_cells;
    _level_0_genome_nb_enhancers           /= (double)_nb_level_0_cells;
    _level_0_genome_nb_operators           /= (double)_nb_level_0_cells;
    _level_0_genome_nb_E_regions           /= (double)_nb_level_0_cells;
    _level_0_genome_nb_TF_regions          /= (double)_nb_level_0_cells;
    _level_0_genome_nb_mixed_regions       /= (double)_nb_level_0_cells;
    _level_0_genome_functional_region_size /= (double)_nb_level_0_cells;
    _level_0_genome_E_region_size          /= (double)_nb_level_0_cells;
    _level_0_genome_TF_region_size         /= (double)_nb_level_0_cells;
    _level_0_genome_mixed_region_size      /= (double)_nb_level_0_cells;
    _level_0_genome_enhancer_size          /= (double)_nb_level_0_cells;
    _level_0_genome_operator_size          /= (double)_nb_level_0_cells;
    _level_0_genome_operon_size            /= (double)_nb_level_0_cells;
    _level_0_genome_E_operon_size          /= (double)_nb_level_0_cells;
    _level_0_genome_TF_operon_size         /= (double)_nb_level_0_cells;
    _level_0_genome_mixed_operon_size      /= (double)_nb_level_0_cells;
    
    /* INHERITED STRUCTURE */
    _level_0_inherited_size             /= (double)_nb_level_0_cells;
    _level_0_inherited_nb_E             /= (double)_nb_level_0_cells;
    _level_0_inherited_nb_TF            /= (double)_nb_level_0_cells;
    _level_0_inherited_nb_inner_enzymes /= (double)_nb_level_0_cells;
    _level_0_inherited_nb_inflow_pumps  /= (double)_nb_level_0_cells;
    _level_0_inherited_nb_outflow_pumps /= (double)_nb_level_0_cells;
  }
  
  if (_nb_level_1_cells > 0)
  {
    /* GENETIC DIVERSITY */
    _level_1_true_diversity              /= (double)_nb_level_1_cells;
    _level_1_max_time_distance           /= (double)_nb_level_1_cells;
    _level_1_max_generation_distance     /= (double)_nb_level_1_cells;
    _level_1_max_point_mutation_distance /= (double)_nb_level_1_cells;
    _level_1_max_hgt_distance            /= (double)_nb_level_1_cells;
    _level_1_max_duplication_distance    /= (double)_nb_level_1_cells;
    _level_1_max_deletion_distance       /= (double)_nb_level_1_cells;
    _level_1_max_inversion_distance      /= (double)_nb_level_1_cells;
    _level_1_max_translocation_distance  /= (double)_nb_level_1_cells;
    
    /* PHENOTYPE */
    _level_1_generations                /= (double)_nb_level_1_cells;
    _level_1_inherited_TF_amount        /= (double)_nb_level_1_cells;
    _level_1_inherited_E_amount         /= (double)_nb_level_1_cells;
    _level_1_TF_amount                  /= (double)_nb_level_1_cells;
    _level_1_E_amount                   /= (double)_nb_level_1_cells;
    _level_1_inherited_metabolic_amount /= (double)_nb_level_1_cells;
    _level_1_metabolic_amount           /= (double)_nb_level_1_cells;
    _level_1_energy                     /= (double)_nb_level_1_cells;
    _level_1_score                      /= (double)_nb_level_1_cells;
    _level_1_lifespan                   /= (double)_nb_level_1_cells;
    _level_1_number_of_divisions        /= (double)_nb_level_1_cells;
    _level_1_toxicity                   /= (double)_nb_level_1_cells;
    _level_1_metabolic_uptake           /= (double)_nb_level_1_cells;
    _level_1_metabolic_release          /= (double)_nb_level_1_cells;
    _level_1_metabolic_growth_rate      /= (double)_nb_level_1_cells;
    _level_1_Dmetabolic_growth_rate     /= (double)_nb_level_1_cells;
    _level_1_grn_nb_nodes               /= (double)_nb_level_1_cells;
    _level_1_grn_nb_edges               /= (double)_nb_level_1_cells;
    _level_1_metabolic_nb_nodes         /= (double)_nb_level_1_cells;
    _level_1_metabolic_nb_edges         /= (double)_nb_level_1_cells;
    _level_1_regulation_redundancy      /= (double)_nb_level_1_cells;
    _level_1_metabolic_redundancy       /= (double)_nb_level_1_cells;
    
    /* GENOME STRUCTURE */
    _level_1_genome_size                   /= (double)_nb_level_1_cells;
    _level_1_functional_size               /= (double)_nb_level_1_cells;
    _level_1_genome_nb_NC                  /= (double)_nb_level_1_cells;
    _level_1_genome_nb_E                   /= (double)_nb_level_1_cells;
    _level_1_genome_nb_TF                  /= (double)_nb_level_1_cells;
    _level_1_genome_nb_BS                  /= (double)_nb_level_1_cells;
    _level_1_genome_nb_P                   /= (double)_nb_level_1_cells;
    _level_1_genome_nb_inner_enzymes       /= (double)_nb_level_1_cells;
    _level_1_genome_nb_inflow_pumps        /= (double)_nb_level_1_cells;
    _level_1_genome_nb_outflow_pumps       /= (double)_nb_level_1_cells;
    _level_1_genome_nb_functional_regions  /= (double)_nb_level_1_cells;
    _level_1_genome_nb_enhancers           /= (double)_nb_level_1_cells;
    _level_1_genome_nb_operators           /= (double)_nb_level_1_cells;
    _level_1_genome_nb_E_regions           /= (double)_nb_level_1_cells;
    _level_1_genome_nb_TF_regions          /= (double)_nb_level_1_cells;
    _level_1_genome_nb_mixed_regions       /= (double)_nb_level_1_cells;
    _level_1_genome_functional_region_size /= (double)_nb_level_1_cells;
    _level_1_genome_E_region_size          /= (double)_nb_level_1_cells;
    _level_1_genome_TF_region_size         /= (double)_nb_level_1_cells;
    _level_1_genome_mixed_region_size      /= (double)_nb_level_1_cells;
    _level_1_genome_enhancer_size          /= (double)_nb_level_1_cells;
    _level_1_genome_operator_size          /= (double)_nb_level_1_cells;
    _level_1_genome_operon_size            /= (double)_nb_level_1_cells;
    _level_1_genome_E_operon_size          /= (double)_nb_level_1_cells;
    _level_1_genome_TF_operon_size         /= (double)_nb_level_1_cells;
    _level_1_genome_mixed_operon_size      /= (double)_nb_level_1_cells;
    
    /* INHERITED STRUCTURE */
    _level_1_inherited_size             /= (double)_nb_level_1_cells;
    _level_1_inherited_nb_E             /= (double)_nb_level_1_cells;
    _level_1_inherited_nb_TF            /= (double)_nb_level_1_cells;
    _level_1_inherited_nb_inner_enzymes /= (double)_nb_level_1_cells;
    _level_1_inherited_nb_inflow_pumps  /= (double)_nb_level_1_cells;
    _level_1_inherited_nb_outflow_pumps /= (double)_nb_level_1_cells;
  }
  
  if (_nb_level_2_cells > 0)
  {
    /* GENETIC DIVERSITY */
    _level_2_true_diversity              /= (double)_nb_level_2_cells;
    _level_2_max_time_distance           /= (double)_nb_level_2_cells;
    _level_2_max_generation_distance     /= (double)_nb_level_2_cells;
    _level_2_max_point_mutation_distance /= (double)_nb_level_2_cells;
    _level_2_max_hgt_distance            /= (double)_nb_level_2_cells;
    _level_2_max_duplication_distance    /= (double)_nb_level_2_cells;
    _level_2_max_deletion_distance       /= (double)_nb_level_2_cells;
    _level_2_max_inversion_distance      /= (double)_nb_level_2_cells;
    _level_2_max_translocation_distance  /= (double)_nb_level_2_cells;
    
    /* PHENOTYPE */
    _level_2_generations                /= (double)_nb_level_2_cells;
    _level_2_inherited_TF_amount        /= (double)_nb_level_2_cells;
    _level_2_inherited_E_amount         /= (double)_nb_level_2_cells;
    _level_2_TF_amount                  /= (double)_nb_level_2_cells;
    _level_2_E_amount                   /= (double)_nb_level_2_cells;
    _level_2_inherited_metabolic_amount /= (double)_nb_level_2_cells;
    _level_2_metabolic_amount           /= (double)_nb_level_2_cells;
    _level_2_energy                     /= (double)_nb_level_2_cells;
    _level_2_score                      /= (double)_nb_level_2_cells;
    _level_2_lifespan                   /= (double)_nb_level_2_cells;
    _level_2_number_of_divisions        /= (double)_nb_level_2_cells;
    _level_2_toxicity                   /= (double)_nb_level_2_cells;
    _level_2_metabolic_uptake           /= (double)_nb_level_2_cells;
    _level_2_metabolic_release          /= (double)_nb_level_2_cells;
    _level_2_metabolic_growth_rate      /= (double)_nb_level_2_cells;
    _level_2_Dmetabolic_growth_rate     /= (double)_nb_level_2_cells;
    _level_2_grn_nb_nodes               /= (double)_nb_level_2_cells;
    _level_2_grn_nb_edges               /= (double)_nb_level_2_cells;
    _level_2_metabolic_nb_nodes         /= (double)_nb_level_2_cells;
    _level_2_metabolic_nb_edges         /= (double)_nb_level_2_cells;
    _level_2_regulation_redundancy      /= (double)_nb_level_2_cells;
    _level_2_metabolic_redundancy       /= (double)_nb_level_2_cells;
    
    /* GENOME STRUCTURE */
    _level_2_genome_size                   /= (double)_nb_level_2_cells;
    _level_2_functional_size               /= (double)_nb_level_2_cells;
    _level_2_genome_nb_NC                  /= (double)_nb_level_2_cells;
    _level_2_genome_nb_E                   /= (double)_nb_level_2_cells;
    _level_2_genome_nb_TF                  /= (double)_nb_level_2_cells;
    _level_2_genome_nb_BS                  /= (double)_nb_level_2_cells;
    _level_2_genome_nb_P                   /= (double)_nb_level_2_cells;
    _level_2_genome_nb_inner_enzymes       /= (double)_nb_level_2_cells;
    _level_2_genome_nb_inflow_pumps        /= (double)_nb_level_2_cells;
    _level_2_genome_nb_outflow_pumps       /= (double)_nb_level_2_cells;
    _level_2_genome_nb_functional_regions  /= (double)_nb_level_2_cells;
    _level_2_genome_nb_enhancers           /= (double)_nb_level_2_cells;
    _level_2_genome_nb_operators           /= (double)_nb_level_2_cells;
    _level_2_genome_nb_E_regions           /= (double)_nb_level_2_cells;
    _level_2_genome_nb_TF_regions          /= (double)_nb_level_2_cells;
    _level_2_genome_nb_mixed_regions       /= (double)_nb_level_2_cells;
    _level_2_genome_functional_region_size /= (double)_nb_level_2_cells;
    _level_2_genome_E_region_size          /= (double)_nb_level_2_cells;
    _level_2_genome_TF_region_size         /= (double)_nb_level_2_cells;
    _level_2_genome_mixed_region_size      /= (double)_nb_level_2_cells;
    _level_2_genome_enhancer_size          /= (double)_nb_level_2_cells;
    _level_2_genome_operator_size          /= (double)_nb_level_2_cells;
    _level_2_genome_operon_size            /= (double)_nb_level_2_cells;
    _level_2_genome_E_operon_size          /= (double)_nb_level_2_cells;
    _level_2_genome_TF_operon_size         /= (double)_nb_level_2_cells;
    _level_2_genome_mixed_operon_size      /= (double)_nb_level_2_cells;
    
    /* INHERITED STRUCTURE */
    _level_2_inherited_size             /= (double)_nb_level_2_cells;
    _level_2_inherited_nb_E             /= (double)_nb_level_2_cells;
    _level_2_inherited_nb_TF            /= (double)_nb_level_2_cells;
    _level_2_inherited_nb_inner_enzymes /= (double)_nb_level_2_cells;
    _level_2_inherited_nb_inflow_pumps  /= (double)_nb_level_2_cells;
    _level_2_inherited_nb_outflow_pumps /= (double)_nb_level_2_cells;
  }

  if (_nb_no_level_cells > 0)
  {
    /* GENETIC DIVERSITY */
    _no_level_true_diversity              /= (double)_nb_no_level_cells;
    _no_level_max_time_distance           /= (double)_nb_no_level_cells;
    _no_level_max_generation_distance     /= (double)_nb_no_level_cells;
    _no_level_max_point_mutation_distance /= (double)_nb_no_level_cells;
    _no_level_max_hgt_distance            /= (double)_nb_no_level_cells;
    _no_level_max_duplication_distance    /= (double)_nb_no_level_cells;
    _no_level_max_deletion_distance       /= (double)_nb_no_level_cells;
    _no_level_max_inversion_distance      /= (double)_nb_no_level_cells;
    _no_level_max_translocation_distance  /= (double)_nb_no_level_cells;
    
    /* PHENOTYPE */
    _no_level_generations                /= (double)_nb_no_level_cells;
    _no_level_inherited_TF_amount        /= (double)_nb_no_level_cells;
    _no_level_inherited_E_amount         /= (double)_nb_no_level_cells;
    _no_level_TF_amount                  /= (double)_nb_no_level_cells;
    _no_level_E_amount                   /= (double)_nb_no_level_cells;
    _no_level_inherited_metabolic_amount /= (double)_nb_no_level_cells;
    _no_level_metabolic_amount           /= (double)_nb_no_level_cells;
    _no_level_energy                     /= (double)_nb_no_level_cells;
    _no_level_score                      /= (double)_nb_no_level_cells;
    _no_level_lifespan                   /= (double)_nb_no_level_cells;
    _no_level_number_of_divisions        /= (double)_nb_no_level_cells;
    _no_level_toxicity                   /= (double)_nb_no_level_cells;
    _no_level_metabolic_uptake           /= (double)_nb_no_level_cells;
    _no_level_metabolic_release          /= (double)_nb_no_level_cells;
    _no_level_metabolic_growth_rate      /= (double)_nb_no_level_cells;
    _no_level_Dmetabolic_growth_rate     /= (double)_nb_no_level_cells;
    _no_level_grn_nb_nodes               /= (double)_nb_no_level_cells;
    _no_level_grn_nb_edges               /= (double)_nb_no_level_cells;
    _no_level_metabolic_nb_nodes         /= (double)_nb_no_level_cells;
    _no_level_metabolic_nb_edges         /= (double)_nb_no_level_cells;
    _no_level_regulation_redundancy      /= (double)_nb_no_level_cells;
    _no_level_metabolic_redundancy       /= (double)_nb_no_level_cells;
    
    /* GENOME STRUCTURE */
    _no_level_genome_size                   /= (double)_nb_no_level_cells;
    _no_level_functional_size               /= (double)_nb_no_level_cells;
    _no_level_genome_nb_NC                  /= (double)_nb_no_level_cells;
    _no_level_genome_nb_E                   /= (double)_nb_no_level_cells;
    _no_level_genome_nb_TF                  /= (double)_nb_no_level_cells;
    _no_level_genome_nb_BS                  /= (double)_nb_no_level_cells;
    _no_level_genome_nb_P                   /= (double)_nb_no_level_cells;
    _no_level_genome_nb_inner_enzymes       /= (double)_nb_no_level_cells;
    _no_level_genome_nb_inflow_pumps        /= (double)_nb_no_level_cells;
    _no_level_genome_nb_outflow_pumps       /= (double)_nb_no_level_cells;
    _no_level_genome_nb_functional_regions  /= (double)_nb_no_level_cells;
    _no_level_genome_nb_enhancers           /= (double)_nb_no_level_cells;
    _no_level_genome_nb_operators           /= (double)_nb_no_level_cells;
    _no_level_genome_nb_E_regions           /= (double)_nb_no_level_cells;
    _no_level_genome_nb_TF_regions          /= (double)_nb_no_level_cells;
    _no_level_genome_nb_mixed_regions       /= (double)_nb_no_level_cells;
    _no_level_genome_functional_region_size /= (double)_nb_no_level_cells;
    _no_level_genome_E_region_size          /= (double)_nb_no_level_cells;
    _no_level_genome_TF_region_size         /= (double)_nb_no_level_cells;
    _no_level_genome_mixed_region_size      /= (double)_nb_no_level_cells;
    _no_level_genome_enhancer_size          /= (double)_nb_no_level_cells;
    _no_level_genome_operator_size          /= (double)_nb_no_level_cells;
    _no_level_genome_operon_size            /= (double)_nb_no_level_cells;
    _no_level_genome_E_operon_size          /= (double)_nb_no_level_cells;
    _no_level_genome_TF_operon_size         /= (double)_nb_no_level_cells;
    _no_level_genome_mixed_operon_size      /= (double)_nb_no_level_cells;
    
    /* INHERITED STRUCTURE */
    _no_level_inherited_size             /= (double)_nb_no_level_cells;
    _no_level_inherited_nb_E             /= (double)_nb_no_level_cells;
    _no_level_inherited_nb_TF            /= (double)_nb_no_level_cells;
    _no_level_inherited_nb_inner_enzymes /= (double)_nb_no_level_cells;
    _no_level_inherited_nb_inflow_pumps  /= (double)_nb_no_level_cells;
    _no_level_inherited_nb_outflow_pumps /= (double)_nb_no_level_cells;
  }
}

/**
 * \brief    Write the trophic network in files
 * \details  --
 * \param    std::string node_filename
 * \param    std::string edge_filename
 * \return   \e void
 */
void TrophicNetwork::write_trophic_network( std::string node_filename, std::string edge_filename )
{
  /*----------------------------*/
  /* 1) Write the list of nodes */
  /*----------------------------*/
  std::ofstream nodes_file(node_filename.c_str(), std::ios::out | std::ios::trunc);
  nodes_file << "id production uptake release level count appearance lifespan\n";
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    TrophicGroup* current_group = _iterator->second;
    nodes_file << current_group->get_identifier() << " " << current_group->get_production_profile() << " " << current_group->get_uptake_profile() << " " << current_group->get_release_profile() << " " << current_group->get_trophic_level() << " " << current_group->get_number_of_cells() << " " << current_group->get_appearance_time() << " " << current_group->get_lifespan() << "\n";
  }
  nodes_file.close();
  
  /*----------------------------*/
  /* 2) Write the list of edges */
  /*----------------------------*/
  std::ofstream edges_file(edge_filename.c_str(), std::ios::out | std::ios::trunc);
  edges_file << "id1 id2 link_type\n";
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    TrophicGroup* current_group = _iterator->second;
    for (size_t j = 0; j < current_group->get_necrophagy_links()->size(); j++)
    {
      edges_file << current_group->get_identifier() << " " << current_group->get_necrophagy_links()->at(j) << " " << "0\n";
    }
    for (size_t j = 0; j < current_group->get_active_release_links()->size(); j++)
    {
      edges_file << current_group->get_identifier() << " " << current_group->get_active_release_links()->at(j) << " " << "1\n";
    }
  }
  edges_file.close();
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Check if the trophic profile is already in the trophic network
 * \details  If the group already exists, its identifier is copied in group_id
 * \param    std::string trophic_profile
 * \return   \e bool
 */
bool TrophicNetwork::groupExists( std::string trophic_profile, unsigned long long int &group_id )
{
  for (_iterator = _group_map.begin(); _iterator != _group_map.end(); ++_iterator)
  {
    if (_iterator->second->get_trophic_profile().compare(trophic_profile) == 0 && _iterator->first != 0)
    {
      group_id = _iterator->second->get_identifier();
      return true;
    }
  }
  return false;
}
