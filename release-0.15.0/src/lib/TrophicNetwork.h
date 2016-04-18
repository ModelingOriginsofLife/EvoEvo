
/**
 * \file      TrophicNetwork.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      14-10-2015
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     TrophicNetwork class declaration
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

#ifndef __EVOEVO__TrophicNetwork__
#define __EVOEVO__TrophicNetwork__

#include <iostream>
#include <unordered_map>
#include <zlib.h>
#include <assert.h>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"
#include "Parameters.h"
#include "TrophicGroup.h"
#include "Population.h"
#include "Environment.h"


class TrophicNetwork
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  TrophicNetwork( void ) = delete;
  TrophicNetwork( Parameters* parameters, Population* population, Environment* environment );
  TrophicNetwork( Parameters* parameters, Population* population, Environment* environment, gzFile backup_file );
  TrophicNetwork( const TrophicNetwork& trophic_network ) = delete;
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~TrophicNetwork( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*------------------------------------------------------------------ Trophic network attributes */
  
  inline unsigned long long int get_current_id( void ) const;
  inline size_t                 get_number_of_groups( void ) const;
  inline TrophicGroup*          get_group( unsigned long long int identifier );
  inline TrophicGroup*          get_first_group( void );
  inline TrophicGroup*          get_next_group( void );
  
  /*------------------------------------------------------------------ Trophic network statistics */
  
  inline size_t get_nb_level_0_groups( void ) const;
  inline size_t get_nb_level_1_groups( void ) const;
  inline size_t get_nb_level_2_groups( void ) const;
  inline size_t get_nb_no_level_groups( void ) const;
  
  inline size_t get_nb_level_0_cells( void ) const;
  inline size_t get_nb_level_1_cells( void ) const;
  inline size_t get_nb_level_2_cells( void ) const;
  inline size_t get_nb_no_level_cells( void ) const;
  
  inline size_t get_nb_group_appearances( void ) const;
  inline size_t get_nb_group_extinctions( void ) const;
  
  inline double get_mean_group_lifespan( void ) const;
  
  inline double get_min_true_diversity( void ) const;
  inline double get_min_max_time_distance( void ) const;
  inline double get_min_max_generation_distance( void ) const;
  inline double get_min_max_point_mutation_distance( void ) const;
  inline double get_min_max_hgt_distance( void ) const;
  inline double get_min_max_duplication_distance( void ) const;
  inline double get_min_max_deletion_distance( void ) const;
  inline double get_min_max_inversion_distance( void ) const;
  inline double get_min_max_translocation_distance( void ) const;
  
  inline double get_mean_true_diversity( void ) const;
  inline double get_mean_max_time_distance( void ) const;
  inline double get_mean_max_generation_distance( void ) const;
  inline double get_mean_max_point_mutation_distance( void ) const;
  inline double get_mean_max_hgt_distance( void ) const;
  inline double get_mean_max_duplication_distance( void ) const;
  inline double get_mean_max_deletion_distance( void ) const;
  inline double get_mean_max_inversion_distance( void ) const;
  inline double get_mean_max_translocation_distance( void ) const;
  
  inline double get_max_true_diversity( void ) const;
  inline double get_max_max_time_distance( void ) const;
  inline double get_max_max_generation_distance( void ) const;
  inline double get_max_max_point_mutation_distance( void ) const;
  inline double get_max_max_hgt_distance( void ) const;
  inline double get_max_max_duplication_distance( void ) const;
  inline double get_max_max_deletion_distance( void ) const;
  inline double get_max_max_inversion_distance( void ) const;
  inline double get_max_max_translocation_distance( void ) const;
  
  /*------------------------------------------------------------------ Identifier management */
  
  inline unsigned long long int get_new_id( void );
  
  /*------------------------------------------------------------------ LEVEL 0 statistics */
  
  /* GENETIC DIVERSITY */
  inline double get_level_0_true_diversity( void ) const;
  inline double get_level_0_max_time_distance( void ) const;
  inline double get_level_0_max_generation_distance( void ) const;
  inline double get_level_0_max_point_mutation_distance( void ) const;
  inline double get_level_0_max_hgt_distance( void ) const;
  inline double get_level_0_max_duplication_distance( void ) const;
  inline double get_level_0_max_deletion_distance( void ) const;
  inline double get_level_0_max_inversion_distance( void ) const;
  inline double get_level_0_max_translocation_distance( void ) const;
  
  /* PHENOTYPE */
  inline double get_level_0_generations( void ) const;
  inline double get_level_0_inherited_TF_amount( void ) const;
  inline double get_level_0_inherited_E_amount( void ) const;
  inline double get_level_0_TF_amount( void ) const;
  inline double get_level_0_E_amount( void ) const;
  inline double get_level_0_inherited_metabolic_amount( void ) const;
  inline double get_level_0_metabolic_amount( void ) const;
  inline double get_level_0_energy( void ) const;
  inline double get_level_0_score( void ) const;
  inline double get_level_0_lifespan( void ) const;
  inline double get_level_0_number_of_divisions( void ) const;
  inline double get_level_0_toxicity( void ) const;
  inline double get_level_0_metabolic_uptake( void ) const;
  inline double get_level_0_metabolic_release( void ) const;
  inline double get_level_0_metabolic_growth_rate( void ) const;
  inline double get_level_0_Dmetabolic_growth_rate( void ) const;
  inline double get_level_0_grn_nb_nodes( void ) const;
  inline double get_level_0_grn_nb_edges( void ) const;
  inline double get_level_0_metabolic_nb_nodes( void ) const;
  inline double get_level_0_metabolic_nb_edges( void ) const;
  inline double get_level_0_regulation_redundancy( void ) const;
  inline double get_level_0_metabolic_redundancy( void ) const;
  
  /* GENOME STRUCTURE */
  inline double get_level_0_genome_size( void ) const;
  inline double get_level_0_functional_size( void ) const;
  inline double get_level_0_genome_nb_NC( void ) const;
  inline double get_level_0_genome_nb_E( void ) const;
  inline double get_level_0_genome_nb_TF( void ) const;
  inline double get_level_0_genome_nb_BS( void ) const;
  inline double get_level_0_genome_nb_P( void ) const;
  inline double get_level_0_genome_nb_inner_enzymes( void ) const;
  inline double get_level_0_genome_nb_inflow_pumps( void ) const;
  inline double get_level_0_genome_nb_outflow_pumps( void ) const;
  inline double get_level_0_genome_nb_functional_regions( void ) const;
  inline double get_level_0_genome_nb_enhancers( void ) const;
  inline double get_level_0_genome_nb_operators( void ) const;
  inline double get_level_0_genome_nb_E_regions( void ) const;
  inline double get_level_0_genome_nb_TF_regions( void ) const;
  inline double get_level_0_genome_nb_mixed_regions( void ) const;
  inline double get_level_0_genome_functional_region_size( void ) const;
  inline double get_level_0_genome_E_region_size( void ) const;
  inline double get_level_0_genome_TF_region_size( void ) const;
  inline double get_level_0_genome_mixed_region_size( void ) const;
  inline double get_level_0_genome_enhancer_size( void ) const;
  inline double get_level_0_genome_operator_size( void ) const;
  inline double get_level_0_genome_operon_size( void ) const;
  inline double get_level_0_genome_E_operon_size( void ) const;
  inline double get_level_0_genome_TF_operon_size( void ) const;
  inline double get_level_0_genome_mixed_operon_size( void ) const;
  
  /* INHERITED STRUCTURE */
  inline double get_level_0_inherited_size( void ) const;
  inline double get_level_0_inherited_nb_E( void ) const;
  inline double get_level_0_inherited_nb_TF( void ) const;
  inline double get_level_0_inherited_nb_inner_enzymes( void ) const;
  inline double get_level_0_inherited_nb_inflow_pumps( void ) const;
  inline double get_level_0_inherited_nb_outflow_pumps( void ) const;
  
  /*------------------------------------------------------------------ LEVEL 1 statistics */
  
  /* GENETIC DIVERSITY */
  inline double get_level_1_true_diversity( void ) const;
  inline double get_level_1_max_time_distance( void ) const;
  inline double get_level_1_max_generation_distance( void ) const;
  inline double get_level_1_max_point_mutation_distance( void ) const;
  inline double get_level_1_max_hgt_distance( void ) const;
  inline double get_level_1_max_duplication_distance( void ) const;
  inline double get_level_1_max_deletion_distance( void ) const;
  inline double get_level_1_max_inversion_distance( void ) const;
  inline double get_level_1_max_translocation_distance( void ) const;
  
  /* PHENOTYPE */
  inline double get_level_1_generations( void ) const;
  inline double get_level_1_inherited_TF_amount( void ) const;
  inline double get_level_1_inherited_E_amount( void ) const;
  inline double get_level_1_TF_amount( void ) const;
  inline double get_level_1_E_amount( void ) const;
  inline double get_level_1_inherited_metabolic_amount( void ) const;
  inline double get_level_1_metabolic_amount( void ) const;
  inline double get_level_1_energy( void ) const;
  inline double get_level_1_score( void ) const;
  inline double get_level_1_lifespan( void ) const;
  inline double get_level_1_number_of_divisions( void ) const;
  inline double get_level_1_toxicity( void ) const;
  inline double get_level_1_metabolic_uptake( void ) const;
  inline double get_level_1_metabolic_release( void ) const;
  inline double get_level_1_metabolic_growth_rate( void ) const;
  inline double get_level_1_Dmetabolic_growth_rate( void ) const;
  inline double get_level_1_grn_nb_nodes( void ) const;
  inline double get_level_1_grn_nb_edges( void ) const;
  inline double get_level_1_metabolic_nb_nodes( void ) const;
  inline double get_level_1_metabolic_nb_edges( void ) const;
  inline double get_level_1_regulation_redundancy( void ) const;
  inline double get_level_1_metabolic_redundancy( void ) const;
  
  /* GENOME STRUCTURE */
  inline double get_level_1_genome_size( void ) const;
  inline double get_level_1_functional_size( void ) const;
  inline double get_level_1_genome_nb_NC( void ) const;
  inline double get_level_1_genome_nb_E( void ) const;
  inline double get_level_1_genome_nb_TF( void ) const;
  inline double get_level_1_genome_nb_BS( void ) const;
  inline double get_level_1_genome_nb_P( void ) const;
  inline double get_level_1_genome_nb_inner_enzymes( void ) const;
  inline double get_level_1_genome_nb_inflow_pumps( void ) const;
  inline double get_level_1_genome_nb_outflow_pumps( void ) const;
  inline double get_level_1_genome_nb_functional_regions( void ) const;
  inline double get_level_1_genome_nb_enhancers( void ) const;
  inline double get_level_1_genome_nb_operators( void ) const;
  inline double get_level_1_genome_nb_E_regions( void ) const;
  inline double get_level_1_genome_nb_TF_regions( void ) const;
  inline double get_level_1_genome_nb_mixed_regions( void ) const;
  inline double get_level_1_genome_functional_region_size( void ) const;
  inline double get_level_1_genome_E_region_size( void ) const;
  inline double get_level_1_genome_TF_region_size( void ) const;
  inline double get_level_1_genome_mixed_region_size( void ) const;
  inline double get_level_1_genome_enhancer_size( void ) const;
  inline double get_level_1_genome_operator_size( void ) const;
  inline double get_level_1_genome_operon_size( void ) const;
  inline double get_level_1_genome_E_operon_size( void ) const;
  inline double get_level_1_genome_TF_operon_size( void ) const;
  inline double get_level_1_genome_mixed_operon_size( void ) const;
  
  /* INHERITED STRUCTURE */
  inline double get_level_1_inherited_size( void ) const;
  inline double get_level_1_inherited_nb_E( void ) const;
  inline double get_level_1_inherited_nb_TF( void ) const;
  inline double get_level_1_inherited_nb_inner_enzymes( void ) const;
  inline double get_level_1_inherited_nb_inflow_pumps( void ) const;
  inline double get_level_1_inherited_nb_outflow_pumps( void ) const;
  
  /*------------------------------------------------------------------ LEVEL 2 statistics */
  
  /* GENETIC DIVERSITY */
  inline double get_level_2_true_diversity( void ) const;
  inline double get_level_2_max_time_distance( void ) const;
  inline double get_level_2_max_generation_distance( void ) const;
  inline double get_level_2_max_point_mutation_distance( void ) const;
  inline double get_level_2_max_hgt_distance( void ) const;
  inline double get_level_2_max_duplication_distance( void ) const;
  inline double get_level_2_max_deletion_distance( void ) const;
  inline double get_level_2_max_inversion_distance( void ) const;
  inline double get_level_2_max_translocation_distance( void ) const;
  
  /* PHENOTYPE */
  inline double get_level_2_generations( void ) const;
  inline double get_level_2_inherited_TF_amount( void ) const;
  inline double get_level_2_inherited_E_amount( void ) const;
  inline double get_level_2_TF_amount( void ) const;
  inline double get_level_2_E_amount( void ) const;
  inline double get_level_2_inherited_metabolic_amount( void ) const;
  inline double get_level_2_metabolic_amount( void ) const;
  inline double get_level_2_energy( void ) const;
  inline double get_level_2_score( void ) const;
  inline double get_level_2_lifespan( void ) const;
  inline double get_level_2_number_of_divisions( void ) const;
  inline double get_level_2_toxicity( void ) const;
  inline double get_level_2_metabolic_uptake( void ) const;
  inline double get_level_2_metabolic_release( void ) const;
  inline double get_level_2_metabolic_growth_rate( void ) const;
  inline double get_level_2_Dmetabolic_growth_rate( void ) const;
  inline double get_level_2_grn_nb_nodes( void ) const;
  inline double get_level_2_grn_nb_edges( void ) const;
  inline double get_level_2_metabolic_nb_nodes( void ) const;
  inline double get_level_2_metabolic_nb_edges( void ) const;
  inline double get_level_2_regulation_redundancy( void ) const;
  inline double get_level_2_metabolic_redundancy( void ) const;
  
  /* GENOME STRUCTURE */
  inline double get_level_2_genome_size( void ) const;
  inline double get_level_2_functional_size( void ) const;
  inline double get_level_2_genome_nb_NC( void ) const;
  inline double get_level_2_genome_nb_E( void ) const;
  inline double get_level_2_genome_nb_TF( void ) const;
  inline double get_level_2_genome_nb_BS( void ) const;
  inline double get_level_2_genome_nb_P( void ) const;
  inline double get_level_2_genome_nb_inner_enzymes( void ) const;
  inline double get_level_2_genome_nb_inflow_pumps( void ) const;
  inline double get_level_2_genome_nb_outflow_pumps( void ) const;
  inline double get_level_2_genome_nb_functional_regions( void ) const;
  inline double get_level_2_genome_nb_enhancers( void ) const;
  inline double get_level_2_genome_nb_operators( void ) const;
  inline double get_level_2_genome_nb_E_regions( void ) const;
  inline double get_level_2_genome_nb_TF_regions( void ) const;
  inline double get_level_2_genome_nb_mixed_regions( void ) const;
  inline double get_level_2_genome_functional_region_size( void ) const;
  inline double get_level_2_genome_E_region_size( void ) const;
  inline double get_level_2_genome_TF_region_size( void ) const;
  inline double get_level_2_genome_mixed_region_size( void ) const;
  inline double get_level_2_genome_enhancer_size( void ) const;
  inline double get_level_2_genome_operator_size( void ) const;
  inline double get_level_2_genome_operon_size( void ) const;
  inline double get_level_2_genome_E_operon_size( void ) const;
  inline double get_level_2_genome_TF_operon_size( void ) const;
  inline double get_level_2_genome_mixed_operon_size( void ) const;
  
  /* INHERITED STRUCTURE */
  inline double get_level_2_inherited_size( void ) const;
  inline double get_level_2_inherited_nb_E( void ) const;
  inline double get_level_2_inherited_nb_TF( void ) const;
  inline double get_level_2_inherited_nb_inner_enzymes( void ) const;
  inline double get_level_2_inherited_nb_inflow_pumps( void ) const;
  inline double get_level_2_inherited_nb_outflow_pumps( void ) const;
  
  /*------------------------------------------------------------------ NO LEVEL statistics */
  
  /* GENETIC DIVERSITY */
  inline double get_no_level_true_diversity( void ) const;
  inline double get_no_level_max_time_distance( void ) const;
  inline double get_no_level_max_generation_distance( void ) const;
  inline double get_no_level_max_point_mutation_distance( void ) const;
  inline double get_no_level_max_hgt_distance( void ) const;
  inline double get_no_level_max_duplication_distance( void ) const;
  inline double get_no_level_max_deletion_distance( void ) const;
  inline double get_no_level_max_inversion_distance( void ) const;
  inline double get_no_level_max_translocation_distance( void ) const;
  
  /* PHENOTYPE */
  inline double get_no_level_generations( void ) const;
  inline double get_no_level_inherited_TF_amount( void ) const;
  inline double get_no_level_inherited_E_amount( void ) const;
  inline double get_no_level_TF_amount( void ) const;
  inline double get_no_level_E_amount( void ) const;
  inline double get_no_level_inherited_metabolic_amount( void ) const;
  inline double get_no_level_metabolic_amount( void ) const;
  inline double get_no_level_energy( void ) const;
  inline double get_no_level_score( void ) const;
  inline double get_no_level_lifespan( void ) const;
  inline double get_no_level_number_of_divisions( void ) const;
  inline double get_no_level_toxicity( void ) const;
  inline double get_no_level_metabolic_uptake( void ) const;
  inline double get_no_level_metabolic_release( void ) const;
  inline double get_no_level_metabolic_growth_rate( void ) const;
  inline double get_no_level_Dmetabolic_growth_rate( void ) const;
  inline double get_no_level_grn_nb_nodes( void ) const;
  inline double get_no_level_grn_nb_edges( void ) const;
  inline double get_no_level_metabolic_nb_nodes( void ) const;
  inline double get_no_level_metabolic_nb_edges( void ) const;
  inline double get_no_level_regulation_redundancy( void ) const;
  inline double get_no_level_metabolic_redundancy( void ) const;
  
  /* GENOME STRUCTURE */
  inline double get_no_level_genome_size( void ) const;
  inline double get_no_level_functional_size( void ) const;
  inline double get_no_level_genome_nb_NC( void ) const;
  inline double get_no_level_genome_nb_E( void ) const;
  inline double get_no_level_genome_nb_TF( void ) const;
  inline double get_no_level_genome_nb_BS( void ) const;
  inline double get_no_level_genome_nb_P( void ) const;
  inline double get_no_level_genome_nb_inner_enzymes( void ) const;
  inline double get_no_level_genome_nb_inflow_pumps( void ) const;
  inline double get_no_level_genome_nb_outflow_pumps( void ) const;
  inline double get_no_level_genome_nb_functional_regions( void ) const;
  inline double get_no_level_genome_nb_enhancers( void ) const;
  inline double get_no_level_genome_nb_operators( void ) const;
  inline double get_no_level_genome_nb_E_regions( void ) const;
  inline double get_no_level_genome_nb_TF_regions( void ) const;
  inline double get_no_level_genome_nb_mixed_regions( void ) const;
  inline double get_no_level_genome_functional_region_size( void ) const;
  inline double get_no_level_genome_E_region_size( void ) const;
  inline double get_no_level_genome_TF_region_size( void ) const;
  inline double get_no_level_genome_mixed_region_size( void ) const;
  inline double get_no_level_genome_enhancer_size( void ) const;
  inline double get_no_level_genome_operator_size( void ) const;
  inline double get_no_level_genome_operon_size( void ) const;
  inline double get_no_level_genome_E_operon_size( void ) const;
  inline double get_no_level_genome_TF_operon_size( void ) const;
  inline double get_no_level_genome_mixed_operon_size( void ) const;
  
  /* INHERITED STRUCTURE */
  inline double get_no_level_inherited_size( void ) const;
  inline double get_no_level_inherited_nb_E( void ) const;
  inline double get_no_level_inherited_nb_TF( void ) const;
  inline double get_no_level_inherited_nb_inner_enzymes( void ) const;
  inline double get_no_level_inherited_nb_inflow_pumps( void ) const;
  inline double get_no_level_inherited_nb_outflow_pumps( void ) const;
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void save( gzFile backup_file );
  void initialize_trophic_network( void );
  void load_population( void );
  void compute_diversity_statistics( void );
  void compute_level_statistics( void );
  void write_trophic_network( std::string node_filename, std::string edge_filename );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  bool groupExists( std::string trophic_profile, unsigned long long int &group_id );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ Trophic network attributes */
  
  Parameters*                                                         _parameters;  /*!< Simulation parameters      */
  Population*                                                         _population;  /*!< Population                 */
  Environment*                                                        _environment; /*!< Environment                */
  unsigned long long int                                              _current_id;  /*!< Current group id           */
  std::unordered_map<unsigned long long int, TrophicGroup*>           _group_map;   /*!< Trophic group map          */
  std::unordered_map<unsigned long long int, TrophicGroup*>::iterator _iterator;    /*!< Trophic group map iterator */
  
  /*------------------------------------------------------------------ Trophic network statistics */
  
  size_t _nb_level_0_groups;  /*!< Number of groups of LEVEL 0 level  */
  size_t _nb_level_1_groups;  /*!< Number of groups of LEVEL 1 level  */
  size_t _nb_level_2_groups;  /*!< Number of groups of LEVEL 2 level  */
  size_t _nb_no_level_groups; /*!< Number of groups of NO LEVEL level */
  
  size_t _nb_level_0_cells;  /*!< Number of cells of LEVEL 0 level  */
  size_t _nb_level_1_cells;  /*!< Number of cells of LEVEL 1 level  */
  size_t _nb_level_2_cells;  /*!< Number of cells of LEVEL 2 level  */
  size_t _nb_no_level_cells; /*!< Number of cells of NO LEVEL level */
  
  size_t _nb_group_appearances; /*!< Number of group appearances per simulation timestep */
  size_t _nb_group_extinctions; /*!< Number of group extinctions per simulation timestep */
  
  double _mean_group_lifespan; /*!< Mean group lifespan */
  
  double _min_true_diversity;              /*!< Minimum true diversity measure          */
  double _min_max_time_distance;           /*!< Minimum maximum time distance           */
  double _min_max_generation_distance;     /*!< Minimum maximum generation distance     */
  double _min_max_point_mutation_distance; /*!< Minimum maximum point mutation distance */
  double _min_max_hgt_distance;            /*!< Minimum maximum HGT distance            */
  double _min_max_duplication_distance;    /*!< Minimum maximum duplication distance    */
  double _min_max_deletion_distance;       /*!< Minimum maximum deletion distance       */
  double _min_max_inversion_distance;      /*!< Minimum maximum inversion distance      */
  double _min_max_translocation_distance;  /*!< Minimum maximum translocation distance  */
  
  double _mean_true_diversity;              /*!< Mean true diversity measure          */
  double _mean_max_time_distance;           /*!< Mean maximum time distance           */
  double _mean_max_generation_distance;     /*!< Mean maximum generation distance     */
  double _mean_max_point_mutation_distance; /*!< Mean maximum point mutation distance */
  double _mean_max_hgt_distance;            /*!< Mean maximum HGT distance            */
  double _mean_max_duplication_distance;    /*!< Mean maximum duplication distance    */
  double _mean_max_deletion_distance;       /*!< Mean maximum deletion distance       */
  double _mean_max_inversion_distance;      /*!< Mean maximum inversion distance      */
  double _mean_max_translocation_distance;  /*!< Mean maximum translocation distance  */
  
  double _max_true_diversity;              /*!< Maximum true diversity measure          */
  double _max_max_time_distance;           /*!< Maximum maximum time distance           */
  double _max_max_generation_distance;     /*!< Maximum maximum generation distance     */
  double _max_max_point_mutation_distance; /*!< Maximum maximum point mutation distance */
  double _max_max_hgt_distance;            /*!< Maximum maximum HGT distance            */
  double _max_max_duplication_distance;    /*!< Maximum maximum duplication distance    */
  double _max_max_deletion_distance;       /*!< Maximum maximum deletion distance       */
  double _max_max_inversion_distance;      /*!< Maximum maximum inversion distance      */
  double _max_max_translocation_distance;  /*!< Maximum maximum translocation distance  */
  
  /*------------------------------------------------------------------ LEVEL 0 statistics */
  
  /* GENETIC DIVERSITY */
  double _level_0_true_diversity;              /*!< True diversity measure          */
  double _level_0_max_time_distance;           /*!< Maximum time distance           */
  double _level_0_max_generation_distance;     /*!< Maximum generation distance     */
  double _level_0_max_point_mutation_distance; /*!< Maximum point mutation distance */
  double _level_0_max_hgt_distance;            /*!< Maximum HGT distance            */
  double _level_0_max_duplication_distance;    /*!< Maximum duplication distance    */
  double _level_0_max_deletion_distance;       /*!< Maximum deletion distance       */
  double _level_0_max_inversion_distance;      /*!< Maximum inversion distance      */
  double _level_0_max_translocation_distance;  /*!< Maximum translocation distance  */
  
  /* PHENOTYPE */
  double _level_0_generations;                /*!< Number of generations                    */
  double _level_0_inherited_TF_amount;        /*!< Inherited TF amount                      */
  double _level_0_inherited_E_amount;         /*!< Inherited E amount                       */
  double _level_0_TF_amount;                  /*!< TF amount                                */
  double _level_0_E_amount;                   /*!< E amount                                 */
  double _level_0_inherited_metabolic_amount; /*!< Inherited metabolic amount               */
  double _level_0_metabolic_amount;           /*!< Metabolic amount                         */
  double _level_0_energy;                     /*!< Energy amount                            */
  double _level_0_score;                      /*!< Score                                    */
  double _level_0_lifespan;                   /*!< Lifespan                                 */
  double _level_0_number_of_divisions;        /*!< Number of divisions                      */
  double _level_0_toxicity;                   /*!< Toxicity accumulation                    */
  double _level_0_metabolic_uptake;           /*!< Amount of metabolic uptake               */
  double _level_0_metabolic_release;          /*!< Amount of metabolic release              */
  double _level_0_metabolic_growth_rate;      /*!< Metabolic growth rate                    */
  double _level_0_Dmetabolic_growth_rate;     /*!< Metabolic growth rate difference         */
  double _level_0_grn_nb_nodes;               /*!< Number of nodes in the GRN               */
  double _level_0_grn_nb_edges;               /*!< Number of edges in the GRN               */
  double _level_0_metabolic_nb_nodes;         /*!< Number of nodes in the metabolic network */
  double _level_0_metabolic_nb_edges;         /*!< Number of edges in the metabolic network */
  double _level_0_regulation_redundancy;      /*!< Regulation redundancy                    */
  double _level_0_metabolic_redundancy;       /*!< Metabolic redundancy                     */
  
  /* GENOME STRUCTURE */
  double _level_0_genome_size;                   /*!< Genome size                                                   */
  double _level_0_functional_size;               /*!< Functional size                                               */
  double _level_0_genome_nb_NC;                  /*!< Number of non coding type (NC) in the genome                  */
  double _level_0_genome_nb_E;                   /*!< Number of enzyme type (E) in the genome                       */
  double _level_0_genome_nb_TF;                  /*!< Number of transcription factor type (TF) in the genome        */
  double _level_0_genome_nb_BS;                  /*!< Number of binding site type (BS) in the genome                */
  double _level_0_genome_nb_P;                   /*!< Number of promoter type (P) in the genome                     */
  double _level_0_genome_nb_inner_enzymes;       /*!< Number of inner enzymes in the genome                         */
  double _level_0_genome_nb_inflow_pumps;        /*!< Number of inflowing pumps in the genome                       */
  double _level_0_genome_nb_outflow_pumps;       /*!< Number of outflowing pumps in the genome                      */
  double _level_0_genome_nb_functional_regions;  /*!< Number of functional regions in the genome                    */
  double _level_0_genome_nb_enhancers;           /*!< Number of enhancers in functional regions in the genome       */
  double _level_0_genome_nb_operators;           /*!< Number of operators in functional regions in the genome       */
  double _level_0_genome_nb_E_regions;           /*!< Number of functional regions containing only TF in the genome */
  double _level_0_genome_nb_TF_regions;          /*!< Number of functional regions containing only TF in the genome */
  double _level_0_genome_nb_mixed_regions;       /*!< Number of functional regions mixing E and TF in the genome    */
  double _level_0_genome_functional_region_size; /*!< Size of functional regions in the genome                      */
  double _level_0_genome_E_region_size;          /*!< Size of E regions in the genome                               */
  double _level_0_genome_TF_region_size;         /*!< Size of TF regions in the genome                              */
  double _level_0_genome_mixed_region_size;      /*!< Size of mixed regions in the genome                           */
  double _level_0_genome_enhancer_size;          /*!< Size of enhancer sites in the genome                          */
  double _level_0_genome_operator_size;          /*!< Size of operator sites in the genome                          */
  double _level_0_genome_operon_size;            /*!< Size of operons in the genome                                 */
  double _level_0_genome_E_operon_size;          /*!< Size of operons containing only E in the genome               */
  double _level_0_genome_TF_operon_size;         /*!< Size of operons containing only TF in the genome              */
  double _level_0_genome_mixed_operon_size;      /*!< Size of operons containing mixing E and TF in the genome      */
  
  /* INHERITED STRUCTURE */
  double _level_0_inherited_size;             /*!< Number of inherited proteins                                   */
  double _level_0_inherited_nb_E;             /*!< Number of enzyme type (E) in inherited proteins                */
  double _level_0_inherited_nb_TF;            /*!< Number of transcription factor type (TF) in inherited proteins */
  double _level_0_inherited_nb_inner_enzymes; /*!< Number of inner enzymes in inherited proteins                  */
  double _level_0_inherited_nb_inflow_pumps;  /*!< Number of inflowing pumps in inherited proteins                */
  double _level_0_inherited_nb_outflow_pumps; /*!< Number of outflowing pumps in inherited proteins               */
  
  /*------------------------------------------------------------------ LEVEL 1 statistics */
  
  /* GENETIC DIVERSITY */
  double _level_1_true_diversity;              /*!< True diversity measure          */
  double _level_1_max_time_distance;           /*!< Maximum time distance           */
  double _level_1_max_generation_distance;     /*!< Maximum generation distance     */
  double _level_1_max_point_mutation_distance; /*!< Maximum point mutation distance */
  double _level_1_max_hgt_distance;            /*!< Maximum HGT distance            */
  double _level_1_max_duplication_distance;    /*!< Maximum duplication distance    */
  double _level_1_max_deletion_distance;       /*!< Maximum deletion distance       */
  double _level_1_max_inversion_distance;      /*!< Maximum inversion distance      */
  double _level_1_max_translocation_distance;  /*!< Maximum translocation distance  */
  
  /* PHENOTYPE */
  double _level_1_generations;                /*!< Number of generations                    */
  double _level_1_inherited_TF_amount;        /*!< Inherited TF amount                      */
  double _level_1_inherited_E_amount;         /*!< Inherited E amount                       */
  double _level_1_TF_amount;                  /*!< TF amount                                */
  double _level_1_E_amount;                   /*!< E amount                                 */
  double _level_1_inherited_metabolic_amount; /*!< Inherited metabolic amount               */
  double _level_1_metabolic_amount;           /*!< Metabolic amount                         */
  double _level_1_energy;                     /*!< Energy amount                            */
  double _level_1_score;                      /*!< Score                                    */
  double _level_1_lifespan;                   /*!< Lifespan                                 */
  double _level_1_number_of_divisions;        /*!< Number of divisions                      */
  double _level_1_toxicity;                   /*!< Toxicity accumulation                    */
  double _level_1_metabolic_uptake;           /*!< Amount of metabolic uptake               */
  double _level_1_metabolic_release;          /*!< Amount of metabolic release              */
  double _level_1_metabolic_growth_rate;      /*!< Metabolic growth rate                    */
  double _level_1_Dmetabolic_growth_rate;     /*!< Metabolic growth rate difference         */
  double _level_1_grn_nb_nodes;               /*!< Number of nodes in the GRN               */
  double _level_1_grn_nb_edges;               /*!< Number of edges in the GRN               */
  double _level_1_metabolic_nb_nodes;         /*!< Number of nodes in the metabolic network */
  double _level_1_metabolic_nb_edges;         /*!< Number of edges in the metabolic network */
  double _level_1_regulation_redundancy;      /*!< Regulation redundancy                    */
  double _level_1_metabolic_redundancy;       /*!< Metabolic redundancy                     */
  
  /* GENOME STRUCTURE */
  double _level_1_genome_size;                   /*!< Genome size                                                   */
  double _level_1_functional_size;               /*!< Functional size                                               */
  double _level_1_genome_nb_NC;                  /*!< Number of non coding type (NC) in the genome                  */
  double _level_1_genome_nb_E;                   /*!< Number of enzyme type (E) in the genome                       */
  double _level_1_genome_nb_TF;                  /*!< Number of transcription factor type (TF) in the genome        */
  double _level_1_genome_nb_BS;                  /*!< Number of binding site type (BS) in the genome                */
  double _level_1_genome_nb_P;                   /*!< Number of promoter type (P) in the genome                     */
  double _level_1_genome_nb_inner_enzymes;       /*!< Number of inner enzymes in the genome                         */
  double _level_1_genome_nb_inflow_pumps;        /*!< Number of inflowing pumps in the genome                       */
  double _level_1_genome_nb_outflow_pumps;       /*!< Number of outflowing pumps in the genome                      */
  double _level_1_genome_nb_functional_regions;  /*!< Number of functional regions in the genome                    */
  double _level_1_genome_nb_enhancers;           /*!< Number of enhancers in functional regions in the genome       */
  double _level_1_genome_nb_operators;           /*!< Number of operators in functional regions in the genome       */
  double _level_1_genome_nb_E_regions;           /*!< Number of functional regions containing only TF in the genome */
  double _level_1_genome_nb_TF_regions;          /*!< Number of functional regions containing only TF in the genome */
  double _level_1_genome_nb_mixed_regions;       /*!< Number of functional regions mixing E and TF in the genome    */
  double _level_1_genome_functional_region_size; /*!< Size of functional regions in the genome                      */
  double _level_1_genome_E_region_size;          /*!< Size of E regions in the genome                               */
  double _level_1_genome_TF_region_size;         /*!< Size of TF regions in the genome                              */
  double _level_1_genome_mixed_region_size;      /*!< Size of mixed regions in the genome                           */
  double _level_1_genome_enhancer_size;          /*!< Size of enhancer sites in the genome                          */
  double _level_1_genome_operator_size;          /*!< Size of operator sites in the genome                          */
  double _level_1_genome_operon_size;            /*!< Size of operons in the genome                                 */
  double _level_1_genome_E_operon_size;          /*!< Size of operons containing only E in the genome               */
  double _level_1_genome_TF_operon_size;         /*!< Size of operons containing only TF in the genome              */
  double _level_1_genome_mixed_operon_size;      /*!< Size of operons containing mixing E and TF in the genome      */
  
  /* INHERITED STRUCTURE */
  double _level_1_inherited_size;             /*!< Number of inherited proteins                                   */
  double _level_1_inherited_nb_E;             /*!< Number of enzyme type (E) in inherited proteins                */
  double _level_1_inherited_nb_TF;            /*!< Number of transcription factor type (TF) in inherited proteins */
  double _level_1_inherited_nb_inner_enzymes; /*!< Number of inner enzymes in inherited proteins                  */
  double _level_1_inherited_nb_inflow_pumps;  /*!< Number of inflowing pumps in inherited proteins                */
  double _level_1_inherited_nb_outflow_pumps; /*!< Number of outflowing pumps in inherited proteins               */
  
  /*------------------------------------------------------------------ LEVEL 2 statistics */
  
  /* GENETIC DIVERSITY */
  double _level_2_true_diversity;              /*!< True diversity measure          */
  double _level_2_max_time_distance;           /*!< Maximum time distance           */
  double _level_2_max_generation_distance;     /*!< Maximum generation distance     */
  double _level_2_max_point_mutation_distance; /*!< Maximum point mutation distance */
  double _level_2_max_hgt_distance;            /*!< Maximum HGT distance            */
  double _level_2_max_duplication_distance;    /*!< Maximum duplication distance    */
  double _level_2_max_deletion_distance;       /*!< Maximum deletion distance       */
  double _level_2_max_inversion_distance;      /*!< Maximum inversion distance      */
  double _level_2_max_translocation_distance;  /*!< Maximum translocation distance  */
  
  /* PHENOTYPE */
  double _level_2_generations;                /*!< Number of generations                    */
  double _level_2_inherited_TF_amount;        /*!< Inherited TF amount                      */
  double _level_2_inherited_E_amount;         /*!< Inherited E amount                       */
  double _level_2_TF_amount;                  /*!< TF amount                                */
  double _level_2_E_amount;                   /*!< E amount                                 */
  double _level_2_inherited_metabolic_amount; /*!< Inherited metabolic amount               */
  double _level_2_metabolic_amount;           /*!< Metabolic amount                         */
  double _level_2_energy;                     /*!< Energy amount                            */
  double _level_2_score;                      /*!< Score                                    */
  double _level_2_lifespan;                   /*!< Lifespan                                 */
  double _level_2_number_of_divisions;        /*!< Number of divisions                      */
  double _level_2_toxicity;                   /*!< Toxicity accumulation                    */
  double _level_2_metabolic_uptake;           /*!< Amount of metabolic uptake               */
  double _level_2_metabolic_release;          /*!< Amount of metabolic release              */
  double _level_2_metabolic_growth_rate;      /*!< Metabolic growth rate                    */
  double _level_2_Dmetabolic_growth_rate;     /*!< Metabolic growth rate difference         */
  double _level_2_grn_nb_nodes;               /*!< Number of nodes in the GRN               */
  double _level_2_grn_nb_edges;               /*!< Number of edges in the GRN               */
  double _level_2_metabolic_nb_nodes;         /*!< Number of nodes in the metabolic network */
  double _level_2_metabolic_nb_edges;         /*!< Number of edges in the metabolic network */
  double _level_2_regulation_redundancy;      /*!< Regulation redundancy                    */
  double _level_2_metabolic_redundancy;       /*!< Metabolic redundancy                     */
  
  /* GENOME STRUCTURE */
  double _level_2_genome_size;                   /*!< Genome size                                                   */
  double _level_2_functional_size;               /*!< Functional size                                               */
  double _level_2_genome_nb_NC;                  /*!< Number of non coding type (NC) in the genome                  */
  double _level_2_genome_nb_E;                   /*!< Number of enzyme type (E) in the genome                       */
  double _level_2_genome_nb_TF;                  /*!< Number of transcription factor type (TF) in the genome        */
  double _level_2_genome_nb_BS;                  /*!< Number of binding site type (BS) in the genome                */
  double _level_2_genome_nb_P;                   /*!< Number of promoter type (P) in the genome                     */
  double _level_2_genome_nb_inner_enzymes;       /*!< Number of inner enzymes in the genome                         */
  double _level_2_genome_nb_inflow_pumps;        /*!< Number of inflowing pumps in the genome                       */
  double _level_2_genome_nb_outflow_pumps;       /*!< Number of outflowing pumps in the genome                      */
  double _level_2_genome_nb_functional_regions;  /*!< Number of functional regions in the genome                    */
  double _level_2_genome_nb_enhancers;           /*!< Number of enhancers in functional regions in the genome       */
  double _level_2_genome_nb_operators;           /*!< Number of operators in functional regions in the genome       */
  double _level_2_genome_nb_E_regions;           /*!< Number of functional regions containing only TF in the genome */
  double _level_2_genome_nb_TF_regions;          /*!< Number of functional regions containing only TF in the genome */
  double _level_2_genome_nb_mixed_regions;       /*!< Number of functional regions mixing E and TF in the genome    */
  double _level_2_genome_functional_region_size; /*!< Size of functional regions in the genome                      */
  double _level_2_genome_E_region_size;          /*!< Size of E regions in the genome                               */
  double _level_2_genome_TF_region_size;         /*!< Size of TF regions in the genome                              */
  double _level_2_genome_mixed_region_size;      /*!< Size of mixed regions in the genome                           */
  double _level_2_genome_enhancer_size;          /*!< Size of enhancer sites in the genome                          */
  double _level_2_genome_operator_size;          /*!< Size of operator sites in the genome                          */
  double _level_2_genome_operon_size;            /*!< Size of operons in the genome                                 */
  double _level_2_genome_E_operon_size;          /*!< Size of operons containing only E in the genome               */
  double _level_2_genome_TF_operon_size;         /*!< Size of operons containing only TF in the genome              */
  double _level_2_genome_mixed_operon_size;      /*!< Size of operons containing mixing E and TF in the genome      */
  
  /* INHERITED STRUCTURE */
  double _level_2_inherited_size;             /*!< Number of inherited proteins                                   */
  double _level_2_inherited_nb_E;             /*!< Number of enzyme type (E) in inherited proteins                */
  double _level_2_inherited_nb_TF;            /*!< Number of transcription factor type (TF) in inherited proteins */
  double _level_2_inherited_nb_inner_enzymes; /*!< Number of inner enzymes in inherited proteins                  */
  double _level_2_inherited_nb_inflow_pumps;  /*!< Number of inflowing pumps in inherited proteins                */
  double _level_2_inherited_nb_outflow_pumps; /*!< Number of outflowing pumps in inherited proteins               */
  
  /*------------------------------------------------------------------ NO LEVEL statistics */
  
  /* GENETIC DIVERSITY */
  double _no_level_true_diversity;              /*!< True diversity measure          */
  double _no_level_max_time_distance;           /*!< Maximum time distance           */
  double _no_level_max_generation_distance;     /*!< Maximum generation distance     */
  double _no_level_max_point_mutation_distance; /*!< Maximum point mutation distance */
  double _no_level_max_hgt_distance;            /*!< Maximum HGT distance            */
  double _no_level_max_duplication_distance;    /*!< Maximum duplication distance    */
  double _no_level_max_deletion_distance;       /*!< Maximum deletion distance       */
  double _no_level_max_inversion_distance;      /*!< Maximum inversion distance      */
  double _no_level_max_translocation_distance;  /*!< Maximum translocation distance  */
  
  /* PHENOTYPE */
  double _no_level_generations;                /*!< Number of generations                    */
  double _no_level_inherited_TF_amount;        /*!< Inherited TF amount                      */
  double _no_level_inherited_E_amount;         /*!< Inherited E amount                       */
  double _no_level_TF_amount;                  /*!< TF amount                                */
  double _no_level_E_amount;                   /*!< E amount                                 */
  double _no_level_inherited_metabolic_amount; /*!< Inherited metabolic amount               */
  double _no_level_metabolic_amount;           /*!< Metabolic amount                         */
  double _no_level_energy;                     /*!< Energy amount                            */
  double _no_level_score;                      /*!< Score                                    */
  double _no_level_lifespan;                   /*!< Lifespan                                 */
  double _no_level_number_of_divisions;        /*!< Number of divisions                      */
  double _no_level_toxicity;                   /*!< Toxicity accumulation                    */
  double _no_level_metabolic_uptake;           /*!< Amount of metabolic uptake               */
  double _no_level_metabolic_release;          /*!< Amount of metabolic release              */
  double _no_level_metabolic_growth_rate;      /*!< Metabolic growth rate                    */
  double _no_level_Dmetabolic_growth_rate;     /*!< Metabolic growth rate difference         */
  double _no_level_grn_nb_nodes;               /*!< Number of nodes in the GRN               */
  double _no_level_grn_nb_edges;               /*!< Number of edges in the GRN               */
  double _no_level_metabolic_nb_nodes;         /*!< Number of nodes in the metabolic network */
  double _no_level_metabolic_nb_edges;         /*!< Number of edges in the metabolic network */
  double _no_level_regulation_redundancy;      /*!< Regulation redundancy                    */
  double _no_level_metabolic_redundancy;       /*!< Metabolic redundancy                     */
  
  /* GENOME STRUCTURE */
  double _no_level_genome_size;                   /*!< Genome size                                                   */
  double _no_level_functional_size;               /*!< Functional size                                               */
  double _no_level_genome_nb_NC;                  /*!< Number of non coding type (NC) in the genome                  */
  double _no_level_genome_nb_E;                   /*!< Number of enzyme type (E) in the genome                       */
  double _no_level_genome_nb_TF;                  /*!< Number of transcription factor type (TF) in the genome        */
  double _no_level_genome_nb_BS;                  /*!< Number of binding site type (BS) in the genome                */
  double _no_level_genome_nb_P;                   /*!< Number of promoter type (P) in the genome                     */
  double _no_level_genome_nb_inner_enzymes;       /*!< Number of inner enzymes in the genome                         */
  double _no_level_genome_nb_inflow_pumps;        /*!< Number of inflowing pumps in the genome                       */
  double _no_level_genome_nb_outflow_pumps;       /*!< Number of outflowing pumps in the genome                      */
  double _no_level_genome_nb_functional_regions;  /*!< Number of functional regions in the genome                    */
  double _no_level_genome_nb_enhancers;           /*!< Number of enhancers in functional regions in the genome       */
  double _no_level_genome_nb_operators;           /*!< Number of operators in functional regions in the genome       */
  double _no_level_genome_nb_E_regions;           /*!< Number of functional regions containing only TF in the genome */
  double _no_level_genome_nb_TF_regions;          /*!< Number of functional regions containing only TF in the genome */
  double _no_level_genome_nb_mixed_regions;       /*!< Number of functional regions mixing E and TF in the genome    */
  double _no_level_genome_functional_region_size; /*!< Size of functional regions in the genome                      */
  double _no_level_genome_E_region_size;          /*!< Size of E regions in the genome                               */
  double _no_level_genome_TF_region_size;         /*!< Size of TF regions in the genome                              */
  double _no_level_genome_mixed_region_size;      /*!< Size of mixed regions in the genome                           */
  double _no_level_genome_enhancer_size;          /*!< Size of enhancer sites in the genome                          */
  double _no_level_genome_operator_size;          /*!< Size of operator sites in the genome                          */
  double _no_level_genome_operon_size;            /*!< Size of operons in the genome                                 */
  double _no_level_genome_E_operon_size;          /*!< Size of operons containing only E in the genome               */
  double _no_level_genome_TF_operon_size;         /*!< Size of operons containing only TF in the genome              */
  double _no_level_genome_mixed_operon_size;      /*!< Size of operons containing mixing E and TF in the genome      */
  
  /* INHERITED STRUCTURE */
  double _no_level_inherited_size;             /*!< Number of inherited proteins                                   */
  double _no_level_inherited_nb_E;             /*!< Number of enzyme type (E) in inherited proteins                */
  double _no_level_inherited_nb_TF;            /*!< Number of transcription factor type (TF) in inherited proteins */
  double _no_level_inherited_nb_inner_enzymes; /*!< Number of inner enzymes in inherited proteins                  */
  double _no_level_inherited_nb_inflow_pumps;  /*!< Number of inflowing pumps in inherited proteins                */
  double _no_level_inherited_nb_outflow_pumps; /*!< Number of outflowing pumps in inherited proteins               */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/
/**
 * \brief    Get current group identifier
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int TrophicNetwork::get_current_id( void ) const
{
  return _current_id;
}

/**
 * \brief    Get the number of groups of the network
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_number_of_groups( void ) const
{
  return _group_map.size();
}

/**
 * \brief    Get the node by its identifier
 * \details  Return NULL if the node do not exist
 * \param    unsigned long long int identifier
 * \return   \e TrophicGroup*
 */
inline TrophicGroup* TrophicNetwork::get_group( unsigned long long int identifier )
{
  if (_group_map.find(identifier) != _group_map.end())
  {
    return _group_map[identifier];
  }
  return NULL;
}

/**
 * \brief    Get the first group of the network
 * \details  Return NULL if the network is empty
 * \param    void
 * \return   \e TrophicGroup*
 */
inline TrophicGroup* TrophicNetwork::get_first_group( void )
{
  _iterator = _group_map.begin();
  if (_iterator != _group_map.end())
  {
    return _iterator->second;
  }
  return NULL;
}

/**
 * \brief    Get the next group
 * \details  Return NULL if the end of the network is reached
 * \param    void
 * \return   \e TrophicGroup*
 */
inline TrophicGroup* TrophicNetwork::get_next_group( void )
{
  _iterator++;
  if (_iterator != _group_map.end())
  {
    return _iterator->second;
  }
  return NULL;
}

/*------------------------------------------------------------------ Trophic network statistics */

/**
 * \brief    Get number of groups of level LEVEL 0
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_0_groups( void ) const
{
  return _nb_level_0_groups;
}

/**
 * \brief    Get number of groups of level LEVEL 1
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_1_groups( void ) const
{
  return _nb_level_1_groups;
}

/**
 * \brief    Get number of groups of level LEVEL 2
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_2_groups( void ) const
{
  return _nb_level_2_groups;
}

/**
 * \brief    Get number of groups of level NO LEVEL
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_no_level_groups( void ) const
{
  return _nb_no_level_groups;
}

/**
 * \brief    Get number of cells of level LEVEL 0
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_0_cells( void ) const
{
  return _nb_level_0_cells;
}

/**
 * \brief    Get number of cells of level LEVEL 1
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_1_cells( void ) const
{
  return _nb_level_1_cells;
}

/**
 * \brief    Get number of cells of level LEVEL 2
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_level_2_cells( void ) const
{
  return _nb_level_2_cells;
}

/**
 * \brief    Get number of cells of level NO LEVEL
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_no_level_cells( void ) const
{
  return _nb_no_level_cells;
}

/**
 * \brief    Get number of group appearances
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_group_appearances( void ) const
{
  return _nb_group_appearances;
}

/**
 * \brief    Get number of group extinctions
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t TrophicNetwork::get_nb_group_extinctions( void ) const
{
  return _nb_group_extinctions;
}

/**
 * \brief    Get mean group lifespan
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline double TrophicNetwork::get_mean_group_lifespan( void ) const
{
  return _mean_group_lifespan;
}

/**
 * \brief    Get minimum true diversity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_min_true_diversity( void ) const
{
  return _min_true_diversity;
}

/**
 * \brief    Get minimum maximum time distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_min_max_time_distance( void ) const
{
  return _min_max_time_distance;
}

/**
 * \brief    Get minimum maximum generation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_min_max_generation_distance( void ) const
{
  return _min_max_generation_distance;
}

/**
 * \brief    Get minimum maximum point mutation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_min_max_point_mutation_distance( void ) const
{
  return _min_max_point_mutation_distance;
}

/**
 * \brief    Get minimum maximum HGT distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_min_max_hgt_distance( void ) const
{
  return _min_max_hgt_distance;
}

/**
 * \brief    Get minimum maximum duplication distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_min_max_duplication_distance( void ) const
{
  return _min_max_duplication_distance;
}

/**
 * \brief    Get minimum maximum deletion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_min_max_deletion_distance( void ) const
{
  return _min_max_deletion_distance;
}

/**
 * \brief    Get minimum maximum inversion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_min_max_inversion_distance( void ) const
{
  return _min_max_inversion_distance;
}

/**
 * \brief    Get minimum maximum translocation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_min_max_translocation_distance( void ) const
{
  return _min_max_translocation_distance;
}

/**
 * \brief    Get mean true diversity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_mean_true_diversity( void ) const
{
  return _mean_true_diversity;
}

/**
 * \brief    Get mean maximum time distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_mean_max_time_distance( void ) const
{
  return _mean_max_time_distance;
}

/**
 * \brief    Get mean maximum generation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_mean_max_generation_distance( void ) const
{
  return _mean_max_generation_distance;
}

/**
 * \brief    Get mean maximum point mutation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_mean_max_point_mutation_distance( void ) const
{
  return _mean_max_point_mutation_distance;
}

/**
 * \brief    Get mean maximum HGT distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_mean_max_hgt_distance( void ) const
{
  return _mean_max_hgt_distance;
}

/**
 * \brief    Get mean maximum duplication distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_mean_max_duplication_distance( void ) const
{
  return _mean_max_duplication_distance;
}

/**
 * \brief    Get mean maximum deletion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_mean_max_deletion_distance( void ) const
{
  return _mean_max_deletion_distance;
}

/**
 * \brief    Get mean maximum inversion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_mean_max_inversion_distance( void ) const
{
  return _mean_max_inversion_distance;
}

/**
 * \brief    Get mean maximum translocation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_mean_max_translocation_distance( void ) const
{
  return _mean_max_translocation_distance;
}

/**
 * \brief    Get maximum true diversity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_max_true_diversity( void ) const
{
  return _max_true_diversity;
}

/**
 * \brief    Get maximum maximum time distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_max_max_time_distance( void ) const
{
  return _max_max_time_distance;
}

/**
 * \brief    Get maximum maximum generation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_max_max_generation_distance( void ) const
{
  return _max_max_generation_distance;
}

/**
 * \brief    Get maximum maximum point mutation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_max_max_point_mutation_distance( void ) const
{
  return _max_max_point_mutation_distance;
}

/**
 * \brief    Get maximum maximum HGT distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_max_max_hgt_distance( void ) const
{
  return _max_max_hgt_distance;
}

/**
 * \brief    Get maximum maximum duplication distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_max_max_duplication_distance( void ) const
{
  return _max_max_duplication_distance;
}

/**
 * \brief    Get maximum maximum deletion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_max_max_deletion_distance( void ) const
{
  return _max_max_deletion_distance;
}

/**
 * \brief    Get maximum maximum inversion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_max_max_inversion_distance( void ) const
{
  return _max_max_inversion_distance;
}

/**
 * \brief    Get maximum maximum translocation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_max_max_translocation_distance( void ) const
{
  return _max_max_translocation_distance;
}

/*------------------------------------------------------------------ Identifier management */

/**
 * \brief    Get a new identifier
 * \details  --
 * \param    void
 * \return   \e unsigned long long int
 */
inline unsigned long long int TrophicNetwork::get_new_id( void )
{
  unsigned long long int new_id = _current_id;
  _current_id++;
  return new_id;
}

/*------------------------------------------------------------------ LEVEL 0 statistics */

/* GENETIC DIVERSITY */

/**
 * \brief    Get level 0 mean true diversity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_true_diversity( void ) const
{
  return _level_0_true_diversity;
}

/**
 * \brief    Get level 0 mean maximum time distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_max_time_distance( void ) const
{
  return _level_0_max_time_distance;
}

/**
 * \brief    Get level 0 mean maximum generation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_max_generation_distance( void ) const
{
  return _level_0_max_generation_distance;
}

/**
 * \brief    Get level 0 mean maximum point mutation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_max_point_mutation_distance( void ) const
{
  return _level_0_max_point_mutation_distance;
}

/**
 * \brief    Get level 0 mean maximum HGT distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_max_hgt_distance( void ) const
{
  return _level_0_max_hgt_distance;
}

/**
 * \brief    Get level 0 mean maximum duplication distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_max_duplication_distance( void ) const
{
  return _level_0_max_duplication_distance;
}

/**
 * \brief    Get level 0 mean maximum deletion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_max_deletion_distance( void ) const
{
  return _level_0_max_deletion_distance;
}

/**
 * \brief    Get level 0 mean maximum inversion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_max_inversion_distance( void ) const
{
  return _level_0_max_inversion_distance;
}

/**
 * \brief    Get level 0 mean maximum translocation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_max_translocation_distance( void ) const
{
  return _level_0_max_translocation_distance;
}


/* PHENOTYPE */

/**
 * \brief    Get level 0 mean generations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_generations( void ) const
{
  return _level_0_generations;
}

/**
 * \brief    Get level 0 mean inherited TF amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_inherited_TF_amount( void ) const
{
  return _level_0_inherited_TF_amount;
}

/**
 * \brief    Get level 0 mean inherited E amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_inherited_E_amount( void ) const
{
  return _level_0_inherited_E_amount;
}

/**
 * \brief    Get level 0 mean TF amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_TF_amount( void ) const
{
  return _level_0_TF_amount;
}

/**
 * \brief    Get level 0 mean E amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_E_amount( void ) const
{
  return _level_0_E_amount;
}

/**
 * \brief    Get level 0 mean inherited metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_inherited_metabolic_amount( void ) const
{
  return _level_0_inherited_metabolic_amount;
}

/**
 * \brief    Get level 0 mean metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_metabolic_amount( void ) const
{
  return _level_0_metabolic_amount;
}

/**
 * \brief    Get level 0 mean energy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_energy( void ) const
{
  return _level_0_energy;
}

/**
 * \brief    Get level 0 mean score
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_score( void ) const
{
  return _level_0_score;
}

/**
 * \brief    Get level 0 mean lifespan
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_lifespan( void ) const
{
  return _level_0_lifespan;
}

/**
 * \brief    Get level 0 mean number of divisions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_number_of_divisions( void ) const
{
  return _level_0_number_of_divisions;
}

/**
 * \brief    Get level 0 mean toxicity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_toxicity( void ) const
{
  return _level_0_toxicity;
}

/**
 * \brief    Get level 0 mean metabolic uptake
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_metabolic_uptake( void ) const
{
  return _level_0_metabolic_uptake;
}

/**
 * \brief    Get level 0 mean metabolic release
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_metabolic_release( void ) const
{
  return _level_0_metabolic_release;
}

/**
 * \brief    Get level 0 mean metabolic growth rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_metabolic_growth_rate( void ) const
{
  return _level_0_metabolic_growth_rate;
}

/**
 * \brief    Get level 0 mean differential metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_Dmetabolic_growth_rate( void ) const
{
  return _level_0_Dmetabolic_growth_rate;
}

/**
 * \brief    Get level 0 mean GRN nodes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_grn_nb_nodes( void ) const
{
  return _level_0_grn_nb_nodes;
}

/**
 * \brief    Get level 0 mean GRN edges
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_grn_nb_edges( void ) const
{
  return _level_0_grn_nb_edges;
}

/**
 * \brief    Get level 0 mean metabolic network nodes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_metabolic_nb_nodes( void ) const
{
  return _level_0_metabolic_nb_nodes;
}

/**
 * \brief    Get level 0 mean metabolic network edges
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_metabolic_nb_edges( void ) const
{
  return _level_0_metabolic_nb_edges;
}

/**
 * \brief    Get level 0 mean regulation redundancy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_regulation_redundancy( void ) const
{
  return _level_0_regulation_redundancy;
}

/**
 * \brief    Get level 0 mean metabolic redundancy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_metabolic_redundancy( void ) const
{
  return _level_0_metabolic_redundancy;
}

/* GENOME STRUCTURE */

/**
 * \brief    Get level 0 genome mean genome size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_size( void ) const
{
  return _level_0_genome_size;
}

/**
 * \brief    Get level 0 genome mean functional size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_functional_size( void ) const
{
  return _level_0_functional_size;
}

/**
 * \brief    Get level 0 genome mean number of genome NC pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_NC( void ) const
{
  return _level_0_genome_nb_NC;
}

/**
 * \brief    Get level 0 genome mean number of genome E pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_E( void ) const
{
  return _level_0_genome_nb_E;
}

/**
 * \brief    Get level 0 genome mean number of genome TF pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_TF( void ) const
{
  return _level_0_genome_nb_TF;
}

/**
 * \brief    Get level 0 genome mean number of genome BS pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_BS( void ) const
{
  return _level_0_genome_nb_BS;
}

/**
 * \brief    Get level 0 genome mean number of genome P pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_P( void ) const
{
  return _level_0_genome_nb_P;
}

/**
 * \brief    Get level 0 genome mean number of inner enzymes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_inner_enzymes( void ) const
{
  return _level_0_genome_nb_inner_enzymes;
}

/**
 * \brief    Get level 0 genome mean number of inflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_inflow_pumps( void ) const
{
  return _level_0_genome_nb_inflow_pumps;
}

/**
 * \brief    Get level 0 genome mean number of outflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_outflow_pumps( void ) const
{
  return _level_0_genome_nb_outflow_pumps;
}

/**
 * \brief    Get level 0 genome mean number of functional regions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_functional_regions( void ) const
{
  return _level_0_genome_nb_functional_regions;
}

/**
 * \brief    Get level 0 genome mean number of enhancers
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_enhancers( void ) const
{
  return _level_0_genome_nb_enhancers;
}

/**
 * \brief    Get level 0 genome mean number of operators
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_operators( void ) const
{
  return _level_0_genome_nb_operators;
}

/**
 * \brief    Get level 0 genome mean E region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_E_regions( void ) const
{
  return _level_0_genome_nb_E_regions;
}

/**
 * \brief    Get level 0 genome mean TF region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_TF_regions( void ) const
{
  return _level_0_genome_nb_TF_regions;
}

/**
 * \brief    Get level 0 genome mean mixed region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_nb_mixed_regions( void ) const
{
  return _level_0_genome_nb_mixed_regions;
}

/**
 * \brief    Get level 0 genome mean TF region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_functional_region_size( void ) const
{
  return _level_0_genome_functional_region_size;
}

/**
 * \brief    Get level 0 genome mean E region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_E_region_size( void ) const
{
  return _level_0_genome_E_region_size;
}

/**
 * \brief    Get level 0 genome mean TF region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_TF_region_size( void ) const
{
  return _level_0_genome_TF_region_size;
}

/**
 * \brief    Get level 0 genome mean mixed region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_mixed_region_size( void ) const
{
  return _level_0_genome_mixed_region_size;
}

/**
 * \brief    Get level 0 genome mean enhancer size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_enhancer_size( void ) const
{
  return _level_0_genome_enhancer_size;
}

/**
 * \brief    Get level 0 genome mean operator size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_operator_size( void ) const
{
  return _level_0_genome_operator_size;
}

/**
 * \brief    Get level 0 genome mean operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_operon_size( void ) const
{
  return _level_0_genome_operon_size;
}

/**
 * \brief    Get level 0 genome mean E operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_E_operon_size( void ) const
{
  return _level_0_genome_E_operon_size;
}

/**
 * \brief    Get level 0 genome mean TF operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_TF_operon_size( void ) const
{
  return _level_0_genome_TF_operon_size;
}

/**
 * \brief    Get level 0 genome mean mixed operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_genome_mixed_operon_size( void ) const
{
  return _level_0_genome_mixed_operon_size;
}


/* INHERITED STRUCTURE */

/**
 * \brief    Get level 0 mean inherited size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_inherited_size( void ) const
{
  return _level_0_inherited_size;
}

/**
 * \brief    Get level 0 mean number of inherited E pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_inherited_nb_E( void ) const
{
  return _level_0_inherited_nb_E;
}

/**
 * \brief    Get level 0 mean number of inherited TF pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_inherited_nb_TF( void ) const
{
  return _level_0_inherited_nb_TF;
}

/**
 * \brief    Get level 0 mean number of inherited inner enzymes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_inherited_nb_inner_enzymes( void ) const
{
  return _level_0_inherited_nb_inner_enzymes;
}

/**
 * \brief    Get level 0 mean number of inherited inflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_inherited_nb_inflow_pumps( void ) const
{
  return _level_0_inherited_nb_inflow_pumps;
}

/**
 * \brief    Get level 0 mean number of inherited outflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_0_inherited_nb_outflow_pumps( void ) const
{
  return _level_0_inherited_nb_outflow_pumps;
}

/*------------------------------------------------------------------ LEVEL 1 statistics */

/* GENETIC DIVERSITY */

/**
 * \brief    Get level 1 mean true diversity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_true_diversity( void ) const
{
  return _level_1_true_diversity;
}

/**
 * \brief    Get level 1 mean maximum time distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_max_time_distance( void ) const
{
  return _level_1_max_time_distance;
}

/**
 * \brief    Get level 1 mean maximum generation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_max_generation_distance( void ) const
{
  return _level_1_max_generation_distance;
}

/**
 * \brief    Get level 1 mean maximum point mutation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_max_point_mutation_distance( void ) const
{
  return _level_1_max_point_mutation_distance;
}

/**
 * \brief    Get level 1 mean maximum HGT distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_max_hgt_distance( void ) const
{
  return _level_1_max_hgt_distance;
}

/**
 * \brief    Get level 1 mean maximum duplication distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_max_duplication_distance( void ) const
{
  return _level_1_max_duplication_distance;
}

/**
 * \brief    Get level 1 mean maximum deletion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_max_deletion_distance( void ) const
{
  return _level_1_max_deletion_distance;
}

/**
 * \brief    Get level 1 mean maximum inversion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_max_inversion_distance( void ) const
{
  return _level_1_max_inversion_distance;
}

/**
 * \brief    Get level 1 mean maximum translocation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_max_translocation_distance( void ) const
{
  return _level_1_max_translocation_distance;
}


/* PHENOTYPE */

/**
 * \brief    Get level 1 mean generations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_generations( void ) const
{
  return _level_1_generations;
}

/**
 * \brief    Get level 1 mean inherited TF amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_inherited_TF_amount( void ) const
{
  return _level_1_inherited_TF_amount;
}

/**
 * \brief    Get level 1 mean inherited E amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_inherited_E_amount( void ) const
{
  return _level_1_inherited_E_amount;
}

/**
 * \brief    Get level 1 mean TF amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_TF_amount( void ) const
{
  return _level_1_TF_amount;
}

/**
 * \brief    Get level 1 mean E amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_E_amount( void ) const
{
  return _level_1_E_amount;
}

/**
 * \brief    Get level 1 mean inherited metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_inherited_metabolic_amount( void ) const
{
  return _level_1_inherited_metabolic_amount;
}

/**
 * \brief    Get level 1 mean metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_metabolic_amount( void ) const
{
  return _level_1_metabolic_amount;
}

/**
 * \brief    Get level 1 mean energy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_energy( void ) const
{
  return _level_1_energy;
}

/**
 * \brief    Get level 1 mean score
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_score( void ) const
{
  return _level_1_score;
}

/**
 * \brief    Get level 1 mean lifespan
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_lifespan( void ) const
{
  return _level_1_lifespan;
}

/**
 * \brief    Get level 1 mean number of divisions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_number_of_divisions( void ) const
{
  return _level_1_number_of_divisions;
}

/**
 * \brief    Get level 1 mean toxicity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_toxicity( void ) const
{
  return _level_1_toxicity;
}

/**
 * \brief    Get level 1 mean metabolic uptake
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_metabolic_uptake( void ) const
{
  return _level_1_metabolic_uptake;
}

/**
 * \brief    Get level 1 mean metabolic release
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_metabolic_release( void ) const
{
  return _level_1_metabolic_release;
}

/**
 * \brief    Get level 1 mean metabolic growth rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_metabolic_growth_rate( void ) const
{
  return _level_1_metabolic_growth_rate;
}

/**
 * \brief    Get level 1 mean differential metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_Dmetabolic_growth_rate( void ) const
{
  return _level_1_Dmetabolic_growth_rate;
}

/**
 * \brief    Get level 1 mean GRN nodes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_grn_nb_nodes( void ) const
{
  return _level_1_grn_nb_nodes;
}

/**
 * \brief    Get level 1 mean GRN edges
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_grn_nb_edges( void ) const
{
  return _level_1_grn_nb_edges;
}

/**
 * \brief    Get level 1 mean metabolic network nodes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_metabolic_nb_nodes( void ) const
{
  return _level_1_metabolic_nb_nodes;
}

/**
 * \brief    Get level 1 mean metabolic network edges
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_metabolic_nb_edges( void ) const
{
  return _level_1_metabolic_nb_edges;
}

/**
 * \brief    Get level 1 mean regulation redundancy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_regulation_redundancy( void ) const
{
  return _level_1_regulation_redundancy;
}

/**
 * \brief    Get level 1 mean metabolic redundancy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_metabolic_redundancy( void ) const
{
  return _level_1_metabolic_redundancy;
}

/* GENOME STRUCTURE */

/**
 * \brief    Get level 1 genome mean genome size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_size( void ) const
{
  return _level_1_genome_size;
}

/**
 * \brief    Get level 1 genome mean functional size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_functional_size( void ) const
{
  return _level_1_functional_size;
}

/**
 * \brief    Get level 1 genome mean number of genome NC pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_NC( void ) const
{
  return _level_1_genome_nb_NC;
}

/**
 * \brief    Get level 1 genome mean number of genome E pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_E( void ) const
{
  return _level_1_genome_nb_E;
}

/**
 * \brief    Get level 1 genome mean number of genome TF pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_TF( void ) const
{
  return _level_1_genome_nb_TF;
}

/**
 * \brief    Get level 1 genome mean number of genome BS pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_BS( void ) const
{
  return _level_1_genome_nb_BS;
}

/**
 * \brief    Get level 1 genome mean number of genome P pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_P( void ) const
{
  return _level_1_genome_nb_P;
}

/**
 * \brief    Get level 1 genome mean number of inner enzymes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_inner_enzymes( void ) const
{
  return _level_1_genome_nb_inner_enzymes;
}

/**
 * \brief    Get level 1 genome mean number of inflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_inflow_pumps( void ) const
{
  return _level_1_genome_nb_inflow_pumps;
}

/**
 * \brief    Get level 1 genome mean number of outflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_outflow_pumps( void ) const
{
  return _level_1_genome_nb_outflow_pumps;
}

/**
 * \brief    Get level 1 genome mean number of functional regions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_functional_regions( void ) const
{
  return _level_1_genome_nb_functional_regions;
}

/**
 * \brief    Get level 1 genome mean number of enhancers
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_enhancers( void ) const
{
  return _level_1_genome_nb_enhancers;
}

/**
 * \brief    Get level 1 genome mean number of operators
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_operators( void ) const
{
  return _level_1_genome_nb_operators;
}

/**
 * \brief    Get level 1 genome mean E region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_E_regions( void ) const
{
  return _level_1_genome_nb_E_regions;
}

/**
 * \brief    Get level 1 genome mean TF region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_TF_regions( void ) const
{
  return _level_1_genome_nb_TF_regions;
}

/**
 * \brief    Get level 1 genome mean mixed region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_nb_mixed_regions( void ) const
{
  return _level_1_genome_nb_mixed_regions;
}

/**
 * \brief    Get level 1 genome mean TF region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_functional_region_size( void ) const
{
  return _level_1_genome_functional_region_size;
}

/**
 * \brief    Get level 1 genome mean E region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_E_region_size( void ) const
{
  return _level_1_genome_E_region_size;
}

/**
 * \brief    Get level 1 genome mean TF region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_TF_region_size( void ) const
{
  return _level_1_genome_TF_region_size;
}

/**
 * \brief    Get level 1 genome mean mixed region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_mixed_region_size( void ) const
{
  return _level_1_genome_mixed_region_size;
}

/**
 * \brief    Get level 1 genome mean enhancer size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_enhancer_size( void ) const
{
  return _level_1_genome_enhancer_size;
}

/**
 * \brief    Get level 1 genome mean operator size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_operator_size( void ) const
{
  return _level_1_genome_operator_size;
}

/**
 * \brief    Get level 1 genome mean operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_operon_size( void ) const
{
  return _level_1_genome_operon_size;
}

/**
 * \brief    Get level 1 genome mean E operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_E_operon_size( void ) const
{
  return _level_1_genome_E_operon_size;
}

/**
 * \brief    Get level 1 genome mean TF operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_TF_operon_size( void ) const
{
  return _level_1_genome_TF_operon_size;
}

/**
 * \brief    Get level 1 genome mean mixed operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_genome_mixed_operon_size( void ) const
{
  return _level_1_genome_mixed_operon_size;
}


/* INHERITED STRUCTURE */

/**
 * \brief    Get level 1 mean inherited size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_inherited_size( void ) const
{
  return _level_1_inherited_size;
}

/**
 * \brief    Get level 1 mean number of inherited E pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_inherited_nb_E( void ) const
{
  return _level_1_inherited_nb_E;
}

/**
 * \brief    Get level 1 mean number of inherited TF pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_inherited_nb_TF( void ) const
{
  return _level_1_inherited_nb_TF;
}

/**
 * \brief    Get level 1 mean number of inherited inner enzymes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_inherited_nb_inner_enzymes( void ) const
{
  return _level_1_inherited_nb_inner_enzymes;
}

/**
 * \brief    Get level 1 mean number of inherited inflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_inherited_nb_inflow_pumps( void ) const
{
  return _level_1_inherited_nb_inflow_pumps;
}

/**
 * \brief    Get level 1 mean number of inherited outflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_1_inherited_nb_outflow_pumps( void ) const
{
  return _level_1_inherited_nb_outflow_pumps;
}

/*------------------------------------------------------------------ LEVEL 2 statistics */

/* GENETIC DIVERSITY */

/**
 * \brief    Get level 2 mean true diversity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_true_diversity( void ) const
{
  return _level_2_true_diversity;
}

/**
 * \brief    Get level 2 mean maximum time distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_max_time_distance( void ) const
{
  return _level_2_max_time_distance;
}

/**
 * \brief    Get level 2 mean maximum generation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_max_generation_distance( void ) const
{
  return _level_2_max_generation_distance;
}

/**
 * \brief    Get level 2 mean maximum point mutation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_max_point_mutation_distance( void ) const
{
  return _level_2_max_point_mutation_distance;
}

/**
 * \brief    Get level 2 mean maximum HGT distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_max_hgt_distance( void ) const
{
  return _level_2_max_hgt_distance;
}

/**
 * \brief    Get level 2 mean maximum duplication distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_max_duplication_distance( void ) const
{
  return _level_2_max_duplication_distance;
}

/**
 * \brief    Get level 2 mean maximum deletion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_max_deletion_distance( void ) const
{
  return _level_2_max_deletion_distance;
}

/**
 * \brief    Get level 2 mean maximum inversion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_max_inversion_distance( void ) const
{
  return _level_2_max_inversion_distance;
}

/**
 * \brief    Get level 2 mean maximum translocation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_max_translocation_distance( void ) const
{
  return _level_2_max_translocation_distance;
}


/* PHENOTYPE */

/**
 * \brief    Get level 2 mean generations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_generations( void ) const
{
  return _level_2_generations;
}

/**
 * \brief    Get level 2 mean inherited TF amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_inherited_TF_amount( void ) const
{
  return _level_2_inherited_TF_amount;
}

/**
 * \brief    Get level 2 mean inherited E amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_inherited_E_amount( void ) const
{
  return _level_2_inherited_E_amount;
}

/**
 * \brief    Get level 2 mean TF amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_TF_amount( void ) const
{
  return _level_2_TF_amount;
}

/**
 * \brief    Get level 2 mean E amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_E_amount( void ) const
{
  return _level_2_E_amount;
}

/**
 * \brief    Get level 2 mean inherited metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_inherited_metabolic_amount( void ) const
{
  return _level_2_inherited_metabolic_amount;
}

/**
 * \brief    Get level 2 mean metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_metabolic_amount( void ) const
{
  return _level_2_metabolic_amount;
}

/**
 * \brief    Get level 2 mean energy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_energy( void ) const
{
  return _level_2_energy;
}

/**
 * \brief    Get level 2 mean score
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_score( void ) const
{
  return _level_2_score;
}

/**
 * \brief    Get level 2 mean lifespan
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_lifespan( void ) const
{
  return _level_2_lifespan;
}

/**
 * \brief    Get level 2 mean number of divisions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_number_of_divisions( void ) const
{
  return _level_2_number_of_divisions;
}

/**
 * \brief    Get level 2 mean toxicity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_toxicity( void ) const
{
  return _level_2_toxicity;
}

/**
 * \brief    Get level 2 mean metabolic uptake
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_metabolic_uptake( void ) const
{
  return _level_2_metabolic_uptake;
}

/**
 * \brief    Get level 2 mean metabolic release
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_metabolic_release( void ) const
{
  return _level_2_metabolic_release;
}

/**
 * \brief    Get level 2 mean metabolic growth rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_metabolic_growth_rate( void ) const
{
  return _level_2_metabolic_growth_rate;
}

/**
 * \brief    Get level 2 mean differential metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_Dmetabolic_growth_rate( void ) const
{
  return _level_2_Dmetabolic_growth_rate;
}

/**
 * \brief    Get level 2 mean GRN nodes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_grn_nb_nodes( void ) const
{
  return _level_2_grn_nb_nodes;
}

/**
 * \brief    Get level 2 mean GRN edges
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_grn_nb_edges( void ) const
{
  return _level_2_grn_nb_edges;
}

/**
 * \brief    Get level 2 mean metabolic network nodes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_metabolic_nb_nodes( void ) const
{
  return _level_2_metabolic_nb_nodes;
}

/**
 * \brief    Get level 2 mean metabolic network edges
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_metabolic_nb_edges( void ) const
{
  return _level_2_metabolic_nb_edges;
}

/**
 * \brief    Get level 2 mean regulation redundancy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_regulation_redundancy( void ) const
{
  return _level_2_regulation_redundancy;
}

/**
 * \brief    Get level 2 mean metabolic redundancy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_metabolic_redundancy( void ) const
{
  return _level_2_metabolic_redundancy;
}

/* GENOME STRUCTURE */

/**
 * \brief    Get level 2 genome mean genome size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_size( void ) const
{
  return _level_2_genome_size;
}

/**
 * \brief    Get level 2 genome mean functional size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_functional_size( void ) const
{
  return _level_2_functional_size;
}

/**
 * \brief    Get level 2 genome mean number of genome NC pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_NC( void ) const
{
  return _level_2_genome_nb_NC;
}

/**
 * \brief    Get level 2 genome mean number of genome E pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_E( void ) const
{
  return _level_2_genome_nb_E;
}

/**
 * \brief    Get level 2 genome mean number of genome TF pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_TF( void ) const
{
  return _level_2_genome_nb_TF;
}

/**
 * \brief    Get level 2 genome mean number of genome BS pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_BS( void ) const
{
  return _level_2_genome_nb_BS;
}

/**
 * \brief    Get level 2 genome mean number of genome P pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_P( void ) const
{
  return _level_2_genome_nb_P;
}

/**
 * \brief    Get level 2 genome mean number of inner enzymes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_inner_enzymes( void ) const
{
  return _level_2_genome_nb_inner_enzymes;
}

/**
 * \brief    Get level 2 genome mean number of inflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_inflow_pumps( void ) const
{
  return _level_2_genome_nb_inflow_pumps;
}

/**
 * \brief    Get level 2 genome mean number of outflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_outflow_pumps( void ) const
{
  return _level_2_genome_nb_outflow_pumps;
}

/**
 * \brief    Get level 2 genome mean number of functional regions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_functional_regions( void ) const
{
  return _level_2_genome_nb_functional_regions;
}

/**
 * \brief    Get level 2 genome mean number of enhancers
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_enhancers( void ) const
{
  return _level_2_genome_nb_enhancers;
}

/**
 * \brief    Get level 2 genome mean number of operators
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_operators( void ) const
{
  return _level_2_genome_nb_operators;
}

/**
 * \brief    Get level 2 genome mean E region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_E_regions( void ) const
{
  return _level_2_genome_nb_E_regions;
}

/**
 * \brief    Get level 2 genome mean TF region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_TF_regions( void ) const
{
  return _level_2_genome_nb_TF_regions;
}

/**
 * \brief    Get level 2 genome mean mixed region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_nb_mixed_regions( void ) const
{
  return _level_2_genome_nb_mixed_regions;
}

/**
 * \brief    Get level 2 genome mean TF region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_functional_region_size( void ) const
{
  return _level_2_genome_functional_region_size;
}

/**
 * \brief    Get level 2 genome mean E region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_E_region_size( void ) const
{
  return _level_2_genome_E_region_size;
}

/**
 * \brief    Get level 2 genome mean TF region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_TF_region_size( void ) const
{
  return _level_2_genome_TF_region_size;
}

/**
 * \brief    Get level 2 genome mean mixed region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_mixed_region_size( void ) const
{
  return _level_2_genome_mixed_region_size;
}

/**
 * \brief    Get level 2 genome mean enhancer size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_enhancer_size( void ) const
{
  return _level_2_genome_enhancer_size;
}

/**
 * \brief    Get level 2 genome mean operator size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_operator_size( void ) const
{
  return _level_2_genome_operator_size;
}

/**
 * \brief    Get level 2 genome mean operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_operon_size( void ) const
{
  return _level_2_genome_operon_size;
}

/**
 * \brief    Get level 2 genome mean E operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_E_operon_size( void ) const
{
  return _level_2_genome_E_operon_size;
}

/**
 * \brief    Get level 2 genome mean TF operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_TF_operon_size( void ) const
{
  return _level_2_genome_TF_operon_size;
}

/**
 * \brief    Get level 2 genome mean mixed operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_genome_mixed_operon_size( void ) const
{
  return _level_2_genome_mixed_operon_size;
}


/* INHERITED STRUCTURE */

/**
 * \brief    Get level 2 mean inherited size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_inherited_size( void ) const
{
  return _level_2_inherited_size;
}

/**
 * \brief    Get level 2 mean number of inherited E pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_inherited_nb_E( void ) const
{
  return _level_2_inherited_nb_E;
}

/**
 * \brief    Get level 2 mean number of inherited TF pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_inherited_nb_TF( void ) const
{
  return _level_2_inherited_nb_TF;
}

/**
 * \brief    Get level 2 mean number of inherited inner enzymes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_inherited_nb_inner_enzymes( void ) const
{
  return _level_2_inherited_nb_inner_enzymes;
}

/**
 * \brief    Get level 2 mean number of inherited inflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_inherited_nb_inflow_pumps( void ) const
{
  return _level_2_inherited_nb_inflow_pumps;
}

/**
 * \brief    Get level 2 mean number of inherited outflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_level_2_inherited_nb_outflow_pumps( void ) const
{
  return _level_2_inherited_nb_outflow_pumps;
}

/*------------------------------------------------------------------ NO LEVEL statistics */

/* GENETIC DIVERSITY */

/**
 * \brief    Get no level mean true diversity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_true_diversity( void ) const
{
  return _no_level_true_diversity;
}

/**
 * \brief    Get no level mean maximum time distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_max_time_distance( void ) const
{
  return _no_level_max_time_distance;
}

/**
 * \brief    Get no level mean maximum generation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_max_generation_distance( void ) const
{
  return _no_level_max_generation_distance;
}

/**
 * \brief    Get no level mean maximum point mutation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_max_point_mutation_distance( void ) const
{
  return _no_level_max_point_mutation_distance;
}

/**
 * \brief    Get no level mean maximum HGT distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_max_hgt_distance( void ) const
{
  return _no_level_max_hgt_distance;
}

/**
 * \brief    Get no level mean maximum duplication distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_max_duplication_distance( void ) const
{
  return _no_level_max_duplication_distance;
}

/**
 * \brief    Get no level mean maximum deletion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_max_deletion_distance( void ) const
{
  return _no_level_max_deletion_distance;
}

/**
 * \brief    Get no level mean maximum inversion distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_max_inversion_distance( void ) const
{
  return _no_level_max_inversion_distance;
}

/**
 * \brief    Get no level mean maximum translocation distance
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_max_translocation_distance( void ) const
{
  return _no_level_max_translocation_distance;
}


/* PHENOTYPE */

/**
 * \brief    Get no level mean generations
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_generations( void ) const
{
  return _no_level_generations;
}

/**
 * \brief    Get no level mean inherited TF amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_inherited_TF_amount( void ) const
{
  return _no_level_inherited_TF_amount;
}

/**
 * \brief    Get no level mean inherited E amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_inherited_E_amount( void ) const
{
  return _no_level_inherited_E_amount;
}

/**
 * \brief    Get no level mean TF amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_TF_amount( void ) const
{
  return _no_level_TF_amount;
}

/**
 * \brief    Get no level mean E amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_E_amount( void ) const
{
  return _no_level_E_amount;
}

/**
 * \brief    Get no level mean inherited metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_inherited_metabolic_amount( void ) const
{
  return _no_level_inherited_metabolic_amount;
}

/**
 * \brief    Get no level mean metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_metabolic_amount( void ) const
{
  return _no_level_metabolic_amount;
}

/**
 * \brief    Get no level mean energy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_energy( void ) const
{
  return _no_level_energy;
}

/**
 * \brief    Get no level mean score
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_score( void ) const
{
  return _no_level_score;
}

/**
 * \brief    Get no level mean lifespan
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_lifespan( void ) const
{
  return _no_level_lifespan;
}

/**
 * \brief    Get no level mean number of divisions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_number_of_divisions( void ) const
{
  return _no_level_number_of_divisions;
}

/**
 * \brief    Get no level mean toxicity
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_toxicity( void ) const
{
  return _no_level_toxicity;
}

/**
 * \brief    Get no level mean metabolic uptake
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_metabolic_uptake( void ) const
{
  return _no_level_metabolic_uptake;
}

/**
 * \brief    Get no level mean metabolic release
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_metabolic_release( void ) const
{
  return _no_level_metabolic_release;
}

/**
 * \brief    Get no level mean metabolic growth rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_metabolic_growth_rate( void ) const
{
  return _no_level_metabolic_growth_rate;
}

/**
 * \brief    Get no level mean differential metabolic amount
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_Dmetabolic_growth_rate( void ) const
{
  return _no_level_Dmetabolic_growth_rate;
}

/**
 * \brief    Get no level mean GRN nodes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_grn_nb_nodes( void ) const
{
  return _no_level_grn_nb_nodes;
}

/**
 * \brief    Get no level mean GRN edges
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_grn_nb_edges( void ) const
{
  return _no_level_grn_nb_edges;
}

/**
 * \brief    Get no level mean metabolic network nodes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_metabolic_nb_nodes( void ) const
{
  return _no_level_metabolic_nb_nodes;
}

/**
 * \brief    Get no level mean metabolic network edges
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_metabolic_nb_edges( void ) const
{
  return _no_level_metabolic_nb_edges;
}

/**
 * \brief    Get no level mean regulation redundancy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_regulation_redundancy( void ) const
{
  return _no_level_regulation_redundancy;
}

/**
 * \brief    Get no level mean metabolic redundancy
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_metabolic_redundancy( void ) const
{
  return _no_level_metabolic_redundancy;
}

/* GENOME STRUCTURE */

/**
 * \brief    Get no level genome mean genome size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_size( void ) const
{
  return _no_level_genome_size;
}

/**
 * \brief    Get no level genome mean functional size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_functional_size( void ) const
{
  return _no_level_functional_size;
}

/**
 * \brief    Get no level genome mean number of genome NC pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_NC( void ) const
{
  return _no_level_genome_nb_NC;
}

/**
 * \brief    Get no level genome mean number of genome E pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_E( void ) const
{
  return _no_level_genome_nb_E;
}

/**
 * \brief    Get no level genome mean number of genome TF pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_TF( void ) const
{
  return _no_level_genome_nb_TF;
}

/**
 * \brief    Get no level genome mean number of genome BS pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_BS( void ) const
{
  return _no_level_genome_nb_BS;
}

/**
 * \brief    Get no level genome mean number of genome P pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_P( void ) const
{
  return _no_level_genome_nb_P;
}

/**
 * \brief    Get no level genome mean number of inner enzymes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_inner_enzymes( void ) const
{
  return _no_level_genome_nb_inner_enzymes;
}

/**
 * \brief    Get no level genome mean number of inflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_inflow_pumps( void ) const
{
  return _no_level_genome_nb_inflow_pumps;
}

/**
 * \brief    Get no level genome mean number of outflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_outflow_pumps( void ) const
{
  return _no_level_genome_nb_outflow_pumps;
}

/**
 * \brief    Get no level genome mean number of functional regions
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_functional_regions( void ) const
{
  return _no_level_genome_nb_functional_regions;
}

/**
 * \brief    Get no level genome mean number of enhancers
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_enhancers( void ) const
{
  return _no_level_genome_nb_enhancers;
}

/**
 * \brief    Get no level genome mean number of operators
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_operators( void ) const
{
  return _no_level_genome_nb_operators;
}

/**
 * \brief    Get no level genome mean E region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_E_regions( void ) const
{
  return _no_level_genome_nb_E_regions;
}

/**
 * \brief    Get no level genome mean TF region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_TF_regions( void ) const
{
  return _no_level_genome_nb_TF_regions;
}

/**
 * \brief    Get no level genome mean mixed region number
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_nb_mixed_regions( void ) const
{
  return _no_level_genome_nb_mixed_regions;
}

/**
 * \brief    Get no level genome mean TF region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_functional_region_size( void ) const
{
  return _no_level_genome_functional_region_size;
}

/**
 * \brief    Get no level genome mean E region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_E_region_size( void ) const
{
  return _no_level_genome_E_region_size;
}

/**
 * \brief    Get no level genome mean TF region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_TF_region_size( void ) const
{
  return _no_level_genome_TF_region_size;
}

/**
 * \brief    Get no level genome mean mixed region size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_mixed_region_size( void ) const
{
  return _no_level_genome_mixed_region_size;
}

/**
 * \brief    Get no level genome mean enhancer size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_enhancer_size( void ) const
{
  return _no_level_genome_enhancer_size;
}

/**
 * \brief    Get no level genome mean operator size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_operator_size( void ) const
{
  return _no_level_genome_operator_size;
}

/**
 * \brief    Get no level genome mean operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_operon_size( void ) const
{
  return _no_level_genome_operon_size;
}

/**
 * \brief    Get no level genome mean E operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_E_operon_size( void ) const
{
  return _no_level_genome_E_operon_size;
}

/**
 * \brief    Get no level genome mean TF operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_TF_operon_size( void ) const
{
  return _no_level_genome_TF_operon_size;
}

/**
 * \brief    Get no level genome mean mixed operon size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_genome_mixed_operon_size( void ) const
{
  return _no_level_genome_mixed_operon_size;
}


/* INHERITED STRUCTURE */

/**
 * \brief    Get no level mean inherited size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_inherited_size( void ) const
{
  return _no_level_inherited_size;
}

/**
 * \brief    Get no level mean number of inherited E pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_inherited_nb_E( void ) const
{
  return _no_level_inherited_nb_E;
}

/**
 * \brief    Get no level mean number of inherited TF pearls
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_inherited_nb_TF( void ) const
{
  return _no_level_inherited_nb_TF;
}

/**
 * \brief    Get no level mean number of inherited inner enzymes
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_inherited_nb_inner_enzymes( void ) const
{
  return _no_level_inherited_nb_inner_enzymes;
}

/**
 * \brief    Get no level mean number of inherited inflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_inherited_nb_inflow_pumps( void ) const
{
  return _no_level_inherited_nb_inflow_pumps;
}

/**
 * \brief    Get no level mean number of inherited outflowing pumps
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double TrophicNetwork::get_no_level_inherited_nb_outflow_pumps( void ) const
{
  return _no_level_inherited_nb_outflow_pumps;
}

/*----------------------------
 * SETTERS
 *----------------------------*/


#endif /* defined(__EVOEVO__TrophicNetwork__) */
