
/**
 * \file      Statistics.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Statistics class declaration
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

#ifndef __EVOEVO__Statistics__
#define __EVOEVO__Statistics__

#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <vector>
#include <assert.h>

#include "Macros.h"
#include "Structs.h"
#include "Enums.h"
#include "Parameters.h"
#include "Population.h"
#include "Environment.h"
#include "TrophicNetwork.h"
#include "Tree.h"


class Statistics
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Statistics( void ) = delete;
  Statistics( Parameters* parameters, Population* population, Environment* environment, TrophicNetwork* trophic_network, Tree* lineage_tree, Tree* phylogenetic_tree, bool clean_statistic_files );
  Statistics( const Statistics& stats );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Statistics( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*----------------------------
   * SETTERS
   *----------------------------*/
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  void clean_files( void );
  void open_files( bool clean_statistic_files );
  void write_headers( void );
  void write_stats( void );
  void flush_files( void );
  void close_files( void );
  
  void init_variables( void );
  void add_individual( size_t pos );
  void add_best_individual( void );
  void compute_mean_and_var( void );
  
  void write_best_genome_file( void );
  void write_best_inherited_proteins_file( void );
  void write_best_genetic_regulation_network_file( void );
  void write_best_metabolic_network_file( void );
  void write_best_metabolic_amounts_file( void );
  void write_last_environment_metabolic_amounts_file( void );
  void write_last_trophic_network_file( void );
  void write_last_lineage_tree_statistics( void );
  void write_last_phylogenetic_tree_statistics( void );
  void write_last_backup_file( size_t last_exp_backup, size_t last_tree_backup );
  
  void plot_population_figures( std::string app_path, std::string executable_name );
  void plot_lineage_tree_figures( std::string app_path, std::string executable_name );
  void plot_phylogenetic_tree_figures( std::string app_path, std::string executable_name );
  void plot_trophic_network_figures( std::string app_path, std::string executable_name );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ MEAN statistical variables */
  
  /* PHENOTYPE */
  double _mean_generations;                /*!< Number of generations                    */
  double _mean_inherited_TF_amount;        /*!< Inherited TF amount                      */
  double _mean_inherited_E_amount;         /*!< Inherited E amount                       */
  double _mean_TF_amount;                  /*!< TF amount                                */
  double _mean_E_amount;                   /*!< E amount                                 */
  double _mean_inherited_metabolic_amount; /*!< Inherited metabolic amount               */
  double _mean_metabolic_amount;           /*!< Metabolic amount                         */
  double _mean_energy;                     /*!< Energy amount                            */
  double _mean_score;                      /*!< Score                                    */
  double _mean_lifespan;                   /*!< Lifespan                                 */
  double _mean_number_of_divisions;        /*!< Number of divisions                      */
  double _mean_toxicity;                   /*!< Toxicity accumulation                    */
  double _mean_metabolic_uptake;           /*!< Amount of metabolic uptake               */
  double _mean_metabolic_release;          /*!< Amount of metabolic release              */
  double _mean_metabolic_growth_rate;      /*!< Metabolic growth rate                    */
  double _mean_Dmetabolic_growth_rate;     /*!< Metabolic growth rate difference         */
  double _mean_grn_nb_nodes;               /*!< Number of nodes in the GRN               */
  double _mean_grn_nb_edges;               /*!< Number of edges in the GRN               */
  double _mean_metabolic_nb_nodes;         /*!< Number of nodes in the metabolic network */
  double _mean_metabolic_nb_edges;         /*!< Number of edges in the metabolic network */
  double _mean_regulation_redundancy;      /*!< Regulation redundancy                    */
  double _mean_metabolic_redundancy;       /*!< Metabolic redundancy                     */
  double _mean_cumulated_error;            /*!< Cumulated error                          */
  
  /* GENOME STRUCTURE */
  double _mean_genome_size;                   /*!< Genome size                                                   */
  double _mean_functional_size;               /*!< Functional size                                               */
  double _mean_genome_nb_NC;                  /*!< Number of non coding type (NC) in the genome                  */
  double _mean_genome_nb_E;                   /*!< Number of enzyme type (E) in the genome                       */
  double _mean_genome_nb_TF;                  /*!< Number of transcription factor type (TF) in the genome        */
  double _mean_genome_nb_BS;                  /*!< Number of binding site type (BS) in the genome                */
  double _mean_genome_nb_P;                   /*!< Number of promoter type (P) in the genome                     */
  double _mean_genome_nb_inner_enzymes;       /*!< Number of inner enzymes in the genome                         */
  double _mean_genome_nb_inflow_pumps;        /*!< Number of inflowing pumps in the genome                       */
  double _mean_genome_nb_outflow_pumps;       /*!< Number of outflowing pumps in the genome                      */
  double _mean_genome_nb_functional_regions;  /*!< Number of functional regions in the genome                    */
  double _mean_genome_nb_enhancers;           /*!< Number of enhancers in functional regions in the genome       */
  double _mean_genome_nb_operators;           /*!< Number of operators in functional regions in the genome       */
  double _mean_genome_nb_E_regions;           /*!< Number of functional regions containing only TF in the genome */
  double _mean_genome_nb_TF_regions;          /*!< Number of functional regions containing only TF in the genome */
  double _mean_genome_nb_mixed_regions;       /*!< Number of functional regions mixing E and TF in the genome    */
  double _mean_genome_functional_region_size; /*!< Size of functional regions in the genome                      */
  double _mean_genome_E_region_size;          /*!< Size of E regions in the genome                               */
  double _mean_genome_TF_region_size;         /*!< Size of TF regions in the genome                              */
  double _mean_genome_mixed_region_size;      /*!< Size of mixed regions in the genome                           */
  double _mean_genome_enhancer_size;          /*!< Size of enhancer sites in the genome                          */
  double _mean_genome_operator_size;          /*!< Size of operator sites in the genome                          */
  double _mean_genome_operon_size;            /*!< Size of operons in the genome                                 */
  double _mean_genome_E_operon_size;          /*!< Size of operons containing only E in the genome               */
  double _mean_genome_TF_operon_size;         /*!< Size of operons containing only TF in the genome              */
  double _mean_genome_mixed_operon_size;      /*!< Size of operons containing mixing E and TF in the genome      */
  
  /* INHERITED STRUCTURE */
  double _mean_inherited_size;             /*!< Number of inherited proteins                                   */
  double _mean_inherited_nb_E;             /*!< Number of enzyme type (E) in inherited proteins                */
  double _mean_inherited_nb_TF;            /*!< Number of transcription factor type (TF) in inherited proteins */
  double _mean_inherited_nb_inner_enzymes; /*!< Number of inner enzymes in inherited proteins                  */
  double _mean_inherited_nb_inflow_pumps;  /*!< Number of inflowing pumps in inherited proteins                */
  double _mean_inherited_nb_outflow_pumps; /*!< Number of outflowing pumps in inherited proteins               */
  
  /*------------------------------------------------------------------ VARIANCE statistical variables */
  
  /* PHENOTYPE */
  double _var_generations;                /*!< Number of generations                    */
  double _var_inherited_TF_amount;        /*!< Inherited TF amount                      */
  double _var_inherited_E_amount;         /*!< Inherited E amount                       */
  double _var_TF_amount;                  /*!< TF amount                                */
  double _var_E_amount;                   /*!< E amount                                 */
  double _var_inherited_metabolic_amount; /*!< Inherited metabolic amount               */
  double _var_metabolic_amount;           /*!< Metabolic amount                         */
  double _var_energy;                     /*!< Energy amount                            */
  double _var_score;                      /*!< Score                                    */
  double _var_lifespan;                   /*!< Lifespan                                 */
  double _var_number_of_divisions;        /*!< Number of divisions                      */
  double _var_toxicity;                   /*!< Toxicity accumulation                    */
  double _var_metabolic_uptake;           /*!< Amount of metabolic uptake               */
  double _var_metabolic_release;          /*!< Amount of metabolic release              */
  double _var_metabolic_growth_rate;      /*!< Metabolic growth rate                    */
  double _var_Dmetabolic_growth_rate;     /*!< Metabolic growth rate difference         */
  double _var_grn_nb_nodes;               /*!< Number of nodes in the GRN               */
  double _var_grn_nb_edges;               /*!< Number of edges in the GRN               */
  double _var_metabolic_nb_nodes;         /*!< Number of nodes in the metabolic network */
  double _var_metabolic_nb_edges;         /*!< Number of edges in the metabolic network */
  double _var_regulation_redundancy;      /*!< Regulation redundancy                    */
  double _var_metabolic_redundancy;       /*!< Metabolic redundancy                     */
  double _var_cumulated_error;            /*!< Cumulated error                          */
  
  /* GENOME STRUCTURE */
  double _var_genome_size;                   /*!< Genome size                                                   */
  double _var_functional_size;               /*!< Functional size                                               */
  double _var_genome_nb_NC;                  /*!< Number of non coding type (NC) in the genome                  */
  double _var_genome_nb_E;                   /*!< Number of enzyme type (E) in the genome                       */
  double _var_genome_nb_TF;                  /*!< Number of transcription factor type (TF) in the genome        */
  double _var_genome_nb_BS;                  /*!< Number of binding site type (BS) in the genome                */
  double _var_genome_nb_P;                   /*!< Number of promoter type (P) in the genome                     */
  double _var_genome_nb_inner_enzymes;       /*!< Number of inner enzymes in the genome                         */
  double _var_genome_nb_inflow_pumps;        /*!< Number of inflowing pumps in the genome                       */
  double _var_genome_nb_outflow_pumps;       /*!< Number of outflowing pumps in the genome                      */
  double _var_genome_nb_functional_regions;  /*!< Number of functional regions in the genome                    */
  double _var_genome_nb_enhancers;           /*!< Number of enhancers in functional regions in the genome       */
  double _var_genome_nb_operators;           /*!< Number of operators in functional regions in the genome       */
  double _var_genome_nb_E_regions;           /*!< Number of functional regions containing only TF in the genome */
  double _var_genome_nb_TF_regions;          /*!< Number of functional regions containing only TF in the genome */
  double _var_genome_nb_mixed_regions;       /*!< Number of functional regions mixing E and TF in the genome    */
  double _var_genome_functional_region_size; /*!< Size of functional regions in the genome                      */
  double _var_genome_E_region_size;          /*!< Size of E regions in the genome                               */
  double _var_genome_TF_region_size;         /*!< Size of TF regions in the genome                              */
  double _var_genome_mixed_region_size;      /*!< Size of mixed regions in the genome                           */
  double _var_genome_enhancer_size;          /*!< Size of enhancer sites in the genome                          */
  double _var_genome_operator_size;          /*!< Size of operator sites in the genome                          */
  double _var_genome_operon_size;            /*!< Size of operons in the genome                                 */
  double _var_genome_E_operon_size;          /*!< Size of operons containing only E in the genome               */
  double _var_genome_TF_operon_size;         /*!< Size of operons containing only TF in the genome              */
  double _var_genome_mixed_operon_size;      /*!< Size of operons containing mixing E and TF in the genome      */
  
  /* INHERITED STRUCTURE */
  double _var_inherited_size;             /*!< Number of inherited proteins                                   */
  double _var_inherited_nb_E;             /*!< Number of enzyme type (E) in inherited proteins                */
  double _var_inherited_nb_TF;            /*!< Number of transcription factor type (TF) in inherited proteins */
  double _var_inherited_nb_inner_enzymes; /*!< Number of inner enzymes in inherited proteins                  */
  double _var_inherited_nb_inflow_pumps;  /*!< Number of inflowing pumps in inherited proteins                */
  double _var_inherited_nb_outflow_pumps; /*!< Number of outflowing pumps in inherited proteins               */
  
  /*------------------------------------------------------------------ BEST statistical variables */
  
  /* PHENOTYPE */
  unsigned long long int _best_id;         /*!< Identifier                               */
  double _best_generations;                /*!< Number of generations                    */
  double _best_inherited_TF_amount;        /*!< Inherited TF amount                      */
  double _best_inherited_E_amount;         /*!< Inherited E amount                       */
  double _best_TF_amount;                  /*!< TF amount                                */
  double _best_E_amount;                   /*!< E amount                                 */
  double _best_inherited_metabolic_amount; /*!< Inherited metabolic amount               */
  double _best_metabolic_amount;           /*!< Metabolic amount                         */
  double _best_energy;                     /*!< Energy amount                            */
  double _best_score;                      /*!< Score                                    */
  double _best_lifespan;                   /*!< Lifespan                                 */
  double _best_number_of_divisions;        /*!< Number of divisions                      */
  double _best_toxicity;                   /*!< Toxicity accumulation                    */
  double _best_metabolic_uptake;           /*!< Amount of metabolic uptake               */
  double _best_metabolic_release;          /*!< Amount of metabolic release              */
  double _best_metabolic_growth_rate;      /*!< Metabolic growth rate                    */
  double _best_Dmetabolic_growth_rate;     /*!< Metabolic growth rate difference         */
  double _best_grn_nb_nodes;               /*!< Number of nodes in the GRN               */
  double _best_grn_nb_edges;               /*!< Number of edges in the GRN               */
  double _best_metabolic_nb_nodes;         /*!< Number of nodes in the metabolic network */
  double _best_metabolic_nb_edges;         /*!< Number of edges in the metabolic network */
  double _best_regulation_redundancy;      /*!< Regulation redundancy                    */
  double _best_metabolic_redundancy;       /*!< Metabolic redundancy                     */
  double _best_trophic_level;              /*!< Trophic network level                    */
  double _best_cumulated_error;            /*!< Cumulated error                          */
  
  /* GENOME STRUCTURE */
  double _best_genome_size;                   /*!< Genome size                                                   */
  double _best_functional_size;               /*!< Functional size                                               */
  double _best_genome_nb_NC;                  /*!< Number of non coding type (NC) in the genome                  */
  double _best_genome_nb_E;                   /*!< Number of enzyme type (E) in the genome                       */
  double _best_genome_nb_TF;                  /*!< Number of transcription factor type (TF) in the genome        */
  double _best_genome_nb_BS;                  /*!< Number of binding site type (BS) in the genome                */
  double _best_genome_nb_P;                   /*!< Number of promoter type (P) in the genome                     */
  double _best_genome_nb_inner_enzymes;       /*!< Number of inner enzymes in the genome                         */
  double _best_genome_nb_inflow_pumps;        /*!< Number of inflowing pumps in the genome                       */
  double _best_genome_nb_outflow_pumps;       /*!< Number of outflowing pumps in the genome                      */
  double _best_genome_nb_functional_regions;  /*!< Number of functional regions in the genome                    */
  double _best_genome_nb_enhancers;           /*!< Number of enhancers in functional regions in the genome       */
  double _best_genome_nb_operators;           /*!< Number of operators in functional regions in the genome       */
  double _best_genome_nb_E_regions;           /*!< Number of functional regions containing only TF in the genome */
  double _best_genome_nb_TF_regions;          /*!< Number of functional regions containing only TF in the genome */
  double _best_genome_nb_mixed_regions;       /*!< Number of functional regions mixing E and TF in the genome    */
  double _best_genome_functional_region_size; /*!< Size of functional regions in the genome                      */
  double _best_genome_E_region_size;          /*!< Size of E regions in the genome                               */
  double _best_genome_TF_region_size;         /*!< Size of TF regions in the genome                              */
  double _best_genome_mixed_region_size;      /*!< Size of mixed regions in the genome                           */
  double _best_genome_enhancer_size;          /*!< Size of enhancer sites in the genome                          */
  double _best_genome_operator_size;          /*!< Size of operator sites in the genome                          */
  double _best_genome_operon_size;            /*!< Size of operons in the genome                                 */
  double _best_genome_E_operon_size;          /*!< Size of operons containing only E in the genome               */
  double _best_genome_TF_operon_size;         /*!< Size of operons containing only TF in the genome              */
  double _best_genome_mixed_operon_size;      /*!< Size of operons containing mixing E and TF in the genome      */
  
  /* INHERITED STRUCTURE */
  double _best_inherited_size;             /*!< Number of inherited proteins                                   */
  double _best_inherited_nb_E;             /*!< Number of enzyme type (E) in inherited proteins                */
  double _best_inherited_nb_TF;            /*!< Number of transcription factor type (TF) in inherited proteins */
  double _best_inherited_nb_inner_enzymes; /*!< Number of inner enzymes in inherited proteins                  */
  double _best_inherited_nb_inflow_pumps;  /*!< Number of inflowing pumps in inherited proteins                */
  double _best_inherited_nb_outflow_pumps; /*!< Number of outflowing pumps in inherited proteins               */
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  void clean_phenotype_mean_file( void );
  void clean_genome_structure_mean_file( void );
  void clean_inherited_proteins_mean_file( void );
  void clean_phenotype_var_file( void );
  void clean_genome_structure_var_file( void );
  void clean_inherited_proteins_var_file( void );
  void clean_phenotype_best_file( void );
  void clean_genome_structure_best_file( void );
  void clean_inherited_proteins_best_file( void );
  void clean_environment_metabolic_amounts_file( void );
  void clean_trophic_network_profile_file( void );
  void clean_level_0_file( void );
  void clean_level_1_file( void );
  void clean_level_2_file( void );
  void clean_no_level_file( void );
  void clean_global_concentrations_file( void );
  void clean_tree_structure_file( void );
  void clean_best_id_file( void );
  
  void write_phenotype_mean_file_header( void );
  void write_genome_structure_mean_file_header( void );
  void write_inherited_proteins_mean_file_header( void );
  void write_phenotype_var_file_header( void );
  void write_genome_structure_var_file_header( void );
  void write_inherited_proteins_var_file_header( void );
  void write_phenotype_best_file_header( void );
  void write_genome_structure_best_file_header( void );
  void write_inherited_proteins_best_file_header( void );
  void write_trophic_network_profile_file_header( void );
  void write_level_0_file_header( void );
  void write_level_1_file_header( void );
  void write_level_2_file_header( void );
  void write_no_level_file_header( void );
  void write_global_concentrations_file_header( void );
  void write_tree_structure_file_header( void );
  void write_best_id_file_header( void );
  
  void write_phenotype_mean_file_stats( void );
  void write_genome_structure_mean_file_stats( void );
  void write_inherited_proteins_mean_file_stats( void );
  void write_phenotype_var_file_stats( void );
  void write_genome_structure_var_file_stats( void );
  void write_inherited_proteins_var_file_stats( void );
  void write_phenotype_best_file_stats( void );
  void write_genome_structure_best_file_stats( void );
  void write_inherited_proteins_best_file_stats( void );
  void write_environment_metabolic_amounts_stats( void );
  void write_trophic_network_profile_stats( void );
  void write_level_0_stats( void );
  void write_level_1_stats( void );
  void write_level_2_stats( void );
  void write_no_level_stats( void );
  void write_global_concentrations_stats( void );
  void write_tree_structure_stats( void );
  void write_best_id_file( void );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ simulation variables */
  
  Parameters*     _parameters;        /*!< Parameters        */
  Population*     _population;        /*!< Population        */
  Environment*    _environment;       /*!< Environment       */
  TrophicNetwork* _trophic_network;   /*!< Trophic network   */
  Tree*           _lineage_tree;      /*!< Lineage tree      */
  Tree*           _phylogenetic_tree; /*!< Phylogenetic tree */
  
  /*------------------------------------------------------------------ files */
  
  std::ofstream _phenotype_mean_file;          /*!< Phenotype mean file          */
  std::ofstream _genome_structure_mean_file;   /*!< Genome structure mean file   */
  std::ofstream _inherited_proteins_mean_file; /*!< Inherited proteins mean file */
  
  std::ofstream _phenotype_var_file;          /*!< Phenotype variance file          */
  std::ofstream _genome_structure_var_file;   /*!< Genome structure variance file   */
  std::ofstream _inherited_proteins_var_file; /*!< Inherited proteins variance file */
  
  std::ofstream _phenotype_best_file;          /*!< Phenotype best file          */
  std::ofstream _genome_structure_best_file;   /*!< Genome structure best file   */
  std::ofstream _inherited_proteins_best_file; /*!< Inherited proteins best file */
  std::ofstream _best_id_file;                 /*!< Best identifier file         */
  
  std::ofstream _best_genome_file;                      /*!< Best genome file                      */
  std::ofstream _best_inherited_proteins_file;          /*!< Best inherited proteins file          */
  std::ofstream _best_grn_nodes_file;                   /*!< Best GRN nodes file                   */
  std::ofstream _best_grn_edges_file;                   /*!< Best GRN edges file                   */
  std::ofstream _best_metabolic_nodes_file;             /*!< Best metabolic nodes file             */
  std::ofstream _best_metabolic_edges_file;             /*!< Best metabolic edges file             */
  std::ofstream _best_metabolic_amounts_file;           /*!< Best metabolic amounts file           */
  std::ofstream _best_inherited_metabolic_amounts_file; /*!< Best inherited metabolic amounts file */
  std::ofstream _environment_metabolic_amounts_file;    /*!< Environment metabolic amounts file    */
  std::ofstream _trophic_network_profile_file;          /*!< Trophic network profile file          */
  std::ofstream _level_0_statistics_file;               /*!< Level 0 statistics file               */
  std::ofstream _level_1_statistics_file;               /*!< Level 1 statistics file               */
  std::ofstream _level_2_statistics_file;               /*!< Level 2 statistics file               */
  std::ofstream _no_level_statistics_file;              /*!< No level statistics file              */
  std::ofstream _global_concentrations_file;            /*!< Evolution of global concentrations    */
  std::ofstream _tree_structure_file;                   /*!< Tree data file                        */
  
  std::ofstream _last_backup_file;                      /*!< Last backup file                      */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/*----------------------------
 * SETTERS
 *----------------------------*/


#endif /* defined(__EVOEVO__Statistics__) */
