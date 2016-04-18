
/**
 * \file      Structs.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Definition of structures
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

#ifndef __EVOEVO__Structs__
#define __EVOEVO__Structs__

#include <iostream>
#include <vector>
#include <cstring>
#include <unordered_map>
#include <gsl/gsl_odeiv2.h>

#include "Macros.h"
#include "Enums.h"
#include "Prng.h"

class Genome;
class InheritedProteins;
class SpeciesList;


/******************************************************************************************/

/**
 * \brief   pearl struct
 * \details Defines the structure of a pearl
 */
typedef struct
{
  
  /*------------------------------------------------------------------ Global attributes */
  
  pearl_type             type;              /*!< Type of pearl            */
  unsigned long long int identifier;        /*!< Pearl identifier         */
  unsigned long long int parent_identifier; /*!< Parental gene identifier */
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  int    s;    /*!< Substrate species tag */
  int    p;    /*!< Product species tag   */
  double km;   /*!< Km constant           */
  double kcat; /*!< Kcat constant         */
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  int    BS_tag;         /*!< Binding site tag      */
  int    coE_tag;        /*!< Co-enzyme species tag */
  bool   free_activity;  /*!< Free activity         */
  bool   bound_activity; /*!< Bound activity        */
  size_t window;         /*!< Window size           */
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  int TF_tag; /*!< Transcription factor tag */
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  double basal_expression_level; /*!< Basal expression level of the promoter */
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  bool functional; /*!< Indicates if he pearl is functional or not */
  
} pearl;

/******************************************************************************************/

/**
 * \brief   String of pearls struct
 * \details Defines the structure of a string of pearls
 */
typedef struct
{
  
  size_t size;        /*!< Size of the array */
  size_t buffer_size; /*!< Buffer size       */
  pearl* x;           /*!< Array of pearls   */

} string_of_pearls;

/******************************************************************************************/

/**
 * \brief   Reaction list struct
 * \details Contains the list of reactions encoded in the genome
 */
typedef struct
{
  
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  Prng* prng; /*!< Pseudorandom numbers generator */
  
  /*------------------------------------------------------------------ list of genetic regulation network equations */
  
  size_t                      grn_N;                         /*!< Number of GRN equations                             */
  std::vector<size_t>         grn_promoter;                  /*!< List of promoter positions                          */
  std::vector<size_t>         grn_Nenhancer;                 /*!< Number of enhancer bindings per promoter            */
  std::vector<size_t>         grn_Noperator;                 /*!< Number of operator bindings per promoter            */
  std::vector<size_t>         grn_Ngenes;                    /*!< Number of genes regulated per promoter              */
  std::vector<size_t>         grn_enhancer_nb_BS;            /*!< Number of binding site per enhancer                 */
  std::vector<size_t>         grn_enhancer_TF_list;          /*!< List of enhancer binding TF indexes per promoter    */
  std::vector<double>         grn_enhancer_affinity_list;    /*!< List of enhancer binding TF affinities per promoter */
  std::vector<int>            grn_enhancer_coe_list;         /*!< List of TF co-enzymes per enhancer                  */
  std::vector<co_enzyme_type> grn_enhancer_coe_type;         /*!< List of co_enzyme types per enhancer                */
  std::vector<size_t>         grn_operator_nb_BS;            /*!< Number of binding site per operator                 */
  std::vector<size_t>         grn_operator_TF_list;          /*!< List of operator binding TF indexes per promoter    */
  std::vector<double>         grn_operator_affinity_list;    /*!< List of operator binding TF affinities per promoter */
  std::vector<int>            grn_operator_coe_list;         /*!< List of TF co-enzymes per operator                  */
  std::vector<co_enzyme_type> grn_operator_coe_type;         /*!< List of co_enzyme types per operator                */
  std::vector<size_t>         grn_regulated_genes;           /*!< List of regulated genes per promoter                */
  std::vector<double>         grn_beta;                      /*!< List of basal expression levels per promoter        */
  double                      grn_hill_theta;                /*!< Parameter theta of the Hill equation                */
  double                      grn_hill_n;                    /*!< Parameter n of the Hill equation                    */
  double                      grn_protein_degradation_rate;  /*!< Protein degradation rate                            */
  double                      grn_energy_transcription_cost; /*!< Energetical cost of transcription                   */
  
  /*------------------------------------------------------------------ list of metabolic equations */
  
  size_t                     metabolic_N;       /*!< Number of metabolic equations */
  std::vector<reaction_type> metabolic_type;    /*!< Type of the equation          */
  std::vector<int>           metabolic_s;       /*!< substrate tag                 */
  std::vector<int>           metabolic_p;       /*!< Product tag                   */
  std::vector<double>        metabolic_km;      /*!< Km constant                   */
  std::vector<double>        metabolic_kcat;    /*!< Kcat constant                 */
  std::vector<double>        metabolic_delta_g; /*!< Energy cost of the reaction   */
  std::vector<size_t>        metabolic_e;       /*!< Enzyme index                  */
  
  /*------------------------------------------------------------------ size of each state vector subspace */
  
  size_t Ngenome;    /*!< Length of the genome state vector       */
  size_t Ninherited; /*!< Length of the inherited proteins vector */
  size_t Ncell;      /*!< Length of the cell state vector         */
  size_t Nenv;       /*!< Length of the environement state vector */
  size_t N;          /*!< Length of the state vector X            */
  
  /*------------------------------------------------------------------ modeling schemes */
  
  bool energy_constraints;    /*!< Indicates if energy constraints are activated   */
  bool membrane_permeability; /*!< Indicates if membrane permeability is activated */
  bool metabolic_inheritance; /*!< Indicates if metabolic inheritance is activated */
  bool enzymatic_inheritance; /*!< Indicates if enzymatic inheritance is activated */
  bool co_enzyme_activity;    /*!< Indicates if co-enzyme activity is activated    */
  
  /*------------------------------------------------------------------ membrane permeability */
  
  double kmp; /*!< Membrane permeability */
  
  /*------------------------------------------------------------------ metabolic exchanges with the environment */
  
  double metabolic_uptake;  /*!< Metabolic uptake  */
  double metabolic_release; /*!< Metabolic release */
  
  /*------------------------------------------------------------------ timestep ratio */
  
  double timestep_ratio; /*!< Timestep ratio between GRN and metabolism timesteps */
  
  /*------------------------------------------------------------------ environment interaction scheme */
  
  bool interacts_with_environment; /*!< Indicates if cells interact with the environment */
  
  /*------------------------------------------------------------------ test variables */
  
#ifdef DEBUG
  Genome*            genome;             /*!< Related cell's genome             */
  InheritedProteins* inherited_proteins; /*!< Related cell's inherited proteins */
  SpeciesList*       cell_species_list;  /*!< Related cell's species list       */
  SpeciesList*       env_species_list;   /*!< Local environment's species list  */
#endif
  
} reaction_list;

/******************************************************************************************/

/**
 * \brief   Distribution law of metabolic range
 * \details Defines the distribution law of metabolic range
 */
typedef struct
{
  
  /*------------------------------------------------------------------ Distribution law */
  
  metabolic_range law; /*!< Distribution law */
  
  /*------------------------------------------------------------------ Uniform law parameters */
  
  double min; /*!< uniform minimum bound */
  double max; /*!< uniform maximum bound */
  
  /*------------------------------------------------------------------ Gaussian law parameters */
  
  double mu;    /*!< Gaussian mean     */
  double sigma; /*!< Gaussian variance */
  
  /*------------------------------------------------------------------ Exponential law parameters */
  
  double lambda; /*!< Exponential mean */
  
} distribution_law;

/******************************************************************************************/

/**
 * \brief   Gompertz law parameters
 * \details Defines the parameters of the Gompertz law used to compute death probability depending on cell age
 */
typedef struct
{
  
  double a; /*!< parameter a (asymptote)                     */
  double b; /*!< parameter b (displacement along the x axis) */
  double c; /*!< parameter c (growth rate)                   */
  
} gompert_law;

/******************************************************************************************/

/**
 * \brief   Environment properties
 * \details --
 */
typedef struct
{
  
  size_t                             number_of_init_cycles;   /*!< Number of environment initialization cycles            */
  distribution_law                   species_tag_range;       /*!< Environment species tag range distribution law         */
  distribution_law                   concentration_range;     /*!< Environment concentration range distribution law       */
  distribution_law                   number_of_species_range; /*!< Environment number of species range distribution law   */
  environment_interaction_scheme     interaction_scheme;      /*!< Environment interaction scheme (interact or fixed)     */
  environment_renewal_scheme         renewal_scheme;          /*!< Environment renewal scheme (keep or clear matter)      */
  environment_variation_scheme       variation_scheme;        /*!< Environment variation mode (random, cyclic, periodic)  */
  environment_variation_localization variation_localization;  /*!< Environment introduction localization (global, center) */
  double                             introduction_rate;       /*!< Environment introduction rate                          */
  double                             diffusion_rate;          /*!< Environment diffusion rate                             */
  double                             degradation_rate;        /*!< Environment degradation rate                           */
  
} environment_properties;

/******************************************************************************************/


#endif /* defined(__EVOEVO__Structs__) */
