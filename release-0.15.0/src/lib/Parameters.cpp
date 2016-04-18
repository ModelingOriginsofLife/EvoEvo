
/**
 * \file      Parameters.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Parameters class definition
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

#include "Parameters.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Default constructor
 * \details  --
 * \param    void
 * \return   \e void
 */
Parameters::Parameters( void )
{
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  _prng = new Prng();
  _seed = 0;
  
  /*------------------------------------------------------------------ parallel computing */
  
  _parallel_computing = false;
  
  /*------------------------------------------------------------------ simulation schemes */
  
  _energy_constraints    = false;
  _membrane_permeability = false;
  _metabolic_inheritance = false;
  _enzymatic_inheritance = false;
  _co_enzyme_activity    = false;
  _score_scheme          = ESSENTIAL_METABOLITES_SUM;
  _selection_threshold   = 0.4;
  
  /*------------------------------------------------------------------ space */
  
  _width  = 32;
  _height = 32;
  
  /*------------------------------------------------------------------ output */
  
  _experiment_backup_step  = 1000;
  _tree_backup_step        = 1000;
  _figures_generation_step = 500;
  
  /*------------------------------------------------------------------ genome level */
  
  _metabolite_tag_initial_range.law    = UNIFORM;
  _metabolite_tag_initial_range.min    = 1;
  _metabolite_tag_initial_range.max    = 5;
  _metabolite_tag_initial_range.mu     = 0.0;
  _metabolite_tag_initial_range.sigma  = 0.0;
  _metabolite_tag_initial_range.lambda = 0.0;
  
  _binding_site_tag_initial_range.law    = UNIFORM;
  _binding_site_tag_initial_range.min    = 1;
  _binding_site_tag_initial_range.max    = 5;
  _binding_site_tag_initial_range.mu     = 0.0;
  _binding_site_tag_initial_range.sigma  = 0.0;
  _binding_site_tag_initial_range.lambda = 0.0;
  
  _co_enzyme_tag_initial_range.law    = UNIFORM;
  _co_enzyme_tag_initial_range.min    = 1;
  _co_enzyme_tag_initial_range.max    = 5;
  _co_enzyme_tag_initial_range.mu     = 0.0;
  _co_enzyme_tag_initial_range.sigma  = 0.0;
  _co_enzyme_tag_initial_range.lambda = 0.0;
  
  _transcription_factor_tag_initial_range.law    = UNIFORM;
  _transcription_factor_tag_initial_range.min    = 1;
  _transcription_factor_tag_initial_range.max    = 5;
  _transcription_factor_tag_initial_range.mu     = 0.0;
  _transcription_factor_tag_initial_range.sigma  = 0.0;
  _transcription_factor_tag_initial_range.lambda = 0.0;
  
  _initial_binding_window                 = 10;
  _initial_number_of_NC_pearls            = 10;
  _initial_number_of_E_pearls             = 10;
  _initial_number_of_TF_pearls            = 10;
  _initial_number_of_BS_pearls            = 10;
  _initial_number_of_P_pearls             = 10;
  _point_mutation_rate                    = 1e-04;
  _duplication_rate                       = 1e-04;
  _deletion_rate                          = 1e-04;
  _translocation_rate                     = 1e-04;
  _inversion_rate                         = 1e-04;
  _transition_rate                        = 1e-05;
  _substrate_tag_mutation_size            = 20;
  _product_tag_mutation_size              = 20;
  _km_mutation_size                       = 0.1;
  _kcat_mutation_size                     = 0.1;
  _binding_size_tag_mutation_size         = 2;
  _co_enzyme_tag_mutation_size            = 2;
  _transcription_factor_tag_mutation_size = 2;
  _basal_expression_level_mutation_size   = 0.1;
  _maximum_genome_size                    = 10000;
  
  /*------------------------------------------------------------------ genetic regulation network level */
  
  _genetic_regulation_network_timestep = 1.0;
  _hill_function_theta                 = 0.5;
  _hill_function_n                     = 4.0;
  _protein_degradation_rate            = 1.0;
  
  /*------------------------------------------------------------------ metabolic network level */
  
  _metabolism_timestep                          = 1.0;
  _essential_metabolites_toxicity_threshold     = 1.0;
  _non_essential_metabolites_toxicity_threshold = 1.0;
  _initial_metabolites_amount_in_cells          = 1.0;
  _energy_reaction_cost_factor                  = 0.1;
  _energy_pumping_cost                          = -0.1;
  _energy_transcription_cost                    = -0.1;
  _energy_toxicity_threshold                    = 10.0;
  
  /*------------------------------------------------------------------ cell level */
  
  _death_probability             = 0.1;
  _variable_lifespan             = false;
  _gompert_law.a                 = 1.0;
  _gompert_law.b                 = 20.0;
  _gompert_law.c                 = 0.1;
  _initial_membrane_permeability = 0.0;
  
  /*------------------------------------------------------------------ population level */
  
  _migration_rate = 0.0;
  _hgt_rate       = 0.0;
  
  /*------------------------------------------------------------------ environment level */
  
  _environment_properties.number_of_init_cycles = 1;
  
  _environment_properties.species_tag_range.law    = UNIFORM;
  _environment_properties.species_tag_range.min    = 1;
  _environment_properties.species_tag_range.max    = 5;
  _environment_properties.species_tag_range.mu     = 0.0;
  _environment_properties.species_tag_range.sigma  = 0.0;
  _environment_properties.species_tag_range.lambda = 0.0;
  
  _environment_properties.concentration_range.law    = UNIFORM;
  _environment_properties.concentration_range.min    = 1;
  _environment_properties.concentration_range.max    = 5;
  _environment_properties.concentration_range.mu     = 0.0;
  _environment_properties.concentration_range.sigma  = 0.0;
  _environment_properties.concentration_range.lambda = 0.0;
  
  _environment_properties.number_of_species_range.law    = UNIFORM;
  _environment_properties.number_of_species_range.min    = 1;
  _environment_properties.number_of_species_range.max    = 5;
  _environment_properties.number_of_species_range.mu     = 0.0;
  _environment_properties.number_of_species_range.sigma  = 0.0;
  _environment_properties.number_of_species_range.lambda = 0.0;
  
  _environment_properties.interaction_scheme     = INTERACTION;
  _environment_properties.renewal_scheme         = KEEP_MATTER;
  _environment_properties.variation_scheme       = RANDOM_SCHEME;
  _environment_properties.variation_localization = RANDOM_LOCALIZATION;
  
  _environment_properties.introduction_rate = 0.0;
  _environment_properties.diffusion_rate    = 0.0;
  _environment_properties.degradation_rate  = 0.0;
  
  /*------------------------------------------------------------------ prime numbers */
  
  build_prime_numbers_list(5000);
}

/**
 * \brief    Contructor from backup file
 * \details  Load Parameter class from backup file
 * \param    size_t backup_time
 * \return   \e void
 */
Parameters::Parameters( size_t backup_time )
{
  std::stringstream backup_file_name;
  backup_file_name << "./parameters/parameters_" << backup_time;
  gzFile backup_file = gzopen(backup_file_name.str().c_str(), "r");
  
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  _prng = new Prng(backup_time);
  gzread( backup_file, &_seed, sizeof(_seed) );
  
  /*------------------------------------------------------------------ parallel computing */
  
  gzread( backup_file, &_parallel_computing, sizeof(_parallel_computing) );
  
  /*------------------------------------------------------------------ simulation schemes */
  
  gzread( backup_file, &_energy_constraints,    sizeof(_energy_constraints) );
  gzread( backup_file, &_membrane_permeability, sizeof(_membrane_permeability) );
  gzread( backup_file, &_metabolic_inheritance, sizeof(_metabolic_inheritance) );
  gzread( backup_file, &_enzymatic_inheritance, sizeof(_enzymatic_inheritance) );
  gzread( backup_file, &_co_enzyme_activity,    sizeof(_co_enzyme_activity) );
  gzread( backup_file, &_score_scheme,          sizeof(_score_scheme) );
  gzread( backup_file, &_selection_threshold,   sizeof(_selection_threshold) );
  
  /*------------------------------------------------------------------ space */
  
  gzread( backup_file, &_width,  sizeof(_width) );
  gzread( backup_file, &_height, sizeof(_height) );
  
  /*------------------------------------------------------------------ output */
  
  gzread( backup_file, &_experiment_backup_step,  sizeof(_experiment_backup_step) );
  gzread( backup_file, &_tree_backup_step,        sizeof(_tree_backup_step) );
  gzread( backup_file, &_figures_generation_step, sizeof(_figures_generation_step) );
  
  /*------------------------------------------------------------------ genome level */
  
  gzread( backup_file, &_metabolite_tag_initial_range.law,    sizeof(_metabolite_tag_initial_range.law) );
  gzread( backup_file, &_metabolite_tag_initial_range.min,    sizeof(_metabolite_tag_initial_range.min) );
  gzread( backup_file, &_metabolite_tag_initial_range.max,    sizeof(_metabolite_tag_initial_range.max) );
  gzread( backup_file, &_metabolite_tag_initial_range.mu,     sizeof(_metabolite_tag_initial_range.mu) );
  gzread( backup_file, &_metabolite_tag_initial_range.sigma,  sizeof(_metabolite_tag_initial_range.sigma) );
  gzread( backup_file, &_metabolite_tag_initial_range.lambda, sizeof(_metabolite_tag_initial_range.lambda) );
  
  gzread( backup_file, &_binding_site_tag_initial_range.law,    sizeof(_binding_site_tag_initial_range.law) );
  gzread( backup_file, &_binding_site_tag_initial_range.min,    sizeof(_binding_site_tag_initial_range.min) );
  gzread( backup_file, &_binding_site_tag_initial_range.max,    sizeof(_binding_site_tag_initial_range.max) );
  gzread( backup_file, &_binding_site_tag_initial_range.mu,     sizeof(_binding_site_tag_initial_range.mu) );
  gzread( backup_file, &_binding_site_tag_initial_range.sigma,  sizeof(_binding_site_tag_initial_range.sigma) );
  gzread( backup_file, &_binding_site_tag_initial_range.lambda, sizeof(_binding_site_tag_initial_range.lambda) );
  
  gzread( backup_file, &_co_enzyme_tag_initial_range.law,    sizeof(_co_enzyme_tag_initial_range.law) );
  gzread( backup_file, &_co_enzyme_tag_initial_range.min,    sizeof(_co_enzyme_tag_initial_range.min) );
  gzread( backup_file, &_co_enzyme_tag_initial_range.max,    sizeof(_co_enzyme_tag_initial_range.max) );
  gzread( backup_file, &_co_enzyme_tag_initial_range.mu,     sizeof(_co_enzyme_tag_initial_range.mu) );
  gzread( backup_file, &_co_enzyme_tag_initial_range.sigma,  sizeof(_co_enzyme_tag_initial_range.sigma) );
  gzread( backup_file, &_co_enzyme_tag_initial_range.lambda, sizeof(_co_enzyme_tag_initial_range.lambda) );
  
  gzread( backup_file, &_transcription_factor_tag_initial_range.law,    sizeof(_transcription_factor_tag_initial_range.law) );
  gzread( backup_file, &_transcription_factor_tag_initial_range.min,    sizeof(_transcription_factor_tag_initial_range.min) );
  gzread( backup_file, &_transcription_factor_tag_initial_range.max,    sizeof(_transcription_factor_tag_initial_range.max) );
  gzread( backup_file, &_transcription_factor_tag_initial_range.mu,     sizeof(_transcription_factor_tag_initial_range.mu) );
  gzread( backup_file, &_transcription_factor_tag_initial_range.sigma,  sizeof(_transcription_factor_tag_initial_range.sigma) );
  gzread( backup_file, &_transcription_factor_tag_initial_range.lambda, sizeof(_transcription_factor_tag_initial_range.lambda) );
  
  gzread( backup_file, &_initial_binding_window,                 sizeof(_initial_binding_window) );
  gzread( backup_file, &_initial_number_of_NC_pearls,            sizeof(_initial_number_of_NC_pearls) );
  gzread( backup_file, &_initial_number_of_E_pearls,             sizeof(_initial_number_of_E_pearls) );
  gzread( backup_file, &_initial_number_of_TF_pearls,            sizeof(_initial_number_of_TF_pearls) );
  gzread( backup_file, &_initial_number_of_BS_pearls,            sizeof(_initial_number_of_BS_pearls) );
  gzread( backup_file, &_initial_number_of_P_pearls,             sizeof(_initial_number_of_P_pearls) );
  gzread( backup_file, &_point_mutation_rate,                    sizeof(_point_mutation_rate) );
  gzread( backup_file, &_duplication_rate,                       sizeof(_duplication_rate) );
  gzread( backup_file, &_deletion_rate,                          sizeof(_deletion_rate) );
  gzread( backup_file, &_translocation_rate,                     sizeof(_translocation_rate) );
  gzread( backup_file, &_inversion_rate,                         sizeof(_inversion_rate) );
  gzread( backup_file, &_transition_rate,                        sizeof(_transition_rate) );
  gzread( backup_file, &_substrate_tag_mutation_size,            sizeof(_substrate_tag_mutation_size) );
  gzread( backup_file, &_product_tag_mutation_size,              sizeof(_product_tag_mutation_size) );
  gzread( backup_file, &_km_mutation_size,                       sizeof(_km_mutation_size) );
  gzread( backup_file, &_kcat_mutation_size,                     sizeof(_kcat_mutation_size) );
  gzread( backup_file, &_binding_size_tag_mutation_size,         sizeof(_binding_size_tag_mutation_size) );
  gzread( backup_file, &_co_enzyme_tag_mutation_size,            sizeof(_co_enzyme_tag_mutation_size) );
  gzread( backup_file, &_transcription_factor_tag_mutation_size, sizeof(_transcription_factor_tag_mutation_size) );
  gzread( backup_file, &_basal_expression_level_mutation_size,   sizeof(_basal_expression_level_mutation_size) );
  gzread( backup_file, &_maximum_genome_size,                    sizeof(_maximum_genome_size) );
  
  /*------------------------------------------------------------------ genetic regulation network level */
  
  gzread( backup_file, &_genetic_regulation_network_timestep, sizeof(_genetic_regulation_network_timestep) );
  gzread( backup_file, &_hill_function_theta,                 sizeof(_hill_function_theta) );
  gzread( backup_file, &_hill_function_n,                     sizeof(_hill_function_n) );
  gzread( backup_file, &_protein_degradation_rate,            sizeof(_protein_degradation_rate) );
  
  /*------------------------------------------------------------------ metabolic network level */
  
  gzread( backup_file, &_metabolism_timestep,                          sizeof(_metabolism_timestep) );
  gzread( backup_file, &_essential_metabolites_toxicity_threshold,     sizeof(_essential_metabolites_toxicity_threshold) );
  gzread( backup_file, &_non_essential_metabolites_toxicity_threshold, sizeof(_non_essential_metabolites_toxicity_threshold) );
  gzread( backup_file, &_initial_metabolites_amount_in_cells,          sizeof(_initial_metabolites_amount_in_cells) );
  gzread( backup_file, &_energy_reaction_cost_factor,                  sizeof(_energy_reaction_cost_factor) );
  gzread( backup_file, &_energy_pumping_cost,                          sizeof(_energy_pumping_cost) );
  gzread( backup_file, &_energy_transcription_cost,                    sizeof(_energy_transcription_cost) );
  gzread( backup_file, &_energy_toxicity_threshold,                    sizeof(_energy_toxicity_threshold) );
  
  /*------------------------------------------------------------------ cell level */
  
  gzread( backup_file, &_death_probability,             sizeof(_death_probability) );
  gzread( backup_file, &_variable_lifespan,             sizeof(_variable_lifespan) );
  gzread( backup_file, &_gompert_law.a,                 sizeof(_gompert_law.a) );
  gzread( backup_file, &_gompert_law.b,                 sizeof(_gompert_law.b) );
  gzread( backup_file, &_gompert_law.c,                 sizeof(_gompert_law.c) );
  gzread( backup_file, &_initial_membrane_permeability, sizeof(_initial_membrane_permeability) );
  
  /*------------------------------------------------------------------ population level */
  
  gzread( backup_file, &_migration_rate, sizeof(_migration_rate) );
  gzread( backup_file, &_hgt_rate,       sizeof(_hgt_rate) );
  
  /*------------------------------------------------------------------ environment level */
  
  gzread( backup_file, &_environment_properties.number_of_init_cycles,    sizeof(_environment_properties.number_of_init_cycles) );
  
  gzread( backup_file, &_environment_properties.species_tag_range.law,    sizeof(_environment_properties.species_tag_range.law) );
  gzread( backup_file, &_environment_properties.species_tag_range.min,    sizeof(_environment_properties.species_tag_range.min) );
  gzread( backup_file, &_environment_properties.species_tag_range.max,    sizeof(_environment_properties.species_tag_range.max) );
  gzread( backup_file, &_environment_properties.species_tag_range.mu,     sizeof(_environment_properties.species_tag_range.mu) );
  gzread( backup_file, &_environment_properties.species_tag_range.sigma,  sizeof(_environment_properties.species_tag_range.sigma) );
  gzread( backup_file, &_environment_properties.species_tag_range.lambda, sizeof(_environment_properties.species_tag_range.lambda) );
  
  gzread( backup_file, &_environment_properties.concentration_range.law,    sizeof(_environment_properties.concentration_range.law) );
  gzread( backup_file, &_environment_properties.concentration_range.min,    sizeof(_environment_properties.concentration_range.min) );
  gzread( backup_file, &_environment_properties.concentration_range.max,    sizeof(_environment_properties.concentration_range.max) );
  gzread( backup_file, &_environment_properties.concentration_range.mu,     sizeof(_environment_properties.concentration_range.mu) );
  gzread( backup_file, &_environment_properties.concentration_range.sigma,  sizeof(_environment_properties.concentration_range.sigma) );
  gzread( backup_file, &_environment_properties.concentration_range.lambda, sizeof(_environment_properties.concentration_range.lambda) );
  
  gzread( backup_file, &_environment_properties.number_of_species_range.law,    sizeof(_environment_properties.number_of_species_range.law) );
  gzread( backup_file, &_environment_properties.number_of_species_range.min,    sizeof(_environment_properties.number_of_species_range.min) );
  gzread( backup_file, &_environment_properties.number_of_species_range.max,    sizeof(_environment_properties.number_of_species_range.max) );
  gzread( backup_file, &_environment_properties.number_of_species_range.mu,     sizeof(_environment_properties.number_of_species_range.mu) );
  gzread( backup_file, &_environment_properties.number_of_species_range.sigma,  sizeof(_environment_properties.number_of_species_range.sigma) );
  gzread( backup_file, &_environment_properties.number_of_species_range.lambda, sizeof(_environment_properties.number_of_species_range.lambda) );
  
  gzread( backup_file, &_environment_properties.interaction_scheme,     sizeof(_environment_properties.interaction_scheme) );
  gzread( backup_file, &_environment_properties.renewal_scheme,         sizeof(_environment_properties.renewal_scheme) );
  gzread( backup_file, &_environment_properties.variation_scheme,       sizeof(_environment_properties.variation_scheme) );
  gzread( backup_file, &_environment_properties.variation_localization, sizeof(_environment_properties.variation_localization) );
  
  gzread( backup_file, &_environment_properties.introduction_rate, sizeof(_environment_properties.introduction_rate) );
  gzread( backup_file, &_environment_properties.diffusion_rate,    sizeof(_environment_properties.diffusion_rate) );
  gzread( backup_file, &_environment_properties.degradation_rate,  sizeof(_environment_properties.degradation_rate) );
  
  gzclose(backup_file);
  
  /*------------------------------------------------------------------ prime numbers */
  
  build_prime_numbers_list(5000);
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const Parameters& parameters
 * \return   \e void
 */
Parameters::Parameters( const Parameters& parameters )
{
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  _prng = new Prng(*parameters._prng);
  _seed = parameters._seed;
  
  /*------------------------------------------------------------------ parallel computing */
  
  _parallel_computing = parameters._parallel_computing;
  
  /*------------------------------------------------------------------ simulation schemes */
  
  _energy_constraints    = parameters._energy_constraints;
  _membrane_permeability = parameters._membrane_permeability;
  _metabolic_inheritance = parameters._metabolic_inheritance;
  _enzymatic_inheritance = parameters._enzymatic_inheritance;
  _co_enzyme_activity    = parameters._co_enzyme_activity;
  _score_scheme          = parameters._score_scheme;
  _selection_threshold   = parameters._selection_threshold;
  
  /*------------------------------------------------------------------ space */
  
  _width  = parameters._width;
  _height = parameters._height;
  
  /*------------------------------------------------------------------ output */
  
  _experiment_backup_step  = parameters._experiment_backup_step;
  _tree_backup_step        = parameters._tree_backup_step;
  _figures_generation_step = parameters._figures_generation_step;
  
  /*------------------------------------------------------------------ genome level */
  
  _metabolite_tag_initial_range.law    = parameters._metabolite_tag_initial_range.law;
  _metabolite_tag_initial_range.min    = parameters._metabolite_tag_initial_range.min;
  _metabolite_tag_initial_range.max    = parameters._metabolite_tag_initial_range.max;
  _metabolite_tag_initial_range.mu     = parameters._metabolite_tag_initial_range.mu;
  _metabolite_tag_initial_range.sigma  = parameters._metabolite_tag_initial_range.sigma;
  _metabolite_tag_initial_range.lambda = parameters._metabolite_tag_initial_range.lambda;
  
  _binding_site_tag_initial_range.law    = parameters._binding_site_tag_initial_range.law;
  _binding_site_tag_initial_range.min    = parameters._binding_site_tag_initial_range.min;
  _binding_site_tag_initial_range.max    = parameters._binding_site_tag_initial_range.max;
  _binding_site_tag_initial_range.mu     = parameters._binding_site_tag_initial_range.mu;
  _binding_site_tag_initial_range.sigma  = parameters._binding_site_tag_initial_range.sigma;
  _binding_site_tag_initial_range.lambda = parameters._binding_site_tag_initial_range.lambda;
  
  _co_enzyme_tag_initial_range.law    = parameters._co_enzyme_tag_initial_range.law;
  _co_enzyme_tag_initial_range.min    = parameters._co_enzyme_tag_initial_range.min;
  _co_enzyme_tag_initial_range.max    = parameters._co_enzyme_tag_initial_range.max;
  _co_enzyme_tag_initial_range.mu     = parameters._co_enzyme_tag_initial_range.mu;
  _co_enzyme_tag_initial_range.sigma  = parameters._co_enzyme_tag_initial_range.sigma;
  _co_enzyme_tag_initial_range.lambda = parameters._co_enzyme_tag_initial_range.lambda;
  
  _transcription_factor_tag_initial_range.law    = parameters._transcription_factor_tag_initial_range.law;
  _transcription_factor_tag_initial_range.min    = parameters._transcription_factor_tag_initial_range.min;
  _transcription_factor_tag_initial_range.max    = parameters._transcription_factor_tag_initial_range.max;
  _transcription_factor_tag_initial_range.mu     = parameters._transcription_factor_tag_initial_range.mu;
  _transcription_factor_tag_initial_range.sigma  = parameters._transcription_factor_tag_initial_range.sigma;
  _transcription_factor_tag_initial_range.lambda = parameters._transcription_factor_tag_initial_range.lambda;
  
  _initial_binding_window                 = parameters._initial_binding_window;
  _initial_number_of_NC_pearls            = parameters._initial_number_of_NC_pearls;
  _initial_number_of_E_pearls             = parameters._initial_number_of_E_pearls;
  _initial_number_of_TF_pearls            = parameters._initial_number_of_TF_pearls;
  _initial_number_of_BS_pearls            = parameters._initial_number_of_BS_pearls;
  _initial_number_of_P_pearls             = parameters._initial_number_of_P_pearls;
  _point_mutation_rate                    = parameters._point_mutation_rate;
  _duplication_rate                       = parameters._duplication_rate;
  _deletion_rate                          = parameters._deletion_rate;
  _translocation_rate                     = parameters._translocation_rate;
  _inversion_rate                         = parameters._inversion_rate;
  _transition_rate                        = parameters._transition_rate;
  _substrate_tag_mutation_size            = parameters._substrate_tag_mutation_size;
  _product_tag_mutation_size              = parameters._product_tag_mutation_size;
  _km_mutation_size                       = parameters._km_mutation_size;
  _kcat_mutation_size                     = parameters._kcat_mutation_size;
  _binding_size_tag_mutation_size         = parameters._binding_size_tag_mutation_size;
  _co_enzyme_tag_mutation_size            = parameters._co_enzyme_tag_mutation_size;
  _transcription_factor_tag_mutation_size = parameters._transcription_factor_tag_mutation_size;
  _basal_expression_level_mutation_size   = parameters._basal_expression_level_mutation_size;
  _maximum_genome_size                    = parameters._maximum_genome_size;
  
  /*------------------------------------------------------------------ genetic regulation network level */
  
  _genetic_regulation_network_timestep = parameters._genetic_regulation_network_timestep;
  _hill_function_theta                 = parameters._hill_function_theta;
  _hill_function_n                     = parameters._hill_function_n;
  _protein_degradation_rate            = parameters._protein_degradation_rate;
  
  /*------------------------------------------------------------------ metabolic network level */
  
  _metabolism_timestep                          = parameters._metabolism_timestep;
  _essential_metabolites_toxicity_threshold     = parameters._essential_metabolites_toxicity_threshold;
  _non_essential_metabolites_toxicity_threshold = parameters._non_essential_metabolites_toxicity_threshold;
  _initial_metabolites_amount_in_cells          = parameters._initial_metabolites_amount_in_cells;
  _energy_reaction_cost_factor                  = parameters._energy_reaction_cost_factor;
  _energy_pumping_cost                          = parameters._energy_pumping_cost;
  _energy_transcription_cost                    = parameters._energy_transcription_cost;
  _energy_toxicity_threshold                    = parameters._energy_toxicity_threshold;
  
  /*------------------------------------------------------------------ cell level */
  
  _death_probability             = parameters._death_probability;
  _variable_lifespan             = parameters._variable_lifespan;
  _gompert_law.a                 = parameters._gompert_law.a;
  _gompert_law.b                 = parameters._gompert_law.b;
  _gompert_law.c                 = parameters._gompert_law.c;
  _initial_membrane_permeability = parameters._initial_membrane_permeability;
  
  /*------------------------------------------------------------------ population level */
  
  _migration_rate = parameters._migration_rate;
  _hgt_rate       = parameters._hgt_rate;
  
  /*------------------------------------------------------------------ environment level */
  
  _environment_properties.number_of_init_cycles = parameters._environment_properties.number_of_init_cycles;
  
  _environment_properties.species_tag_range.law    = parameters._environment_properties.species_tag_range.law;
  _environment_properties.species_tag_range.min    = parameters._environment_properties.species_tag_range.min;
  _environment_properties.species_tag_range.max    = parameters._environment_properties.species_tag_range.max;
  _environment_properties.species_tag_range.mu     = parameters._environment_properties.species_tag_range.mu;
  _environment_properties.species_tag_range.sigma  = parameters._environment_properties.species_tag_range.sigma;
  _environment_properties.species_tag_range.lambda = parameters._environment_properties.species_tag_range.lambda;
  
  _environment_properties.concentration_range.law    = parameters._environment_properties.concentration_range.law;
  _environment_properties.concentration_range.min    = parameters._environment_properties.concentration_range.min;
  _environment_properties.concentration_range.max    = parameters._environment_properties.concentration_range.max;
  _environment_properties.concentration_range.mu     = parameters._environment_properties.concentration_range.mu;
  _environment_properties.concentration_range.sigma  = parameters._environment_properties.concentration_range.sigma;
  _environment_properties.concentration_range.lambda = parameters._environment_properties.concentration_range.lambda;
  
  _environment_properties.number_of_species_range.law    = parameters._environment_properties.number_of_species_range.law;
  _environment_properties.number_of_species_range.min    = parameters._environment_properties.number_of_species_range.min;
  _environment_properties.number_of_species_range.max    = parameters._environment_properties.number_of_species_range.max;
  _environment_properties.number_of_species_range.mu     = parameters._environment_properties.number_of_species_range.mu;
  _environment_properties.number_of_species_range.sigma  = parameters._environment_properties.number_of_species_range.sigma;
  _environment_properties.number_of_species_range.lambda = parameters._environment_properties.number_of_species_range.lambda;
  
  _environment_properties.interaction_scheme     = parameters._environment_properties.interaction_scheme;
  _environment_properties.renewal_scheme         = parameters._environment_properties.renewal_scheme;
  _environment_properties.variation_scheme       = parameters._environment_properties.variation_scheme;
  _environment_properties.variation_localization = parameters._environment_properties.variation_localization;
  
  _environment_properties.introduction_rate = parameters._environment_properties.introduction_rate;
  _environment_properties.diffusion_rate    = parameters._environment_properties.diffusion_rate;
  _environment_properties.degradation_rate  = parameters._environment_properties.degradation_rate;
  
  /*------------------------------------------------------------------ prime numbers */
  
  build_prime_numbers_list(5000);
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
Parameters::~Parameters( void )
{
  delete _prng;
  _prng = NULL;
  delete[] _prime_numbers;
  _prime_numbers = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Load parameters from a parameters file
 * \details  Parameters are loaded from a text file (default name: parameters.txt)
 * \param    std::string filename
 * \return   \e bool (succeed to load parameters)
 */
bool Parameters::load_parameters_from_file( std::string filename )
{
  std::ifstream file(filename.c_str(), std::ios::in);
  if(!file)
  {
    std::cout << "Error: " << filename << " file not found.\n";
    return false;
  }
  else
  {
    std::string line;
    std::string param_name;
    while(getline(file, line))
    {
      std::vector<std::string> words;
      if(parse_line(&words, line))
      {
        /*------------------------------------------------------------------ pseudorandom numbers generator */
        
        if ( strcmp(words[0].c_str(), "SEED") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _seed;
        }
        
        /*------------------------------------------------------------------ parallel computing */
        
        else if ( strcmp(words[0].c_str(), "PARALLEL_COMPUTING") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _parallel_computing = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _parallel_computing = true;
          }
          else
          {
            std::cout << "Error : PARALLEL_COMPUTING wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        
        /*------------------------------------------------------------------ simulation schemes */
        
        else if ( strcmp(words[0].c_str(), "ENERGY_CONSTRAINTS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _energy_constraints = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _energy_constraints = true;
          }
          else
          {
            std::cout << "Error : ENERGY_CONSTRAINTS wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "MEMBRANE_PERMEABILITY") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _membrane_permeability = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _membrane_permeability = true;
          }
          else
          {
            std::cout << "Error : MEMBRANE_PERMEABILITY wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "METABOLIC_INHERITANCE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _metabolic_inheritance = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _metabolic_inheritance = true;
          }
          else
          {
            std::cout << "Error : METABOLIC_INHERITANCE wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENZYMATIC_INHERITANCE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _enzymatic_inheritance = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _enzymatic_inheritance = true;
          }
          else
          {
            std::cout << "Error : ENZYMATIC_INHERITANCE wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "CO_ENZYME_ACTIVITY") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string state;
          flux >> param_name >> state;
          if (strcmp(state.c_str(), "NO") == 0 || strcmp(state.c_str(), "no") == 0)
          {
            _co_enzyme_activity = false;
          }
          else if (strcmp(state.c_str(), "YES") == 0 || strcmp(state.c_str(), "yes") == 0)
          {
            _co_enzyme_activity = true;
          }
          else
          {
            std::cout << "Error : CO_ENZYME_ACTIVITY wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "SCORE_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string law;
          flux >> param_name >> law;
          if (strcmp(law.c_str(), "SUM") == 0)
          {
            _score_scheme = ESSENTIAL_METABOLITES_SUM;
          }
          else if (strcmp(law.c_str(), "SUM_MINUS_DEV") == 0)
          {
            _score_scheme = ESSENTIAL_METABOLITES_SUM_MINUS_DEVIATION;
          }
          else if (strcmp(law.c_str(), "COMBINATORIAL") == 0)
          {
            _score_scheme = ESSENTIAL_METABOLITES_COMBINATORIAL_CONTRIBUTION;
          }
          else
          {
            std::cout << "Error : SCORE_SCHEME wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "SELECTION_THRESHOLD") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _selection_threshold;
          assert(_selection_threshold >= 0.0);
          assert(_selection_threshold <= 1.0);
        }
        
        /*------------------------------------------------------------------ space */
        
        else if ( strcmp(words[0].c_str(), "WIDTH") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _width;
          assert(_width > 0);
        }
        else if ( strcmp(words[0].c_str(), "HEIGHT") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _height;
          assert(_height > 0);
        }
        
        /*------------------------------------------------------------------ output */
        
        else if ( strcmp(words[0].c_str(), "EXPERIMENT_BACKUP_STEP") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _experiment_backup_step;
        }
        else if ( strcmp(words[0].c_str(), "TREE_BACKUP_STEP") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _tree_backup_step;
        }
        else if ( strcmp(words[0].c_str(), "FIGURES_GENERATION_STEP") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _figures_generation_step;
        }
        
        /*------------------------------------------------------------------ genome level */
        
        else if ( strcmp(words[0].c_str(), "METABOLITE_TAG_INITIAL_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string law;
          flux >> param_name >> law;
          if (strcmp(law.c_str(), "UNIFORM") == 0)
          {
            _metabolite_tag_initial_range.law = UNIFORM;
            flux >> _metabolite_tag_initial_range.min >> _metabolite_tag_initial_range.max;
            _metabolite_tag_initial_range.mu     = 0.0;
            _metabolite_tag_initial_range.sigma  = 0.0;
            _metabolite_tag_initial_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "GAUSSIAN") == 0)
          {
            _metabolite_tag_initial_range.law = GAUSSIAN;
            flux >> _metabolite_tag_initial_range.mu >> _metabolite_tag_initial_range.sigma;
            _metabolite_tag_initial_range.min    = 0;
            _metabolite_tag_initial_range.max    = 0;
            _metabolite_tag_initial_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "EXPONENTIAL") == 0)
          {
            _metabolite_tag_initial_range.law = EXPONENTIAL;
            flux >> _metabolite_tag_initial_range.lambda;
            _metabolite_tag_initial_range.min   = 0;
            _metabolite_tag_initial_range.max   = 0;
            _metabolite_tag_initial_range.mu    = 0.0;
            _metabolite_tag_initial_range.sigma = 0.0;
          }
          else
          {
            std::cout << "Error : METABOLITE_TAG_INITIAL_RANGE wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "BINDING_SIZE_TAG_INITIAL_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string law;
          flux >> param_name >> law;
          if (strcmp(law.c_str(), "UNIFORM") == 0)
          {
            _binding_site_tag_initial_range.law = UNIFORM;
            flux >> _binding_site_tag_initial_range.min >> _binding_site_tag_initial_range.max;
            _binding_site_tag_initial_range.mu     = 0.0;
            _binding_site_tag_initial_range.sigma  = 0.0;
            _binding_site_tag_initial_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "GAUSSIAN") == 0)
          {
            _binding_site_tag_initial_range.law = GAUSSIAN;
            flux >> _binding_site_tag_initial_range.mu >> _binding_site_tag_initial_range.sigma;
            _binding_site_tag_initial_range.min    = 0;
            _binding_site_tag_initial_range.max    = 0;
            _binding_site_tag_initial_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "EXPONENTIAL") == 0)
          {
            _binding_site_tag_initial_range.law = EXPONENTIAL;
            flux >> _binding_site_tag_initial_range.lambda;
            _binding_site_tag_initial_range.min   = 0;
            _binding_site_tag_initial_range.max   = 0;
            _binding_site_tag_initial_range.mu    = 0.0;
            _binding_site_tag_initial_range.sigma = 0.0;
          }
          else
          {
            std::cout << "Error : BINDING_SIZE_TAG_INITIAL_RANGE wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "CO_ENZYME_TAG_INITIAL_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string law;
          flux >> param_name >> law;
          if (strcmp(law.c_str(), "UNIFORM") == 0)
          {
            _co_enzyme_tag_initial_range.law = UNIFORM;
            flux >> _co_enzyme_tag_initial_range.min >> _co_enzyme_tag_initial_range.max;
            _co_enzyme_tag_initial_range.mu     = 0.0;
            _co_enzyme_tag_initial_range.sigma  = 0.0;
            _co_enzyme_tag_initial_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "GAUSSIAN") == 0)
          {
            _co_enzyme_tag_initial_range.law = GAUSSIAN;
            flux >> _co_enzyme_tag_initial_range.mu >> _co_enzyme_tag_initial_range.sigma;
            _co_enzyme_tag_initial_range.min    = 0;
            _co_enzyme_tag_initial_range.max    = 0;
            _co_enzyme_tag_initial_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "EXPONENTIAL") == 0)
          {
            _co_enzyme_tag_initial_range.law = EXPONENTIAL;
            flux >> _co_enzyme_tag_initial_range.lambda;
            _co_enzyme_tag_initial_range.min   = 0;
            _co_enzyme_tag_initial_range.max   = 0;
            _co_enzyme_tag_initial_range.mu    = 0.0;
            _co_enzyme_tag_initial_range.sigma = 0.0;
          }
          else
          {
            std::cout << "Error : CO_ENZYME_TAG_INITIAL_RANGE wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "TRANSCRIPTION_FACTOR_TAG_INITIAL_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string law;
          flux >> param_name >> law;
          if (strcmp(law.c_str(), "UNIFORM") == 0)
          {
            _transcription_factor_tag_initial_range.law = UNIFORM;
            flux >> _transcription_factor_tag_initial_range.min >> _transcription_factor_tag_initial_range.max;
            _transcription_factor_tag_initial_range.mu     = 0.0;
            _transcription_factor_tag_initial_range.sigma  = 0.0;
            _transcription_factor_tag_initial_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "GAUSSIAN") == 0)
          {
            _transcription_factor_tag_initial_range.law = GAUSSIAN;
            flux >> _transcription_factor_tag_initial_range.mu >> _transcription_factor_tag_initial_range.sigma;
            _transcription_factor_tag_initial_range.min    = 0;
            _transcription_factor_tag_initial_range.max    = 0;
            _transcription_factor_tag_initial_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "EXPONENTIAL") == 0)
          {
            _transcription_factor_tag_initial_range.law = EXPONENTIAL;
            flux >> _transcription_factor_tag_initial_range.lambda;
            _transcription_factor_tag_initial_range.min   = 0;
            _transcription_factor_tag_initial_range.max   = 0;
            _transcription_factor_tag_initial_range.mu    = 0.0;
            _transcription_factor_tag_initial_range.sigma = 0.0;
          }
          else
          {
            std::cout << "Error : TRANSCRIPTION_FACTOR_TAG_INITIAL_RANGE wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_BINDING_WINDOW") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_binding_window;
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_NUMBER_OF_NON_CODING_PEARLS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_number_of_NC_pearls;
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_NUMBER_OF_ENZYME_PEARLS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_number_of_E_pearls;
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_NUMBER_OF_TRANSCRIPTION_FACTOR_PEARLS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_number_of_TF_pearls;
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_NUMBER_OF_BINDING_SITE_PEARLS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_number_of_BS_pearls;
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_NUMBER_OF_PROMOTER_PEARLS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_number_of_P_pearls;
        }
        else if ( strcmp(words[0].c_str(), "POINT_MUTATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _point_mutation_rate;
          assert(_point_mutation_rate >= 0.0);
          assert(_point_mutation_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "DUPLICATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _duplication_rate;
          assert(_duplication_rate >= 0.0);
          assert(_duplication_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "DELETION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _deletion_rate;
          assert(_deletion_rate >= 0.0);
          assert(_deletion_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "TRANSLOCATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _translocation_rate;
          assert(_translocation_rate >= 0.0);
          assert(_translocation_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "INVERSION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _inversion_rate;
          assert(_inversion_rate >= 0.0);
          assert(_inversion_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "TRANSITION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _transition_rate;
          assert(_transition_rate >= 0.0);
          assert(_transition_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "SUBSTRATE_TAG_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _substrate_tag_mutation_size;
          assert(_substrate_tag_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "PRODUCT_TAG_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _product_tag_mutation_size;
          assert(_product_tag_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "KM_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _km_mutation_size;
          assert(_km_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "KCAT_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _kcat_mutation_size;
          assert(_kcat_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "BINDING_SIZE_TAG_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _binding_size_tag_mutation_size;
          assert(_binding_size_tag_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "CO_ENZYME_TAG_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _co_enzyme_tag_mutation_size;
          assert(_co_enzyme_tag_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _transcription_factor_tag_mutation_size;
          assert(_transcription_factor_tag_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "BASAL_EXPRESSION_LEVEL_MUTATION_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _basal_expression_level_mutation_size;
          assert(_basal_expression_level_mutation_size >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "MAXIMUM_GENOME_SIZE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _maximum_genome_size;
          assert(_maximum_genome_size > 0);
        }
        
        /*------------------------------------------------------------------ genetic regulation network level */
        
        else if ( strcmp(words[0].c_str(), "GENETIC_REGULATION_NETWORK_TIMESTEP") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _genetic_regulation_network_timestep;
          assert(_genetic_regulation_network_timestep > 0.0);
        }
        else if ( strcmp(words[0].c_str(), "HILL_FUNCTION_THETA") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _hill_function_theta;
          assert(_hill_function_theta >= 0.0);
          assert(_hill_function_theta <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "HILL_FUNCTION_N") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _hill_function_n;
          assert(_hill_function_n >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "PROTEIN_DEGRADATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _protein_degradation_rate;
          assert(_protein_degradation_rate >= 0.0);
          assert(_protein_degradation_rate <= 1.0);
        }
        
        /*------------------------------------------------------------------ metabolic network level */
        
        else if ( strcmp(words[0].c_str(), "METABOLISM_TIMESTEP") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _metabolism_timestep;
          assert(_metabolism_timestep > 0.0);
        }
        else if ( strcmp(words[0].c_str(), "ESSENTIAL_METABOLITES_TOXICITY_THRESHOLD") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _essential_metabolites_toxicity_threshold;
          assert(_essential_metabolites_toxicity_threshold > 0.0);
        }
        else if ( strcmp(words[0].c_str(), "NON_ESSENTIAL_METABOLITES_TOXICITY_THRESHOLD") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _non_essential_metabolites_toxicity_threshold;
          assert(_non_essential_metabolites_toxicity_threshold > 0.0);
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_METABOLITES_AMOUNT_IN_CELLS") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_metabolites_amount_in_cells;
          assert(_initial_metabolites_amount_in_cells >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "ENERGY_REACTION_COST_FACTOR") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _energy_reaction_cost_factor;
          assert(_energy_reaction_cost_factor >= 0.0);
          assert(_energy_reaction_cost_factor <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "ENERGY_PUMPING_COST") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _energy_pumping_cost;
        }
        else if ( strcmp(words[0].c_str(), "ENERGY_TRANSCRIPTION_COST") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _energy_transcription_cost;
        }
        else if ( strcmp(words[0].c_str(), "ENERGY_TOXICITY_THRESHOLD") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _energy_toxicity_threshold;
        }
        
        /*------------------------------------------------------------------ cell level */
        
        else if ( strcmp(words[0].c_str(), "DEATH_PROBABILITY") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _death_probability;
          assert(_death_probability >= 0.0);
          assert(_death_probability <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "VARIABLE_LIFESPAN") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string choice;
          flux >> param_name >> choice;
          if (strcmp(choice.c_str(), "NO") == 0 || strcmp(choice.c_str(), "no") == 0)
          {
            _variable_lifespan = false;
          }
          else if (strcmp(choice.c_str(), "YES") == 0 || strcmp(choice.c_str(), "yes") == 0)
          {
            _variable_lifespan = true;
          }
          else
          {
            std::cout << "Error : VARIABLE_LIFESPAN wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "GOMPERTZ_LAW") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _gompert_law.a >> _gompert_law.b >> _gompert_law.c;
        }
        else if ( strcmp(words[0].c_str(), "INITIAL_MEMBRANE_PERMEABILITY") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _initial_membrane_permeability;
          assert(_initial_membrane_permeability >= 0.0);
          assert(_initial_membrane_permeability <= 1.0);
        }
        
        /*------------------------------------------------------------------ population level */
        
        else if ( strcmp(words[0].c_str(), "MIGRATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _migration_rate;
          assert(_migration_rate >= 0.0);
          assert(_migration_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "HGT_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _hgt_rate;
          assert(_hgt_rate >= 0.0);
          assert(_hgt_rate <= 1.0);
        }
        
        /*------------------------------------------------------------------ environment level */
        
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_INITIALIZATION_CYCLES") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _environment_properties.number_of_init_cycles;
          assert(_environment_properties.number_of_init_cycles >= 1);
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_SPECIES_TAG_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string law;
          flux >> param_name >> law;
          if (strcmp(law.c_str(), "UNIFORM") == 0)
          {
            _environment_properties.species_tag_range.law = UNIFORM;
            flux >> _environment_properties.species_tag_range.min >> _environment_properties.species_tag_range.max;
            _environment_properties.species_tag_range.mu     = 0.0;
            _environment_properties.species_tag_range.sigma  = 0.0;
            _environment_properties.species_tag_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "GAUSSIAN") == 0)
          {
            _environment_properties.species_tag_range.law = GAUSSIAN;
            flux >> _environment_properties.species_tag_range.mu >> _environment_properties.species_tag_range.sigma;
            _environment_properties.species_tag_range.min    = 0;
            _environment_properties.species_tag_range.max    = 0;
            _environment_properties.species_tag_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "EXPONENTIAL") == 0)
          {
            _environment_properties.species_tag_range.law = EXPONENTIAL;
            flux >> _environment_properties.species_tag_range.lambda;
            _environment_properties.species_tag_range.min   = 0;
            _environment_properties.species_tag_range.max   = 0;
            _environment_properties.species_tag_range.mu    = 0.0;
            _environment_properties.species_tag_range.sigma = 0.0;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_SPECIES_TAG_RANGE wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_CONCENTRATION_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string law;
          flux >> param_name >> law;
          if (strcmp(law.c_str(), "UNIFORM") == 0)
          {
            _environment_properties.concentration_range.law = UNIFORM;
            flux >> _environment_properties.concentration_range.min >> _environment_properties.concentration_range.max;
            _environment_properties.concentration_range.mu     = 0.0;
            _environment_properties.concentration_range.sigma  = 0.0;
            _environment_properties.concentration_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "GAUSSIAN") == 0)
          {
            _environment_properties.concentration_range.law = GAUSSIAN;
            flux >> _environment_properties.concentration_range.mu >> _environment_properties.concentration_range.sigma;
            _environment_properties.concentration_range.min    = 0;
            _environment_properties.concentration_range.max    = 0;
            _environment_properties.concentration_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "EXPONENTIAL") == 0)
          {
            _environment_properties.concentration_range.law = EXPONENTIAL;
            flux >> _environment_properties.concentration_range.lambda;
            _environment_properties.concentration_range.min   = 0;
            _environment_properties.concentration_range.max   = 0;
            _environment_properties.concentration_range.mu    = 0.0;
            _environment_properties.concentration_range.sigma = 0.0;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_CONCENTRATION_RANGE wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_NUMBER_OF_SPECIES_RANGE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string law;
          flux >> param_name >> law;
          if (strcmp(law.c_str(), "UNIFORM") == 0)
          {
            _environment_properties.number_of_species_range.law = UNIFORM;
            flux >> _environment_properties.number_of_species_range.min >> _environment_properties.number_of_species_range.max;
            _environment_properties.number_of_species_range.mu     = 0.0;
            _environment_properties.number_of_species_range.sigma  = 0.0;
            _environment_properties.number_of_species_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "GAUSSIAN") == 0)
          {
            _environment_properties.number_of_species_range.law = GAUSSIAN;
            flux >> _environment_properties.number_of_species_range.mu >> _environment_properties.number_of_species_range.sigma;
            _environment_properties.number_of_species_range.min    = 0;
            _environment_properties.number_of_species_range.max    = 0;
            _environment_properties.number_of_species_range.lambda = 0.0;
          }
          else if (strcmp(law.c_str(), "EXPONENTIAL") == 0)
          {
            _environment_properties.number_of_species_range.law = EXPONENTIAL;
            flux >> _environment_properties.number_of_species_range.lambda;
            _environment_properties.number_of_species_range.min   = 0;
            _environment_properties.number_of_species_range.max   = 0;
            _environment_properties.number_of_species_range.mu    = 0.0;
            _environment_properties.number_of_species_range.sigma = 0.0;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_NUMBER_OF_SPECIES_RANGE wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_INTERACTION_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string scheme;
          flux >> param_name >> scheme;
          if (strcmp(scheme.c_str(), "NO_INTERACTION") == 0)
          {
            _environment_properties.interaction_scheme = NO_INTERACTION;
          }
          else if (strcmp(scheme.c_str(), "INTERACTION") == 0)
          {
            _environment_properties.interaction_scheme = INTERACTION;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_INTERACTION_SCHEME wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_RENEWAL_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string scheme;
          flux >> param_name >> scheme;
          if (strcmp(scheme.c_str(), "KEEP_MATTER") == 0)
          {
            _environment_properties.renewal_scheme = KEEP_MATTER;
          }
          else if (strcmp(scheme.c_str(), "CLEAR_MATTER") == 0)
          {
            _environment_properties.renewal_scheme = CLEAR_MATTER;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_RENEWAL_SCHEME wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_VARIATION_SCHEME") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string scheme;
          flux >> param_name >> scheme;
          if (strcmp(scheme.c_str(), "RANDOM") == 0)
          {
            _environment_properties.variation_scheme = RANDOM_SCHEME;
          }
          else if (strcmp(scheme.c_str(), "PERIODIC") == 0)
          {
            _environment_properties.variation_scheme = PERIODIC_SCHEME;
          }
          else if (strcmp(scheme.c_str(), "CYCLIC") == 0)
          {
            _environment_properties.variation_scheme = CYCLIC_SCHEME;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_VARIATION_SCHEME wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_VARIATION_LOCALIZATION") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          std::string localization;
          flux >> param_name >> localization;
          if (strcmp(localization.c_str(), "RANDOM") == 0)
          {
            _environment_properties.variation_localization = RANDOM_LOCALIZATION;
          }
          else if (strcmp(localization.c_str(), "GLOBAL") == 0)
          {
            _environment_properties.variation_localization = GLOBAL_LOCALIZATION;
          }
          else if (strcmp(localization.c_str(), "CENTER") == 0)
          {
            _environment_properties.variation_localization = CENTER_LOCALIZATION;
          }
          else
          {
            std::cout << "Error : ENVIRONMENT_VARIATION_LOCALIZATION wrong value at line:\n " << line.c_str() << ".\n";
            exit(EXIT_FAILURE);
          }
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_INTRODUCTION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _environment_properties.introduction_rate;
          assert(_environment_properties.introduction_rate >= 0.0);
          assert(_environment_properties.introduction_rate <= 1.0);
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_DIFFUSION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _environment_properties.diffusion_rate;
          assert(_environment_properties.diffusion_rate >= 0.0);
        }
        else if ( strcmp(words[0].c_str(), "ENVIRONMENT_DEGRADATION_RATE") == 0 )
        {
          std::stringstream flux;
          flux.str(line.c_str());
          flux >> param_name >> _environment_properties.degradation_rate;
          assert(_environment_properties.degradation_rate >= 0.0);
        }
        else
        {
          std::cout << "Unknown parameter '" << words[0] << "'. Exit.\n";
          exit(EXIT_FAILURE);
        }
      }
    }
    file.close();
  }
  return true;
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    size_t backup_time
 * \return   \e void
 */
void Parameters::save( size_t backup_time )
{
  std::stringstream backup_file_name;
  backup_file_name << "./parameters/parameters_" << backup_time;
  gzFile backup_file = gzopen(backup_file_name.str().c_str(), "w");
  
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  _prng->save( backup_time );
  gzwrite( backup_file, &_seed, sizeof(_seed) );
  
  /*------------------------------------------------------------------ parallel computing */
  
  gzwrite( backup_file, &_parallel_computing, sizeof(_parallel_computing) );
  
  /*------------------------------------------------------------------ simulation schemes */
  
  gzwrite( backup_file, &_energy_constraints,    sizeof(_energy_constraints) );
  gzwrite( backup_file, &_membrane_permeability, sizeof(_membrane_permeability) );
  gzwrite( backup_file, &_metabolic_inheritance, sizeof(_metabolic_inheritance) );
  gzwrite( backup_file, &_enzymatic_inheritance, sizeof(_enzymatic_inheritance) );
  gzwrite( backup_file, &_co_enzyme_activity,    sizeof(_co_enzyme_activity) );
  gzwrite( backup_file, &_score_scheme,          sizeof(_score_scheme) );
  gzwrite( backup_file, &_selection_threshold,   sizeof(_selection_threshold) );
  
  /*------------------------------------------------------------------ space */
  
  gzwrite( backup_file, &_width,  sizeof(_width) );
  gzwrite( backup_file, &_height, sizeof(_height) );
  
  /*------------------------------------------------------------------ output */
  
  gzwrite( backup_file, &_experiment_backup_step,  sizeof(_experiment_backup_step) );
  gzwrite( backup_file, &_tree_backup_step,        sizeof(_tree_backup_step) );
  gzwrite( backup_file, &_figures_generation_step, sizeof(_figures_generation_step) );
  
  /*------------------------------------------------------------------ genome level */
  
  gzwrite( backup_file, &_metabolite_tag_initial_range.law,    sizeof(_metabolite_tag_initial_range.law) );
  gzwrite( backup_file, &_metabolite_tag_initial_range.min,    sizeof(_metabolite_tag_initial_range.min) );
  gzwrite( backup_file, &_metabolite_tag_initial_range.max,    sizeof(_metabolite_tag_initial_range.max) );
  gzwrite( backup_file, &_metabolite_tag_initial_range.mu,     sizeof(_metabolite_tag_initial_range.mu) );
  gzwrite( backup_file, &_metabolite_tag_initial_range.sigma,  sizeof(_metabolite_tag_initial_range.sigma) );
  gzwrite( backup_file, &_metabolite_tag_initial_range.lambda, sizeof(_metabolite_tag_initial_range.lambda) );
  
  gzwrite( backup_file, &_binding_site_tag_initial_range.law,    sizeof(_binding_site_tag_initial_range.law) );
  gzwrite( backup_file, &_binding_site_tag_initial_range.min,    sizeof(_binding_site_tag_initial_range.min) );
  gzwrite( backup_file, &_binding_site_tag_initial_range.max,    sizeof(_binding_site_tag_initial_range.max) );
  gzwrite( backup_file, &_binding_site_tag_initial_range.mu,     sizeof(_binding_site_tag_initial_range.mu) );
  gzwrite( backup_file, &_binding_site_tag_initial_range.sigma,  sizeof(_binding_site_tag_initial_range.sigma) );
  gzwrite( backup_file, &_binding_site_tag_initial_range.lambda, sizeof(_binding_site_tag_initial_range.lambda) );
  
  gzwrite( backup_file, &_co_enzyme_tag_initial_range.law,    sizeof(_co_enzyme_tag_initial_range.law) );
  gzwrite( backup_file, &_co_enzyme_tag_initial_range.min,    sizeof(_co_enzyme_tag_initial_range.min) );
  gzwrite( backup_file, &_co_enzyme_tag_initial_range.max,    sizeof(_co_enzyme_tag_initial_range.max) );
  gzwrite( backup_file, &_co_enzyme_tag_initial_range.mu,     sizeof(_co_enzyme_tag_initial_range.mu) );
  gzwrite( backup_file, &_co_enzyme_tag_initial_range.sigma,  sizeof(_co_enzyme_tag_initial_range.sigma) );
  gzwrite( backup_file, &_co_enzyme_tag_initial_range.lambda, sizeof(_co_enzyme_tag_initial_range.lambda) );
  
  gzwrite( backup_file, &_transcription_factor_tag_initial_range.law,    sizeof(_transcription_factor_tag_initial_range.law) );
  gzwrite( backup_file, &_transcription_factor_tag_initial_range.min,    sizeof(_transcription_factor_tag_initial_range.min) );
  gzwrite( backup_file, &_transcription_factor_tag_initial_range.max,    sizeof(_transcription_factor_tag_initial_range.max) );
  gzwrite( backup_file, &_transcription_factor_tag_initial_range.mu,     sizeof(_transcription_factor_tag_initial_range.mu) );
  gzwrite( backup_file, &_transcription_factor_tag_initial_range.sigma,  sizeof(_transcription_factor_tag_initial_range.sigma) );
  gzwrite( backup_file, &_transcription_factor_tag_initial_range.lambda, sizeof(_transcription_factor_tag_initial_range.lambda) );
  
  gzwrite( backup_file, &_initial_binding_window,                 sizeof(_initial_binding_window) );
  gzwrite( backup_file, &_initial_number_of_NC_pearls,            sizeof(_initial_number_of_NC_pearls) );
  gzwrite( backup_file, &_initial_number_of_E_pearls,             sizeof(_initial_number_of_E_pearls) );
  gzwrite( backup_file, &_initial_number_of_TF_pearls,            sizeof(_initial_number_of_TF_pearls) );
  gzwrite( backup_file, &_initial_number_of_BS_pearls,            sizeof(_initial_number_of_BS_pearls) );
  gzwrite( backup_file, &_initial_number_of_P_pearls,             sizeof(_initial_number_of_P_pearls) );
  gzwrite( backup_file, &_point_mutation_rate,                    sizeof(_point_mutation_rate) );
  gzwrite( backup_file, &_duplication_rate,                       sizeof(_duplication_rate) );
  gzwrite( backup_file, &_deletion_rate,                          sizeof(_deletion_rate) );
  gzwrite( backup_file, &_translocation_rate,                     sizeof(_translocation_rate) );
  gzwrite( backup_file, &_inversion_rate,                         sizeof(_inversion_rate) );
  gzwrite( backup_file, &_transition_rate,                        sizeof(_transition_rate) );
  gzwrite( backup_file, &_substrate_tag_mutation_size,            sizeof(_substrate_tag_mutation_size) );
  gzwrite( backup_file, &_product_tag_mutation_size,              sizeof(_product_tag_mutation_size) );
  gzwrite( backup_file, &_km_mutation_size,                       sizeof(_km_mutation_size) );
  gzwrite( backup_file, &_kcat_mutation_size,                     sizeof(_kcat_mutation_size) );
  gzwrite( backup_file, &_binding_size_tag_mutation_size,         sizeof(_binding_size_tag_mutation_size) );
  gzwrite( backup_file, &_co_enzyme_tag_mutation_size,            sizeof(_co_enzyme_tag_mutation_size) );
  gzwrite( backup_file, &_transcription_factor_tag_mutation_size, sizeof(_transcription_factor_tag_mutation_size) );
  gzwrite( backup_file, &_basal_expression_level_mutation_size,   sizeof(_basal_expression_level_mutation_size) );
  gzwrite( backup_file, &_maximum_genome_size,                    sizeof(_maximum_genome_size) );
  
  /*------------------------------------------------------------------ genetic regulation network level */
  
  gzwrite( backup_file, &_genetic_regulation_network_timestep, sizeof(_genetic_regulation_network_timestep) );
  gzwrite( backup_file, &_hill_function_theta,                 sizeof(_hill_function_theta) );
  gzwrite( backup_file, &_hill_function_n,                     sizeof(_hill_function_n) );
  gzwrite( backup_file, &_protein_degradation_rate,            sizeof(_protein_degradation_rate) );
  
  /*------------------------------------------------------------------ metabolic network level */
  
  gzwrite( backup_file, &_metabolism_timestep,                          sizeof(_metabolism_timestep) );
  gzwrite( backup_file, &_essential_metabolites_toxicity_threshold,     sizeof(_essential_metabolites_toxicity_threshold) );
  gzwrite( backup_file, &_non_essential_metabolites_toxicity_threshold, sizeof(_non_essential_metabolites_toxicity_threshold) );
  
  gzwrite( backup_file, &_initial_metabolites_amount_in_cells,          sizeof(_initial_metabolites_amount_in_cells) );
  
  gzwrite( backup_file, &_energy_reaction_cost_factor,                  sizeof(_energy_reaction_cost_factor) );
  gzwrite( backup_file, &_energy_pumping_cost,                          sizeof(_energy_pumping_cost) );
  gzwrite( backup_file, &_energy_transcription_cost,                    sizeof(_energy_transcription_cost) );
  gzwrite( backup_file, &_energy_toxicity_threshold,                    sizeof(_energy_toxicity_threshold) );
  
  /*------------------------------------------------------------------ cell level */
  
  gzwrite( backup_file, &_death_probability,             sizeof(_death_probability) );
  gzwrite( backup_file, &_variable_lifespan,             sizeof(_variable_lifespan) );
  gzwrite( backup_file, &_gompert_law.a,                 sizeof(_gompert_law.a) );
  gzwrite( backup_file, &_gompert_law.b,                 sizeof(_gompert_law.b) );
  gzwrite( backup_file, &_gompert_law.c,                 sizeof(_gompert_law.c) );
  gzwrite( backup_file, &_initial_membrane_permeability, sizeof(_initial_membrane_permeability) );
  
  /*------------------------------------------------------------------ population level */
  
  gzwrite( backup_file, &_migration_rate, sizeof(_migration_rate) );
  gzwrite( backup_file, &_hgt_rate,       sizeof(_hgt_rate) );
  
  /*------------------------------------------------------------------ environment level */
  
  gzwrite( backup_file, &_environment_properties.number_of_init_cycles,    sizeof(_environment_properties.number_of_init_cycles) );
  
  gzwrite( backup_file, &_environment_properties.species_tag_range.law,    sizeof(_environment_properties.species_tag_range.law) );
  gzwrite( backup_file, &_environment_properties.species_tag_range.min,    sizeof(_environment_properties.species_tag_range.min) );
  gzwrite( backup_file, &_environment_properties.species_tag_range.max,    sizeof(_environment_properties.species_tag_range.max) );
  gzwrite( backup_file, &_environment_properties.species_tag_range.mu,     sizeof(_environment_properties.species_tag_range.mu) );
  gzwrite( backup_file, &_environment_properties.species_tag_range.sigma,  sizeof(_environment_properties.species_tag_range.sigma) );
  gzwrite( backup_file, &_environment_properties.species_tag_range.lambda, sizeof(_environment_properties.species_tag_range.lambda) );
  
  gzwrite( backup_file, &_environment_properties.concentration_range.law,    sizeof(_environment_properties.concentration_range.law) );
  gzwrite( backup_file, &_environment_properties.concentration_range.min,    sizeof(_environment_properties.concentration_range.min) );
  gzwrite( backup_file, &_environment_properties.concentration_range.max,    sizeof(_environment_properties.concentration_range.max) );
  gzwrite( backup_file, &_environment_properties.concentration_range.mu,     sizeof(_environment_properties.concentration_range.mu) );
  gzwrite( backup_file, &_environment_properties.concentration_range.sigma,  sizeof(_environment_properties.concentration_range.sigma) );
  gzwrite( backup_file, &_environment_properties.concentration_range.lambda, sizeof(_environment_properties.concentration_range.lambda) );
  
  gzwrite( backup_file, &_environment_properties.number_of_species_range.law,    sizeof(_environment_properties.number_of_species_range.law) );
  gzwrite( backup_file, &_environment_properties.number_of_species_range.min,    sizeof(_environment_properties.number_of_species_range.min) );
  gzwrite( backup_file, &_environment_properties.number_of_species_range.max,    sizeof(_environment_properties.number_of_species_range.max) );
  gzwrite( backup_file, &_environment_properties.number_of_species_range.mu,     sizeof(_environment_properties.number_of_species_range.mu) );
  gzwrite( backup_file, &_environment_properties.number_of_species_range.sigma,  sizeof(_environment_properties.number_of_species_range.sigma) );
  gzwrite( backup_file, &_environment_properties.number_of_species_range.lambda, sizeof(_environment_properties.number_of_species_range.lambda) );
  
  gzwrite( backup_file, &_environment_properties.interaction_scheme,     sizeof(_environment_properties.interaction_scheme) );
  gzwrite( backup_file, &_environment_properties.renewal_scheme,         sizeof(_environment_properties.renewal_scheme) );
  gzwrite( backup_file, &_environment_properties.variation_scheme,       sizeof(_environment_properties.variation_scheme) );
  gzwrite( backup_file, &_environment_properties.variation_localization, sizeof(_environment_properties.variation_localization) );
  
  gzwrite( backup_file, &_environment_properties.introduction_rate, sizeof(_environment_properties.introduction_rate) );
  gzwrite( backup_file, &_environment_properties.diffusion_rate,    sizeof(_environment_properties.diffusion_rate) );
  gzwrite( backup_file, &_environment_properties.degradation_rate,  sizeof(_environment_properties.degradation_rate) );
  
  gzclose(backup_file);
  
}

/**
 * \brief    Write parameters in a parameter file
 * \details  --
 * \param    std::string filename
 * \return   \e void
 */
void Parameters::write( std::string filename )
{
  std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
  file << "\n";
  file << "########################################################\n";
  file << "# PSEUDORANDOM NUMBERS GENERATOR\n";
  file << "########################################################\n";
  file << "\n";
  file << "SEED  " << _seed << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# PARALLEL COMPUTING\n";
  file << "########################################################\n";
  file << "\n";
  file << "PARALLEL_COMPUTING  ";
  if (!_parallel_computing)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# SIMULATION SCHEMES\n";
  file << "########################################################\n";
  file << "\n";
  
  file << "ENERGY_CONSTRAINTS          ";
  if (!_energy_constraints)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  
  file << "MEMBRANE_PERMEABILITY       ";
  if (!_membrane_permeability)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  
  file << "METABOLIC_INHERITANCE       ";
  if (!_metabolic_inheritance)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  
  file << "ENZYMATIC_INHERITANCE       ";
  if (!_enzymatic_inheritance)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  
  file << "CO_ENZYME_ACTIVITY          ";
  if (!_co_enzyme_activity)
  {
    file << "NO ";
  }
  else
  {
    file << "YES";
  }
  file << " (YES, NO)\n";
  
  file << "\n";
  file << "SCORE_SCHEME                ";
  switch (_score_scheme)
  {
    case ESSENTIAL_METABOLITES_SUM:
      file << "SUM (SUM, SUM_MINUS_DEV, COMBINATORIAL)\n";
      break;
    case ESSENTIAL_METABOLITES_SUM_MINUS_DEVIATION:
      file << "SUM_MINUS_DEV (SUM, SUM_MINUS_DEV, COMBINATORIAL)\n";
      break;
    case ESSENTIAL_METABOLITES_COMBINATORIAL_CONTRIBUTION:
      file << "COMBINATORIAL (SUM, SUM_MINUS_DEV, COMBINATORIAL)\n";
      break;
  }
  file << "SELECTION_THRESHOLD         " << _selection_threshold << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# SPACE\n";
  file << "########################################################\n";
  file << "\n";
  file << "WIDTH   " << _width << "\n";
  file << "HEIGHT  " << _height << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# OUTPUT\n";
  file << "########################################################\n";
  file << "\n";
  file << "EXPERIMENT_BACKUP_STEP   " << _experiment_backup_step << "\n";
  file << "TREE_BACKUP_STEP         " << _tree_backup_step << "\n";
  file << "FIGURES_GENERATION_STEP  " << _figures_generation_step << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# GENOME LEVEL\n";
  file << "########################################################\n";
  file << "\n";
  file << "METABOLITE_TAG_INITIAL_RANGE            ";
  switch (_metabolite_tag_initial_range.law)
  {
    case UNIFORM:
      file << "UNIFORM " << _metabolite_tag_initial_range.min << " " << _metabolite_tag_initial_range.max << "\n";
      break;
    case GAUSSIAN:
      file << "GAUSSIAN " << _metabolite_tag_initial_range.mu << " " << _metabolite_tag_initial_range.sigma << "\n";
      break;
    case EXPONENTIAL:
      file << "EXPONENTIAL " << _metabolite_tag_initial_range.lambda << "\n";
      break;
  }
  file << "BINDING_SIZE_TAG_INITIAL_RANGE          ";
  switch (_binding_site_tag_initial_range.law)
  {
    case UNIFORM:
      file << "UNIFORM " << _binding_site_tag_initial_range.min << " " << _binding_site_tag_initial_range.max << "\n";
      break;
    case GAUSSIAN:
      file << "GAUSSIAN " << _binding_site_tag_initial_range.mu << " " << _binding_site_tag_initial_range.sigma << "\n";
      break;
    case EXPONENTIAL:
      file << "EXPONENTIAL " << _binding_site_tag_initial_range.lambda << "\n";
      break;
  }
  file << "CO_ENZYME_TAG_INITIAL_RANGE             ";
  switch (_co_enzyme_tag_initial_range.law)
  {
    case UNIFORM:
      file << "UNIFORM " << _co_enzyme_tag_initial_range.min << " " << _co_enzyme_tag_initial_range.max << "\n";
      break;
    case GAUSSIAN:
      file << "GAUSSIAN " << _co_enzyme_tag_initial_range.mu << " " << _co_enzyme_tag_initial_range.sigma << "\n";
      break;
    case EXPONENTIAL:
      file << "EXPONENTIAL " << _co_enzyme_tag_initial_range.lambda << "\n";
      break;
  }
  file << "TRANSCRIPTION_FACTOR_TAG_INITIAL_RANGE  ";
  switch (_transcription_factor_tag_initial_range.law)
  {
    case UNIFORM:
      file << "UNIFORM " << _transcription_factor_tag_initial_range.min << " " << _transcription_factor_tag_initial_range.max << "\n";
      break;
    case GAUSSIAN:
      file << "GAUSSIAN " << _transcription_factor_tag_initial_range.mu << " " << _transcription_factor_tag_initial_range.sigma << "\n";
      break;
    case EXPONENTIAL:
      file << "EXPONENTIAL " << _transcription_factor_tag_initial_range.lambda << "\n";
      break;
  }
  file << "\n";
  file << "INITIAL_BINDING_WINDOW  " << _initial_binding_window << "\n";
  file << "\n";
  file << "INITIAL_NUMBER_OF_NON_CODING_PEARLS            " << _initial_number_of_NC_pearls << "\n";
  file << "INITIAL_NUMBER_OF_ENZYME_PEARLS                " << _initial_number_of_E_pearls << "\n";
  file << "INITIAL_NUMBER_OF_TRANSCRIPTION_FACTOR_PEARLS  " << _initial_number_of_TF_pearls << "\n";
  file << "INITIAL_NUMBER_OF_BINDING_SITE_PEARLS          " << _initial_number_of_BS_pearls << "\n";
  file << "INITIAL_NUMBER_OF_PROMOTER_PEARLS              " << _initial_number_of_P_pearls << "\n";
  file << "\n";
  file << "POINT_MUTATION_RATE  " << _point_mutation_rate << "\n";
  file << "DUPLICATION_RATE     " << _duplication_rate << "\n";
  file << "DELETION_RATE        " << _deletion_rate << "\n";
  file << "TRANSLOCATION_RATE   " << _translocation_rate << "\n";
  file << "INVERSION_RATE       " << _inversion_rate << "\n";
  file << "TRANSITION_RATE      " << _transition_rate << "\n";
  file << "\n";
  file << "SUBSTRATE_TAG_MUTATION_SIZE             " << _substrate_tag_mutation_size << "\n";
  file << "PRODUCT_TAG_MUTATION_SIZE               " << _product_tag_mutation_size << "\n";
  file << "KM_MUTATION_SIZE                        " << _km_mutation_size << "\n";
  file << "KCAT_MUTATION_SIZE                      " << _kcat_mutation_size << "\n";
  file << "BINDING_SIZE_TAG_MUTATION_SIZE          " << _binding_size_tag_mutation_size << "\n";
  file << "CO_ENZYME_TAG_MUTATION_SIZE             " << _co_enzyme_tag_mutation_size << "\n";
  file << "TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE  " << _transcription_factor_tag_mutation_size << "\n";
  file << "BASAL_EXPRESSION_LEVEL_MUTATION_SIZE    " << _basal_expression_level_mutation_size << "\n";
  file << "\n";
  file << "MAXIMUM_GENOME_SIZE  " << _maximum_genome_size << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# GENETIC REGULATION NETWORK LEVEL\n";
  file << "########################################################\n";
  file << "\n";
  file << "GENETIC_REGULATION_NETWORK_TIMESTEP  " << _genetic_regulation_network_timestep << "\n";
  file << "\n";
  file << "HILL_FUNCTION_THETA       " << _hill_function_theta << "\n";
  file << "HILL_FUNCTION_N           " << _hill_function_n << "\n";
  file << "PROTEIN_DEGRADATION_RATE  " << _protein_degradation_rate << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# METABOLIC NETWORK LEVEL\n";
  file << "########################################################\n";
  file << "\n";
  file << "METABOLISM_TIMESTEP  " << _metabolism_timestep << "\n";
  file << "\n";
  file << "ESSENTIAL_METABOLITES_TOXICITY_THRESHOLD      " << _essential_metabolites_toxicity_threshold << "\n";
  file << "NON_ESSENTIAL_METABOLITES_TOXICITY_THRESHOLD  " << _non_essential_metabolites_toxicity_threshold << "\n";
  file << "\n";
  file << "INITIAL_METABOLITES_AMOUNT_IN_CELLS  " << _initial_metabolites_amount_in_cells << "\n";
  file << "\n";
  file << "ENERGY_REACTION_COST_FACTOR  " << _energy_reaction_cost_factor << "\n";
  file << "ENERGY_PUMPING_COST          " << _energy_pumping_cost << "\n";
  file << "ENERGY_TRANSCRIPTION_COST    " << _energy_transcription_cost << "\n";
  file << "ENERGY_TOXICITY_THRESHOLD    " << _energy_toxicity_threshold << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# CELL LEVEL\n";
  file << "########################################################\n";
  file << "\n";
  file << "DEATH_PROBABILITY  " << _death_probability << "\n";
  file << "\n";
  file << "VARIABLE_LIFESPAN  ";
  if (!_variable_lifespan)
  {
    file << "NO\n";
  }
  else
  {
    file << "YES\n";
  }
  file << "GOMPERTZ_LAW       " << _gompert_law.a << " " << _gompert_law.b << " " << _gompert_law.c << "\n";
  file << "\n";
  file << "INITIAL_MEMBRANE_PERMEABILITY  " << _initial_membrane_permeability << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# POPULATION LEVEL\n";
  file << "########################################################\n";
  file << "\n";
  file << "MIGRATION_RATE  " << _migration_rate << "\n";
  file << "HGT_RATE        " << _hgt_rate << "\n";
  file << "\n";
  file << "\n";
  file << "########################################################\n";
  file << "# ENVIRONMENT LEVEL\n";
  file << "########################################################\n";
  file << "\n";
  file << "ENVIRONMENT_INITIALIZATION_CYCLES    " << _environment_properties.number_of_init_cycles << "\n";;
  
  file << "ENVIRONMENT_SPECIES_TAG_RANGE        ";
  switch (_environment_properties.species_tag_range.law)
  {
    case UNIFORM:
      file << "UNIFORM " << _environment_properties.species_tag_range.min << " " << _environment_properties.species_tag_range.max << "\n";
      break;
    case GAUSSIAN:
      file << "GAUSSIAN " << _environment_properties.species_tag_range.mu << " " << _environment_properties.species_tag_range.sigma << "\n";
      break;
    case EXPONENTIAL:
      file << "EXPONENTIAL " << _environment_properties.species_tag_range.lambda << "\n";
      break;
  }
  file << "ENVIRONMENT_CONCENTRATION_RANGE      ";
  switch (_environment_properties.concentration_range.law)
  {
    case UNIFORM:
      file << "UNIFORM " << _environment_properties.concentration_range.min << " " << _environment_properties.concentration_range.max << "\n";
      break;
    case GAUSSIAN:
      file << "GAUSSIAN " << _environment_properties.concentration_range.mu << " " << _environment_properties.concentration_range.sigma << "\n";
      break;
    case EXPONENTIAL:
      file << "EXPONENTIAL " << _environment_properties.concentration_range.lambda << "\n";
      break;
  }
  file << "ENVIRONMENT_NUMBER_OF_SPECIES_RANGE  ";
  switch (_environment_properties.number_of_species_range.law)
  {
    case UNIFORM:
      file << "UNIFORM " << _environment_properties.number_of_species_range.min << " " << _environment_properties.number_of_species_range.max << "\n";
      break;
    case GAUSSIAN:
      file << "GAUSSIAN " << _environment_properties.number_of_species_range.mu << " " << _environment_properties.number_of_species_range.sigma << "\n";
      break;
    case EXPONENTIAL:
      file << "EXPONENTIAL " << _environment_properties.number_of_species_range.lambda << "\n";
      break;
  }
  file << "\n";
  file << "ENVIRONMENT_INTERACTION_SCHEME      ";
  switch (_environment_properties.interaction_scheme)
  {
    case NO_INTERACTION:
      file << "NO_INTERACTION (NO_INTERACTION/INTERACTION)\n";
      break;
    case INTERACTION:
      file << "INTERACTION (NO_INTERACTION/INTERACTION)\n";
      break;
  }
  file << "ENVIRONMENT_RENEWAL_SCHEME          ";
  switch (_environment_properties.renewal_scheme)
  {
    case KEEP_MATTER:
      file << "KEEP_MATTER (KEEP_MATTER/CLEAR_MATTER)\n";
      break;
    case CLEAR_MATTER:
      file << "CLEAR_MATTER (KEEP_MATTER/CLEAR_MATTER)\n";
      break;
  }
  file << "ENVIRONMENT_VARIATION_SCHEME        ";
  switch (_environment_properties.variation_scheme)
  {
    case RANDOM_SCHEME:
      file << "RANDOM (RANDOM/PERIODIC/CYCLIC)\n";
      break;
    case PERIODIC_SCHEME:
      file << "PERIODIC (RANDOM/PERIODIC/CYCLIC)\n";
      break;
    case CYCLIC_SCHEME:
      file << "CYCLIC (RANDOM/PERIODIC/CYCLIC)\n";
      break;
  }
  file << "ENVIRONMENT_VARIATION_LOCALIZATION  ";
  switch (_environment_properties.variation_localization)
  {
    case RANDOM_LOCALIZATION:
      file << "RANDOM (RANDOM/GLOBAL/CENTER)\n";
      break;
    case GLOBAL_LOCALIZATION:
      file << "GLOBAL (RANDOM/GLOBAL/CENTER)\n";
      break;
    case CENTER_LOCALIZATION:
      file << "CENTER (RANDOM/GLOBAL/CENTER)\n";
      break;
  }
  file << "\n";
  file << "ENVIRONMENT_INTRODUCTION_RATE  " << _environment_properties.introduction_rate << "\n";
  file << "ENVIRONMENT_DIFFUSION_RATE     " << _environment_properties.diffusion_rate << "\n";
  file << "ENVIRONMENT_DEGRADATION_RATE   " << _environment_properties.degradation_rate << "\n";
  file << "\n";
  file << "\n";
  file.close();
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Parse a line in words
 * \details  --
 * \param    char** words
 * \param    char* line
 * \return   \e void
 */
bool Parameters::parse_line( std::vector<std::string>* words, std::string line )
{
  if (line.size() == 0)
  {
    return false;
  }
  if (line[0] == '#' || line[0] == '\n')
  {
    return false;
  }
  words->clear();
  bool new_word   = false;
  bool word_saved = false;
  int pos       = -1;
  int length    = -1;
  for (int i = 0; i < (int)line.length(); i++)
  {
    /* if find a new character, start word saving */
    if (line[i] != ' ' && line[i] != '\n' && !new_word)
    {
      new_word   = true;
      word_saved = false;
      pos      = i;
      length   = 0;
    }
    /* else stop it */
    else if (line[i] == ' ' || line[i] == '\n')
    {
      new_word = false;
    }
    /* if a new word is found, save word length */
    if (new_word)
    {
      length++;
    }
    if (!new_word && !word_saved && pos >= 0 && length > 0)
    {
      words->push_back(line.substr(pos, length));
      word_saved = true;
    }
  }
  if (new_word && !word_saved && pos >= 0 && length > 0)
  {
    words->push_back(line.substr(pos, length));
  }
  return true;
}

/**
 * \brief    Build the list of the prime numbers until maximum
 * \details  Uses 'crible' method
 * \param    size_t maximum
 * \return   \e void
 */
void Parameters::build_prime_numbers_list( size_t maximum )
{
  _prime_numbers = NULL;
  _prime_numbers = new int[maximum];
  for (size_t i = 1; i <= maximum; i++)
  {
    _prime_numbers[i-1] = 1;
  }
  _prime_numbers[0] = 0;
  for (size_t i = 2; i <= maximum; i++)
  {
    size_t multiple = 2*i;
    while (multiple <= maximum)
    {
      _prime_numbers[multiple-1] = 0;
      multiple += i;
    }
  }
}
