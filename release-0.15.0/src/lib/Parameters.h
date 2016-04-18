
/**
 * \file      Parameters.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Parameters class declaration
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

#ifndef __EVOEVO__Parameters__
#define __EVOEVO__Parameters__

#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <zlib.h>
#include <assert.h>

#include "Macros.h"
#include "Enums.h"
#include "Structs.h"
#include "Prng.h"


class Parameters
{
  
public:
  
  /*----------------------------
   * CONSTRUCTORS
   *----------------------------*/
  Parameters( void );
  Parameters( size_t backup_time );
  Parameters( const Parameters& parameters );
  
  /*----------------------------
   * DESTRUCTORS
   *----------------------------*/
  ~Parameters( void );
  
  /*----------------------------
   * GETTERS
   *----------------------------*/
  
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  inline Prng*             get_prng( void );
  inline unsigned long int get_seed( void ) const;
  
  /*------------------------------------------------------------------ parallel computing */
  
  inline bool get_parallel_computing( void ) const;
  
  /*------------------------------------------------------------------ simulation schemes */
  
  inline bool         get_energy_constraints( void ) const;
  inline bool         get_membrane_permeability( void ) const;
  inline bool         get_metabolic_inheritance( void ) const;
  inline bool         get_enzymatic_inheritance( void ) const;
  inline bool         get_co_enzyme_activity( void ) const;
  inline score_scheme get_score_scheme( void ) const;
  inline double       get_selection_threshold( void ) const;
  
  /*------------------------------------------------------------------ space */
  
  inline size_t get_width( void ) const;
  inline size_t get_height( void ) const;
  
  /*------------------------------------------------------------------ output */
  
  inline size_t get_experiment_backup_step( void ) const;
  inline size_t get_tree_backup_step( void ) const;
  inline size_t get_figures_generation_step( void ) const;
  
  /*------------------------------------------------------------------ genome level */
  
  inline distribution_law* get_metabolite_tag_initial_range( void );
  inline distribution_law* get_binding_site_tag_initial_range( void );
  inline distribution_law* get_co_enzyme_tag_initial_range( void );
  inline distribution_law* get_transcription_factor_tag_initial_range( void );
  
  inline size_t get_initial_binding_window( void ) const;
  
  inline size_t get_initial_number_of_NC_pearls( void ) const;
  inline size_t get_initial_number_of_E_pearls( void ) const;
  inline size_t get_initial_number_of_TF_pearls( void ) const;
  inline size_t get_initial_number_of_BS_pearls( void ) const;
  inline size_t get_initial_number_of_P_pearls( void ) const;
  
  inline double get_point_mutation_rate( void ) const;
  inline double get_duplication_rate( void ) const;
  inline double get_deletion_rate( void ) const;
  inline double get_translocation_rate( void ) const;
  inline double get_inversion_rate( void ) const;
  inline double get_transition_rate( void ) const;
  
  inline double get_substrate_tag_mutation_size( void ) const;
  inline double get_product_tag_mutation_size( void ) const;
  inline double get_km_mutation_size( void ) const;
  inline double get_kcat_mutation_size( void ) const;
  inline double get_binding_size_tag_mutation_size( void ) const;
  inline double get_co_enzyme_tag_mutation_size( void ) const;
  inline double get_transcription_factor_tag_mutation_size( void ) const;
  inline double get_basal_expression_level_mutation_size( void ) const;
  
  inline size_t get_maximum_genome_size( void ) const;
  
  /*------------------------------------------------------------------ genetic regulation network level */
  
  inline double get_genetic_regulation_network_timestep( void ) const;
  
  inline double get_hill_function_theta( void ) const;
  inline double get_hill_function_n( void ) const;
  inline double get_protein_degradation_rate( void ) const;
  
  /*------------------------------------------------------------------ metabolic network level */
  
  inline double get_metabolism_timestep( void ) const;
  
  inline double get_essential_metabolites_toxicity_threshold( void );
  inline double get_non_essential_metabolites_toxicity_threshold( void );
  
  inline double get_initial_metabolites_amount_in_cells( void ) const;
  
  inline double get_energy_reaction_cost_factor( void ) const;
  inline double get_energy_pumping_cost( void ) const;
  inline double get_energy_transcription_cost( void ) const;
  inline double get_energy_toxicity_threshold( void ) const;
  
  /*------------------------------------------------------------------ cell level */
  
  inline double get_death_probability( void ) const;
  
  inline bool         get_variable_lifespan( void ) const;
  inline gompert_law* get_gompert_law( void );
  inline double       get_initial_membrane_permeability( void ) const;
  
  /*------------------------------------------------------------------ population level */
  
  inline double get_migration_rate( void ) const;
  inline double get_hgt_rate( void ) const;
  
  /*------------------------------------------------------------------ environment level */
  
  inline environment_properties* get_environment_properties( void );
  
  /*------------------------------------------------------------------ prime numbers */
  
  inline int* get_prime_numbers( void );

  /*----------------------------
   * SETTERS
   *----------------------------*/
  
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  inline void set_prng( Prng* prng );
  inline void set_seed( unsigned long int seed );
  
  /*------------------------------------------------------------------ parallel computing */
  
  inline void set_parallel_computing( bool parallel_computing );
  
  /*------------------------------------------------------------------ simulation schemes */
  
  inline void set_energy_constraints( bool energy_constraints );
  inline void set_membrane_permeability( bool membrane_permeability );
  inline void set_metabolic_inheritance( bool metabolic_inheritance );
  inline void set_enzymatic_inheritance( bool enzymatic_inheritance );
  inline void set_co_enzyme_activity( bool co_enzyme_activity );
  inline void set_score_scheme( score_scheme scheme );
  inline void set_selection_threshold( double selection_threshold );
  
  /*------------------------------------------------------------------ space */
  
  inline void set_width( size_t width );
  inline void set_height( size_t height );
  
  /*------------------------------------------------------------------ output */
  
  inline void set_experiment_backup_step( size_t step );
  inline void set_tree_backup_step( size_t step );
  inline void set_figures_generation_step( size_t step );
  
  /*------------------------------------------------------------------ genome level */
  
  inline void set_metabolite_tag_initial_range( const distribution_law* range_law );
  inline void set_binding_site_tag_initial_range( const distribution_law* range_law );
  inline void set_co_enzyme_tag_initial_range( const distribution_law* range_law );
  inline void set_transcription_factor_tag_initial_range( const distribution_law* range_law );
  
  inline void set_initial_binding_window( size_t window );
  
  inline void set_initial_number_of_NC_pearls( size_t number_of_NC_pearls );
  inline void set_initial_number_of_E_pearls( size_t number_of_E_pearls );
  inline void set_initial_number_of_TF_pearls( size_t number_of_TF_pearls );
  inline void set_initial_number_of_BS_pearls( size_t number_of_BS_pearls );
  inline void set_initial_number_of_P_pearls( size_t number_of_P_pearls );
  
  inline void set_point_mutation_rate( double rate );
  inline void set_duplication_rate( double rate );
  inline void set_deletion_rate( double rate );
  inline void set_translocation_rate( double rate );
  inline void set_inversion_rate( double rate );
  inline void set_transition_rate( double rate );
  
  inline void set_substrate_tag_mutation_size( double size );
  inline void set_product_tag_mutation_size( double size );
  inline void set_km_mutation_size( double size );
  inline void set_kcat_mutation_size( double size );
  inline void set_binding_size_tag_mutation_size( double size );
  inline void set_co_enzyme_tag_mutation_size( double size );
  inline void set_transcription_factor_tag_mutation_size( double size );
  inline void set_basal_expression_level_mutation_size( double size );
  
  inline void set_maximum_genome_size( size_t maximum_genome_size );
  
  /*------------------------------------------------------------------ genetic regulation network level */
  
  inline void set_genetic_regulation_network_timestep( double genetic_regulation_network_timestep );
  
  inline void set_hill_function_theta( double hill_function_theta );
  inline void set_hill_function_n( double hill_function_n );
  inline void set_protein_degradation_rate( double protein_degradation_rate );
  
  /*------------------------------------------------------------------ metabolic network level */
  
  inline void set_metabolism_timestep( double metabolism_timestep );
  
  inline void set_essential_metabolites_toxicity_threshold( double toxicity_threshold );
  inline void set_non_essential_metabolites_toxicity_threshold( double toxicity_threshold );
  
  inline void set_initial_metabolites_amount_in_cells( double metabolites_amount_in_cells );
  
  inline void set_energy_reaction_cost_factor( double energy_reaction_cost_factor );
  inline void set_energy_pumping_cost( double energy_pumping_cost );
  inline void set_energy_transcription_cost( double energy_transcription_cost );
  inline void set_energy_toxicity_threshold( double energy_toxicity_threshold );
  
  /*------------------------------------------------------------------ cell level */
  
  inline void set_death_probability( double death_probability );
  
  inline void set_variable_lifespan( bool variable_lifespan );
  inline void set_gompert_law( double a, double b, double c );
  
  inline void set_initial_membrane_permeability( double permeability );
  
  /*------------------------------------------------------------------ population level */
  
  inline void set_migration_rate( double migration_rate );
  inline void set_hgt_rate( double hgt_rate );
  
  /*------------------------------------------------------------------ environment level */
  
  inline void set_environment_properties( const environment_properties* properties );
  
  /*----------------------------
   * PUBLIC METHODS
   *----------------------------*/
  bool load_parameters_from_file( std::string filename );
  void save( size_t backup_time );
  void write( std::string filename );
  
  /*----------------------------
   * PUBLIC ATTRIBUTES
   *----------------------------*/
  
protected:
  
  /*----------------------------
   * PROTECTED METHODS
   *----------------------------*/
  bool parse_line( std::vector<std::string>* words, std::string line );
  void build_prime_numbers_list( size_t maximum );
  
  /*----------------------------
   * PROTECTED ATTRIBUTES
   *----------------------------*/
  
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  Prng*             _prng; /*!< Pseudorandom numbers generator */
  unsigned long int _seed; /*!< Seed of the prng               */
  
  /*------------------------------------------------------------------ parallel computing */
  
  bool _parallel_computing; /*!< Indicates if the parallel computing is activated */
  
  /*------------------------------------------------------------------ simulation schemes */
  
  bool         _energy_constraints;    /*!< Indicates if energy constraints are activated        */
  bool         _membrane_permeability; /*!< Indicates if membrane permeability is activated      */
  bool         _metabolic_inheritance; /*!< Indicates if metabolic inheritance is activated      */
  bool         _enzymatic_inheritance; /*!< Indicates if enzymatic inheritance is activated      */
  bool         _co_enzyme_activity;    /*!< Indicates if co-enzyme activity is activated         */
  score_scheme _score_scheme;          /*!< Score scheme used to compute the score               */
  double       _selection_threshold;   /*!< Selection threshold when selection shecme is fitprop */
  
  /*------------------------------------------------------------------ space */
  
  size_t _width;  /*!< Spatial grid width  */
  size_t _height; /*!< Spatial grid height */
  
  /*------------------------------------------------------------------ output */
  
  size_t _experiment_backup_step;  /*!< Step of experiment backups */
  size_t _tree_backup_step;        /*!< Step of tree backups       */
  size_t _figures_generation_step; /*!< Step of figures generation */
  
  /*------------------------------------------------------------------ genome level */
  
  distribution_law _metabolite_tag_initial_range;           /*!< Metabolite tag initial range           */
  distribution_law _binding_site_tag_initial_range;         /*!< Binding site tag initial range         */
  distribution_law _co_enzyme_tag_initial_range;            /*!< Co-enzyme tag initial range            */
  distribution_law _transcription_factor_tag_initial_range; /*!< Transcription factor tag initial range */
  
  size_t           _initial_binding_window;                 /*!< Initial binding window of a TF on a BS */
  
  size_t           _initial_number_of_NC_pearls;            /*!< Initial number of NC pearls            */
  size_t           _initial_number_of_E_pearls;             /*!< Initial number of E pearls             */
  size_t           _initial_number_of_TF_pearls;            /*!< Initial number of TF pearls            */
  size_t           _initial_number_of_BS_pearls;            /*!< Initial number of BS pearls            */
  size_t           _initial_number_of_P_pearls;             /*!< Initial number of P pearls             */
  
  double           _point_mutation_rate;                    /*!< Point mutation rate                    */
  double           _duplication_rate;                       /*!< Duplication rate                       */
  double           _deletion_rate;                          /*!< Deletion rate                          */
  double           _translocation_rate;                     /*!< Translocation rate                     */
  double           _inversion_rate;                         /*!< Inversion rate                         */
  double           _transition_rate;                        /*!< Transition rate                        */
  
  double           _substrate_tag_mutation_size;            /*!< Substrate tag mutation size            */
  double           _product_tag_mutation_size;              /*!< Product tag mutation size              */
  double           _km_mutation_size;                       /*!< Km constant mutation size              */
  double           _kcat_mutation_size;                     /*!< Kcat constant mutation size            */
  double           _binding_size_tag_mutation_size;         /*!< Binding site tag mutation size         */
  double           _co_enzyme_tag_mutation_size;            /*!< Co-enzyme tag mutation size            */
  double           _transcription_factor_tag_mutation_size; /*!< Transcription factor tag mutation size */
  double           _basal_expression_level_mutation_size;   /*!< Basal expression level mutation size   */
  
  size_t           _maximum_genome_size;                    /*!< Maximum genome size                    */
  
  /*------------------------------------------------------------------ genetic regulation network level */
  
  double _genetic_regulation_network_timestep; /*!< Genetic regulation network ODE timestep */
  
  double _hill_function_theta;                 /*!< Parameter theta of the Hill function    */
  double _hill_function_n;                     /*!< Parameter n of the Hill function        */
  double _protein_degradation_rate;            /*!< Protein degradation rate                */
  
  /*------------------------------------------------------------------ metabolic network level */
  
  double _metabolism_timestep;                          /*!< Metabolism ODE timestep                          */
  
  double _essential_metabolites_toxicity_threshold;     /*!< Essential metabolites toxicity threshold         */
  double _non_essential_metabolites_toxicity_threshold; /*!< Non essential metabolites toxicity threshold     */
  
  double _initial_metabolites_amount_in_cells;          /*!< Initial metabolites amount in cells              */
  
  double _energy_reaction_cost_factor;                  /*!< Differential factor for enzymatic reaction costs */
  double _energy_pumping_cost;                          /*!< Energetical cost of pumps                        */
  double _energy_transcription_cost;                    /*!< Energetical cost of transcription                */
  double _energy_toxicity_threshold;                    /*!< Energy toxicity threshold                        */
  
  /*------------------------------------------------------------------ cell level */
  
  double      _death_probability;             /*!< Death probability                                     */
  bool        _variable_lifespan;             /*!< Variable Lifespan boolean                             */
  gompert_law _gompert_law;                   /*!< Gompertz law (death probability law)                  */
  double      _initial_membrane_permeability; /*!< Initial membrane permeability (diffusion coefficient) */
  
  /*------------------------------------------------------------------ population level */
  
  double _migration_rate; /*!< Migration rate */
  double _hgt_rate;       /*!< HGT rate       */
  
  /*------------------------------------------------------------------ environment level */
  
  environment_properties _environment_properties; /*!< Environment properties */
  
  /*------------------------------------------------------------------ prime numbers */
  
  int* _prime_numbers; /*!< List of prime numbers */
  
};


/*----------------------------
 * GETTERS
 *----------------------------*/

/*------------------------------------------------------------------ pseudorandom numbers generator */

/**
 * \brief    Get prng
 * \details  Return the pseudorandom numbers generator (Mersenne twister 19937 generator from gsl)
 * \param    void
 * \return   \e Prng*
 */
inline Prng* Parameters::get_prng( void )
{
  return _prng;
}

/**
 * \brief    Get PRNG seed
 * \details  Return the seed of the prng
 * \param    void
 * \return   \e unsigned long int
 */
inline unsigned long int Parameters::get_seed( void ) const
{
  return _seed;
}

/*------------------------------------------------------------------ parallel computing */

/**
 * \brief    Get parallel computing
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_parallel_computing( void ) const
{
  return _parallel_computing;
}

/*------------------------------------------------------------------ simulation schemes */

/**
 * \brief    Get the energy constraints boolean
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_energy_constraints( void ) const
{
  return _energy_constraints;
}

/**
 * \brief    Get the membrane permeability boolean
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_membrane_permeability( void ) const
{
  return _membrane_permeability;
}

/**
 * \brief    Get the metabolic inheritance boolean
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_metabolic_inheritance( void ) const
{
  return _metabolic_inheritance;
}

/**
 * \brief    Get the enzymatic inheritance boolean
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_enzymatic_inheritance( void ) const
{
  return _enzymatic_inheritance;
}

/**
 * \brief    Get the co-enzyme activity boolean
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_co_enzyme_activity( void ) const
{
  return _co_enzyme_activity;
}

/**
 * \brief    Get the score scheme
 * \details  --
 * \param    void
 * \return   \e score_scheme
 */
inline score_scheme Parameters::get_score_scheme( void ) const
{
  return _score_scheme;
}

/**
 * \brief    Get the selection threshold
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_selection_threshold( void ) const
{
  return _selection_threshold;
}

/*------------------------------------------------------------------ space */

/**
 * \brief    Get population's grid width
 * \details  Return the population's grid width
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_width( void ) const
{
  return _width;
}

/**
 * \brief    Get population's grid height
 * \details  Return the population's grid height
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_height( void ) const
{
  return _height;
}

/*------------------------------------------------------------------ output */

/**
 * \brief    Get experiment backup step
 * \details  Return the step of each experiment backup
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_experiment_backup_step( void ) const
{
  return _experiment_backup_step;
}

/**
 * \brief    Get tree backup step
 * \details  Return the step of each tree backup
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_tree_backup_step( void ) const
{
  return _tree_backup_step;
}

/**
 * \brief    Get figures generation step
 * \details  Return the step of each figures generation
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_figures_generation_step( void ) const
{
  return _figures_generation_step;
}

/*------------------------------------------------------------------ genome level */

/**
 * \brief    Get metabolite tag initial range
 * \details  --
 * \param    void
 * \return   \e distribution_law*
 */
inline distribution_law* Parameters::get_metabolite_tag_initial_range( void )
{
  return &_metabolite_tag_initial_range;
}

/**
 * \brief    Binding site tag initial range
 * \details  --
 * \param    void
 * \return   \e distribution_law*
 */
inline distribution_law* Parameters::get_binding_site_tag_initial_range( void )
{
  return &_binding_site_tag_initial_range;
}

/**
 * \brief    Co-enzyme tag initial range
 * \details  --
 * \param    void
 * \return   \e distribution_law*
 */
inline distribution_law* Parameters::get_co_enzyme_tag_initial_range( void )
{
  return &_co_enzyme_tag_initial_range;
}

/**
 * \brief    Get transcription factor tag initial range
 * \details  --
 * \param    void
 * \return   \e distribution_law*
 */
inline distribution_law* Parameters::get_transcription_factor_tag_initial_range( void )
{
  return &_transcription_factor_tag_initial_range;
}

/**
 * \brief    Get initial binding window of a TF on a BS
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_initial_binding_window( void ) const
{
  return _initial_binding_window;
}

/**
 * \brief    Get initial number of non coding pearls (NC)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_initial_number_of_NC_pearls( void ) const
{
  return _initial_number_of_NC_pearls;
}

/**
 * \brief    Get initial number of enzyme pearls (E)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_initial_number_of_E_pearls( void ) const
{
  return _initial_number_of_E_pearls;
}

/**
 * \brief    Get initial number of transcription factor pearls (TF)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_initial_number_of_TF_pearls( void ) const
{
  return _initial_number_of_TF_pearls;
}

/**
 * \brief    Get initial number of binding site pearls (BS)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_initial_number_of_BS_pearls( void ) const
{
  return _initial_number_of_BS_pearls;
}

/**
 * \brief    Get initial number of promoter pearls (P)
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_initial_number_of_P_pearls( void ) const
{
  return _initial_number_of_P_pearls;
}

/**
 * \brief    Get point mutation rate
 * \details  Return point mutation rate used to initialize individuals
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_point_mutation_rate( void ) const
{
  return _point_mutation_rate;
}

/**
 * \brief    Get duplication rate
 * \details  Return duplication rate used to initialize individuals
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_duplication_rate( void ) const
{
  return _duplication_rate;
}

/**
 * \brief    Get deletion rate
 * \details  Return deletion rate used to initialize individuals
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_deletion_rate( void ) const
{
  return _deletion_rate;
}

/**
 * \brief    Get translocation rate
 * \details  Return translocation rate used to initialize individuals
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_translocation_rate( void ) const
{
  return _translocation_rate;
}

/**
 * \brief    Get inversion rate
 * \details  Return inversion rate used to initialize individuals
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_inversion_rate( void ) const
{
  return _inversion_rate;
}

/**
 * \brief    Get initial transition rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_transition_rate( void ) const
{
  return _transition_rate;
}

/**
 * \brief    Get substrate tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_substrate_tag_mutation_size( void ) const
{
  return _substrate_tag_mutation_size;
}

/**
 * \brief    Get product tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_product_tag_mutation_size( void ) const
{
  return _product_tag_mutation_size;
}

/**
 * \brief    Get km constant mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_km_mutation_size( void ) const
{
  return _km_mutation_size;
}

/**
 * \brief    Get kcat constant mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_kcat_mutation_size( void ) const
{
  return _kcat_mutation_size;
}

/**
 * \brief    Get binding site tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_binding_size_tag_mutation_size( void ) const
{
  return _binding_size_tag_mutation_size;
}

/**
 * \brief    Get co-enzyme tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_co_enzyme_tag_mutation_size( void ) const
{
  return _co_enzyme_tag_mutation_size;
}

/**
 * \brief    Get transcription factor tag mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_transcription_factor_tag_mutation_size( void ) const
{
  return _transcription_factor_tag_mutation_size;
}

/**
 * \brief    Get basal expression level mutation size
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_basal_expression_level_mutation_size( void ) const
{
  return _basal_expression_level_mutation_size;
}

/**
 * \brief    Get the maximum genome size
 * \details  --
 * \param    void
 * \return   \e size_t
 */
inline size_t Parameters::get_maximum_genome_size( void ) const
{
  return _maximum_genome_size;
}

/*------------------------------------------------------------------ genetic regulation network level */

/**
 * \brief    Get the genetic regulation network ODE timestep
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_genetic_regulation_network_timestep( void ) const
{
  return _genetic_regulation_network_timestep;
}

/**
 * \brief    Get the parameter theta of the Hill function
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_hill_function_theta( void ) const
{
  return _hill_function_theta;
}

/**
 * \brief    Get the parameter n of the Hill function
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_hill_function_n( void ) const
{
  return _hill_function_n;
}

/**
 * \brief    Get the protein degradation rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_protein_degradation_rate( void ) const
{
  return _protein_degradation_rate;
}

/*------------------------------------------------------------------ metabolic network level */

/**
 * \brief    Get metabolism ODE timestep
 * \details  Return the timestep of ODE solver per simulation timestep
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_metabolism_timestep( void ) const
{
  return _metabolism_timestep;
}

/**
 * \brief    Get essential metabolites toxicity threshold
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_essential_metabolites_toxicity_threshold( void )
{
  return _essential_metabolites_toxicity_threshold;
}

/**
 * \brief    Get non essential metabolites toxicity threshold
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_non_essential_metabolites_toxicity_threshold( void )
{
  return _non_essential_metabolites_toxicity_threshold;
}

/**
 * \brief    Get initial metabolites amount in cells
 * \details  Return the initial metabolite concentration in cells
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_initial_metabolites_amount_in_cells( void ) const
{
  return _initial_metabolites_amount_in_cells;
}

/**
 * \brief    Get energy differential factor for enzymatic reaction costs
 * \details  Return the differential cost when producing or consuming energy during enzymatic reactions
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_energy_reaction_cost_factor( void ) const
{
  return _energy_reaction_cost_factor;
}

/**
 * \brief    Get energy pumping cost
 * \details  Return the energy cost of pump activity
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_energy_pumping_cost( void ) const
{
  return _energy_pumping_cost;
}

/**
 * \brief    Get energy transcription cost
 * \details  Return the energy cost of transcription
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_energy_transcription_cost( void ) const
{
  return _energy_transcription_cost;
}

/**
 * \brief    Get energy toxicity threshold
 * \details  Return the energy toxicity threshold
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_energy_toxicity_threshold( void ) const
{
  return _energy_toxicity_threshold;
}

/*------------------------------------------------------------------ cell level */

/**
 * \brief    Get death probability
 * \details  Return individual death probability
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_death_probability( void ) const
{
  return _death_probability;
}

/**
 * \brief    Get variable lifespan boolean
 * \details  --
 * \param    void
 * \return   \e bool
 */
inline bool Parameters::get_variable_lifespan( void ) const
{
  return _variable_lifespan;
}

/**
 * \brief    Get Gompertz law
 * \details  --
 * \param    void
 * \return   \e gompert_law*
 */
inline gompert_law* Parameters::get_gompert_law( void )
{
  return &_gompert_law;
}

/**
 * \brief    Get initial membrane permeability
 * \details  Defines the initial permeability of the cell membrane
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_initial_membrane_permeability( void ) const
{
  return _initial_membrane_permeability;
}

/*------------------------------------------------------------------ population level */

/**
 * \brief    Get the migration rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_migration_rate( void ) const
{
  return _migration_rate;
}

/**
 * \brief    Get the HGT rate
 * \details  --
 * \param    void
 * \return   \e double
 */
inline double Parameters::get_hgt_rate( void ) const
{
  return _hgt_rate;
}

/*------------------------------------------------------------------ environment level */

/**
 * \brief    Get the environment properties
 * \details  --
 * \param    void
 * \return   \e environment_properties*
 */
inline environment_properties* Parameters::get_environment_properties( void )
{
  return &_environment_properties;
}

/*------------------------------------------------------------------ prime numbers */

/**
 * \brief    Get list of prime numbers
 * \details  Return the first 25 prime numbers. This list is used to evaluate individuals if fitness evaluation is 'PRIME_NUMBERS'
 * \param    void
 * \return   \e int*
 */
inline int* Parameters::get_prime_numbers( void )
{
  return _prime_numbers;
}

/*----------------------------
 * SETTERS
 *----------------------------*/

/*------------------------------------------------------------------ pseudorandom numbers generator */

/**
 * \brief    Set prng
 * \details  --
 * \param    Prng* prng
 * \return   \e void
 */
inline void Parameters::set_prng( Prng* prng )
{
  delete _prng;
  _prng = NULL;
  _prng = new Prng(*prng);
}

/**
 * \brief    Set PRNG seed
 * \details  --
 * \param    unsigned long int seed
 * \return   \e void
 */
inline void Parameters::set_seed( unsigned long int seed )
{
  assert(seed > 0);
  _seed = seed;
}

/*------------------------------------------------------------------ parallel computing */

/**
 * \brief    Set parallel computing
 * \details  --
 * \param    bool parallel_computing
 * \return   \e void
 */
inline void Parameters::set_parallel_computing( bool parallel_computing )
{
  _parallel_computing = parallel_computing;
}

/*------------------------------------------------------------------ simulation schemes */

/**
 * \brief    Set the energy constraints boolean
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Parameters::set_energy_constraints( bool energy_constraints )
{
  _energy_constraints = energy_constraints;
}

/**
 * \brief    Set the membrane permeability boolean
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Parameters::set_membrane_permeability( bool membrane_permeability )
{
  _membrane_permeability = membrane_permeability;
}

/**
 * \brief    Set the metabolic inheritance boolean
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Parameters::set_metabolic_inheritance( bool metabolic_inheritance )
{
  _metabolic_inheritance = metabolic_inheritance;
}

/**
 * \brief    Set the enzymatic inheritance boolean
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Parameters::set_enzymatic_inheritance( bool enzymatic_inheritance )
{
  _enzymatic_inheritance = enzymatic_inheritance;
}

/**
 * \brief    Set the co-enzyme activity boolean
 * \details  --
 * \param    void
 * \return   \e void
 */
inline void Parameters::set_co_enzyme_activity( bool co_enzyme_activity )
{
  _co_enzyme_activity = co_enzyme_activity;
}

/**
 * \brief    Set the score scheme
 * \details  --
 * \param    score_scheme scheme
 * \return   \e void
 */
inline void Parameters::set_score_scheme( score_scheme scheme )
{
  _score_scheme = scheme;
}

/**
 * \brief    Set the selection threshold
 * \details  --
 * \param    double selection_threshold
 * \return   \e void
 */
inline void Parameters::set_selection_threshold( double selection_threshold )
{
  assert(selection_threshold >= 0.0);
  assert(selection_threshold <= 1.0);
  _selection_threshold = selection_threshold;
}

/*------------------------------------------------------------------ space */

/**
 * \brief    Set population's grid width
 * \details  --
 * \param    size_t width
 * \return   \e void
 */
inline void Parameters::set_width( size_t width )
{
  _width = width;
}

/**
 * \brief    Set population's grid height
 * \details  --
 * \param    size_t height
 * \return   \e void
 */
inline void Parameters::set_height( size_t height )
{
  _height = height;
}

/*------------------------------------------------------------------ output */

/**
 * \brief    Set experiment backup step
 * \details  --
 * \param    size_t step
 * \return   \e void
 */
inline void Parameters::set_experiment_backup_step( size_t step )
{
  _experiment_backup_step = step;
}

/**
 * \brief    Set tree backup step
 * \details  --
 * \param    size_t step
 * \return   \e void
 */
inline void Parameters::set_tree_backup_step( size_t step )
{
  _tree_backup_step = step;
}

/**
 * \brief    Set figures generation step
 * \details  --
 * \param    size_t step
 * \return   \e void
 */
inline void Parameters::set_figures_generation_step( size_t step )
{
  _figures_generation_step = step;
}

/*------------------------------------------------------------------ genome level */

/**
 * \brief    Set metabolite tag initial range
 * \details  --
 * \param    distribution_law* law
 * \return   \e void
 */
inline void Parameters::set_metabolite_tag_initial_range( const distribution_law* range_law )
{
  _metabolite_tag_initial_range.law    = range_law->law;
  _metabolite_tag_initial_range.min    = range_law->min;
  _metabolite_tag_initial_range.max    = range_law->max;
  _metabolite_tag_initial_range.mu     = range_law->mu;
  _metabolite_tag_initial_range.sigma  = range_law->sigma;
  _metabolite_tag_initial_range.lambda = range_law->lambda;
}

/**
 * \brief    Set binding site tag initial range
 * \details  --
 * \param    distribution_law* law
 * \return   \e void
 */
inline void Parameters::set_binding_site_tag_initial_range( const distribution_law* range_law )
{
  _binding_site_tag_initial_range.law    = range_law->law;
  _binding_site_tag_initial_range.min    = range_law->min;
  _binding_site_tag_initial_range.max    = range_law->max;
  _binding_site_tag_initial_range.mu     = range_law->mu;
  _binding_site_tag_initial_range.sigma  = range_law->sigma;
  _binding_site_tag_initial_range.lambda = range_law->lambda;
}

/**
 * \brief    Set co-enzyme tag initial range
 * \details  --
 * \param    distribution_law* law
 * \return   \e void
 */
inline void Parameters::set_co_enzyme_tag_initial_range( const distribution_law* range_law )
{
  _co_enzyme_tag_initial_range.law    = range_law->law;
  _co_enzyme_tag_initial_range.min    = range_law->min;
  _co_enzyme_tag_initial_range.max    = range_law->max;
  _co_enzyme_tag_initial_range.mu     = range_law->mu;
  _co_enzyme_tag_initial_range.sigma  = range_law->sigma;
  _co_enzyme_tag_initial_range.lambda = range_law->lambda;
}

/**
 * \brief    Set transcription factor tag initial range
 * \details  --
 * \param    distribution_law* law
 * \return   \e void
 */
inline void Parameters::set_transcription_factor_tag_initial_range( const distribution_law* range_law )
{
  _transcription_factor_tag_initial_range.law    = range_law->law;
  _transcription_factor_tag_initial_range.min    = range_law->min;
  _transcription_factor_tag_initial_range.max    = range_law->max;
  _transcription_factor_tag_initial_range.mu     = range_law->mu;
  _transcription_factor_tag_initial_range.sigma  = range_law->sigma;
  _transcription_factor_tag_initial_range.lambda = range_law->lambda;
}

/**
 * \brief    Set initial binding window of a TF on a BS
 * \details  --
 * \param    size_t window
 * \return   \e void
 */
inline void Parameters::set_initial_binding_window( size_t window )
{
  _initial_binding_window = window;
}

/**
 * \brief    Set initial number of non coding pearls (NC)
 * \details  --
 * \param    size_t number_of_NC_pearls
 * \return   \e void
 */
inline void Parameters::set_initial_number_of_NC_pearls( size_t number_of_NC_pearls )
{
  _initial_number_of_NC_pearls = number_of_NC_pearls;
}

/**
 * \brief    Set initial number of enzyme pearls (E)
 * \details  --
 * \param    size_t number_of_E_pearls
 * \return   \e void
 */
inline void Parameters::set_initial_number_of_E_pearls( size_t number_of_E_pearls )
{
  _initial_number_of_E_pearls = number_of_E_pearls;
}

/**
 * \brief    Set initial number of transcription factor pearls (TF)
 * \details  --
 * \param    size_t number_of_TF_pearls
 * \return   \e void
 */
inline void Parameters::set_initial_number_of_TF_pearls( size_t number_of_TF_pearls )
{
  _initial_number_of_TF_pearls = number_of_TF_pearls;
}

/**
 * \brief    Set initial number of binding site pearls (BS)
 * \details  --
 * \param    size_t number_of_BS_pearls
 * \return   \e void
 */
inline void Parameters::set_initial_number_of_BS_pearls( size_t number_of_BS_pearls )
{
  _initial_number_of_BS_pearls = number_of_BS_pearls;
}

/**
 * \brief    Set initial number of promoter pearls (P)
 * \details  --
 * \param    size_t number_of_P_pearls
 * \return   \e void
 */
inline void Parameters::set_initial_number_of_P_pearls( size_t number_of_P_pearls )
{
  _initial_number_of_P_pearls = number_of_P_pearls;
}

/**
 * \brief    Set point mutation rates
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_point_mutation_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _point_mutation_rate = rate;
}

/**
 * \brief    Set duplication rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_duplication_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _duplication_rate = rate;
}

/**
 * \brief    Set deletion rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_deletion_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _deletion_rate = rate;
}

/**
 * \brief    Set translocation rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_translocation_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _translocation_rate = rate;
}

/**
 * \brief    Set inversion rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_inversion_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _inversion_rate = rate;
}

/**
 * \brief    Set transition rate
 * \details  --
 * \param    double rate
 * \return   \e void
 */
inline void Parameters::set_transition_rate( double rate )
{
  assert(rate >= 0.0);
  assert(rate <= 1.0);
  _transition_rate = rate;
}

/**
 * \brief    Set substrate tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_substrate_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _substrate_tag_mutation_size = size;
}

/**
 * \brief    Set product tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_product_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _product_tag_mutation_size = size;
}

/**
 * \brief    Set km constant mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_km_mutation_size( double size )
{
  assert(size >= 0.0);
  _km_mutation_size = size;
}

/**
 * \brief    Set kcat constant mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_kcat_mutation_size( double size )
{
  assert(size >= 0.0);
  _kcat_mutation_size = size;
}

/**
 * \brief    Set binding site tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_binding_size_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _binding_size_tag_mutation_size = size;
}

/**
 * \brief    Set co-enzyme tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_co_enzyme_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _co_enzyme_tag_mutation_size = size;
}

/**
 * \brief    Set transcription factor tag mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_transcription_factor_tag_mutation_size( double size )
{
  assert(size >= 0.0);
  _transcription_factor_tag_mutation_size = size;
}

/**
 * \brief    Set basal expression level mutation size
 * \details  --
 * \param    double size
 * \return   \e void
 */
inline void Parameters::set_basal_expression_level_mutation_size( double size )
{
  assert(size >= 0.0);
  _basal_expression_level_mutation_size = size;
}

/**
 * \brief    Set the maximum genome size
 * \details  --
 * \param    size_t maximum_genome_size
 * \return   \e void
 */
inline void Parameters::set_maximum_genome_size( size_t maximum_genome_size )
{
  _maximum_genome_size = maximum_genome_size;
}

/*------------------------------------------------------------------ genetic regulation network level */

/**
 * \brief    Set the genetic regulation network ODE timestep
 * \details  --
 * \param    double genetic_regulation_network_timestep
 * \return   \e void
 */
inline void Parameters::set_genetic_regulation_network_timestep( double genetic_regulation_network_timestep )
{
  assert(genetic_regulation_network_timestep > 0.0);
  _genetic_regulation_network_timestep = genetic_regulation_network_timestep;
}

/**
 * \brief    Set the parameter theta of the Hill function
 * \details  --
 * \param    double hill_function_theta
 * \return   \e void
 */
inline void Parameters::set_hill_function_theta( double hill_function_theta )
{
  assert(hill_function_theta >= 0.0);
  assert(hill_function_theta <= 1.0);
  _hill_function_theta = hill_function_theta;
}

/**
 * \brief    Set the parameter n of the Hill function
 * \details  --
 * \param    double hill_function_n
 * \return   \e void
 */
inline void Parameters::set_hill_function_n( double hill_function_n )
{
  assert(hill_function_n >= 0.0);
  _hill_function_n = hill_function_n;
}

/**
 * \brief    Set the protein degradation rate
 * \details  --
 * \param    double protein_degradation_rate
 * \return   \e void
 */
inline void Parameters::set_protein_degradation_rate( double protein_degradation_rate )
{
  assert(protein_degradation_rate >= 0.0);
  assert(protein_degradation_rate <= 1.0);
  _protein_degradation_rate = protein_degradation_rate;
}

/*------------------------------------------------------------------ metabolic network level */

/**
 * \brief    Set metabolism ODE timestep
 * \details  --
 * \param    double metabolism_timestep
 * \return   \e void
 */
inline void Parameters::set_metabolism_timestep( double metabolism_timestep )
{
  assert(metabolism_timestep > 0.0);
  _metabolism_timestep = metabolism_timestep;
}

/**
 * \brief    Set essential metabolites toxicity threshold
 * \details  --
 * \param    double toxicity_threshold
 * \return   \e void
 */
inline void Parameters::set_essential_metabolites_toxicity_threshold( double toxicity_threshold )
{
  assert(toxicity_threshold > 0.0);
  _essential_metabolites_toxicity_threshold = toxicity_threshold;
}

/**
 * \brief    Set non essential metabolites toxicity threshold
 * \details  --
 * \param    double toxic_threshold
 * \return   \e void
 */
inline void Parameters::set_non_essential_metabolites_toxicity_threshold( double toxicity_threshold )
{
  assert(toxicity_threshold > 0.0);
  _non_essential_metabolites_toxicity_threshold = toxicity_threshold;
}

/**
 * \brief    Set initial metabolites amount in cells
 * \details  --
 * \param    double amount
 * \return   \e void
 */
inline void Parameters::set_initial_metabolites_amount_in_cells( double amount )
{
  assert(amount >= 0.0);
  _initial_metabolites_amount_in_cells = amount;
}

/**
 * \brief    Set energy differential factor for enzymatic reation costs
 * \details  --
 * \param    double energy_reaction_cost_factor
 * \return   \e void
 */
inline void Parameters::set_energy_reaction_cost_factor( double energy_reaction_cost_factor )
{
  assert(energy_reaction_cost_factor >= 0.0);
  assert(energy_reaction_cost_factor <= 1.0);
  _energy_reaction_cost_factor = energy_reaction_cost_factor;
}

/**
 * \brief    Set energy pumping cost
 * \details  --
 * \param    double energy_pumping_cost
 * \return   \e void
 */
inline void Parameters::set_energy_pumping_cost( double energy_pumping_cost )
{
  _energy_pumping_cost = energy_pumping_cost;
}

/**
 * \brief    Set energy transcription cost
 * \details  --
 * \param    double energy_transcription_cost
 * \return   \e void
 */
inline void Parameters::set_energy_transcription_cost( double energy_transcription_cost )
{
  _energy_transcription_cost = energy_transcription_cost;
}
/**
 * \brief    Set energy toxicity threshold
 * \details  --
 * \param    double energy_toxicity_threshold
 * \return   \e void
 */
inline void Parameters::set_energy_toxicity_threshold( double energy_toxicity_threshold )
{
  _energy_toxicity_threshold = energy_toxicity_threshold;
}

/*------------------------------------------------------------------ cell level */

/**
 * \brief    Set death probability
 * \details  --
 * \param    double death_probability
 * \return   \e void
 */
inline void Parameters::set_death_probability( double death_probability )
{
  assert(death_probability >= 0.0);
  assert(death_probability <= 1.0);
  _death_probability = death_probability;
}

/**
 * \brief    Set variable lifespan boolean
 * \details  --
 * \param    bool variable_lifespan
 * \return   \e void
 */
inline void Parameters::set_variable_lifespan( bool variable_lifespan )
{
  _variable_lifespan = variable_lifespan;
}

/**
 * \brief    Set Gompertz law parameters
 * \details  --
 * \param    double a
 * \param    double b
 * \param    double c
 * \return   \e void
 */
inline void Parameters::set_gompert_law( double a, double b, double c )
{
  _gompert_law.a = a;
  _gompert_law.b = b;
  _gompert_law.c = c;
}

/**
 * \brief    Set initial membrane permeability
 * \details  Defines the permeability of the cell membrane
 * \param    double permeability
 * \return   \e void
 */
inline void Parameters::set_initial_membrane_permeability( double permeability )
{
  assert(permeability >= 0.0);
  assert(permeability <= 1.0);
  _initial_membrane_permeability = permeability;
}

/*------------------------------------------------------------------ population level */

/**
 * \brief    Set the migration rate
 * \details  --
 * \param    double migration_rate
 * \return   \e void
 */
inline void Parameters::set_migration_rate( double migration_rate )
{
  assert(migration_rate >= 0.0);
  assert(migration_rate <= 1.0);
  _migration_rate = migration_rate;
}

/**
 * \brief    Set the hgt rate
 * \details  --
 * \param    double hgt_rate
 * \return   \e void
 */
inline void Parameters::set_hgt_rate( double hgt_rate )
{
  assert(hgt_rate >= 0.0);
  assert(hgt_rate <= 1.0);
  _hgt_rate = hgt_rate;
}

/*------------------------------------------------------------------ environment level */

/**
 * \brief    Set the environment properties
 * \details  --
 * \param    environment_properties* properties
 * \return   \e void
 */
inline void Parameters::set_environment_properties( const environment_properties* properties )
{
  _environment_properties.number_of_init_cycles = properties->number_of_init_cycles;
  
  _environment_properties.species_tag_range.law    = properties->species_tag_range.law;
  _environment_properties.species_tag_range.min    = properties->species_tag_range.min;
  _environment_properties.species_tag_range.max    = properties->species_tag_range.max;
  _environment_properties.species_tag_range.mu     = properties->species_tag_range.mu;
  _environment_properties.species_tag_range.sigma  = properties->species_tag_range.sigma;
  _environment_properties.species_tag_range.lambda = properties->species_tag_range.lambda;
  
  _environment_properties.concentration_range.law    = properties->concentration_range.law;
  _environment_properties.concentration_range.min    = properties->concentration_range.min;
  _environment_properties.concentration_range.max    = properties->concentration_range.max;
  _environment_properties.concentration_range.mu     = properties->concentration_range.mu;
  _environment_properties.concentration_range.sigma  = properties->concentration_range.sigma;
  _environment_properties.concentration_range.lambda = properties->concentration_range.lambda;
  
  _environment_properties.number_of_species_range.law    = properties->number_of_species_range.law;
  _environment_properties.number_of_species_range.min    = properties->number_of_species_range.min;
  _environment_properties.number_of_species_range.max    = properties->number_of_species_range.max;
  _environment_properties.number_of_species_range.mu     = properties->number_of_species_range.mu;
  _environment_properties.number_of_species_range.sigma  = properties->number_of_species_range.sigma;
  _environment_properties.number_of_species_range.lambda = properties->number_of_species_range.lambda;
  
  _environment_properties.interaction_scheme     = properties->interaction_scheme;
  _environment_properties.renewal_scheme         = properties->renewal_scheme;
  _environment_properties.variation_scheme       = properties->variation_scheme;
  _environment_properties.variation_localization = properties->variation_localization;
  
  _environment_properties.introduction_rate = properties->introduction_rate;
  _environment_properties.diffusion_rate    = properties->diffusion_rate;
  _environment_properties.degradation_rate  = properties->degradation_rate;
}


#endif /* defined(__EVOEVO__Parameters__) */
