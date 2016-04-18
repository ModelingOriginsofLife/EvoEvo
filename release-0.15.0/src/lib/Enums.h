
/**
 * \file      Enums.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Definition of enumerations
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

#ifndef __EVOEVO__Enums__
#define __EVOEVO__Enums__


/******************************************************************************************/

/**
 * \brief   State of the experiment
 * \details Experiment can be freshly created (using create executable) or loaded from backup.
 */
enum experiment_state
{
  NEW_EXPERIMENT = 0, /*!< Freshly created experiment    */
  FROM_BACKUP    = 1  /*!< Loaded from backup experiment */
};

/******************************************************************************************/

/**
 * \brief   Distribution law describing the metabolic range
 * \details --
 */
enum metabolic_range
{
  UNIFORM     = 0, /*!< Uniform distribution                */
  GAUSSIAN    = 1, /*!< Gaussian distribution               */
  EXPONENTIAL = 2  /*!< Exponential decreasing distribution */
};

/******************************************************************************************/

/**
 * \brief   Mutation rates
 * \details Defines mutation rate types index
 */
enum mutation_rate
{
  POINT_MUTATION_RATE                    = 0,  /*!< Point mutation rate                    */
  DUPLICATION_RATE                       = 1,  /*!< Duplication rate                       */
  DELETION_RATE                          = 2,  /*!< Deletion rate                          */
  TRANSLOCATION_RATE                     = 3,  /*!< Translocation rate                     */
  INVERSION_RATE                         = 4,  /*!< Inversion rate                         */
  TRANSITION_RATE                        = 5,  /*!< Transition rate                        */
  SUBSTRATE_TAG_MUTATION_SIZE            = 6,  /*!< Substrate tag mutation size            */
  PRODUCT_TAG_MUTATION_SIZE              = 7,  /*!< Product tag mutation size              */
  KM_MUTATION_SIZE                       = 8,  /*!< Km constant mutation size              */
  KCAT_MUTATION_SIZE                     = 9,  /*!< Kcat constant mutation size            */
  BINDING_SIZE_TAG_MUTATION_SIZE         = 10, /*!< Binding site tag mutation size         */
  CO_ENZYME_TAG_MUTATION_SIZE            = 11, /*!< Co-enzyme tag mutation size            */
  TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE = 12, /*!< Transcription factor tag mutation size */
  BASAL_EXPRESSION_LEVEL_MUTATION_SIZE   = 13  /*!< Basal expression level mutation size   */
};

/******************************************************************************************/

/**
 * \brief   Mutation type
 * \details Defines the type of each mutation event.
 */
enum mutation_type
{
  NONE           = 0, /*!< No event                 */
  POINT_MUTATION = 1, /*!< Point mutation event     */
  DUPLICATION    = 2, /*!< Duplication event        */
  DELETION       = 3, /*!< Deletion event           */
  TRANSLOCATION  = 4, /*!< Translocation event      */
  INVERSION      = 5, /*!< Inversion event          */
  HGT            = 6  /*!< Horizontal gene transfer */
};

/******************************************************************************************/

/**
 * \brief   Pearl type
 * \details Defines the type of pearl.
 */
enum pearl_type
{
  NON_CODING           = 0, /*!< Non coding type (NC)           */
  ENZYME               = 1, /*!< Enzyme type (E)                */
  TRANSCRIPTION_FACTOR = 2, /*!< Transcription factor type (TF) */
  BINDING_SITE         = 3, /*!< Binding site type (BS)         */
  PROMOTER             = 4, /*!< Promoter type (P)              */

  /* Following types are only used in MutationVector class */
  
  NC_TO_E_TRANSITION  = 5,  /*!< Point mutation leading to a transition from NC to E  */
  NC_TO_TF_TRANSITION = 6,  /*!< Point mutation leading to a transition from NC to TF */
  NC_TO_BS_TRANSITION = 7,  /*!< Point mutation leading to a transition from NC to BS */
  NC_TO_P_TRANSITION  = 8,  /*!< Point mutation leading to a transition from NC to P  */
  
  E_TO_NC_TRANSITION  = 9,  /*!< Point mutation leading to a transition from E to NC  */
  E_TO_TF_TRANSITION  = 10, /*!< Point mutation leading to a transition from E to TF  */
  E_TO_BS_TRANSITION  = 11, /*!< Point mutation leading to a transition from E to BS  */
  E_TO_P_TRANSITION   = 12, /*!< Point mutation leading to a transition from E to P   */
  
  TF_TO_NC_TRANSITION = 13, /*!< Point mutation leading to a transition from TF to NC */
  TF_TO_E_TRANSITION  = 14, /*!< Point mutation leading to a transition from TF to E  */
  TF_TO_BS_TRANSITION = 15, /*!< Point mutation leading to a transition from TF to BS */
  TF_TO_P_TRANSITION  = 16, /*!< Point mutation leading to a transition from TF to P  */
  
  BS_TO_NC_TRANSITION = 17, /*!< Point mutation leading to a transition from BS to NC */
  BS_TO_E_TRANSITION  = 18, /*!< Point mutation leading to a transition from BS to E  */
  BS_TO_TF_TRANSITION = 19, /*!< Point mutation leading to a transition from BS to TF */
  BS_TO_P_TRANSITION  = 20, /*!< Point mutation leading to a transition from BS to P  */
  
  P_TO_NC_TRANSITION  = 21, /*!< Point mutation leading to a transition from P to NC  */
  P_TO_E_TRANSITION   = 22, /*!< Point mutation leading to a transition from P to E   */
  P_TO_TF_TRANSITION  = 23, /*!< Point mutation leading to a transition from P to TF  */
  P_TO_BS_TRANSITION  = 24  /*!< Point mutation leading to a transition from P to BS  */
  
};

/******************************************************************************************/

/**
 * \brief   Co-enzyme activity type
 * \details Defines the type of co-enzyme activity
 */
enum co_enzyme_type
{
  REPRESSOR        = 0, /*!< The co-enzyme represses the TF */
  ACTIVATOR        = 1, /*!< The co-enzyme activates the TF */
  ALWAYS_ACTIVE    = 2, /*!< The TF is always active        */
  ALWAYS_REPRESSED = 3  /*!< The TF is always repressed     */
};

/******************************************************************************************/

/**
 * \brief   Reaction type
 * \details Defines the type of the reaction
 */
enum reaction_type
{
  CATALYTIC_ACTIVITY       = 0, /*!< Catalytic activity (transition in the metabolic space)                         */
  INFLOWING_PUMP_ACTIVITY  = 1, /*!< Inflowing pump activity (exchange between the cell and its local environment)  */
  OUTFLOWING_PUMP_ACTIVITY = 2  /*!< Outflowing pump activity (exchange between the cell and its local environment) */
};

/******************************************************************************************/

/**
 * \brief   Event type
 * \details Defines 
 */
enum event_type
{
  BIRTH = 0, /* Cell birth */
  DEATH = 1  /* Cell death */
};

/******************************************************************************************/

/**
 * \brief   Node class
 * \details Defines the class of a node in the tree (master root, root or normal).
 */
enum node_class
{
  MASTER_ROOT = 1, /*!< The node is the master root */
  ROOT        = 2, /*!< The node is a root          */
  NORMAL      = 3  /*!< The node is normal          */
};

/******************************************************************************************/

/**
 * \brief   Node state
 * \details Defines the state of a node in the tree, depending on cell's status (dead or alive).
 */
enum node_state
{
  DEAD  = 1, /*!< The cell is dead  */
  ALIVE = 2  /*!< The cell is alive */
};

/******************************************************************************************/

/**
 * \brief   Cell's state
 * \details Defines the state of a the cell used to update it at each simulation time step
 */
enum cell_state
{
  NOTHING   = 0, /*!< No action to lead        */
  TO_KILL   = 1, /*!< The cell must be killed  */
  TO_UPDATE = 2  /*!< The cell must be updated */
};

/******************************************************************************************/

/**
 * \brief   Score scheme
 * \details Defines the score computation method
 */
enum score_scheme
{
  ESSENTIAL_METABOLITES_SUM                        = 0, /*!< Score is the sum of the essential metabolites           */
  ESSENTIAL_METABOLITES_SUM_MINUS_DEVIATION        = 1, /*!< Score is the sum minus the standard deviation           */
  ESSENTIAL_METABOLITES_COMBINATORIAL_CONTRIBUTION = 2  /*!< Score is the sum of the essential metabolites complexes */
};

/******************************************************************************************/

/**
 * \brief   Trophic level
 * \details Defines the trophic level (0, 1, 2)
 */
enum trophic_level
{
  LEVEL_0  = 0, /*!< Level 0 (exogenous nutrients consumers)          */
  LEVEL_1  = 1, /*!< Level 1 (exogenous and cells products consumers) */
  LEVEL_2  = 2, /*!< Level 2 (cells products consumers only)          */
  NO_LEVEL = 3  /*!< No level (the cell has no inflowing pumps)       */
};

/******************************************************************************************/

/**
 * \brief   Environment interaction scheme
 * \details Defines if the environment is fixed or in interaction with the population
 */
enum environment_interaction_scheme
{
  NO_INTERACTION = 0, /*!< The environment is fixed                      */
  INTERACTION    = 1  /*!< The environment interacts with the population */
};

/******************************************************************************************/

/**
 * \brief   Environment renewal scheme
 * \details Defines if matter is keeped or cleaned at each variation
 */
enum environment_renewal_scheme
{
  KEEP_MATTER  = 0, /*!< Keep matter in the environment */
  CLEAR_MATTER = 1  /*!< Clear the environment          */
};

/******************************************************************************************/

/**
 * \brief   Environment variation scheme
 * \details Defines the type of variation undergone by the environment
 */
enum environment_variation_scheme
{
  RANDOM_SCHEME   = 0, /*!< Variation is at random depending on introduction probability                               */
  PERIODIC_SCHEME = 1, /*!< Variation is periodic, with the 1/(introduction probability) being the period              */
  CYCLIC_SCHEME   = 2  /*!< Variation is cyclic (cosinus like), with the 1/(introduction probability) being the period */
};

/******************************************************************************************/

/**
 * \brief   Environment variation localization
 * \details Defines the localization scheme used to vary the environment
 */
enum environment_variation_localization
{
  RANDOM_LOCALIZATION = 0, /*!< The variation occurs at a random place      */
  GLOBAL_LOCALIZATION = 1, /*!< The variation affects the whole environment */
  CENTER_LOCALIZATION = 2  /*!< The variation occurs only at the center     */
};

/******************************************************************************************/

enum grn_node_type
{
  GRN_TF  = 0, /*!< Genetic regulation network TF node  */
  GRN_E   = 1, /*!< Genetic regulation network E node   */
  GRN_COE = 2  /*!< Genetic regulation network CoE node */
};

/******************************************************************************************/


#endif /* defined(__EVOEVO__Enums__) */
