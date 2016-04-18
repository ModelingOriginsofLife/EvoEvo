
/**
 * \file      Simulation.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Simulation class definition
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

#include "Simulation.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \return   \e void
 */
Simulation::Simulation( Parameters* parameters )
{
  _parameters        = parameters;
  _population        = new Population(_parameters);
  _environment       = new Environment(_parameters);
  _trophic_network   = new TrophicNetwork(_parameters, _population, _environment);
  _lineage_tree      = new Tree(_parameters);
  _phylogenetic_tree = new Tree(_parameters);
  _statistics        = new Statistics(_parameters, _population, _environment, _trophic_network, _lineage_tree, _phylogenetic_tree, true);
  _min_score         = 0.0;
  _max_score         = 0.0;
  _width             = _parameters->get_width();
  _height            = _parameters->get_height();
}

/**
 * \brief    Constructor from backup files
 * \details  Load the simulation from a specified backup file
 * \param    Parameters* parameters
 * \param    size_t backup_time
 * \param    bool clean_statistic_files
 * \return   \e void
 */
Simulation::Simulation( Parameters* parameters, size_t backup_time, bool clean_statistic_files )
{
  /*---------------------------*/
  /* 1) set parameters         */
  /*---------------------------*/
  _parameters = parameters;
  
  /*---------------------------*/
  /* 2) load population        */
  /*---------------------------*/
  std::stringstream pop_file_name;
  pop_file_name << "./population/population_" << backup_time;
  gzFile pop_file = gzopen(pop_file_name.str().c_str(), "r");
  _population = new Population(_parameters, pop_file);
  gzclose(pop_file);
  
  /*---------------------------*/
  /* 3) load environment       */
  /*---------------------------*/
  std::stringstream env_file_name;
  env_file_name << "./environment/environment_" << backup_time;
  gzFile env_file = gzopen(env_file_name.str().c_str(), "r");
  _environment = new Environment(_parameters, env_file);
  gzclose(env_file);
  
  /*---------------------------*/
  /* 4) load trophic network   */
  /*---------------------------*/
  std::stringstream trophic_file_name;
  trophic_file_name << "./trophic_network/trophic_network_" << backup_time;
  gzFile trophic_file = gzopen(trophic_file_name.str().c_str(), "r");
  _trophic_network = new TrophicNetwork(_parameters, _population, _environment, trophic_file);
  gzclose(trophic_file);
  
  /*---------------------------*/
  /* 5) load lineage tree      */
  /*---------------------------*/
  std::stringstream tree_file_name;
  tree_file_name << "./tree/lineage_tree_" << backup_time;
  gzFile tree_file = gzopen(tree_file_name.str().c_str(), "r");
  _lineage_tree = new Tree(_parameters, _population, tree_file);
  gzclose(tree_file);
  
  /*---------------------------*/
  /* 6) load phylogenetic tree */
  /*---------------------------*/
  tree_file_name.str("");
  tree_file_name << "./tree/phylogenetic_tree_" << backup_time;
  tree_file = gzopen(tree_file_name.str().c_str(), "r");
  _phylogenetic_tree = new Tree(_parameters, _population, tree_file);
  gzclose(tree_file);
  
  /*---------------------------*/
  /* 7) initialize statistics  */
  /*---------------------------*/
  _statistics = new Statistics(_parameters, _population, _environment, _trophic_network, _lineage_tree, _phylogenetic_tree, clean_statistic_files);
  _statistics->flush_files();
  
  /*---------------------------*/
  /* 8) initialize others      */
  /*---------------------------*/
  _min_score = 0.0;
  _max_score = 0.0;
  _width     = _parameters->get_width();
  _height    = _parameters->get_height();
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
Simulation::~Simulation( void )
{
  delete _population;
  _population = NULL;
  delete _environment;
  _environment = NULL;
  delete _trophic_network;
  _trophic_network = NULL;
  delete _lineage_tree;
  _lineage_tree = NULL;
  delete _phylogenetic_tree;
  _phylogenetic_tree = NULL;
  delete _statistics;
  _statistics = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Initialize the simulation
 * \details  --
 * \param    experiment_state state
 * \return   \e void
 */
void Simulation::initialize( experiment_state state )
{
  /*-----------------------------------*/
  /* 1) initialize the environment     */
  /*-----------------------------------*/
  if (state == NEW_EXPERIMENT)
  {
    _environment->set_inflowing_amount(0.0);
    _environment->set_outflowing_amount(0.0);
    initialize_environment();
    for (size_t i = 1; i < _parameters->get_environment_properties()->number_of_init_cycles; i++)
    {
      update_environment();
    }
  }
  
  /*-----------------------------------*/
  /* 2) initialize the population      */
  /*-----------------------------------*/
  initialize_population(state);
  if (_parameters->get_migration_rate() > 0.0)
  {
    mix();
  }
  
  /*-----------------------------------*/
  /* 3) initialize the trophic network */
  /*-----------------------------------*/
  if (state == NEW_EXPERIMENT)
  {
    initialize_trophic_network();
    _trophic_network->load_population();
    compute_trophic_groups_diversity();
    _trophic_network->compute_diversity_statistics();
    _trophic_network->compute_level_statistics();
  }
  
  /*-----------------------------------*/
  /* 4) update environment state       */
  /*-----------------------------------*/
  _environment->update_environment_amount();
  
  /*-----------------------------------*/
  /* 5) apply tests                    */
  /*-----------------------------------*/
#ifdef DEBUG
  test_lineage_tree_structure();
  test_phylogenetic_tree_structure();
#endif
}

/**
 * \brief    Initialize the simulation from an evolved one (cells have already been modified)
 * \details  This method should only be called at experiment creation
 * \return   \e void
 */
void Simulation::initialize_from_evolved_population( void )
{
  /*-----------------------------------*/
  /* 1) initialize the environment     */
  /*-----------------------------------*/
  _environment->set_inflowing_amount(0.0);
  _environment->set_outflowing_amount(0.0);
  initialize_environment();
  for (size_t i = 1; i < _parameters->get_environment_properties()->number_of_init_cycles; i++)
  {
    update_environment();
  }
  
  /*-----------------------------------*/
  /* 2) initialize the population      */
  /*-----------------------------------*/
  initialize_population_from_evolved_population();
  if (_parameters->get_migration_rate() > 0.0)
  {
    mix();
  }
  
  /*-----------------------------------*/
  /* 3) initialize the trophic network */
  /*-----------------------------------*/
  initialize_trophic_network();
  _trophic_network->load_population();
  compute_trophic_groups_diversity();
  _trophic_network->compute_diversity_statistics();
  _trophic_network->compute_level_statistics();
  
  /*-----------------------------------*/
  /* 4) update environment state       */
  /*-----------------------------------*/
  _environment->update_environment_amount();
  
  /*-----------------------------------*/
  /* 5) apply tests                    */
  /*-----------------------------------*/
#ifdef DEBUG
  test_lineage_tree_structure();
  test_phylogenetic_tree_structure();
#endif
}

/**
 * \brief    Update the simulation
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::update( void )
{
  /*-------------------------------*/
  /* 1) update the environment     */
  /*-------------------------------*/
  if (_population->get_time() > 0)
  {
    _environment->set_inflowing_amount(0.0);
    _environment->set_outflowing_amount(0.0);
    update_environment();
  }
  
  /*-------------------------------*/
  /* 2) update the population      */
  /*-------------------------------*/
  _population->set_previous_size();
  update_population();
  _population->compute_growth_rate();
  if (_parameters->get_migration_rate() > 0.0)
  {
    mix();
  }
  
  /*-------------------------------*/
  /* 3) update the trophic network */
  /*-------------------------------*/
  _trophic_network->load_population();
  compute_trophic_groups_diversity();
  _trophic_network->compute_diversity_statistics();
  _trophic_network->compute_level_statistics();
  
  /*-------------------------------*/
  /* 4) update environment state   */
  /*-------------------------------*/
  _environment->update_environment_amount();
  /*
  if (_pop->get_time()%500 == 0)
  {
    evaluate_species_lists_size();
  }
  */
  
  /*-------------------------------*/
  /* 5) apply tests                */
  /*-------------------------------*/
#ifdef DEBUG
  test_lineage_tree_structure();
  test_phylogenetic_tree_structure();
#endif
  
}

/**
 * \brief    Save experiment in backup
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::save_experiment( void )
{
  /*-----------------------------*/
  /* 1) save the parameters      */
  /*-----------------------------*/
  _parameters->save(_population->get_time());
  
  /*-----------------------------*/
  /* 2) save the population      */
  /*-----------------------------*/
  std::stringstream pop_file_name;
  pop_file_name << "./population/population_" << _population->get_time();
  gzFile pop_file = gzopen(pop_file_name.str().c_str(), "w");
  _population->save(pop_file);
  gzclose(pop_file);
  
  /*-----------------------------*/
  /* 3) save the environment     */
  /*-----------------------------*/
  std::stringstream env_file_name;
  env_file_name << "./environment/environment_" << _population->get_time();
  gzFile env_file = gzopen(env_file_name.str().c_str(), "w");
  _environment->save(env_file);
  gzclose(env_file);
  
  /*-----------------------------*/
  /* 4) save the trophic network */
  /*-----------------------------*/
  std::stringstream trophic_file_name;
  trophic_file_name << "./trophic_network/trophic_network_" << _population->get_time();
  gzFile trophic_file = gzopen(trophic_file_name.str().c_str(), "w");
  _trophic_network->save(trophic_file);
  gzclose(trophic_file);
}

/**
 * \brief    Save trees in backup
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::save_trees( void )
{
  /*-------------------------------*/
  /* 1) save the lineage tree      */
  /*-------------------------------*/
  std::stringstream tree_file_name;
  tree_file_name << "./tree/lineage_tree_" << _population->get_time();
  gzFile tree_file = gzopen(tree_file_name.str().c_str(), "w");
  _lineage_tree->save(tree_file);
  gzclose(tree_file);
  
  /*-------------------------------*/
  /* 2) save the phylogenetic tree */
  /*-------------------------------*/
  tree_file_name.str("");
  tree_file_name << "./tree/phylogenetic_tree_" << _population->get_time();
  tree_file = gzopen(tree_file_name.str().c_str(), "w");
  _phylogenetic_tree->save(tree_file);
  gzclose(tree_file);
}

/**
 * \brief    Add random pearls to each individual, and shuffle it if shuffle is true
 * \details  --
 * \param    size_t NC_type
 * \param    size_t E_type
 * \param    size_t TF_type
 * \param    size_t BS_type
 * \param    size_t P_type
 * \param    bool shuffle
 * \return   \e void
 */
void Simulation::add_random_pearls( size_t NC_type, size_t E_type, size_t TF_type, size_t BS_type, size_t P_type, bool shuffle )
{
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    Cell* cell = _population->get_cell(pos);
    for (size_t i = 0; i < NC_type; i++)
    {
      pearl p;
      initialize_pearl(p);
      draw_random_pearl(p, NON_CODING);
      cell->get_genome()->add_pearl(&p);
    }
    for (size_t i = 0; i < E_type; i++)
    {
      pearl p;
      initialize_pearl(p);
      draw_random_pearl(p, ENZYME);
      cell->get_genome()->add_pearl(&p);
    }
    for (size_t i = 0; i < TF_type; i++)
    {
      pearl p;
      initialize_pearl(p);
      draw_random_pearl(p, TRANSCRIPTION_FACTOR);
      cell->get_genome()->add_pearl(&p);
    }
    for (size_t i = 0; i < BS_type; i++)
    {
      pearl p;
      initialize_pearl(p);
      draw_random_pearl(p, BINDING_SITE);
      cell->get_genome()->add_pearl(&p);
    }
    for (size_t i = 0; i < P_type; i++)
    {
      pearl p;
      initialize_pearl(p);
      draw_random_pearl(p, PROMOTER);
      cell->get_genome()->add_pearl(&p);
    }
    if (shuffle)
    {
      cell->get_genome()->shuffle();
    }
  }
}

/**
 * \brief    Write the clustering data
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::write_clustering_data( void )
{
  /*---------------------*/
  /* 1) Open data file   */
  /*---------------------*/
  std::stringstream filename;
  filename << "clustering_data_" << _population->get_time() << ".txt";
  std::ofstream file(filename.str(), std::ios::out | std::ios::trunc);
  
  /*---------------------*/
  /* 2) Write the header */
  /*---------------------*/
  file << "t" << " " << "id" << " " <<
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
  "trophic_group" << " " <<
  "trophic_level" << " " <<
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
  "mixed_operon_size" << " " <<
  "inherited_size" << " " <<
  "inherited_nb_E" << " " <<
  "inherited_nb_TF" << " " <<
  "inherited_nb_inner_enzymes" << " " <<
  "inherited_nb_inflow_pumps" << " " <<
  "inherited_nb_outflow_pumps";
  size_t N = _environment->get_species_lists_size();
  for (size_t i = 0; i < N; i++)
  {
    file << " p" << i+1;
  }
  for (size_t i = 0; i < N; i++)
  {
    file << " u" << i+1;
  }
  for (size_t i = 0; i < N; i++)
  {
    file << " r" << i+1;
  }
  file << "\n";
  
  /*---------------------*/
  /* 3) Write the data   */
  /*---------------------*/
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    Cell* cell = _population->get_cell(pos);
    if (cell->isAlive())
    {
      std::stringstream line;
      
      line << _population->get_time() << " " <<
      cell->get_id() << " " <<
      
      /* PHENOTYPE */
      
      cell->get_generation() << " " <<
      cell->get_inherited_TF_amount() << " " <<
      cell->get_inherited_E_amount() << " " <<
      cell->get_TF_amount() << " " <<
      cell->get_E_amount() << " " <<
      cell->get_inherited_metabolic_amount() << " " <<
      cell->get_species_list()->get_amount() << " " <<
      cell->get_energy() << " " <<
      cell->get_score() << " " <<
      cell->get_lifespan() << " " <<
      cell->get_number_of_divisions() << " " <<
      cell->get_toxicity() << " " <<
      cell->get_metabolic_uptake() << " " <<
      cell->get_metabolic_release() << " " <<
      cell->get_metabolic_growth_rate() << " " <<
      cell->get_Dmetabolic_growth_rate() << " " <<
      cell->get_grn_nb_nodes() << " " <<
      cell->get_grn_nb_edges() << " " <<
      cell->get_metabolic_nb_nodes()<< " " <<
      cell->get_metabolic_nb_edges() << " " <<
      cell->get_replication_report()->get_mean_regulation_redundancy() << " " <<
      cell->get_replication_report()->get_mean_metabolic_redundancy() << " " <<
      cell->get_trophic_group() << " " <<
      cell->get_trophic_level() << " " <<
      
      /* GENOME STRUCTURE */
      
      cell->get_genome()->get_size() << " " <<
      cell->get_replication_report()->get_genome_functional_size() << " " <<
      cell->get_genome()->get_nb_NC() << " " <<
      cell->get_genome()->get_nb_E() << " " <<
      cell->get_genome()->get_nb_TF() << " " <<
      cell->get_genome()->get_nb_BS() << " " <<
      cell->get_genome()->get_nb_P() << " " <<
      cell->get_genome()->get_nb_inner_enzymes() << " " <<
      cell->get_genome()->get_nb_inflow_pumps() << " " <<
      cell->get_genome()->get_nb_outflow_pumps() << " " <<
      cell->get_replication_report()->get_nb_functional_regions() << " " <<
      cell->get_replication_report()->get_nb_enhancers() << " " <<
      cell->get_replication_report()->get_nb_operators() << " " <<
      cell->get_replication_report()->get_nb_E_regions() << " " <<
      cell->get_replication_report()->get_nb_TF_regions() << " " <<
      cell->get_replication_report()->get_nb_mixed_regions() << " " <<
      cell->get_replication_report()->get_mean_functional_region_size() << " " <<
      cell->get_replication_report()->get_mean_E_region_size() << " " <<
      cell->get_replication_report()->get_mean_TF_region_size() << " " <<
      cell->get_replication_report()->get_mean_mixed_region_size() << " " <<
      cell->get_replication_report()->get_mean_enhancer_size() << " " <<
      cell->get_replication_report()->get_mean_operator_size() << " " <<
      cell->get_replication_report()->get_mean_operon_size() << " " <<
      cell->get_replication_report()->get_mean_E_operon_size() << " " <<
      cell->get_replication_report()->get_mean_TF_operon_size() << " " <<
      cell->get_replication_report()->get_mean_mixed_operon_size() << " " <<
      
      /* INHERITED STRUCTURE */
      
      cell->get_inherited_proteins()->get_size() << " " <<
      cell->get_inherited_proteins()->get_nb_E() << " " <<
      cell->get_inherited_proteins()->get_nb_TF() << " " <<
      cell->get_inherited_proteins()->get_nb_inner_enzymes() << " " <<
      cell->get_inherited_proteins()->get_nb_inflow_pumps() << " " <<
      cell->get_inherited_proteins()->get_nb_outflow_pumps();
      
      /* TROPHIC PROFILE */
      
      TrophicGroup* group = _trophic_network->get_group(cell->get_trophic_group());
      for (size_t i = 0; i < group->get_trophic_profile().size(); i++)
      {
        line << " " << group->get_trophic_profile()[i];
      }
      line << "\n";
      
      /* WRITE LINE */
      
      file << line.str();
    }
  }
  
  file.close();
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Initialize environment
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::initialize_environment( void )
{
  bool initialize = true;
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) If the introduction is at a random time (RANDOM_SCHEME) */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  if (_parameters->get_environment_properties()->variation_scheme == RANDOM_SCHEME)
  {
    /*--------------------------------------*/
    /* 1.A) Draw in the event probability   */
    /*--------------------------------------*/
    bool apply_change = (_parameters->get_prng()->uniform() < _parameters->get_environment_properties()->introduction_rate);
    
    /*--------------------------------------*/
    /* 1.B) Clear the environment if needed */
    /*--------------------------------------*/
    if (apply_change && _parameters->get_environment_properties()->renewal_scheme == CLEAR_MATTER)
    {
      _environment->clear_environment();
    }
    
    /*--------------------------------------*/
    /* 1.C) Apply environmental change      */
    /*--------------------------------------*/
    if (apply_change && _parameters->get_environment_properties()->variation_localization == RANDOM_LOCALIZATION)
    {
      update_environment_random_localization(initialize, 1.0);
    }
    else if (apply_change && _parameters->get_environment_properties()->variation_localization == GLOBAL_LOCALIZATION)
    {
      update_environment_global_localization(initialize, 1.0);
    }
    else if (apply_change && _parameters->get_environment_properties()->variation_localization == CENTER_LOCALIZATION)
    {
      update_environment_center_localization(initialize, 1.0);
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) If the introduction is periodic (PERIODIC_SCHEME)       */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  else if (_parameters->get_environment_properties()->variation_scheme == PERIODIC_SCHEME)
  {
    /*--------------------------------------*/
    /* 1.A) Evaluate the period of change   */
    /*--------------------------------------*/
    int  period       = (int)1.0/_parameters->get_environment_properties()->introduction_rate;
    bool apply_change = (_population->get_time()%period == 0);
    
    /*--------------------------------------*/
    /* 1.B) Clear the environment if needed */
    /*--------------------------------------*/
    if (apply_change && _parameters->get_environment_properties()->renewal_scheme == CLEAR_MATTER)
    {
      _environment->clear_environment();
    }
    
    /*--------------------------------------*/
    /* 1.C) Apply environmental change      */
    /*--------------------------------------*/
    if (apply_change && _parameters->get_environment_properties()->variation_localization == RANDOM_LOCALIZATION)
    {
      update_environment_random_localization(initialize, 1.0);
    }
    else if (apply_change && _parameters->get_environment_properties()->variation_localization == GLOBAL_LOCALIZATION)
    {
      update_environment_global_localization(initialize, 1.0);
    }
    else if (apply_change && _parameters->get_environment_properties()->variation_localization == CENTER_LOCALIZATION)
    {
      update_environment_center_localization(initialize, 1.0);
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) If the introduction is cyclic (CYCLIC_SCHEME)           */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  else if (_parameters->get_environment_properties()->variation_scheme == CYCLIC_SCHEME)
  {
    /*--------------------------------------*/
    /* 1.A) Evaluate the cyclic factor      */
    /*--------------------------------------*/
    double period = 1.0/_parameters->get_environment_properties()->introduction_rate;
    double factor = (sin(_population->get_time()*2.0*3.14159265359/period)+1.0)/2.0;
    
    /*--------------------------------------*/
    /* 1.B) Clear the environment if needed */
    /*--------------------------------------*/
    if (_parameters->get_environment_properties()->renewal_scheme == CLEAR_MATTER)
    {
      _environment->clear_environment();
    }
    
    /*--------------------------------------*/
    /* 1.C) Apply environmental change      */
    /*--------------------------------------*/
    if (_parameters->get_environment_properties()->variation_localization == RANDOM_LOCALIZATION)
    {
      update_environment_random_localization(initialize, factor);
    }
    else if (_parameters->get_environment_properties()->variation_localization == GLOBAL_LOCALIZATION)
    {
      update_environment_global_localization(initialize, factor);
    }
    else if (_parameters->get_environment_properties()->variation_localization == CENTER_LOCALIZATION)
    {
      update_environment_center_localization(initialize, factor);
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4) Compute diffusion and degradation                       */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  _environment->compute_diffusion_and_degradation();
}

/**
 * \brief    Initialize population
 * \details  --
 * \param    experiment_state state
 * \return   \e void
 */
void Simulation::initialize_population( experiment_state state )
{
  _min_score = 1e+6;
  _max_score = 0.0;
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    Cell* cell = _population->get_cell(pos);
    /*-------------------------------------*/
    /* A) if experiment is freshly created */
    /*-------------------------------------*/
    if (state == NEW_EXPERIMENT)
    {
      /*--------------------------------------------------*/
      /* A.1) set coordinates                             */
      /*--------------------------------------------------*/
      size_t y = pos%_height;
      size_t x = (pos-y)/_height;
      cell->set_x(x);
      cell->set_y(y);
      
      /*--------------------------------------------------*/
      /* A.2) build the genome                            */
      /*--------------------------------------------------*/
      for (size_t i = 0; i < _parameters->get_initial_number_of_NC_pearls(); i++)
      {
        pearl p;
        initialize_pearl(p);
        draw_random_pearl(p, NON_CODING);
        cell->get_genome()->add_pearl(&p);
      }
      for (size_t i = 0; i < _parameters->get_initial_number_of_E_pearls(); i++)
      {
        pearl p;
        initialize_pearl(p);
        draw_random_pearl(p, ENZYME);
        cell->get_genome()->add_pearl(&p);
      }
      for (size_t i = 0; i < _parameters->get_initial_number_of_TF_pearls(); i++)
      {
        pearl p;
        initialize_pearl(p);
        draw_random_pearl(p, TRANSCRIPTION_FACTOR);
        cell->get_genome()->add_pearl(&p);
      }
      for (size_t i = 0; i < _parameters->get_initial_number_of_BS_pearls(); i++)
      {
        pearl p;
        initialize_pearl(p);
        draw_random_pearl(p, BINDING_SITE);
        cell->get_genome()->add_pearl(&p);
      }
      for (size_t i = 0; i < _parameters->get_initial_number_of_P_pearls(); i++)
      {
        pearl p;
        initialize_pearl(p);
        draw_random_pearl(p, PROMOTER);
        cell->get_genome()->add_pearl(&p);
      }
      cell->get_genome()->shuffle();
      
      /*--------------------------------------------------*/
      /* A.3) mutate the genome                           */
      /*--------------------------------------------------*/
      cell->mutate();
      
      /*--------------------------------------------------*/
      /* A.4) initialize metabolic amounts                */
      /*--------------------------------------------------*/
      int    tag  = draw_species_tag(_parameters->get_metabolite_tag_initial_range());
      double conc = _parameters->get_initial_metabolites_amount_in_cells();
      cell->get_species_list()->add(tag, conc);
      cell->initialize_inherited_species_list();
      
      /*--------------------------------------------------*/
      /* A.5) load genome in species list and environment */
      /*--------------------------------------------------*/
      cell->load_genome_in_species_lists();
      
      /*--------------------------------------------------*/
      /* A.6) initialize energy                           */
      /*--------------------------------------------------*/
      cell->set_energy(0.0);
      
      /*--------------------------------------------------*/
      /* A.7) activate the cell                           */
      /*--------------------------------------------------*/
      cell->activate();
      cell->set_alive(true);
      
      /*--------------------------------------------------*/
      /* A.8) initialize time variables                   */
      /*--------------------------------------------------*/
      cell->set_birth_time(_population->get_time());
      
      /*--------------------------------------------------*/
      /* A.9) set cell id                                 */
      /*--------------------------------------------------*/
      cell->set_parent_id(0);
      cell->set_id(_population->get_new_id());
      
      /*--------------------------------------------------*/
      /* A.10) compute score                              */
      /*--------------------------------------------------*/
      cell->synchronize_state_vectors(_environment);
      cell->load_genome_in_ODE_system(_environment, false, true);
      cell->update(_population->get_time());
      if (_min_score > cell->get_score())
      {
        _min_score = cell->get_score();
      }
      if (_max_score < cell->get_score())
      {
        _max_score = cell->get_score();
        _population->set_best(cell->get_id(), pos);
      }
      
      /*--------------------------------------------------*/
      /* A.11) Add the cell as a root in trees            */
      /*--------------------------------------------------*/
      if (_parameters->get_tree_backup_step() > 0)
      {
        _lineage_tree->add_root(cell);
        _phylogenetic_tree->add_root(cell);
      }
    }
    /*------------------------------------------------*/
    /* B) else if experiment is loaded from backup    */
    /*------------------------------------------------*/
    else if (state == FROM_BACKUP)
    {
      /*-----------------------------------------------------*/
      /* B.1) Synchronize cell and environment state vectors */
      /*-----------------------------------------------------*/
      cell->synchronize_state_vectors(_environment);
      /*-----------------------------------------------------*/
      /* B.2) load genome in ODE system                      */
      /*-----------------------------------------------------*/
      cell->load_genome_in_ODE_system(_environment, true, false);
    }
  }
}

/**
 * \brief    Initialize the population from an evolved population (cells have already been modified)
 * \details  This method should only be called at experiment creation
 * \return   \e void
 */
void Simulation::initialize_population_from_evolved_population( void )
{
  _min_score = 1e+6;
  _max_score = 0.0;
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    Cell* cell = _population->get_cell(pos);
    
    /*------------------------------------------------*/
    /* 1) set coordinates                             */
    /*------------------------------------------------*/
    size_t y = pos%_height;
    size_t x = (pos-y)/_height;
    cell->set_x(x);
    cell->set_y(y);
    
    /*------------------------------------------------*/
    /* 2) load genome in species list and environment */
    /*------------------------------------------------*/
    cell->load_genome_in_species_lists();
    
    /*------------------------------------------------*/
    /* 3) activate the cell                           */
    /*------------------------------------------------*/
    cell->activate();
    cell->set_alive(true);
    
    /*------------------------------------------------*/
    /* 4) initialize time variables                   */
    /*------------------------------------------------*/
    cell->set_birth_time(_population->get_time());
    
    /*------------------------------------------------*/
    /* 5) set cell id                                 */
    /*------------------------------------------------*/
    cell->set_parent_id(0);
    cell->set_id(_population->get_new_id());
    
    /*------------------------------------------------*/
    /* 6) compute score                               */
    /*------------------------------------------------*/
    cell->synchronize_state_vectors(_environment);
    cell->load_genome_in_ODE_system(_environment, false, true);
    cell->update(_population->get_time());
    if (_min_score > cell->get_score())
    {
      _min_score = cell->get_score();
    }
    if (_max_score < cell->get_score())
    {
      _max_score = cell->get_score();
      _population->set_best(cell->get_id(), pos);
    }
    
    /*------------------------------------------------*/
    /* 7) Add the cell as a root in trees             */
    /*------------------------------------------------*/
    if (_parameters->get_tree_backup_step() > 0)
    {
      _lineage_tree->add_root(cell);
      _phylogenetic_tree->add_root(cell);
    }
  }
}

/**
 * \brief    Initialize trophic network
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::initialize_trophic_network( void )
{
  _trophic_network->initialize_trophic_network();
}

/**
 * \brief    Initialize a default pearl
 * \details  --
 * \param    pearl& p
 * \return   \e void
 */
void Simulation::initialize_pearl( pearl& p )
{
  /*------------------------------------------------------------------ Global attributes */
  
  p.type              = NON_CODING;
  p.identifier        = 0;
  p.parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  p.s    = 1;
  p.p    = 1;
  p.km   = pow(10.0, KM_MIN_LOG);
  p.kcat = pow(10.0, KCAT_MIN_LOG);
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  p.BS_tag         = 0;
  p.coE_tag        = 1;
  p.free_activity  = false;
  p.bound_activity = false;
  p.window         = 0;
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  p.TF_tag = 0;
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  p.basal_expression_level = 0.0;
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  p.functional = false;
}

/**
 * \brief    Draw pearl's attributes at random
 * \details  --
 * \param    pearl& p
 * \param    pearl_type type
 * \return   \e void
 */
void Simulation::draw_random_pearl( pearl& p, pearl_type type )
{
  /*------------------------------------------------------------------ Global attributes */
  
  p.type              = type;
  p.identifier        = 0;
  p.parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  p.s    = draw_species_tag(_parameters->get_metabolite_tag_initial_range());
  p.s    = (p.s > 0 ? p.s : 1);
  p.p    = draw_species_tag(_parameters->get_metabolite_tag_initial_range());
  p.p    = (p.p > 0 ? p.p : 1);
  p.km   = pow(10.0, _parameters->get_prng()->uniform()*(KM_MAX_LOG-KM_MIN_LOG)+KM_MIN_LOG);
  p.kcat = pow(10.0, _parameters->get_prng()->uniform()*(KCAT_MAX_LOG-KCAT_MIN_LOG)+KCAT_MIN_LOG);
  p.kcat = p.kcat*(_parameters->get_prng()->uniform() < 0.5 ? -1.0 : 1.0);
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  p.BS_tag         = draw_species_tag(_parameters->get_binding_site_tag_initial_range());
  p.coE_tag        = draw_species_tag(_parameters->get_co_enzyme_tag_initial_range());
  p.coE_tag        = (p.coE_tag > 0 ? p.coE_tag : 1);
  p.free_activity  = (_parameters->get_prng()->uniform() < 0.5 ? true : false);
  p.bound_activity = (_parameters->get_prng()->uniform() < 0.5 ? true : false);
  p.window         = _parameters->get_initial_binding_window();
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  p.TF_tag = draw_species_tag(_parameters->get_transcription_factor_tag_initial_range());
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  p.basal_expression_level = _parameters->get_prng()->uniform();
}

/**
 * \brief    Update environment
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::update_environment( void )
{
  bool initialize = false;
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) If the introduction is at a random time (RANDOM_SCHEME) */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  if (_parameters->get_environment_properties()->variation_scheme == RANDOM_SCHEME)
  {
    /*--------------------------------------*/
    /* 1.A) Draw in the event probability   */
    /*--------------------------------------*/
    bool apply_change = (_parameters->get_prng()->uniform() < _parameters->get_environment_properties()->introduction_rate);
    
    /*--------------------------------------*/
    /* 1.B) Clear the environment if needed */
    /*--------------------------------------*/
    if (apply_change && _parameters->get_environment_properties()->renewal_scheme == CLEAR_MATTER)
    {
      _environment->clear_environment();
    }
    
    /*--------------------------------------*/
    /* 1.C) Apply environmental change      */
    /*--------------------------------------*/
    if (apply_change && _parameters->get_environment_properties()->variation_localization == RANDOM_LOCALIZATION)
    {
      update_environment_random_localization(initialize, 1.0);
    }
    else if (apply_change && _parameters->get_environment_properties()->variation_localization == GLOBAL_LOCALIZATION)
    {
      update_environment_global_localization(initialize, 1.0);
    }
    else if (apply_change && _parameters->get_environment_properties()->variation_localization == CENTER_LOCALIZATION)
    {
      update_environment_center_localization(initialize, 1.0);
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) If the introduction is periodic (PERIODIC_SCHEME)       */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  else if (_parameters->get_environment_properties()->variation_scheme == PERIODIC_SCHEME)
  {
    /*--------------------------------------*/
    /* 1.A) Evaluate the period of change   */
    /*--------------------------------------*/
    int  period       = (int)1.0/_parameters->get_environment_properties()->introduction_rate;
    bool apply_change = (_population->get_time()%period == 0);
    
    /*--------------------------------------*/
    /* 1.B) Clear the environment if needed */
    /*--------------------------------------*/
    if (apply_change && _parameters->get_environment_properties()->renewal_scheme == CLEAR_MATTER)
    {
      _environment->clear_environment();
    }
    
    /*--------------------------------------*/
    /* 1.C) Apply environmental change      */
    /*--------------------------------------*/
    if (apply_change && _parameters->get_environment_properties()->variation_localization == RANDOM_LOCALIZATION)
    {
      update_environment_random_localization(initialize, 1.0);
    }
    else if (apply_change && _parameters->get_environment_properties()->variation_localization == GLOBAL_LOCALIZATION)
    {
      update_environment_global_localization(initialize, 1.0);
    }
    else if (apply_change && _parameters->get_environment_properties()->variation_localization == CENTER_LOCALIZATION)
    {
      update_environment_center_localization(initialize, 1.0);
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) If the introduction is cyclic (CYCLIC_SCHEME)           */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  else if (_parameters->get_environment_properties()->variation_scheme == CYCLIC_SCHEME)
  {
    /*--------------------------------------*/
    /* 1.A) Evaluate the cyclic factor      */
    /*--------------------------------------*/
    double period = 1.0/_parameters->get_environment_properties()->introduction_rate;
    double factor = (sin(_population->get_time()*2.0*3.14159265359/period)+1.0)/2.0;
    
    /*--------------------------------------*/
    /* 1.B) Clear the environment if needed */
    /*--------------------------------------*/
    if (_parameters->get_environment_properties()->renewal_scheme == CLEAR_MATTER)
    {
      _environment->clear_environment();
    }
    
    /*--------------------------------------*/
    /* 1.C) Apply environmental change      */
    /*--------------------------------------*/
    if (_parameters->get_environment_properties()->variation_localization == RANDOM_LOCALIZATION)
    {
      update_environment_random_localization(initialize, factor);
    }
    else if (_parameters->get_environment_properties()->variation_localization == GLOBAL_LOCALIZATION)
    {
      update_environment_global_localization(initialize, factor);
    }
    else if (_parameters->get_environment_properties()->variation_localization == CENTER_LOCALIZATION)
    {
      update_environment_center_localization(initialize, factor);
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 4) Compute diffusion and degradation                       */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  _environment->compute_diffusion_and_degradation();
}

/**
 * \brief    Update population
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::update_population( void )
{
  /*-------------------------------------------------------------------------*/
  /* 1) initialize population size and statistics                            */
  /*-------------------------------------------------------------------------*/
  _population->set_population_size(0);
  _statistics->init_variables();
  _min_score              = 1e+6;
  _max_score              = 0.0;
  cell_state* cell_states = new cell_state[_width*_height];
  
  /*-------------------------------------------------------------------------*/
  /* 2) explore population: save gaps, kill or update cells unable to divide */
  /*-------------------------------------------------------------------------*/
  std::vector<size_t> gaps;
  
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    Cell* cell = _population->get_cell(pos);
    cell->untag();
    cell_states[pos] = NOTHING;
    
    /*---------------------------------------------------------*/
    /* 2.1) if the cell is alive                               */
    /*---------------------------------------------------------*/
    if (cell->isAlive())
    {
      /* A) increase population size --------*/
      _population->increase_population_size(1);
      
      /* B) compute statistics --------------*/
      _statistics->add_individual(pos);
      if (_min_score > cell->get_score())
      {
        _min_score = cell->get_score();
      }
      if (_max_score < cell->get_score())
      {
        _max_score = cell->get_score();
        _population->set_best(cell->get_id(), pos);
        _statistics->add_best_individual();
      }
      
      /* C) evaluate the individual ---------*/
      size_t child_position = 0;
      bool isActive         = cell->isActive();
      bool amount           = (cell->get_amount() > MINIMUM_CONCENTRATION);
      bool score            = (cell->get_score() > MINIMUM_SCORE);
      bool energy           = false;
      if ((_parameters->get_energy_constraints() && cell->get_energy() > 0.0) || !_parameters->get_energy_constraints())
      {
        energy = true;
      }
      bool dying = false;
      if (_parameters->get_variable_lifespan())
      {
        double X       = _population->get_time()-cell->get_birth_time()+cell->get_number_of_divisions()+cell->get_toxicity();
        double p_death = _parameters->get_gompert_law()->a*exp(-_parameters->get_gompert_law()->b*exp(-_parameters->get_gompert_law()->c*X));
        dying = (_parameters->get_prng()->uniform() < p_death || !energy || !amount || !isActive || !score);
      }
      else
      {
        dying = (_parameters->get_prng()->uniform() < _parameters->get_death_probability() || !energy || !amount || !isActive || !score);
      }
      bool updating = !find_empty_cases(pos, child_position);
      
      /* D) update or kill individual -------*/
      if (dying)
      {
        cell_states[pos] = TO_KILL;
        cell->tag();
      }
      else if (!dying && updating)
      {
        cell_states[pos] = TO_UPDATE;
        cell->tag();
      }
      else if (!dying && !updating)
      {
        cell_states[pos] = TO_UPDATE;
      }
    }
    /*---------------------------------------------------------*/
    /* 2.2) else if the cell is dead, add gap to the gaps list */
    /*---------------------------------------------------------*/
    else
    {
      gaps.push_back(pos);
    }
  }
  
  /*-------------------------------------------------------------------------*/
  /* 3) write best individual genome and metabolic network                   */
  /*-------------------------------------------------------------------------*/
  _statistics->compute_mean_and_var();
  
  /*-------------------------------------------------------------------------*/
  /* 4) shuffle gaps vector                                                  */
  /*-------------------------------------------------------------------------*/
  int size = (int)gaps.size();
  for (int i = 0; i < size; i++)
  {
    size_t pos1 = _parameters->get_prng()->uniform(0, size-1);
    size_t pos2 = _parameters->get_prng()->uniform(0, size-1);
    size_t tmp  = gaps[pos1];
    gaps[pos1]  = gaps[pos2];
    gaps[pos2]  = tmp;
  }
  
  /*-------------------------------------------------------------------------*/
  /* 5) explore gaps and compute cell divisions                              */
  /*-------------------------------------------------------------------------*/
  double th = _parameters->get_selection_threshold();
  
  /* A) for each gap -------------------------------------------------------*/
  for (size_t gap_index = 0; gap_index < gaps.size(); gap_index++)
  {
    size_t y_gap = gaps[gap_index]%_height;
    size_t x_gap = (gaps[gap_index]-y_gap)/_height;
    size_t selected_position = 0;
    double propensity = 0.0, best_propensity = 0.0;
    
    /* B) explore gap moore neighbourhood ----------------------------------*/
    for (int i = -1; i < 2; i++)
    {
      for (int j = -1; j < 2; j++)
      {
        if ( !(i == 0 && j == 0) )
        {
          size_t x     = (x_gap + i + _width) % _width;
          size_t y     = (y_gap + j + _height) % _height;
          Cell*  cell  = _population->get_cell(x, y);
          double score = cell->get_score();
          
          /* C) if the cell is able to divide itself, save its coordinates -*/
          if (cell->isActive() && !cell->isTagged() && score > th*_max_score)
          {
            propensity = score;
            if (best_propensity < propensity)
            {
              best_propensity   = propensity;
              selected_position = x*_height+y;
            }
          }
        }
      }
    }
    
    /* D) if best propensity is above threshold, divide the cell ------------*/
    if (best_propensity > 0.0)
    {
      /* Compute cell replication */
      divide(selected_position, gaps[gap_index]);
      
      /* Update parent cell ------*/
      _population->get_cell(selected_position)->update_number_of_divisions();
      cell_states[selected_position] = TO_UPDATE;
      _population->get_cell(selected_position)->tag();
      
      /* Update daughter cell ----*/
      cell_states[gaps[gap_index]] = TO_UPDATE;
      _population->get_cell(gaps[gap_index])->tag();
    }
  }
  
  /*-------------------------------------------------------------------------*/
  /* 6) update cell states                                                   */
  /*-------------------------------------------------------------------------*/
  if (!_parameters->get_parallel_computing())
  {
    for (size_t pos = 0; pos < _parameters->get_width()*_parameters->get_height(); pos++)
    {
      if (cell_states[pos] == TO_KILL)
      {
        kill(pos);
      }
      else if (cell_states[pos] == TO_UPDATE)
      {
        update_cell(pos);
      }
    }
  }
  else if (_parameters->get_parallel_computing())
  {
    tbb::task_group tasks;
    for (size_t pos = 0; pos < _width*_height; pos++)
    {
      if (cell_states[pos] == TO_KILL)
      {
        tasks.run([=]{kill(pos);});
      }
      else if (cell_states[pos] == TO_UPDATE)
      {
        tasks.run([=]{update_cell(pos);});
      }
    }
    tasks.wait();
  }
  delete cell_states;
  cell_states = NULL;
  
  /*-------------------------------------------------------------------------*/
  /* 7) update trees                                                         */
  /*-------------------------------------------------------------------------*/
  if (_parameters->get_tree_backup_step() > 0)
  {
    _lineage_tree->clean_cell_map();
    _lineage_tree->prune();
    _phylogenetic_tree->clean_cell_map();
    _phylogenetic_tree->prune();
    _phylogenetic_tree->shorten();
  }
  
  /*-------------------------------------------------------------------------*/
  /* 8) update population time                                               */
  /*-------------------------------------------------------------------------*/
  _population->update_time();
}

/**
 * \brief    For each trophic group, compute diversity measures
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::compute_trophic_groups_diversity( void )
{
  /*-------------------------------------------------------------------*/
  /* For each trophic group, get the cells list and compute statistics */
  /*-------------------------------------------------------------------*/
  TrophicGroup* group = _trophic_network->get_first_group();
  while (group != NULL)
  {
    unsigned long long int group_id = group->get_identifier();
    if (group_id != 0)
    {
      std::vector<unsigned long long int> cell_list;
      for (size_t i = 0; i < _width*_height; i++)
      {
        if (_population->get_cell(i)->isAlive() && _population->get_cell(i)->get_trophic_group() == group_id)
        {
          cell_list.push_back(_population->get_cell(i)->get_id());
        }
      }
      _phylogenetic_tree->compute_trophic_group_statistics(&cell_list, group);
    }
    group = _trophic_network->get_next_group();
  }
}

/**
 * \brief    Update the environment with random localization
 * \details  --
 * \param    bool initialize
 * \param    double factor
 * \return   \e void
 */
void Simulation::update_environment_random_localization( bool initialize, double factor )
{
  double inflowing_amount = 0.0;
  size_t x_draw            = (size_t)_parameters->get_prng()->uniform(0, (int)_width);
  size_t y_draw            = (size_t)_parameters->get_prng()->uniform(0, (int)_height);
  size_t number_of_species = draw_number_of_species(&(_parameters->get_environment_properties()->number_of_species_range));
  for (size_t nb = 0; nb < number_of_species; nb++)
  {
    /*---------------------------------*/
    /* Set the new environmental state */
    /*---------------------------------*/
    int species_tag      = draw_species_tag(&(_parameters->get_environment_properties()->species_tag_range));
    double concentration = draw_concentration(&(_parameters->get_environment_properties()->concentration_range));
    if (initialize)
    {
      _environment->set(x_draw, y_draw, species_tag, concentration*factor);
    }
    else
    {
      _environment->add(x_draw, y_draw, species_tag, concentration*factor);
    }
    
    /*---------------------------------*/
    /* Save matter flows               */
    /*---------------------------------*/
    inflowing_amount += concentration;
  }
  _environment->set_inflowing_amount(inflowing_amount);
}

/**
 * \brief    Update the environment with global localization
 * \details  --
 * \param    bool initialize
 * \param    double factor
 * \return   \e void
 */
void Simulation::update_environment_global_localization( bool initialize, double factor )
{
  double inflowing_amount = 0.0;
  for (size_t x = 0; x < _width; x++)
  {
    for (size_t y = 0; y < _height; y++)
    {
      size_t number_of_species = draw_number_of_species(&(_parameters->get_environment_properties()->number_of_species_range));
      for (size_t nb = 0; nb < number_of_species; nb++)
      {
        /*---------------------------------*/
        /* Set the new environmental state */
        /*---------------------------------*/
        int species_tag      = draw_species_tag(&(_parameters->get_environment_properties()->species_tag_range));
        double concentration = draw_concentration(&(_parameters->get_environment_properties()->concentration_range));
        if (initialize)
        {
          _environment->set(x, y, species_tag, concentration*factor);
        }
        else
        {
          _environment->add(x, y, species_tag, concentration*factor);
        }
        
        /*---------------------------------*/
        /* Save matter flows               */
        /*---------------------------------*/
        inflowing_amount += concentration;
      }
    }
  }
  _environment->set_inflowing_amount(inflowing_amount);
}

/**
 * \brief    Update the environment with centered localization
 * \details  --
 * \param    bool initialize
 * \param    double factor
 * \return   \e void
 */
void Simulation::update_environment_center_localization( bool initialize, double factor )
{
  double inflowing_amount = 0.0;
  size_t number_of_species = draw_number_of_species(&(_parameters->get_environment_properties()->number_of_species_range));
  for (size_t nb = 0; nb < number_of_species; nb++)
  {
    /*---------------------------------*/
    /* Set the new environmental state */
    /*---------------------------------*/
    int species_tag      = draw_species_tag(&(_parameters->get_environment_properties()->species_tag_range));
    double concentration = draw_concentration(&(_parameters->get_environment_properties()->concentration_range));
    if (initialize)
    {
      _environment->set(_width/2, _height/2, species_tag, concentration*factor);
    }
    else
    {
      _environment->add(_width/2, _height/2, species_tag, concentration*factor);
    }
    
    /*---------------------------------*/
    /* Save matter flows               */
    /*---------------------------------*/
    inflowing_amount += concentration;
  }
  _environment->set_inflowing_amount(inflowing_amount);
}

/**
 * \brief    Mix the population
 * \details  --
 * \param    void
 * \return   \e void
 */
void Simulation::mix( void )
{
  size_t n = _parameters->get_prng()->binomial(_width*_height, _parameters->get_migration_rate());
  while (n > 0)
  {
    size_t x1 = _parameters->get_prng()->uniform(0, (int)_width-1);
    size_t y1 = _parameters->get_prng()->uniform(0, (int)_height-1);
    size_t x2 = _parameters->get_prng()->uniform(0, (int)_width-1);
    size_t y2 = _parameters->get_prng()->uniform(0, (int)_height-1);
    Cell* tmp1 = _population->get_cell(x1, y1);
    tmp1->set_x(x2);
    tmp1->set_y(y2);
    Cell* tmp2 = _population->get_cell(x2, y2);
    tmp2->set_x(x1);
    tmp2->set_y(y1);
    _population->set_cell(x1*_height+y1, tmp2);
    _population->set_cell(x2*_height+y2, tmp1);
    n--;
  }
}

/**
 * \brief    Kill the cell at position i
 * \details  --
 * \param    size_t i
 * \return   \e void
 */
void Simulation::kill( size_t i )
{
  assert(i < _width*_height);
  Cell* cell = _population->get_cell(i);
  size_t y = i%_height;
  size_t x = (i-y)/_height;
  
  /*--------------------------------------------------*/
  /* 1) copy metabolite concentrations in environment */
  /*--------------------------------------------------*/
  if (_parameters->get_environment_properties()->interaction_scheme == INTERACTION)
  {
    SpeciesList *cell_list = cell->get_species_list();
    for (int tag = 1; tag <= (int)cell_list->get_size(); tag++)
    {
      _environment->add(x, y, tag, cell_list->get(tag));
    }
  }
  
  /*--------------------------------------------------*/
  /* 2) kill the cell                                 */
  /*--------------------------------------------------*/
  if (_parameters->get_tree_backup_step() > 0)
  {
    _lineage_tree->freeze_node(cell->get_id(), _population->get_time());
    _phylogenetic_tree->freeze_node(cell->get_id(), _population->get_time());
  }
  cell->kill(_population->get_time());
}

/**
 * \brief    Divide cell at position i in position child_position
 * \details  --
 * \param    size_t i
 * \param    size_t child_position
 * \return   \e void
 */
void Simulation::divide( size_t i, size_t child_position )
{
  assert(i < _width*_height);
  assert(child_position < _width*_height);
  Cell* parent = _population->get_cell(i);
  if (_parameters->get_tree_backup_step() > 0)
  {
    _lineage_tree->freeze_node(parent->get_id(), _population->get_time());
    _phylogenetic_tree->freeze_node(parent->get_id(), _population->get_time());
  }
  
  /*---------------------------------------------------*/
  /* 1) create child                                   */
  /*---------------------------------------------------*/
  _population->new_cell(i, child_position);
  Cell* child = _population->get_cell(child_position);
  
  /*---------------------------------------------------*/
  /* 2) mutate genomes                                 */
  /*---------------------------------------------------*/
  parent->mutate();
  child->mutate();
  
  /*---------------------------------------------------*/
  /* 3) load new genomes in species list               */
  /*---------------------------------------------------*/
  parent->load_genome_in_species_lists();
  child->load_genome_in_species_lists();
  
  /*---------------------------------------------------*/
  /* 4) synchronize state vectors                      */
  /*---------------------------------------------------*/
  parent->synchronize_state_vectors(_environment);
  child->synchronize_state_vectors(_environment);
  
  /*---------------------------------------------------*/
  /* 5) load genome in ODE system                      */
  /*---------------------------------------------------*/
  parent->load_genome_in_ODE_system(_environment, false, false);
  child->load_genome_in_ODE_system(_environment, false, false);
  
  /*---------------------------------------------------*/
  /* 6) update tree                                    */
  /*---------------------------------------------------*/
  if (_parameters->get_tree_backup_step() > 0)
  {
    _lineage_tree->add_division(parent, parent, child);
    _phylogenetic_tree->add_division(parent, parent, child);
  }
}

/**
 * \brief    Update the cell at position i
 * \details  --
 * \param    size_t i
 * \return   \e void
 */
void Simulation::update_cell( size_t i )
{
  /*----------------------------------------------------------------*/
  /* 1) check cell position                                         */
  /*----------------------------------------------------------------*/
  assert(i < _width*_height);
  
  /*----------------------------------------------------------------*/
  /* 2) get cell                                                    */
  /*----------------------------------------------------------------*/
  Cell* cell = _population->get_cell(i);
  
  /*----------------------------------------------------------------*/
  /* 3) synchronize species lists size between cell and environment */
  /*----------------------------------------------------------------*/
  cell->synchronize_state_vectors(_environment);
  
  /*----------------------------------------------------------------*/
  /* 4) solve ODEs, update species lists and compute score          */
  /*----------------------------------------------------------------*/
  cell->update(_population->get_time());
}

/**
 * \brief    Find empty cases around cell at position i
 * \details  --
 * \param    size_t i
 * \param    size_t& child_position
 * \return   \e bool
 */
bool Simulation::find_empty_cases( size_t i, size_t& child_position )
{
  assert(i < _width*_height);
  size_t x_coord[8];
  size_t y_coord[8];
  int    N = 0;
  size_t y = i%_height;
  size_t x = (i-y)/_height;
  for (int i = -1; i < 2; i++)
  {
    for (int j = -1; j < 2; j++)
    {
      if ( !(i == 0 && j == 0) )
      {
        size_t new_x = (x + i + _width) % _width;
        size_t new_y = (y + j + _height) % _height;
        if ( !_population->get_cell(new_x, new_y)->isAlive() )
        {
          x_coord[N] = new_x;
          y_coord[N] = new_y;
          N++;
        }
      }
    }
  }
  if (N == 0)
  {
    return false;
  }
  else
  {
    size_t draw = _parameters->get_prng()->uniform(0, N-1);
    child_position = x_coord[draw]*_height + y_coord[draw];
    return true;
  }
}

/**
 * \brief    Evaluate the size of species lists
 * \details  This method evaluates cells and environment species lists size and then reduce it if needed.
 * \param    void
 * \return   \e void
 */
void Simulation::evaluate_species_lists_size( void )
{
  size_t maximum_size = 0;
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 1) Evaluate the cells                  */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    Cell* cell = _population->get_cell(pos);
    if (cell->isAlive())
    {
      /*-------------------------------------------*/
      /* 1.1) Evaluate the reaction list           */
      /*-------------------------------------------*/
      reaction_list* rlist = cell->get_ode()->get_reaction_list();
      for (size_t i = 0; i < rlist->metabolic_N; i++)
      {
        if (maximum_size < (size_t)rlist->metabolic_s[i]-1)
        {
          maximum_size = (size_t)rlist->metabolic_s[i]-1;
        }
        if (maximum_size < (size_t)rlist->metabolic_p[i]-1)
        {
          maximum_size = (size_t)rlist->metabolic_p[i]-1;
        }
      }
      
      /*-------------------------------------------*/
      /* 1.2) Evaluate the inherited species list  */
      /*-------------------------------------------*/
      double* X = cell->get_inherited_species_list()->get_X();
      for (size_t i = 0; i < cell->get_inherited_species_list()->get_size(); i++)
      {
        if (X[i] > MINIMUM_CONCENTRATION && maximum_size < i)
        {
          maximum_size = i;
        }
      }
      
      /*-------------------------------------------*/
      /* 1.3) Evaluate the species list            */
      /*-------------------------------------------*/
      X = cell->get_species_list()->get_X();
      for (size_t i = 0; i < cell->get_species_list()->get_size(); i++)
      {
        if (X[i] > MINIMUM_CONCENTRATION && maximum_size < i)
        {
          maximum_size = i;
        }
      }
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 2) Evaluate the environment            */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  double* X = _environment->get_X()->get_X();
  for (size_t i = 0; i < _environment->get_species_lists_size(); i++)
  {
    if (X[i] > MINIMUM_CONCENTRATION && maximum_size < i)
    {
      maximum_size = i;
    }
  }
  
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  /* 3) If the maximum size is lower than   */
  /*    the environment species lists size, */
  /*    reduce every species lists to       */
  /*    maximum size                        */
  /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
  
  if (maximum_size < _environment->get_species_lists_size())
  {
    std::cout << "> Lower vectors size from " << _environment->get_species_lists_size() << " to " << maximum_size << "\n";
    for (size_t pos = 0; pos < _width*_height; pos++)
    {
      Cell* cell = _population->get_cell(pos);
      if (cell->isAlive())
      {
        cell->get_inherited_species_list()->decrease_size(maximum_size);
        cell->get_species_list()->decrease_size(maximum_size);
      }
    }
    _environment->decrease_species_lists_size(maximum_size);
  }
}

/**
 * \brief    Test the structure of the lineage tree
 * \details  --
 * \param    void
 * \return   \e void
 */
#ifdef DEBUG
void Simulation::test_lineage_tree_structure( void )
{
  /*--------------------------*/
  /* 1) Check tree structure  */
  /*--------------------------*/
  size_t master_root_count = 0;
  size_t root_count        = 0;
  size_t normal_count      = 0;
  size_t dead_count        = 0;
  size_t alive_count       = 0;
  Node* node = _lineage_tree->get_first_node();
  while (node != NULL)
  {
    /*----------------------------*/
    /* 1.1) Test master root node */
    /*----------------------------*/
    if (node->isMasterRoot())
    {
      master_root_count++;
      dead_count++;
      assert(node->get_parent() == NULL);
      assert(!node->isRoot());
      assert(!node->isNormal());
      assert(node->isDead());
      assert(!node->isAlive());
      assert(node->get_alive_cell() == NULL);
      assert(node->get_replication_report() == NULL);
    }
    /*----------------------------*/
    /* 1.2) Test root nodes       */
    /*----------------------------*/
    else if (node->isRoot())
    {
      root_count++;
      if (node->isDead())
      {
        assert(!node->isAlive());
        assert(node->get_alive_cell() == NULL);
        dead_count++;
      }
      else if (node->isAlive())
      {
        assert(!node->isDead());
        assert(_population->get_cell_by_id(node->get_alive_cell()->get_id()) != NULL);
        alive_count++;
      }
      assert(node->get_replication_report() != NULL);
      assert(!node->isMasterRoot());
      assert(!node->isNormal());
      assert(node->get_parent()->isMasterRoot());
      for (size_t i = 0; i < node->get_number_of_children(); i++)
      {
        assert(!node->get_child(i)->isMasterRoot());
        assert(!node->get_child(i)->isRoot());
        assert(node->get_child(i)->isNormal());
        assert(node->get_child(i)->get_parent()->get_id() == node->get_id());
      }
    }
    /*----------------------------*/
    /* 1.3) Test normal nodes     */
    /*----------------------------*/
    else if (node->isNormal())
    {
      normal_count++;
      if (node->isDead())
      {
        assert(!node->isAlive());
        assert(node->get_alive_cell() == NULL);
        dead_count++;
      }
      else if (node->isAlive())
      {
        assert(!node->isDead());
        assert(_population->get_cell_by_id(node->get_alive_cell()->get_id()) != NULL);
        alive_count++;
      }
      assert(node->get_replication_report() != NULL);
      assert(!node->isMasterRoot());
      assert(!node->isRoot());
      assert(node->get_parent()->isRoot() || node->get_parent()->isNormal());
      for (size_t i = 0; i < node->get_number_of_children(); i++)
      {
        assert(!node->get_child(i)->isMasterRoot());
        assert(!node->get_child(i)->isRoot());
        assert(node->get_child(i)->isNormal());
        assert(node->get_child(i)->get_parent()->get_id() == node->get_id());
      }
    }
    
    node = _lineage_tree->get_next_node();
  }
  
  /*--------------------------*/
  /* 2) Test nodes counts     */
  /*--------------------------*/
  assert(master_root_count == 1);
  assert(master_root_count+root_count+normal_count == _lineage_tree->get_number_of_nodes());
  assert(dead_count+alive_count == _lineage_tree->get_number_of_nodes());
  std::vector<unsigned long long int> alive_nodes;
  _lineage_tree->get_alive_nodes(&alive_nodes);
  assert(alive_nodes.size() == alive_count);
  for (size_t i = 0; i < alive_nodes.size(); i++)
  {
    Node* node = _lineage_tree->get_node(alive_nodes[i]);
    assert(!node->isMasterRoot());
    assert(node->isRoot() || node->isNormal());
    assert(node->isAlive());
    assert(node->get_alive_cell()->isAlive());
    Cell* cell = _population->get_cell_by_id(node->get_alive_cell()->get_id());
    assert(cell->isAlive());
    assert(node->get_replication_report()->get_number_of_events() == cell->get_replication_report()->get_number_of_events());
    (void)cell;
  }
  
  /*--------------------------*/
  /* 3) Check the master root */
  /*--------------------------*/
  Node* master_root = _lineage_tree->get_node(0);
  assert(master_root->isMasterRoot());
  assert(master_root->get_number_of_children() == root_count);
  (void)master_root;
}
#endif

/**
 * \brief    Test the structure of the phylogenetic tree
 * \details  --
 * \param    void
 * \return   \e void
 */
#ifdef DEBUG
void Simulation::test_phylogenetic_tree_structure( void )
{
  /*--------------------------*/
  /* 1) Check tree structure  */
  /*--------------------------*/
  size_t master_root_count = 0;
  size_t root_count        = 0;
  size_t normal_count      = 0;
  size_t dead_count        = 0;
  size_t alive_count       = 0;
  Node* node = _phylogenetic_tree->get_first_node();
  while (node != NULL)
  {
    /*----------------------------*/
    /* 1.1) Test master root node */
    /*----------------------------*/
    if (node->isMasterRoot())
    {
      master_root_count++;
      dead_count++;
      assert(node->get_parent() == NULL);
      assert(!node->isRoot());
      assert(!node->isNormal());
      assert(node->isDead());
      assert(!node->isAlive());
      assert(node->get_alive_cell() == NULL);
      assert(node->get_replication_report() == NULL);
    }
    /*----------------------------*/
    /* 1.2) Test root nodes       */
    /*----------------------------*/
    else if (node->isRoot())
    {
      root_count++;
      if (node->isDead())
      {
        assert(!node->isAlive());
        assert(node->get_alive_cell() == NULL);
        dead_count++;
      }
      else if (node->isAlive())
      {
        assert(!node->isDead());
        assert(_population->get_cell_by_id(node->get_alive_cell()->get_id()) != NULL);
        alive_count++;
      }
      assert(node->get_replication_report() != NULL);
      assert(!node->isMasterRoot());
      assert(!node->isNormal());
      assert(node->get_parent()->isMasterRoot());
      bool alive_child = false;
      for (size_t i = 0; i < node->get_number_of_children(); i++)
      {
        assert(!node->get_child(i)->isMasterRoot());
        assert(!node->get_child(i)->isRoot());
        assert(node->get_child(i)->isNormal());
        assert(node->get_child(i)->get_parent()->get_id() == node->get_id());
        if (node->get_child(i)->isAlive())
        {
          assert(_population->get_cell_by_id(node->get_child(i)->get_alive_cell()->get_id()) != NULL);
          alive_child = true;
        }
      }
      if (node->isDead() && !alive_child)
      {
        assert(node->get_number_of_children() == 2);
      }
    }
    /*----------------------------*/
    /* 1.3) Test normal nodes     */
    /*----------------------------*/
    else if (node->isNormal())
    {
      normal_count++;
      if (node->isDead())
      {
        assert(node->get_number_of_children() == 2);
        assert(!node->isAlive());
        assert(node->get_alive_cell() == NULL);
        dead_count++;
      }
      else if (node->isAlive())
      {
        assert(!node->isDead());
        assert(_population->get_cell_by_id(node->get_alive_cell()->get_id()) != NULL);
        alive_count++;
      }
      assert(node->get_replication_report() != NULL);
      assert(!node->isMasterRoot());
      assert(!node->isRoot());
      assert(node->get_parent()->isRoot() || node->get_parent()->isNormal());
      bool alive_child = false;
      for (size_t i = 0; i < node->get_number_of_children(); i++)
      {
        assert(!node->get_child(i)->isMasterRoot());
        assert(!node->get_child(i)->isRoot());
        assert(node->get_child(i)->isNormal());
        assert(node->get_child(i)->get_parent()->get_id() == node->get_id());
        if (node->get_child(i)->isAlive())
        {
          assert(_population->get_cell_by_id(node->get_child(i)->get_alive_cell()->get_id()) != NULL);
          alive_child = true;
        }
      }
      if (node->isDead() && !alive_child)
      {
        assert(node->get_number_of_children() == 2);
      }
    }
    
    node = _phylogenetic_tree->get_next_node();
  }
  
  /*--------------------------*/
  /* 2) Test nodes counts     */
  /*--------------------------*/
  assert(master_root_count == 1);
  assert(master_root_count+root_count+normal_count == _phylogenetic_tree->get_number_of_nodes());
  assert(dead_count+alive_count == _phylogenetic_tree->get_number_of_nodes());
  std::vector<unsigned long long int> alive_nodes;
  _phylogenetic_tree->get_alive_nodes(&alive_nodes);
  assert(alive_nodes.size() == alive_count);
  for (size_t i = 0; i < alive_nodes.size(); i++)
  {
    Node* node = _phylogenetic_tree->get_node(alive_nodes[i]);
    assert(!node->isMasterRoot());
    assert(node->isRoot() || node->isNormal());
    assert(node->isAlive());
    assert(node->get_alive_cell()->isAlive());
    Cell* cell = _population->get_cell_by_id(node->get_alive_cell()->get_id());
    assert(cell->isAlive());
    assert(node->get_replication_report()->get_number_of_events() == cell->get_replication_report()->get_number_of_events());
    (void)cell;
  }
  
  /*--------------------------*/
  /* 3) Check the master root */
  /*--------------------------*/
  Node* master_root = _phylogenetic_tree->get_node(0);
  assert(master_root->isMasterRoot());
  assert(master_root->get_number_of_children() == root_count);
  (void)master_root;
}
#endif
