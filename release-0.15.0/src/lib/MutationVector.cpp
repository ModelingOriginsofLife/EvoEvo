
/**
 * \file      MutationVector.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     MutationVector class definition
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

#include "MutationVector.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Default constructor
 * \details  --
 * \param    void
 * \return   \e void
 */
MutationVector::MutationVector( void )
{
  _dX = new pearl;
  
  /*------------------------------------------------------------------ Global attributes */
  
  _dX->type              = NON_CODING;
  _dX->identifier        = 0;
  _dX->parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  _dX->s    = 0;
  _dX->p    = 0;
  _dX->km   = 0.0;
  _dX->kcat = 0.0;
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  _dX->BS_tag         = 0;
  _dX->coE_tag        = 0;
  _dX->free_activity  = false;
  _dX->bound_activity = false;
  _dX->window         = 0;
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  _dX->TF_tag = 0;
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  _dX->basal_expression_level = 0.0;
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  _dX->functional = false;
}

/**
 * \brief    Constructor from backup file
 * \details  Load MutationVector class from backup file
 * \param    gzFile backup_file
 * \return   \e void
 */
MutationVector::MutationVector( gzFile backup_file )
{
  _dX = new pearl;
  load_pearl(backup_file, _dX);
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const MutationVector& vector
 * \return   \e void
 */
MutationVector::MutationVector( const MutationVector& vector )
{
  _dX = new pearl;
  
  /*------------------------------------------------------------------ Global attributes */
  
  _dX->type              = vector._dX->type;
  _dX->identifier        = vector._dX->identifier;
  _dX->parent_identifier = vector._dX->parent_identifier;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  _dX->s    = vector._dX->s;
  _dX->p    = vector._dX->p;
  _dX->km   = vector._dX->km;
  _dX->kcat = vector._dX->kcat;
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  _dX->BS_tag         = vector._dX->BS_tag;
  _dX->coE_tag        = vector._dX->coE_tag;
  _dX->free_activity  = vector._dX->free_activity;
  _dX->bound_activity = vector._dX->bound_activity;
  _dX->window         = vector._dX->window;
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  _dX->TF_tag = vector._dX->TF_tag;
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  _dX->basal_expression_level = vector._dX->basal_expression_level;
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  _dX->functional = vector._dX->functional;
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
MutationVector::~MutationVector( void )
{
  delete _dX;
  _dX = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/*
 * \brief    Draw vector
 * \details  --
 * \param    Prng* prng
 * \param    const double* mutation_rates
 * \return   \e bool
 */
bool MutationVector::draw( Prng* prng, const double* mutation_rates )
{
  bool mutate = false;
  
  /*-----------------------------------------------------*/
  /* 1) mutate pearl type                                */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[TRANSITION_RATE])
  {
    /* Draw uniformly between the five types of pearls */
    double probas[5]    = {1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0};
    pearl_type new_type = (pearl_type)prng->roulette_wheel(probas, 1.0, 5);
    if (new_type != _dX->type)
    {
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* A) the pearl was NON_CODING type           */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      if (_dX->type == NON_CODING && new_type == ENZYME)
      {
        _dX->type = NC_TO_E_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = NC_TO_TF_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == BINDING_SITE)
      {
        _dX->type = NC_TO_BS_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == PROMOTER)
      {
        _dX->type = NC_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* B) the pearl was ENZYME type               */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == ENZYME && new_type == NON_CODING)
      {
        _dX->type = E_TO_NC_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = E_TO_TF_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == BINDING_SITE)
      {
        _dX->type = E_TO_BS_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == PROMOTER)
      {
        _dX->type = E_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* C) the pearl was TRANSCRIPTION_FACTOR type */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == NON_CODING)
      {
        _dX->type = TF_TO_NC_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == ENZYME)
      {
        _dX->type = TF_TO_E_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == BINDING_SITE)
      {
        _dX->type = TF_TO_BS_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == PROMOTER)
      {
        _dX->type = TF_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* D) the pearl was BINDING_SITE type         */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == BINDING_SITE && new_type == NON_CODING)
      {
        _dX->type = BS_TO_NC_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == ENZYME)
      {
        _dX->type = BS_TO_E_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = BS_TO_TF_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == PROMOTER)
      {
        _dX->type = BS_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* E) the pearl was PROMOTER type             */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == PROMOTER && new_type == NON_CODING)
      {
        _dX->type = P_TO_NC_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == ENZYME)
      {
        _dX->type = P_TO_E_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = P_TO_TF_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == BINDING_SITE)
      {
        _dX->type = P_TO_BS_TRANSITION;
      }
      mutate = true;
    }
  }
  
  /*-----------------------------------------------------*/
  /* 2) mutate enzyme type (E) attributes                */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->s = prng->uniform(-mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE], mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->p = prng->uniform(-mutation_rates[PRODUCT_TAG_MUTATION_SIZE], mutation_rates[PRODUCT_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->km = prng->gaussian(0.0, mutation_rates[KM_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->kcat = prng->gaussian(0.0, mutation_rates[KCAT_MUTATION_SIZE]);
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 3) mutate transcription factor type (TF) attributes */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->BS_tag = prng->uniform(-mutation_rates[BINDING_SIZE_TAG_MUTATION_SIZE], mutation_rates[BINDING_SIZE_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->coE_tag = prng->uniform(-mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE], mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->free_activity = !_dX->free_activity;
    mutate = true;
  }
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->bound_activity = !_dX->bound_activity;
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 4) mutate binding site type (BS) attributes         */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->TF_tag = prng->uniform(-mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE], mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 5) mutate promoter type (P) attributes              */
  /*-----------------------------------------------------*/
  if (prng->uniform() < mutation_rates[POINT_MUTATION_RATE])
  {
    _dX->basal_expression_level = prng->gaussian(0.0, mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE]);
    mutate = true;
  }
  return mutate;
}

/*
 * \brief    Draw vector at breakpoints
 * \details  --
 * \param    Prng* prng
 * \param    const double* mutation_rates
 * \return   \e bool
 */
bool MutationVector::breakpoint_draw( Prng* prng, const double* mutation_rates )
{
  bool mutate = false;
  
  /*-----------------------------------------------------*/
  /* 1) mutate pearl type                                */
  /*-----------------------------------------------------*/
  if (prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    /* Draw uniformly between the five types of pearls */
    double probas[5]    = {1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0};
    pearl_type new_type = (pearl_type)prng->roulette_wheel(probas, 1.0, 5);
    if (new_type != _dX->type)
    {
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* A) the pearl was NON_CODING type           */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      if (_dX->type == NON_CODING && new_type == ENZYME)
      {
        _dX->type = NC_TO_E_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = NC_TO_TF_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == BINDING_SITE)
      {
        _dX->type = NC_TO_BS_TRANSITION;
      }
      else if (_dX->type == NON_CODING && new_type == PROMOTER)
      {
        _dX->type = NC_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* B) the pearl was ENZYME type               */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == ENZYME && new_type == NON_CODING)
      {
        _dX->type = E_TO_NC_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = E_TO_TF_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == BINDING_SITE)
      {
        _dX->type = E_TO_BS_TRANSITION;
      }
      else if (_dX->type == ENZYME && new_type == PROMOTER)
      {
        _dX->type = E_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* C) the pearl was TRANSCRIPTION_FACTOR type */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == NON_CODING)
      {
        _dX->type = TF_TO_NC_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == ENZYME)
      {
        _dX->type = TF_TO_E_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == BINDING_SITE)
      {
        _dX->type = TF_TO_BS_TRANSITION;
      }
      else if (_dX->type == TRANSCRIPTION_FACTOR && new_type == PROMOTER)
      {
        _dX->type = TF_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* D) the pearl was BINDING_SITE type         */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == BINDING_SITE && new_type == NON_CODING)
      {
        _dX->type = BS_TO_NC_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == ENZYME)
      {
        _dX->type = BS_TO_E_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = BS_TO_TF_TRANSITION;
      }
      else if (_dX->type == BINDING_SITE && new_type == PROMOTER)
      {
        _dX->type = BS_TO_P_TRANSITION;
      }
      
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      /* E) the pearl was PROMOTER type             */
      /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
      else if (_dX->type == PROMOTER && new_type == NON_CODING)
      {
        _dX->type = P_TO_NC_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == ENZYME)
      {
        _dX->type = P_TO_E_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == TRANSCRIPTION_FACTOR)
      {
        _dX->type = P_TO_TF_TRANSITION;
      }
      else if (_dX->type == PROMOTER && new_type == BINDING_SITE)
      {
        _dX->type = P_TO_BS_TRANSITION;
      }
      mutate = true;
    }
  }
  
  /*-----------------------------------------------------*/
  /* 2) mutate enzyme type (E) attributes                */
  /*-----------------------------------------------------*/
  if (prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    _dX->s = prng->uniform(-mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE], mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    _dX->p = prng->uniform(-mutation_rates[PRODUCT_TAG_MUTATION_SIZE], mutation_rates[PRODUCT_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    _dX->km = prng->gaussian(0.0, mutation_rates[KM_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    _dX->kcat = prng->gaussian(0.0, mutation_rates[KCAT_MUTATION_SIZE]);
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 3) mutate transcription factor type (TF) attributes */
  /*-----------------------------------------------------*/
  if (prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    _dX->BS_tag = prng->uniform(-mutation_rates[BINDING_SIZE_TAG_MUTATION_SIZE], mutation_rates[BINDING_SIZE_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    _dX->coE_tag = prng->uniform(-mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE], mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  if (prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    _dX->free_activity = !_dX->free_activity;
    mutate = true;
  }
  if (prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    _dX->bound_activity = !_dX->bound_activity;
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 4) mutate binding site type (BS) attributes         */
  /*-----------------------------------------------------*/
  if (prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    _dX->TF_tag = prng->uniform(-mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE], mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE]);
    mutate = true;
  }
  
  /*-----------------------------------------------------*/
  /* 5) mutate promoter type (P) attributes              */
  /*-----------------------------------------------------*/
  if (prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    _dX->basal_expression_level = prng->gaussian(0.0, mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE]);
    mutate = true;
  }
  return mutate;
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void MutationVector::save( gzFile backup_file )
{
  save_pearl(backup_file, _dX);
}

/**
 * \brief    Clear random vector
 * \details  --
 * \param    void
 * \return   \e void
 */
void MutationVector::clear( void )
{
  /*------------------------------------------------------------------ Global attributes */
  
  _dX->type              = NON_CODING;
  _dX->identifier        = 0;
  _dX->parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  _dX->s    = 0;
  _dX->p    = 0;
  _dX->km   = 0.0;
  _dX->kcat = 0.0;
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  _dX->BS_tag         = 0;
  _dX->coE_tag        = 0;
  _dX->free_activity  = false;
  _dX->bound_activity = false;
  _dX->window         = 0;
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  _dX->TF_tag = 0;
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  _dX->basal_expression_level = 0.0;
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  _dX->functional = false;
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Load a pearl in backup file
 * \details  --
 * \param    gzFile backup_file
 * \param    pearl* p
 * \return   \e void
 */
void MutationVector::load_pearl( gzFile backup_file, pearl* p )
{
  /*------------------------------------------------------------------ global attributes */
  
  gzread( backup_file, &p->type,              sizeof(p->type) );
  gzread( backup_file, &p->identifier,        sizeof(p->identifier) );
  gzread( backup_file, &p->parent_identifier, sizeof(p->parent_identifier) );
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  gzread( backup_file, &p->s,    sizeof(p->s) );
  gzread( backup_file, &p->p,    sizeof(p->p) );
  gzread( backup_file, &p->km,   sizeof(p->km) );
  gzread( backup_file, &p->kcat, sizeof(p->kcat) );
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  gzread( backup_file, &p->BS_tag,         sizeof(p->BS_tag) );
  gzread( backup_file, &p->coE_tag,        sizeof(p->coE_tag) );
  gzread( backup_file, &p->free_activity,  sizeof(p->free_activity) );
  gzread( backup_file, &p->bound_activity, sizeof(p->bound_activity) );
  gzread( backup_file, &p->window,         sizeof(p->window) );
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzread( backup_file, &p->TF_tag, sizeof(p->TF_tag) );
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  gzread( backup_file, &p->basal_expression_level, sizeof(p->basal_expression_level) );
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  gzread( backup_file, &p->functional, sizeof(p->functional) );
}

/**
 * \brief    Save a pearl in backup file
 * \details  --
 * \param    gzFile backup_file
 * \param    pearl* p
 * \return   \e void
 */
void MutationVector::save_pearl( gzFile backup_file, pearl* p )
{
  /*------------------------------------------------------------------ global attributes */
  
  gzwrite( backup_file, &p->type,              sizeof(p->type) );
  gzwrite( backup_file, &p->identifier,        sizeof(p->identifier) );
  gzwrite( backup_file, &p->parent_identifier, sizeof(p->parent_identifier) );
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  gzwrite( backup_file, &p->s,    sizeof(p->s) );
  gzwrite( backup_file, &p->p,    sizeof(p->p) );
  gzwrite( backup_file, &p->km,   sizeof(p->km) );
  gzwrite( backup_file, &p->kcat, sizeof(p->kcat) );
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  gzwrite( backup_file, &p->BS_tag,         sizeof(p->BS_tag) );
  gzwrite( backup_file, &p->coE_tag,        sizeof(p->coE_tag) );
  gzwrite( backup_file, &p->free_activity,  sizeof(p->free_activity) );
  gzwrite( backup_file, &p->bound_activity, sizeof(p->bound_activity) );
  gzwrite( backup_file, &p->window,         sizeof(p->window) );
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzwrite( backup_file, &p->TF_tag, sizeof(p->TF_tag) );
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  gzwrite( backup_file, &p->basal_expression_level, sizeof(p->basal_expression_level) );
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  gzwrite( backup_file, &p->functional, sizeof(p->functional) );
}
