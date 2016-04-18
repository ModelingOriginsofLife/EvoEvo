
/**
 * \file      InheritedProteins.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     InheritedProteins class definition
 */

/****************************************************************************
 * Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * E-mail: charles.rocabert@inria.fr
 * Web: http://evoevo.eu/
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

#include "InheritedProteins.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \return   \e void
 */
InheritedProteins::InheritedProteins( Parameters* parameters )
{
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = parameters;
  
  /*------------------------------------------------------------------ string of pearls */
  
  create_string_of_pearls();
  
  /*------------------------------------------------------------------ concentration vector */
  
  _concentration_vector = NULL;
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_E             = 0;
  _nb_TF            = 0;
  _nb_inner_enzymes = 0;
  _nb_inflow_pumps  = 0;
  _nb_outflow_pumps = 0;
  
  /*------------------------------------------------------------------ pearl positions */
  
  _Ei  = NULL;
  _TFi = NULL;
}

/**
 * \brief    Constructor from backup file
 * \details  Load InheritedProteins class from backup file
 * \param    Parameters* parameters
 * \param    gzFile backup_file
 * \return   \e void
 */
InheritedProteins::InheritedProteins( Parameters* parameters, gzFile backup_file )
{
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = parameters;
  
  /*------------------------------------------------------------------ string of pearls */
  
  load_string_of_pearls(backup_file);
  
  /*------------------------------------------------------------------ concentration vector */
  
  _concentration_vector = NULL;
  if (_string_of_pearls->size > 0)
  {
    _concentration_vector = new double[_string_of_pearls->size];
    for (size_t i = 0; i < _string_of_pearls->size; i++)
    {
      gzread( backup_file, &_concentration_vector[i], sizeof(_concentration_vector[i]) );
    }
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  gzread( backup_file, &_nb_E,             sizeof(_nb_E) );
  gzread( backup_file, &_nb_TF,            sizeof(_nb_TF) );
  gzread( backup_file, &_nb_inner_enzymes, sizeof(_nb_inner_enzymes) );
  gzread( backup_file, &_nb_inflow_pumps,  sizeof(_nb_inflow_pumps) );
  gzread( backup_file, &_nb_outflow_pumps, sizeof(_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ pearl positions */
  
  _Ei = NULL;
  if (_nb_E > 0)
  {
    _Ei = new size_t[_nb_E];
    for (size_t i = 0; i < _nb_E; i++)
    {
      gzread( backup_file, &_Ei[i], sizeof(_Ei[i]) );
    }
  }
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    for (size_t i = 0; i < _nb_TF; i++)
    {
      gzread( backup_file, &_TFi[i], sizeof(_TFi[i]) );
    }
  }
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const InheritedProteins& inherited_proteins
 * \return   \e void
 */
InheritedProteins::InheritedProteins( const InheritedProteins& inherited_proteins )
{
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = inherited_proteins._parameters;
  
  /*------------------------------------------------------------------ string of pearls */
  
  copy_string_of_pearls(inherited_proteins._string_of_pearls);
  
  /*------------------------------------------------------------------ concentration vector */
  
  _concentration_vector = NULL;
  if (_string_of_pearls->size > 0)
  {
    _concentration_vector = new double[_string_of_pearls->size];
    memcpy(_concentration_vector, inherited_proteins._concentration_vector, sizeof(double)*_string_of_pearls->size);
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_E             = inherited_proteins._nb_E;
  _nb_TF            = inherited_proteins._nb_TF;
  _nb_inner_enzymes = inherited_proteins._nb_inner_enzymes;
  _nb_inflow_pumps  = inherited_proteins._nb_inflow_pumps;
  _nb_outflow_pumps = inherited_proteins._nb_outflow_pumps;
  
  /*------------------------------------------------------------------ pearl positions */
  
  _Ei = NULL;
  if (_nb_E > 0)
  {
    _Ei = new size_t[_nb_E];
    memcpy(_Ei,  inherited_proteins._Ei,  sizeof(size_t)*_nb_E);
  }
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    memcpy(_TFi, inherited_proteins._TFi, sizeof(size_t)*_nb_TF);
  }
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
InheritedProteins::~InheritedProteins( void )
{
  delete_string_of_pearls();
  delete[] _concentration_vector;
  _concentration_vector = NULL;
  delete[] _Ei;
  _Ei  = NULL;
  delete[] _TFi;
  _TFi = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Initialize the concentration vector at zero (after mutation)
 * \details  --
 * \param    void
 * \return   \e void
 */
void InheritedProteins::initialize_concentration_vector( void )
{
  delete[] _concentration_vector;
  _concentration_vector = new double[_string_of_pearls->size];
  for (size_t i = 0; i < _string_of_pearls->size; i++)
  {
    _concentration_vector[i] = 0.0;
  }
}

/**
 * \brief    Build the indexes list
 * \details  --
 * \param    void
 * \return   \e void
 */
void InheritedProteins::build_index_list( void )
{
  /*--------------------------------------------------*/
  /* 1) count pearl types                             */
  /*--------------------------------------------------*/
  _nb_E             = 0;
  _nb_TF            = 0;
  _nb_inner_enzymes = 0;
  _nb_inflow_pumps  = 0;
  _nb_outflow_pumps = 0;
  pearl* str        = _string_of_pearls->x;
  for (size_t pos = 0; pos < _string_of_pearls->size; pos++)
  {
    if ( str[pos].type == ENZYME )
    {
      _nb_E++;
      if (str[pos].kcat > 0.0 && str[pos].s == str[pos].p)
      {
        _nb_inflow_pumps++;
      }
      else if (str[pos].kcat < 0.0 && str[pos].s == str[pos].p)
      {
        _nb_outflow_pumps++;
      }
      else if (str[pos].kcat != 0.0 && str[pos].s != str[pos].p)
      {
        _nb_inner_enzymes++;
      }
    }
    else if (str[pos].type == TRANSCRIPTION_FACTOR)
    {
      _nb_TF++;
    }
  }
  
  /*--------------------------------------------------*/
  /* 2) compute the indexes list for each pearl types */
  /*--------------------------------------------------*/
  delete[] _Ei;
  _Ei = NULL;
  if (_nb_E > 0)
  {
    _Ei = new size_t[_nb_E];
  }
  delete[] _TFi;
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
  }
  size_t Ecount  = 0;
  size_t TFcount = 0;
  for (size_t pos = 0; pos < _string_of_pearls->size; pos++)
  {
    if (str[pos].type == ENZYME)
    {
      _Ei[Ecount] = pos;
      Ecount++;
    }
    else if (str[pos].type == TRANSCRIPTION_FACTOR)
    {
      _TFi[TFcount] = pos;
      TFcount++;
    }
  }
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void InheritedProteins::save( gzFile backup_file )
{
  /*------------------------------------------------------------------ string of pearls */
  
  save_string_of_pearls(backup_file);
  
  /*------------------------------------------------------------------ concentration vector */
  
  if (_string_of_pearls->size > 0)
  {
    for (size_t i = 0; i < _string_of_pearls->size; i++)
    {
      gzwrite( backup_file, &_concentration_vector[i], sizeof(_concentration_vector[i]) );
    }
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  gzwrite( backup_file, &_nb_E,             sizeof(_nb_E) );
  gzwrite( backup_file, &_nb_TF,            sizeof(_nb_TF) );
  gzwrite( backup_file, &_nb_inner_enzymes, sizeof(_nb_inner_enzymes) );
  gzwrite( backup_file, &_nb_inflow_pumps,  sizeof(_nb_inflow_pumps) );
  gzwrite( backup_file, &_nb_outflow_pumps, sizeof(_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ pearl positions */
  
  if (_nb_E > 0)
  {
    for (size_t i = 0; i < _nb_E; i++)
    {
      gzwrite( backup_file, &_Ei[i], sizeof(_Ei[i]) );
    }
  }
  if (_nb_TF > 0)
  {
    for (size_t i = 0; i < _nb_TF; i++)
    {
      gzwrite( backup_file, &_TFi[i], sizeof(_TFi[i]) );
    }
  }
}

/**
 * \brief    Replace inherited proteins data
 * \details  --
 * \param    InheritedProteins* inherited_proteins
 * \return   \e void
 */
void InheritedProteins::replace_data( InheritedProteins* inherited_proteins )
{
  /*------------------------------------------------------------------ string of pearls */
  
  delete_string_of_pearls();
  copy_string_of_pearls(inherited_proteins->get_string_of_pearl());
  
  /*------------------------------------------------------------------ concentration vector */
  
  delete[] _concentration_vector;
  _concentration_vector = NULL;
  if (_string_of_pearls->size > 0)
  {
    _concentration_vector = new double[_string_of_pearls->size];
    memcpy(_concentration_vector, inherited_proteins->get_concentration_vector(), sizeof(double)*_string_of_pearls->size);
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_E             = inherited_proteins->get_nb_E();
  _nb_TF            = inherited_proteins->get_nb_TF();
  _nb_inner_enzymes = inherited_proteins->get_nb_inner_enzymes();
  _nb_inflow_pumps  = inherited_proteins->get_nb_inflow_pumps();
  _nb_outflow_pumps = inherited_proteins->get_nb_outflow_pumps();
  
  /*------------------------------------------------------------------ pearl positions */
  
  delete[] _Ei;
  _Ei = NULL;
  if (_nb_E > 0)
  {
    _Ei = new size_t[_nb_E];
    memcpy(_Ei,  inherited_proteins->get_Ei(),  sizeof(size_t)*_nb_E);
  }
  delete[] _TFi;
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    memcpy(_TFi, inherited_proteins->get_TFi(), sizeof(size_t)*_nb_TF);
  }
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Create a default string of pearls
 * \details  --
 * \param    void
 * \return   \e void
 */
void InheritedProteins::create_string_of_pearls( void )
{
  _string_of_pearls              = new string_of_pearls;
  _string_of_pearls->x           = new pearl[INHERITED_PROTEINS_BUFFER];
  _string_of_pearls->size        = 0;
  _string_of_pearls->buffer_size = INHERITED_PROTEINS_BUFFER;
}

/**
 * \brief    Copy a string of pearls
 * \details  --
 * \param    const string_of_pearls* model
 * \return   \e void
 */
void InheritedProteins::copy_string_of_pearls( const string_of_pearls* model )
{
  _string_of_pearls              = new string_of_pearls;
  _string_of_pearls->x           = new pearl[model->buffer_size];
  memcpy(_string_of_pearls->x, model->x, sizeof(pearl)*model->buffer_size);
  _string_of_pearls->size        = model->size;
  _string_of_pearls->buffer_size = model->buffer_size;
}

/**
 * \brief    Delete the string of pearls
 * \details  --
 * \param    void
 * \return   \e void
 */
void InheritedProteins::delete_string_of_pearls( void )
{
  delete[] _string_of_pearls->x;
  _string_of_pearls->x = NULL;
  delete _string_of_pearls;
  _string_of_pearls = NULL;
}

/**
 * \brief    Load the string of pearls from backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void InheritedProteins::load_string_of_pearls( gzFile backup_file )
{
  _string_of_pearls = new string_of_pearls;
  gzread( backup_file, &_string_of_pearls->size,        sizeof(_string_of_pearls->size) );
  gzread( backup_file, &_string_of_pearls->buffer_size, sizeof(_string_of_pearls->buffer_size) );
  _string_of_pearls->x = new pearl[_string_of_pearls->buffer_size];
  for (size_t i = 0; i < _string_of_pearls->size; i++)
  {
    load_pearl(backup_file, _string_of_pearls->x[i]);
  }
}

/**
 * \brief    Save the string of pearls in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void InheritedProteins::save_string_of_pearls( gzFile backup_file )
{
  gzwrite( backup_file, &_string_of_pearls->size,        sizeof(_string_of_pearls->size) );
  gzwrite( backup_file, &_string_of_pearls->buffer_size, sizeof(_string_of_pearls->buffer_size) );
  for (size_t i = 0; i < _string_of_pearls->size; i++)
  {
    save_pearl(backup_file, _string_of_pearls->x[i]);
  }
}

/**
 * \brief    Load a pearl in backup file
 * \details  --
 * \param    gzFile backup_file
 * \param    pearl& p
 * \return   \e void
 */
void InheritedProteins::load_pearl( gzFile backup_file, pearl& p )
{
  /*------------------------------------------------------------------ global attributes */
  
  gzread( backup_file, &p.type,              sizeof(p.type) );
  gzread( backup_file, &p.identifier,        sizeof(p.identifier) );
  gzread( backup_file, &p.parent_identifier, sizeof(p.parent_identifier) );
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  gzread( backup_file, &p.s,    sizeof(p.s) );
  gzread( backup_file, &p.p,    sizeof(p.p) );
  gzread( backup_file, &p.km,   sizeof(p.km) );
  gzread( backup_file, &p.kcat, sizeof(p.kcat) );
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  gzread( backup_file, &p.BS_tag,         sizeof(p.BS_tag) );
  gzread( backup_file, &p.coE_tag,        sizeof(p.coE_tag) );
  gzread( backup_file, &p.free_activity,  sizeof(p.free_activity) );
  gzread( backup_file, &p.bound_activity, sizeof(p.bound_activity) );
  gzread( backup_file, &p.window,         sizeof(p.window));
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzread( backup_file, &p.TF_tag, sizeof(p.TF_tag));
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  gzread( backup_file, &p.basal_expression_level, sizeof(p.basal_expression_level) );
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  gzread( backup_file, &p.functional, sizeof(p.functional) );
}

/**
 * \brief    Save a pearl in backup file
 * \details  --
 * \param    gzFile backup_file
 * \param    pearl& p
 * \return   \e void
 */
void InheritedProteins::save_pearl( gzFile backup_file, pearl& p )
{
  /*------------------------------------------------------------------ global attributes */
  
  gzwrite( backup_file, &p.type,              sizeof(p.type) );
  gzwrite( backup_file, &p.identifier,        sizeof(p.identifier) );
  gzwrite( backup_file, &p.parent_identifier, sizeof(p.parent_identifier) );
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  gzwrite( backup_file, &p.s,    sizeof(p.s) );
  gzwrite( backup_file, &p.p,    sizeof(p.p) );
  gzwrite( backup_file, &p.km,   sizeof(p.km) );
  gzwrite( backup_file, &p.kcat, sizeof(p.kcat) );
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  gzwrite( backup_file, &p.BS_tag,         sizeof(p.BS_tag) );
  gzwrite( backup_file, &p.coE_tag,        sizeof(p.coE_tag) );
  gzwrite( backup_file, &p.free_activity,  sizeof(p.free_activity) );
  gzwrite( backup_file, &p.bound_activity, sizeof(p.bound_activity) );
  gzwrite( backup_file, &p.window,         sizeof(p.window));
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzwrite( backup_file, &p.TF_tag, sizeof(p.TF_tag));
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  gzwrite( backup_file, &p.basal_expression_level, sizeof(p.basal_expression_level) );
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  gzwrite( backup_file, &p.functional, sizeof(p.functional) );
}

/**
 * \brief    Increase buffer size relatively to the new inherited proteins size
 * \details  --
 * \param    size_t new_size
 * \return   \e void
 */
void InheritedProteins::increase_buffer_size( size_t new_size )
{
  assert(new_size <= _string_of_pearls->size*2);
  _string_of_pearls->buffer_size = (new_size/INHERITED_PROTEINS_BUFFER+1)*INHERITED_PROTEINS_BUFFER;
  pearl* new_x = new pearl[_string_of_pearls->buffer_size];
  memcpy(new_x, _string_of_pearls->x, sizeof(pearl)*_string_of_pearls->size);
  delete[] _string_of_pearls->x;
  _string_of_pearls->x = new_x;
}

/**
 * \brief    Decrease buffer size relatively to the inherited proteins size
 * \details  --
 * \param    void
 * \return   \e void
 */
void InheritedProteins::decrease_buffer_size( void )
{
  if (_string_of_pearls->buffer_size > INHERITED_PROTEINS_BUFFER)
  {
    _string_of_pearls->buffer_size = (_string_of_pearls->size/INHERITED_PROTEINS_BUFFER+1)*INHERITED_PROTEINS_BUFFER;
    assert(_string_of_pearls->buffer_size >= _string_of_pearls->size);
    pearl* new_x = new pearl[_string_of_pearls->buffer_size];
    memcpy(new_x, _string_of_pearls->x, sizeof(pearl)*_string_of_pearls->size);
    delete[] _string_of_pearls->x;
    _string_of_pearls->x = new_x;
  }
}
