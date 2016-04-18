
/**
 * \file      Environment.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Environment class definition
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

#include "Environment.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \return   \e void
 */
Environment::Environment( Parameters* parameters )
{
  _parameters = parameters;
  _width      = parameters->get_width();
  _height     = parameters->get_height();
  _grid       = new SpeciesList*[_width*_height];
  for (size_t i = 0; i < _width*_height; i++)
  {
    _grid[i] = new SpeciesList();
  }
  _X                  = new SpeciesList();
  _species_lists_size = SPECIES_LIST_BUFFER;
  _total_amount       = 0.0;
  _min_amount         = 0.0;
  _max_amount         = 0.0;
  _inflowing_amount   = 0.0;
  _outflowing_amount  = 0.0;
#ifdef DEBUG
  for (size_t i = 0; i < _width*_height; i++)
  {
    assert(_grid[i]->get_size() == _species_lists_size);
  }
  assert(_X->get_size() == _species_lists_size);
#endif
}

/**
 * \brief    Constructor from backup file
 * \details  Load Environment class from backup file
 * \param    Parameters* parameters
 * \param    gzFile backup_file
 * \return   \e void
 */
Environment::Environment( Parameters* parameters, gzFile backup_file )
{
  _parameters = parameters;
  gzread( backup_file, &_width,  sizeof(_width) );
  gzread( backup_file, &_height, sizeof(_height) );
  _grid = new SpeciesList*[_width*_height];
  for (size_t i = 0; i < _width*_height; i++)
  {
    _grid[i] = new SpeciesList(backup_file);
  }
  _X = new SpeciesList(backup_file);
  gzread( backup_file, &_species_lists_size, sizeof(_species_lists_size) );
  gzread( backup_file, &_total_amount,       sizeof(_total_amount) );
  gzread( backup_file, &_min_amount,         sizeof(_min_amount) );
  gzread( backup_file, &_max_amount,         sizeof(_max_amount) );
  gzread( backup_file, &_inflowing_amount,   sizeof(_inflowing_amount) );
  gzread( backup_file, &_outflowing_amount,  sizeof(_outflowing_amount) );
#ifdef DEBUG
  for (size_t i = 0; i < _width*_height; i++)
  {
    assert(_grid[i]->get_size() == _species_lists_size);
  }
  assert(_X->get_size() == _species_lists_size);
#endif
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const Environment& environment 
 * \return   \e void
 */
Environment::Environment( Environment const &environment )
{
  _parameters = environment._parameters;
  _width      = environment._width;
  _height     = environment._height;
  _grid       = new SpeciesList*[_width*_height];
  for (size_t i = 0; i < _width*_height; i++)
  {
    _grid[i] = new SpeciesList(*environment._grid[i]);
  }
  _X                  = new SpeciesList(*environment._X);
  _species_lists_size = environment._species_lists_size;
  _total_amount       = environment._total_amount;
  _min_amount         = environment._min_amount;
  _max_amount         = environment._max_amount;
  _inflowing_amount   = environment._inflowing_amount;
  _outflowing_amount  = environment._outflowing_amount;
#ifdef DEBUG
  for (size_t i = 0; i < _width*_height; i++)
  {
    assert(_grid[i]->get_size() == _species_lists_size);
  }
  assert(_X->get_size() == _species_lists_size);
#endif
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
Environment::~Environment( void )
{
  for (size_t i = 0; i < _width*_height; i++)
  {
    delete _grid[i];
    _grid[i] = NULL;
  }
  delete[] _grid;
  _grid = NULL;
  delete _X;
  _X = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Update environment amount
 * \details  Compute total amount and X vector
 * \param    void
 * \return   \e void
 */
void Environment::update_environment_amount( void )
{
  _total_amount = 0.0;
  _min_amount   = 1e+06;
  _max_amount   = 0.0;
  _X->reset(false);
  assert(_X->get_size() == _species_lists_size);
  for (size_t pos = 0; pos < _width*_height; pos++)
  {
    assert(_grid[pos]->get_size() == _species_lists_size);
    for (int met = 0; met < (int)_species_lists_size; met++)
    {
      _X->add(met+1, _grid[pos]->get(met+1));
    }
    _total_amount  += _grid[pos]->get_amount();
    if (_min_amount > _grid[pos]->get_amount())
    {
      _min_amount = _grid[pos]->get_amount();
    }
    if (_max_amount < _grid[pos]->get_amount())
    {
      _max_amount = _grid[pos]->get_amount();
    }
  }
}

/**
 * \brief    Compute diffusion and degradation
 * \details  Compute diffusion and degradation of metabolites in environment
 * \param    void
 * \return   \e void
 */
void Environment::compute_diffusion_and_degradation( void )
{
  /*-----------------------------------------------------------------*/
  /* 1) Copy the grid content and compute diffusion first step       */
  /*-----------------------------------------------------------------*/
  double* new_grid = new double[_width*_height*_species_lists_size];
  double diff_rate = _parameters->get_environment_properties()->diffusion_rate;
  double degr_rate = _parameters->get_environment_properties()->degradation_rate;
  for (size_t x = 0; x < _width; x++)
  {
    for (size_t y = 0; y < _height; y++)
    {
      memcpy(&new_grid[x*_height*_species_lists_size+y*_species_lists_size], _grid[x*_height+y]->get_X(), sizeof(double)*_species_lists_size);
      for (int i = -1; i < 2; i++)
      {
        for (int j = -1; j < 2; j++)
        {
          size_t  cur_x = (x + i + _width)  % _width;
          size_t  cur_y = (y + j + _height) % _height;
          double* cur_v = _grid[cur_x*_height+cur_y]->get_X();
          for (size_t index = 0; index < _species_lists_size; index++)
          {
            new_grid[x*_height*_species_lists_size+y*_species_lists_size+index] += cur_v[index]*diff_rate;
          }
        }
      }
    }
  }
  /*-----------------------------------------------------------------*/
  /* 2) Compute the second step of diffusion and compute degradation */
  /*-----------------------------------------------------------------*/
  for (size_t x = 0; x < _width; x++)
  {
    for (size_t y = 0; y < _height; y++)
    {
      for (int index = 0; index < (int)_species_lists_size; index++)
      {
        double fresh_amount = new_grid[x*_height*_species_lists_size+y*_species_lists_size+index] - 9.0*_grid[x*_height+y]->get_X()[index]*diff_rate;
        double new_amount   = fresh_amount*(1.0-degr_rate);
        _outflowing_amount += fresh_amount*degr_rate;
        _grid[x*_height+y]->set(index+1, new_amount);
      }
    }
  }
  delete[] new_grid;
  new_grid = NULL;
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void Environment::save( gzFile backup_file )
{
  gzwrite( backup_file, &_width,  sizeof(_width) );
  gzwrite( backup_file, &_height, sizeof(_height) );
  for (size_t i = 0; i < _width*_height; i++)
  {
    _grid[i]->save(backup_file);
  }
  _X->save(backup_file);
  gzwrite( backup_file, &_species_lists_size, sizeof(_species_lists_size) );
  gzwrite( backup_file, &_total_amount,       sizeof(_total_amount) );
  gzwrite( backup_file, &_min_amount,         sizeof(_min_amount) );
  gzwrite( backup_file, &_max_amount,         sizeof(_max_amount) );
  gzwrite( backup_file, &_inflowing_amount,   sizeof(_inflowing_amount) );
  gzwrite( backup_file, &_outflowing_amount,  sizeof(_outflowing_amount) );
}

/**
 * \brief    Write the 1D metabolic state vector of the environment
 * \details  --
 * \param    std::string filename
 * \return   \e void
 */
void Environment::write_metabolic_state_vector( std::string filename )
{
  if (_species_lists_size > 0)
  {
    std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
    for (int i = 0; i < (int)_species_lists_size; i++)
    {
      file << i+1 << " " << _X->get(i+1) << "\n";
    }
    file.close();
  }
}

/**
 * \brief    Write the 1D metabolic state vector of one local environment
 * \details  --
 * \param    std::string filename
 * \param    size_t pos
 * \return   \e void
 */
void Environment::write_local_metabolic_state_vector( std::string filename, size_t pos )
{
  assert(pos < _width*_height);
  if (_species_lists_size > 0)
  {
    std::ofstream file(filename.c_str(), std::ios::out | std::ios::trunc);
    for (int i = 0; i < (int)_species_lists_size; i++)
    {
      file << i+1 << " " << _grid[pos]->get(i+1) << "\n";
    }
    file.close();
  }
}

/**
 * \brief    Clear the environment
 * \details  Remove all metabolites from the environment
 * \param    void
 * \return   \e void
 */
void Environment::clear_environment( void )
{
  for (size_t i = 0; i < _width*_height; i++)
  {
    _outflowing_amount += _grid[i]->get_amount();
    _grid[i]->reset(false);
  }
  _X->reset(false);
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/
