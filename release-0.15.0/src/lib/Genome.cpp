
/**
 * \file      Genome.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Genome class definition
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

#include "Genome.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    Parameters* parameters
 * \param    ReplicationReport* replication_report
 * \return   \e void
 */
Genome::Genome( Parameters* parameters, ReplicationReport* replication_report )
{
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  _prng = NULL;
  if (parameters != NULL)
  {
    _prng = parameters->get_prng();
  }
  
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = parameters;
  
  /*------------------------------------------------------------------ string of pearls */
  
  create_string_of_pearls();
  
  /*------------------------------------------------------------------ concentration vector */
  
  _concentration_vector = NULL;
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_NC            = 0;
  _nb_E             = 0;
  _nb_TF            = 0;
  _nb_BS            = 0;
  _nb_P             = 0;
  _nb_inner_enzymes = 0;
  _nb_inflow_pumps  = 0;
  _nb_outflow_pumps = 0;
  
  /*------------------------------------------------------------------ replication report */
  
  _replication_report = replication_report;
  
  /*------------------------------------------------------------------ pearl positions */
  
  _TFi = NULL;
  _Pi  = NULL;
}

/**
 * \brief    Constructor from backup file
 * \details  Load Genome class from backup file
 * \param    Parameters* parameters
 * \param    ReplicationReport* replication_report
 * \param    gzFile backup_file
 * \return   \e void
 */
Genome::Genome( Parameters* parameters, ReplicationReport* replication_report, gzFile backup_file )
{
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  _prng = NULL;
  if (parameters != NULL)
  {
    _prng = parameters->get_prng();
  }
  
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
  
  gzread( backup_file, &_nb_NC,            sizeof(_nb_NC) );
  gzread( backup_file, &_nb_E,             sizeof(_nb_E) );
  gzread( backup_file, &_nb_TF,            sizeof(_nb_TF) );
  gzread( backup_file, &_nb_BS,            sizeof(_nb_BS) );
  gzread( backup_file, &_nb_P,             sizeof(_nb_P) );
  gzread( backup_file, &_nb_inner_enzymes, sizeof(_nb_inner_enzymes) );
  gzread( backup_file, &_nb_inflow_pumps,  sizeof(_nb_inflow_pumps) );
  gzread( backup_file, &_nb_outflow_pumps, sizeof(_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ replication report */
  
  _replication_report = replication_report;
  
  /*------------------------------------------------------------------ pearl positions */
  
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    for (size_t i = 0; i < _nb_TF; i++)
    {
      gzread( backup_file, &_TFi[i], sizeof(_TFi[i]) );
    }
  }
  _Pi = NULL;
  if (_nb_P > 0)
  {
    _Pi = new size_t[_nb_P];
    for (size_t i = 0; i < _nb_P; i++)
    {
      gzread( backup_file, &_Pi[i], sizeof(_Pi[i]) );
    }
  }
}

/**
 * \brief    Copy constructor
 * \details  --
 * \param    const Genome& genome
 * \param    ReplicationReport* replication_report
 * \return   \e void
 */
Genome::Genome( const Genome& genome, ReplicationReport* replication_report )
{
  /*------------------------------------------------------------------ pseudorandom numbers generator */
  
  _prng = genome._prng;
  
  /*------------------------------------------------------------------ simulation parameters */
  
  _parameters = genome._parameters;
  
  /*------------------------------------------------------------------ string of pearls */
  
  copy_string_of_pearls(genome._string_of_pearls);
  
  /*------------------------------------------------------------------ concentration vector */
  
  _concentration_vector = NULL;
  if (_string_of_pearls->size)
  {
    _concentration_vector = new double[_string_of_pearls->size];
    memcpy(_concentration_vector, genome._concentration_vector, sizeof(double)*_string_of_pearls->size);
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_NC            = genome._nb_NC;
  _nb_E             = genome._nb_E;
  _nb_TF            = genome._nb_TF;
  _nb_BS            = genome._nb_BS;
  _nb_P             = genome._nb_P;
  _nb_inner_enzymes = genome._nb_inner_enzymes;
  _nb_inflow_pumps  = genome._nb_inflow_pumps;
  _nb_outflow_pumps = genome._nb_outflow_pumps;
  
  /*------------------------------------------------------------------ replication report */
  
  _replication_report = replication_report;
  
  /*------------------------------------------------------------------ pearl positions */
  
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    memcpy(_TFi, genome._TFi, sizeof(size_t)*_nb_TF);
  }
  _Pi = NULL;
  if (_nb_P > 0)
  {
    _Pi = new size_t[_nb_P];
    memcpy(_Pi, genome._Pi, sizeof(size_t)*_nb_P);
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
Genome::~Genome( void )
{
  delete_string_of_pearls();
  delete[] _concentration_vector;
  _concentration_vector = NULL;
  delete[] _TFi;
  _TFi = NULL;
  delete[] _Pi;
  _Pi = NULL;
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
void Genome::initialize_concentration_vector( void )
{
  delete[] _concentration_vector;
  _concentration_vector = new double[_string_of_pearls->size];
  for (size_t i = 0; i < _string_of_pearls->size; i++)
  {
    _concentration_vector[i] = 0.0;
  }
}

/**
 * \brief    Mutate the genome
 * \details  Apply rearrangements, then apply point mutations
 * \param    const double* mutation_rates
 * \return   \e void
 */
void Genome::mutate( const double* mutation_rates )
{
  /*--------------------------------------------------*/
  /* 1) clear the replication report                  */
  /*--------------------------------------------------*/
  //_replication_report->clear();
  
  /*--------------------------------------------------*/
  /* 2) do rearrangements                             */
  /*--------------------------------------------------*/
  do_rearrangements(mutation_rates);
  
  /*--------------------------------------------------*/
  /* 3) do horizontal gene transfer                   */
  /*--------------------------------------------------*/
  if (_parameters->get_hgt_rate() > 0.0)
  {
    do_hgt();
  }
  
  /*--------------------------------------------------*/
  /* 4) do point mutations                            */
  /*--------------------------------------------------*/
  do_point_mutations(mutation_rates);
  
  /*--------------------------------------------------*/
  /* 5) compute the indexes list for each pearl types */
  /*--------------------------------------------------*/
  pearl* str = _string_of_pearls->x;
  delete[] _TFi;
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
  }
  delete[] _Pi;
  _Pi = NULL;
  if (_nb_P > 0)
  {
    _Pi  = new size_t[_nb_P];
  }
  size_t TFcount = 0;
  size_t Pcount  = 0;
  for (size_t pos = 0; pos < _string_of_pearls->size; pos++)
  {
    if (str[pos].type == TRANSCRIPTION_FACTOR)
    {
      _TFi[TFcount] = pos;
      TFcount++;
    }
    else if (str[pos].type == PROMOTER)
    {
      _Pi[Pcount] = pos;
      Pcount++;
    }
    str[pos].functional = false;
  }
  
  /*--------------------------------------------------*/
  /* 6) initialize concentration vector               */
  /*--------------------------------------------------*/
  initialize_concentration_vector();
}

/**
 * \brief    Shuffle the genome
 * \details  Shuffle genome order
 * \param    void
 * \return   \e void
 */
void Genome::shuffle( void )
{
  for (size_t i = 0; i < _string_of_pearls->size; i++)
  {
    size_t pos1 = (size_t)_prng->uniform(0, (int)_string_of_pearls->size-1);
    size_t pos2 = (size_t)_prng->uniform(0, (int)_string_of_pearls->size-1);
    pearl tmp   = _string_of_pearls->x[pos1];
    _string_of_pearls->x[pos1] = _string_of_pearls->x[pos2];
    _string_of_pearls->x[pos2] = tmp;
  }
}

/**
 * \brief    Save in backup file
 * \details  --
 * \param    gzFile backup_file
 * \return   \e void
 */
void Genome::save( gzFile backup_file )
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
  
  gzwrite( backup_file, &_nb_NC,            sizeof(_nb_NC) );
  gzwrite( backup_file, &_nb_E,             sizeof(_nb_E) );
  gzwrite( backup_file, &_nb_TF,            sizeof(_nb_TF) );
  gzwrite( backup_file, &_nb_BS,            sizeof(_nb_BS) );
  gzwrite( backup_file, &_nb_P,             sizeof(_nb_P) );
  gzwrite( backup_file, &_nb_inner_enzymes, sizeof(_nb_inner_enzymes) );
  gzwrite( backup_file, &_nb_inflow_pumps,  sizeof(_nb_inflow_pumps) );
  gzwrite( backup_file, &_nb_outflow_pumps, sizeof(_nb_outflow_pumps) );
  
  /*------------------------------------------------------------------ pearl positions */
  
  if (_nb_TF > 0)
  {
    for (size_t i = 0; i < _nb_TF; i++)
    {
      gzwrite( backup_file, &_TFi[i], sizeof(_TFi[i]) );
    }
  }
  if (_nb_P > 0)
  {
    for (size_t i = 0; i < _nb_P; i++)
    {
      gzwrite( backup_file, &_Pi[i], sizeof(_Pi[i]) );
    }
  }
}

/**
 * \brief    Replace genome data
 * \details  --
 * \param    Genome* genome
 * \return   \e void
 */
void Genome::replace_data( Genome* genome )
{
  /*------------------------------------------------------------------ string of pearls */
  
  delete_string_of_pearls();
  copy_string_of_pearls(genome->get_string_of_pearl());
  
  /*------------------------------------------------------------------ concentration vector */
  
  delete[] _concentration_vector;
  _concentration_vector = NULL;
  if (_string_of_pearls->size > 0)
  {
    _concentration_vector = new double[_string_of_pearls->size];
    memcpy(_concentration_vector, genome->get_concentration_vector(), sizeof(double)*_string_of_pearls->size);
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  _nb_NC            = genome->get_nb_NC();
  _nb_E             = genome->get_nb_E();
  _nb_TF            = genome->get_nb_TF();
  _nb_BS            = genome->get_nb_BS();
  _nb_P             = genome->get_nb_P();
  _nb_inner_enzymes = genome->get_nb_inner_enzymes();
  _nb_inflow_pumps  = genome->get_nb_inflow_pumps();
  _nb_outflow_pumps = genome->get_nb_outflow_pumps();
  
  /*------------------------------------------------------------------ pearl positions */
  
  delete[] _TFi;
  _TFi = NULL;
  if (_nb_TF > 0)
  {
    _TFi = new size_t[_nb_TF];
    memcpy(_TFi, genome->get_TFi(), sizeof(size_t)*_nb_TF);
  }
  delete[] _Pi;
  _Pi = NULL;
  if (_nb_P > 0)
  {
    _Pi = new size_t[_nb_P];
    memcpy(_Pi, genome->get_Pi(), sizeof(size_t)*_nb_P);
  }
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Do point mutations
 * \details  Browse the genome and apply Gene::mutate() on each gene
 * \param    const double* mutation_rates
 * \return   \e void
 */
void Genome::do_point_mutations( const double* mutation_rates )
{
  /*--------------------------------------------------*/
  /* 1) mutate pearls and count pearl types           */
  /*--------------------------------------------------*/
  _nb_NC            = 0;
  _nb_E             = 0;
  _nb_TF            = 0;
  _nb_BS            = 0;
  _nb_P             = 0;
  _nb_inner_enzymes = 0;
  _nb_inflow_pumps  = 0;
  _nb_outflow_pumps = 0;
  pearl* str        = _string_of_pearls->x;
  for (size_t pos = 0; pos < _string_of_pearls->size; pos++)
  {
    /*-------------------------*/
    /* 1.1) mutate the gene    */
    /*-------------------------*/
    do_pearl_mutation(mutation_rates, pos);
    /*-------------------------*/
    /* 1.2) compute statistics */
    /*-------------------------*/
    if (str[pos].type == NON_CODING)
    {
      _nb_NC++;
    }
    else if (str[pos].type == ENZYME)
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
    else if (str[pos].type == BINDING_SITE)
    {
      _nb_BS++;
    }
    else if (str[pos].type == PROMOTER)
    {
      _nb_P++;
    }
  }
}

/*
 * \brief    Do rearrangements (duplications + deletions + translocations + inversions)
 * \details  Each breakpoint can mutate one of the two adjacent genes
 * \param    const double* mutation_rates
 * \return   \e void
 */
void Genome::do_rearrangements( const double* mutation_rates )
{
  if (_string_of_pearls->size == 0)
  {
    return;
  }
  size_t ndup = _prng->binomial(_string_of_pearls->size, mutation_rates[DUPLICATION_RATE]);
  size_t ndel = _prng->binomial(_string_of_pearls->size, mutation_rates[DELETION_RATE]);
  size_t ntra = _prng->binomial(_string_of_pearls->size, mutation_rates[TRANSLOCATION_RATE]);
  size_t ninv = _prng->binomial(_string_of_pearls->size, mutation_rates[INVERSION_RATE]);
  size_t nrea = ndup + ndel + ntra + ninv;
  for (size_t i = nrea; i >= 1; i--)
  {
    size_t draw = (size_t)_prng->uniform(0, (int)i);
    /*------------------*/
    /* do duplication   */
    /*------------------*/
    if ( draw < ndup )
    {
      do_duplication();
      ndup--;
    }
    /*------------------*/
    /* do deletion      */
    /*------------------*/
    else if ( draw < ndup + ndel )
    {
      do_deletion();
      ndel--;
    }
    /*------------------*/
    /* do translocation */
    /*------------------*/
    else if ( draw < ndup + ndel + ntra )
    {
      do_translocation();
      ntra--;
    }
    /*------------------*/
    /* do inversion     */
    /*------------------*/
    else
    {
      do_inversion();
      ninv--;
    }
    check_genome_size();
  }
}

/**
 * \brief    Do a duplication
 * \details  Draw uniformly two breakpoints in the genome (start and end), copy the sequence at a third breakpoint (insert). Pasted pearls at insertion point undergo BLX-a crossover
 * \param    void
 * \return   \e void
 */
void Genome::do_duplication( void )
{
  if (_string_of_pearls->size == 0)
  {
    return;
  }
  size_t start     = (size_t)_prng->uniform(0, (int)_string_of_pearls->size-1);
  size_t end       = (size_t)_prng->uniform(0, (int)_string_of_pearls->size-1);
  size_t insert    = (size_t)_prng->uniform(0, (int)_string_of_pearls->size-1);
  pearl* duplicata = NULL;
  pearl* shift     = NULL;
  size_t size      = 0;
  size_t nb_NC     = 0;
  size_t nb_E      = 0;
  size_t nb_TF     = 0;
  size_t nb_BS     = 0;
  size_t nb_P      = 0;
  
  /*----------------------------------------------------------------------*/
  /* 1) if start <= end                                                   */
  /*----------------------------------------------------------------------*/
  if (start <= end)
  {
    /* 1.1) compute the duplicata size -------*/
    size = end-start+1;
    
    /* 1.2) evaluate buffer size -------------*/
    if (_string_of_pearls->size+size > _string_of_pearls->buffer_size)
    {
      increase_buffer_size(_string_of_pearls->size+size);
    }
    
    /* 1.3) copy the duplicata ---------------*/
    duplicata = new pearl[size];
    memcpy(duplicata, &_string_of_pearls->x[start], sizeof(pearl)*size);
    
    /* 1.4) copy the piece of genome to shift */
    shift = new pearl[_string_of_pearls->size-insert];
    memcpy(shift, &_string_of_pearls->x[insert], sizeof(pearl)*(_string_of_pearls->size-insert));
    
    /* 1.5) count the number of each type ----*/
    for (size_t i = 0; i < size; i++)
    {
      if (duplicata[i].type == NON_CODING)
      {
        nb_NC++;
      }
      else if (duplicata[i].type == ENZYME)
      {
        nb_E++;
      }
      else if (duplicata[i].type == TRANSCRIPTION_FACTOR)
      {
        nb_TF++;
      }
      else if (duplicata[i].type == BINDING_SITE)
      {
        nb_BS++;
      }
      else if (duplicata[i].type == PROMOTER)
      {
        nb_P++;
      }
    }
  }
  
  /*----------------------------------------------------------------------*/
  /* 2) else if start > end                                               */
  /*----------------------------------------------------------------------*/
  else if (start > end)
  {
    /* 2.1) compute the duplicata size -------*/
    size = _string_of_pearls->size-start+end+1;
    
    /* 2.2) evaluate buffer size -------------*/
    if (_string_of_pearls->size+size > _string_of_pearls->buffer_size)
    {
      increase_buffer_size(_string_of_pearls->size+size);
    }
    
    /* 2.3) copy the duplicata ---------------*/
    duplicata = new pearl[size];
    memcpy(duplicata, &_string_of_pearls->x[start], sizeof(pearl)*(_string_of_pearls->size-start));
    memcpy(&duplicata[_string_of_pearls->size-start], _string_of_pearls->x, sizeof(pearl)*(end+1));
    
    /* 2.4) copy the piece of genome to shift */
    shift = new pearl[_string_of_pearls->size-insert];
    memcpy(shift, &_string_of_pearls->x[insert], sizeof(pearl)*(_string_of_pearls->size-insert));
    
    /* 2.5) count the number of each type ----*/
    for (size_t i = 0; i < size; i++)
    {
      if (duplicata[i].type == NON_CODING)
      {
        nb_NC++;
      }
      else if (duplicata[i].type == ENZYME)
      {
        nb_E++;
      }
      else if (duplicata[i].type == TRANSCRIPTION_FACTOR)
      {
        nb_TF++;
      }
      else if (duplicata[i].type == BINDING_SITE)
      {
        nb_BS++;
      }
      else if (duplicata[i].type == PROMOTER)
      {
        nb_P++;
      }
    }
  }
  
  /*----------------------------------------------------------------------*/
  /* 3) rebuild the sequence from both the duplicata and the shift pieces */
  /*----------------------------------------------------------------------*/
  
  /* 3.1) insert the duplicata at the right place ----------------*/
  memcpy(&_string_of_pearls->x[insert], duplicata, sizeof(pearl)*size);
  
  /* 3.2) insert the piece of genome to shift after the duplicata */
  memcpy(&_string_of_pearls->x[insert+size], shift, sizeof(pearl)*(_string_of_pearls->size-insert));
  
  /* 3.3) update size --------------------------------------------*/
  _string_of_pearls->size += size;
  
  /*----------------------------------------------------------------------*/
  /* 4) delete sequences                                                  */
  /*----------------------------------------------------------------------*/
  delete[] duplicata;
  duplicata = NULL;
  delete[] shift;
  shift = NULL;
  
  /*----------------------------------------------------------------------*/
  /* 5) perform crossover                                                 */
  /*----------------------------------------------------------------------*/
  size_t pos1 = (insert-1+_string_of_pearls->size)%_string_of_pearls->size;
  size_t pos2 = (insert+_string_of_pearls->size)%_string_of_pearls->size;
  //do_BLX_alpha_crossover(pos1, pos2);
  do_permutation(pos1, pos2);
  pos1 = (insert+size-1+_string_of_pearls->size)%_string_of_pearls->size;
  pos2 = (insert+size+_string_of_pearls->size)%_string_of_pearls->size;
  //do_BLX_alpha_crossover(pos1, pos2);
  do_permutation(pos1, pos2);
  
  /*----------------------------------------------------------------------*/
  /* 6) save mutation events                                              */
  /*----------------------------------------------------------------------*/
  MutationEvent* new_event = new MutationEvent(DUPLICATION, start, end, insert, size, nb_NC, nb_E, nb_TF, nb_BS, nb_P);
  _replication_report->add_mutation_event(new_event);
}

/**
 * \brief    Do a deletion
 * \details  Draw uniformly two breakpoints in the genome (start and end), delete the sequence. Pasted pearls at insertion point undergo BLX-a crossover
 * \param    void
 * \return   \e void
 */
void Genome::do_deletion( void )
{
  if (_string_of_pearls->size == 0)
  {
    return;
  }
  size_t start = (size_t)_prng->uniform(0, (int)_string_of_pearls->size-1);
  size_t end   = (size_t)_prng->uniform(0, (int)_string_of_pearls->size-1);
  size_t size  = 0;
  size_t nb_NC = 0;
  size_t nb_E  = 0;
  size_t nb_TF = 0;
  size_t nb_BS = 0;
  size_t nb_P  = 0;
  
  size_t index = start;
  while (index != end)
  {
    if (_string_of_pearls->x[index].type == NON_CODING)
    {
      nb_NC++;
    }
    else if (_string_of_pearls->x[index].type == ENZYME)
    {
      nb_E++;
    }
    else if (_string_of_pearls->x[index].type == TRANSCRIPTION_FACTOR)
    {
      nb_TF++;
    }
    else if (_string_of_pearls->x[index].type == BINDING_SITE)
    {
      nb_BS++;
    }
    else if (_string_of_pearls->x[index].type == PROMOTER)
    {
      nb_P++;
    }
    index = (index+1+_string_of_pearls->size)%_string_of_pearls->size;
  }
  if (_string_of_pearls->x[index].type == NON_CODING)
  {
    nb_NC++;
  }
  else if (_string_of_pearls->x[index].type == ENZYME)
  {
    nb_E++;
  }
  else if (_string_of_pearls->x[index].type == TRANSCRIPTION_FACTOR)
  {
    nb_TF++;
  }
  else if (_string_of_pearls->x[index].type == BINDING_SITE)
  {
    nb_BS++;
  }
  else if (_string_of_pearls->x[index].type == PROMOTER)
  {
    nb_P++;
  }
  
  /*----------------------------------------------------------------------*/
  /* 1) if start <= end                                                   */
  /*----------------------------------------------------------------------*/
  if (start <= end)
  {
    /* 2.1) perform cross-over -----------*/
    size_t pos1 = (start-1+_string_of_pearls->size)%_string_of_pearls->size;
    size_t pos2 = (end+1+_string_of_pearls->size)%_string_of_pearls->size;
    //do_BLX_alpha_crossover(pos1, pos2);
    do_permutation(pos1, pos2);
    
    /* 2.2) compute the deletion size ----*/
    size = end-start+1;
    
    /* 2.3) delete the sequence ----------*/
    if (_string_of_pearls->size-end-1 > 0)
    {
      memcpy(&_string_of_pearls->x[start], &_string_of_pearls->x[end+1], sizeof(pearl)*(_string_of_pearls->size-end-1));
    }
    _string_of_pearls->size -= size;
    
    /* 2.4) evaluate buffer size ---------*/
    if (_string_of_pearls->size <= _string_of_pearls->buffer_size-3*GENOME_BUFFER/2)
    {
      decrease_buffer_size();
    }
  }
  
  /*----------------------------------------------------------------------*/
  /* 2) else if start > end                                               */
  /*----------------------------------------------------------------------*/
  else if (start > end)
  {
    /* 3.1) perform cross-over -----------*/
    size_t pos1 = (end-1+_string_of_pearls->size)%_string_of_pearls->size;
    size_t pos2 = (start+1+_string_of_pearls->size)%_string_of_pearls->size;
    //do_BLX_alpha_crossover(pos1, pos2);
    do_permutation(pos1, pos2);
    
    /* 3.2) compute the duplicata size ---*/
    size = _string_of_pearls->size-start+end+1;
    
    /* 3.3) delete the sequence ----------*/
    if (start-end-1 > 0)
    {
      memcpy(_string_of_pearls->x, &_string_of_pearls->x[end+1], sizeof(pearl)*(start-end-1));
    }
    _string_of_pearls->size = size;
    
    /* 3.4) evaluate buffer size ---------*/
    if (_string_of_pearls->size < _string_of_pearls->buffer_size-3*GENOME_BUFFER/2)
    {
      decrease_buffer_size();
    }
  }
  
  /*----------------------------------------------------------------------*/
  /* 3) save mutation event                                               */
  /*----------------------------------------------------------------------*/
  MutationEvent* new_event = new MutationEvent(DELETION, start, end, 0, size, nb_NC, nb_E, nb_TF, nb_BS, nb_P);
  _replication_report->add_mutation_event(new_event);
}

/**
 * \brief    Do a translocation
 * \details  Draw uniformly two breakpoints in the genome (start and end), translocate the sequence at a third breakpoint. Pasted pearls at insertion point undergo BLX-a crossover
 * \param    void
 * \return   \e void
 */
void Genome::do_translocation( void )
{
  if (_string_of_pearls->size == 0)
  {
    return;
  }
  size_t start     = (size_t)_prng->uniform(0, (int)_string_of_pearls->size-1);
  size_t end       = (size_t)_prng->uniform(0, (int)_string_of_pearls->size-1);
  size_t insert    = 0;
  pearl* duplicata = NULL;
  pearl* shift     = NULL;
  size_t size      = 0;
  size_t nb_NC     = 0;
  size_t nb_E      = 0;
  size_t nb_TF     = 0;
  size_t nb_BS     = 0;
  size_t nb_P      = 0;
  
  /*----------------------------------------------------------------------*/
  /* 1) if start <= end                                                   */
  /*----------------------------------------------------------------------*/
  if (start <= end)
  {
    /* 1.1) compute the duplicata size -------*/
    size = end-start+1;
    
    /* 1.2) copy the duplicata ---------------*/
    duplicata = new pearl[size];
    memcpy(duplicata, &_string_of_pearls->x[start], sizeof(pearl)*size);
    
    /* 1.3) delete the sequence --------------*/
    memcpy(&_string_of_pearls->x[start], &_string_of_pearls->x[end+1], sizeof(pearl)*(_string_of_pearls->size-end-1));
    _string_of_pearls->size -= size;
    
    /* 1.4) draw insert index ----------------*/
    if (_string_of_pearls->size == 0)
    {
      insert = 0;
    }
    else
    {
      insert = _prng->uniform(0, (int)_string_of_pearls->size-1);
    }
    
    /* 1.5) copy the piece of genome to shift */
    shift = new pearl[_string_of_pearls->size-insert];
    memcpy(shift, &_string_of_pearls->x[insert], sizeof(pearl)*(_string_of_pearls->size-insert));
  }
  
  /*----------------------------------------------------------------------*/
  /* 2) else if start > end                                               */
  /*----------------------------------------------------------------------*/
  else if (start > end)
  {
    /* 2.1) compute the duplicata size -------*/
    size = _string_of_pearls->size-start+end+1;
    
    /* 2.2) copy the duplicata ---------------*/
    duplicata = new pearl[size];
    memcpy(duplicata, &_string_of_pearls->x[start], sizeof(pearl)*(_string_of_pearls->size-start));
    memcpy(&duplicata[_string_of_pearls->size-start], _string_of_pearls->x, sizeof(pearl)*(end+1));
    
    /* 2.3) delete the sequence --------------*/
    memcpy(_string_of_pearls->x, &_string_of_pearls->x[end+1], sizeof(pearl)*(start-end-1));
    _string_of_pearls->size -= size;
    
    /* 2.4) draw insert index ----------------*/
    if (_string_of_pearls->size == 0)
    {
      insert = 0;
    }
    else
    {
      insert = _prng->uniform(0, (int)_string_of_pearls->size-1);
    }
    
    /* 2.5) copy the piece of genome to shift */
    shift = new pearl[_string_of_pearls->size-insert];
    memcpy(shift, &_string_of_pearls->x[insert], sizeof(pearl)*(_string_of_pearls->size-insert));
  }
  
  /*----------------------------------------------------------------------*/
  /* 3) rebuild the sequence from both the duplicata and the shift pieces */
  /*----------------------------------------------------------------------*/
  
  /* 3.1) insert the duplicata at the right place ----------------*/
  memcpy(&_string_of_pearls->x[insert], duplicata, sizeof(pearl)*size);
  
  /* 3.2) insert the piece of genome to shift after the duplicata */
  memcpy(&_string_of_pearls->x[insert+size], shift, sizeof(pearl)*(_string_of_pearls->size-insert));
  
  /* 3.3) update size --------------------------------------------*/
  _string_of_pearls->size += size;
  
  /*----------------------------------------------------------------------*/
  /* 4) delete sequences                                                  */
  /*----------------------------------------------------------------------*/
  delete[] duplicata;
  duplicata = NULL;
  delete[] shift;
  shift = NULL;
  
  /*----------------------------------------------------------------------*/
  /* 5) perform crossover                                                 */
  /*----------------------------------------------------------------------*/
  size_t pos1 = (insert-1+_string_of_pearls->size)%_string_of_pearls->size;
  size_t pos2 = (insert+_string_of_pearls->size)%_string_of_pearls->size;
  //do_BLX_alpha_crossover(pos1, pos2);
  do_permutation(pos1, pos2);
  pos1 = (insert+size-1+_string_of_pearls->size)%_string_of_pearls->size;
  pos2 = (insert+size+_string_of_pearls->size)%_string_of_pearls->size;
  //do_BLX_alpha_crossover(pos1, pos2);
  do_permutation(pos1, pos2);
  
  /*----------------------------------------------------------------------*/
  /* 6) save mutation event                                               */
  /*----------------------------------------------------------------------*/
  MutationEvent* new_event = new MutationEvent(TRANSLOCATION, start, end, insert, size, nb_NC, nb_E, nb_TF, nb_BS, nb_P);
  _replication_report->add_mutation_event(new_event);
}

/**
 * \brief    Do an inversion
 * \details  Draw uniformly two breakpoints in the genome (start end end), revert the sequence. Pasted pearls at insertion point undergo BLX-a crossover
 * \param    void
 * \return   \e void
 */
void Genome::do_inversion( void )
{
  if (_string_of_pearls->size == 0)
  {
    return;
  }
  size_t start         = (size_t)_prng->uniform(0, (int)_string_of_pearls->size-1);
  size_t end           = (size_t)_prng->uniform(0, (int)_string_of_pearls->size-1);
  size_t current_start = start;
  size_t current_end   = end;
  size_t size          = 0;
  size_t nb_NC         = 0;
  size_t nb_E          = 0;
  size_t nb_TF         = 0;
  size_t nb_BS         = 0;
  size_t nb_P          = 0;
  pearl* tmp           = new pearl;
  
  /*----------------------------------------------------------------------*/
  /* 1) if start <= end                                                   */
  /*----------------------------------------------------------------------*/
  if (start <= end)
  {
    size = end-start+1;
    while (current_start < current_end)
    {
      memcpy(tmp, &_string_of_pearls->x[current_start], sizeof(pearl));
      memcpy(&_string_of_pearls->x[current_start], &_string_of_pearls->x[current_end], sizeof(pearl));
      memcpy(&_string_of_pearls->x[current_end], tmp, sizeof(pearl));
      current_start = (current_start+1+_string_of_pearls->size)%_string_of_pearls->size;
      current_end   = (current_end-1+_string_of_pearls->size)%_string_of_pearls->size;
    }
    size_t pos1 = (start-1+_string_of_pearls->size)%_string_of_pearls->size;
    size_t pos2 = (start+_string_of_pearls->size)%_string_of_pearls->size;
    //do_BLX_alpha_crossover(pos1, pos2);
    do_permutation(pos1, pos2);
    pos1 = (end+_string_of_pearls->size)%_string_of_pearls->size;
    pos2 = (end+1+_string_of_pearls->size)%_string_of_pearls->size;
    //do_BLX_alpha_crossover(pos1, pos2);
    do_permutation(pos1, pos2);
  }
  
  /*----------------------------------------------------------------------*/
  /* 2) else if start > end                                               */
  /*----------------------------------------------------------------------*/
  else if (start > end)
  {
    size = _string_of_pearls->size-start+end+1;
    while (current_start > current_end)
    {
      memcpy(tmp, &_string_of_pearls->x[current_start], sizeof(pearl));
      memcpy(&_string_of_pearls->x[current_start], &_string_of_pearls->x[current_end], sizeof(pearl));
      memcpy(&_string_of_pearls->x[current_end], tmp, sizeof(pearl));
      current_start = (current_start+1+_string_of_pearls->size)%_string_of_pearls->size;
      current_end   = (current_end-1+_string_of_pearls->size)%_string_of_pearls->size;
    }
    size_t pos1 = (end-1+_string_of_pearls->size)%_string_of_pearls->size;
    size_t pos2 = (end+_string_of_pearls->size)%_string_of_pearls->size;
    //do_BLX_alpha_crossover(pos1, pos2);
    do_permutation(pos1, pos2);
    pos1 = (start+_string_of_pearls->size)%_string_of_pearls->size;
    pos2 = (start+1+_string_of_pearls->size)%_string_of_pearls->size;
    //do_BLX_alpha_crossover(pos1, pos2);
    do_permutation(pos1, pos2);
  }
  
  /*----------------------------------------------------------------------*/
  /* 3) delete sequences                                                  */
  /*----------------------------------------------------------------------*/
  delete[] tmp;
  tmp = NULL;
  
  /*----------------------------------------------------------------------*/
  /* 4) save mutation event                                               */
  /*----------------------------------------------------------------------*/
  MutationEvent* new_event = new MutationEvent(INVERSION, start, end, 0, size, nb_NC, nb_E, nb_TF, nb_BS, nb_P);
  _replication_report->add_mutation_event(new_event);
}

/**
 * \brief    Perform horizontal gene transfer (HGT)
 * \details  The number of new random sequences introduced in the genome depends on a Poisson law of parameter 'HGT_RATE'. Pasted pearls at insertion point undergo BLX-a crossover
 * \param    void
 * \return   \e void
 */
void Genome::do_hgt( void )
{
  /*--------------------------------------*/
  /* 1) Draw the number of HGT to perform */
  /*--------------------------------------*/
  int nb_draws = _parameters->get_prng()->poisson(_parameters->get_hgt_rate());
  
  /*--------------------------------------*/
  /* 2) For each new HGT                  */
  /*--------------------------------------*/
  for (int i = 0; i < nb_draws; i++)
  {
    /* 2.1) Draw the HGT sequence */
    size_t insert = 0;
    pearl* shift  = NULL;
    size_t nb_NC  = 0;
    size_t nb_E   = 0;
    size_t nb_TF  = 0;
    size_t nb_BS  = 0;
    size_t nb_P   = 0;
    size_t size   = (size_t)_prng->uniform(HGT_MIN_SIZE, HGT_MAX_SIZE);
    pearl* hgt    = new pearl[size];
    for (size_t j = 0; j < size; j++)
    {
      draw_random_pearl(hgt[j]);
      if (hgt[j].type == NON_CODING)
      {
        nb_NC++;
      }
      else if (hgt[j].type == ENZYME)
      {
        nb_E++;
      }
      else if (hgt[j].type == TRANSCRIPTION_FACTOR)
      {
        nb_TF++;
      }
      else if (hgt[j].type == BINDING_SITE)
      {
        nb_BS++;
      }
      else if (hgt[j].type == PROMOTER)
      {
        nb_P++;
      }
    }
    
    /* 2.2) Insert the HGT sequence in the genome */
    if (_string_of_pearls->size == 0)
    {
      insert = 0;
    }
    else
    {
      insert = _prng->uniform(0, (int)_string_of_pearls->size-1);
    }
    
    shift = new pearl[_string_of_pearls->size-insert];
    memcpy(shift, &_string_of_pearls->x[insert], sizeof(pearl)*(_string_of_pearls->size-insert));
    memcpy(&_string_of_pearls->x[insert], hgt, sizeof(pearl)*size);
    memcpy(&_string_of_pearls->x[insert+size], shift, sizeof(pearl)*(_string_of_pearls->size-insert));
    _string_of_pearls->size += size;
    
    /* 2.3) Delete sequences */
    delete[] hgt;
    hgt = NULL;
    delete[] shift;
    shift = NULL;
    
    /* 2.4) perform crossover */
    size_t pos1 = (insert-1+_string_of_pearls->size)%_string_of_pearls->size;
    size_t pos2 = (insert+_string_of_pearls->size)%_string_of_pearls->size;
    //do_BLX_alpha_crossover(pos1, pos2);
    do_permutation(pos1, pos2);
    pos1 = (insert+size-1+_string_of_pearls->size)%_string_of_pearls->size;
    pos2 = (insert+size+_string_of_pearls->size)%_string_of_pearls->size;
    //do_BLX_alpha_crossover(pos1, pos2);
    do_permutation(pos1, pos2);
    
    /* 2.5) Save the mutation event */
    MutationEvent* new_event = new MutationEvent(HGT, insert, size, nb_NC, nb_E, nb_TF, nb_BS, nb_P);
    _replication_report->add_mutation_event(new_event);
  }
}

/**
 * \brief    Do a breakpoint
 * \details  A breakpoint mutate one of the two adjacent genes
 * \param    const double* mutation_rates
 * \param    size_t breakpoint
 * \return   \e void
 */
void Genome::do_breakpoint( const double* mutation_rates, size_t breakpoint )
{
  if (_string_of_pearls->size == 0)
  {
    return;
  }
  if (_prng->uniform() < 0.5)
  {
    do_pearl_mutation_at_breakpoints(mutation_rates, breakpoint);
  }
  else
  {
    size_t previous = (breakpoint-1+_string_of_pearls->size)%_string_of_pearls->size;
    do_pearl_mutation_at_breakpoints(mutation_rates, previous);
  }
}

/**
 * \brief    Mutate a pearl
 * \details  Apply point mutations on a pearl
 * \param    const double* mutation_rates
 * \param    size_t pos
 * \return   \e bool
 */
void Genome::do_pearl_mutation( const double* mutation_rates, size_t pos )
{
  MutationVector* dX           = new MutationVector();
  pearl& p                     = _string_of_pearls->x[pos];
  dX->get_dX()->type           = p.type;
  dX->get_dX()->free_activity  = p.free_activity;
  dX->get_dX()->bound_activity = p.bound_activity;
  dX->get_dX()->window         = p.window;
  bool mutate                  = dX->draw(_prng, mutation_rates);
  
  /*----------------------------------------*/
  /* if at least one point mutation occured */
  /*----------------------------------------*/
  if (mutate)
  {
    pearl* dx = dX->get_dX();
    
    /*------------------------------------------------------------------ Global attributes */
    
    if (dx->type == NON_CODING ||
        dx->type == E_TO_NC_TRANSITION ||
        dx->type == TF_TO_NC_TRANSITION ||
        dx->type == BS_TO_NC_TRANSITION ||
        dx->type == P_TO_NC_TRANSITION)
    {
      p.type = NON_CODING;
    }
    else if (dx->type == ENZYME ||
             dx->type == NC_TO_E_TRANSITION ||
             dx->type == TF_TO_E_TRANSITION ||
             dx->type == BS_TO_E_TRANSITION ||
             dx->type == P_TO_E_TRANSITION)
    {
      p.type = ENZYME;
    }
    else if (dx->type == TRANSCRIPTION_FACTOR ||
             dx->type == NC_TO_TF_TRANSITION ||
             dx->type == E_TO_TF_TRANSITION ||
             dx->type == BS_TO_TF_TRANSITION ||
             dx->type == P_TO_TF_TRANSITION)
    {
      p.type = TRANSCRIPTION_FACTOR;
    }
    else if (dx->type == BINDING_SITE ||
             dx->type == NC_TO_BS_TRANSITION ||
             dx->type == E_TO_BS_TRANSITION ||
             dx->type == TF_TO_BS_TRANSITION ||
             dx->type == P_TO_BS_TRANSITION)
    {
      p.type = BINDING_SITE;
    }
    else if (dx->type == PROMOTER ||
             dx->type == NC_TO_P_TRANSITION ||
             dx->type == E_TO_P_TRANSITION ||
             dx->type == TF_TO_P_TRANSITION ||
             dx->type == BS_TO_P_TRANSITION)
    {
      p.type = PROMOTER;
    }
    
    /*------------------------------------------------------------------ Enzyme type (E) attributes */
    
    p.s += dx->s;
    p.s  = (p.s > 0 ? p.s : 1);
    p.p += dx->p;
    p.p  = (p.p > 0 ? p.p : 1);
    
    p.km = pow(10.0, log10(p.km)+dx->km);
    if (p.km < pow(10.0, KM_MIN_LOG))
    {
      p.km = pow(10.0, KM_MIN_LOG);
    }
    else if (p.km > pow(10.0, KM_MAX_LOG))
    {
      p.km = pow(10.0, KM_MAX_LOG);
    }
    
    double kcat = fabs(p.kcat);
    kcat = pow(10.0, log10(kcat)+dx->kcat);
    if (p.kcat < 0.0)
    {
      if (kcat < pow(10.0, KCAT_MIN_LOG))
      {
        p.kcat = pow(10.0, KCAT_MIN_LOG);
      }
      else if (kcat > pow(10.0, KCAT_MAX_LOG))
      {
        p.kcat = -pow(10.0, KCAT_MAX_LOG);
      }
      else
      {
        p.kcat = -kcat;
      }
    }
    else if (p.kcat > 0.0)
    {
      if (kcat < pow(10.0, KCAT_MIN_LOG))
      {
        p.kcat = -pow(10.0, KCAT_MIN_LOG);
      }
      else if (kcat > pow(10.0, KCAT_MAX_LOG))
      {
        p.kcat = pow(10.0, KCAT_MAX_LOG);
      }
      else
      {
        p.kcat = kcat;
      }
    }
    
    /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
    
    p.BS_tag         += dx->BS_tag;
    p.coE_tag        += dx->coE_tag;
    p.coE_tag         = (p.coE_tag > 0 ? p.coE_tag : 1);
    p.free_activity   = dx->free_activity;
    p.bound_activity  = dx->bound_activity;
    p.window          = dx->window;
    
    /*------------------------------------------------------------------ Binding site type (BS) attributes */
    
    p.TF_tag += dx->TF_tag;
    
    /*------------------------------------------------------------------ Promoter type (P) attributes */
    
    p.basal_expression_level += dx->basal_expression_level;
    p.basal_expression_level  = (p.basal_expression_level > 0.0 ? p.basal_expression_level : 0.0);
    p.basal_expression_level  = (p.basal_expression_level < 1.0 ? p.basal_expression_level : 1.0);
    
    /*------------------------------------------------------------------ Save point mutation event */
    
    MutationEvent* new_event = new MutationEvent(POINT_MUTATION, pos, dX);
    _replication_report->add_mutation_event(new_event);
  }
  else
  {
    delete dX;
    dX = NULL;
  }
}

/**
 * \brief    Mutate a pearl at breakpoints
 * \details  Apply point mutations on a pearl
 * \param    const double* mutation_rates
 * \param    size_t pos
 * \return   \e void
 */
void Genome::do_pearl_mutation_at_breakpoints( const double* mutation_rates, size_t pos )
{
  MutationVector* dX           = new MutationVector();
  pearl& p                     = _string_of_pearls->x[pos];
  dX->get_dX()->type           = p.type;
  dX->get_dX()->free_activity  = p.free_activity;
  dX->get_dX()->bound_activity = p.bound_activity;
  dX->get_dX()->window         = p.window;
  bool mutate                  = dX->breakpoint_draw(_prng, mutation_rates);
  
  /*----------------------------------------*/
  /* if at least one point mutation occured */
  /*----------------------------------------*/
  if (mutate)
  {
    pearl* dx = dX->get_dX();
    
    /*------------------------------------------------------------------ Global attributes */
    
    if (dx->type == NON_CODING ||
        dx->type == E_TO_NC_TRANSITION ||
        dx->type == TF_TO_NC_TRANSITION ||
        dx->type == BS_TO_NC_TRANSITION ||
        dx->type == P_TO_NC_TRANSITION)
    {
      p.type = NON_CODING;
    }
    else if (dx->type == ENZYME ||
             dx->type == NC_TO_E_TRANSITION ||
             dx->type == TF_TO_E_TRANSITION ||
             dx->type == BS_TO_E_TRANSITION ||
             dx->type == P_TO_E_TRANSITION)
    {
      p.type = ENZYME;
    }
    else if (dx->type == TRANSCRIPTION_FACTOR ||
             dx->type == NC_TO_TF_TRANSITION ||
             dx->type == E_TO_TF_TRANSITION ||
             dx->type == BS_TO_TF_TRANSITION ||
             dx->type == P_TO_TF_TRANSITION)
    {
      p.type = TRANSCRIPTION_FACTOR;
    }
    else if (dx->type == BINDING_SITE ||
             dx->type == NC_TO_BS_TRANSITION ||
             dx->type == E_TO_BS_TRANSITION ||
             dx->type == TF_TO_BS_TRANSITION ||
             dx->type == P_TO_BS_TRANSITION)
    {
      p.type = BINDING_SITE;
    }
    else if (dx->type == PROMOTER ||
             dx->type == NC_TO_P_TRANSITION ||
             dx->type == E_TO_P_TRANSITION ||
             dx->type == TF_TO_P_TRANSITION ||
             dx->type == BS_TO_P_TRANSITION)
    {
      p.type = PROMOTER;
    }
    
    /*------------------------------------------------------------------ Enzyme type (E) attributes */
    
    p.s += dx->s;
    p.s  = (p.s > 0 ? p.s : 1);
    p.p += dx->p;
    p.p  = (p.p > 0 ? p.p : 1);
    p.s += dx->s;
    p.s  = (p.s > 0 ? p.s : 1);
    p.p += dx->p;
    p.p  = (p.p > 0 ? p.p : 1);
    
    p.km = pow(10.0, log10(p.km)+dx->km);
    if (p.km < pow(10.0, KM_MIN_LOG))
    {
      p.km = pow(10.0, KM_MIN_LOG);
    }
    else if (p.km > pow(10.0, KM_MAX_LOG))
    {
      p.km = pow(10.0, KM_MAX_LOG);
    }
    
    double kcat = fabs(p.kcat);
    kcat = pow(10.0, log10(kcat)+dx->kcat);
    if (p.kcat < 0.0)
    {
      if (kcat < pow(10.0, KCAT_MIN_LOG))
      {
        p.kcat = pow(10.0, KCAT_MIN_LOG);
      }
      else if (kcat > pow(10.0, KCAT_MAX_LOG))
      {
        p.kcat = -pow(10.0, KCAT_MAX_LOG);
      }
      else
      {
        p.kcat = -kcat;
      }
    }
    else if (p.kcat > 0.0)
    {
      if (kcat < pow(10.0, KCAT_MIN_LOG))
      {
        p.kcat = -pow(10.0, KCAT_MIN_LOG);
      }
      else if (kcat > pow(10.0, KCAT_MAX_LOG))
      {
        p.kcat = pow(10.0, KCAT_MAX_LOG);
      }
      else
      {
        p.kcat = kcat;
      }
    }
    
    /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
    
    p.BS_tag         += dx->BS_tag;
    p.coE_tag        += dx->coE_tag;
    p.coE_tag         = (p.coE_tag > 0 ? p.coE_tag : 1);
    p.free_activity   = dx->free_activity;
    p.bound_activity  = dx->bound_activity;
    p.window          = dx->window;
    
    /*------------------------------------------------------------------ Binding site type (BS) attributes */
    
    p.TF_tag += dx->TF_tag;
    
    /*------------------------------------------------------------------ Promoter type (P) attributes */
    
    p.basal_expression_level += dx->basal_expression_level;
    p.basal_expression_level  = (p.basal_expression_level > 0.0 ? p.basal_expression_level : 0.0);
    p.basal_expression_level  = (p.basal_expression_level < 1.0 ? p.basal_expression_level : 1.0);
    
    /*------------------------------------------------------------------ Save point mutation event */
    
    MutationEvent* new_event = new MutationEvent(POINT_MUTATION, pos, dX);
    _replication_report->add_mutation_event(new_event);
  }
  else
  {
    delete dX;
    dX = NULL;
  }
}

/**
 * \brief    Do a permutation
 * \details  Perform a permutation between pos1 and pos2 tuples
 * \param    size_t pos1
 * \param    size_t pos2
 * \return   \e void
 */
void Genome::do_permutation( size_t pos1, size_t pos2 )
{
  pearl& p1 = _string_of_pearls->x[pos1];
  pearl& p2 = _string_of_pearls->x[pos2];
  
  /*------------------------------------------------------------------ Global attributes */
  
  if (_prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    pearl_type tmp = p1.type;
    p1.type = p2.type;
    p2.type = tmp;
  }
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  if (_prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    int tmp = p1.s;
    p1.s = p2.s;
    p2.s = tmp;
  }
  if (_prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    int tmp = p1.p;
    p1.p = p2.p;
    p2.p = tmp;
  }
  if (_prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    double tmp = p1.km;
    p1.km = p2.km;
    p2.km = tmp;
  }
  if (_prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    double tmp = p1.kcat;
    p1.kcat = p2.kcat;
    p2.kcat = tmp;
  }
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  if (_prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    int tmp = p1.BS_tag;
    p1.BS_tag = p2.BS_tag;
    p2.BS_tag = tmp;
  }
  if (_prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    int tmp = p1.coE_tag;
    p1.coE_tag = p2.coE_tag;
    p2.coE_tag = tmp;
  }
  if (_prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    bool tmp = p1.free_activity;
    p1.free_activity = p2.free_activity;
    p2.free_activity = tmp;
  }
  if (_prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    bool tmp = p1.bound_activity;
    p1.bound_activity = p2.bound_activity;
    p2.bound_activity = tmp;
  }
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  if (_prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    int tmp = p1.TF_tag;
    p1.TF_tag = p2.TF_tag;
    p2.TF_tag = tmp;
  }
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  if (_prng->uniform() < BREAKPOINT_MUTATION_RATE)
  {
    double tmp = p1.basal_expression_level;
    p1.basal_expression_level = p2.basal_expression_level;
    p2.basal_expression_level = tmp;
  }
}

/**
 * \brief    Create a default string of pearls
 * \details  --
 * \param    void
 * \return   \e void
 */
void Genome::create_string_of_pearls( void )
{
  _string_of_pearls              = new string_of_pearls;
  _string_of_pearls->x           = new pearl[GENOME_BUFFER];
  _string_of_pearls->size        = 0;
  _string_of_pearls->buffer_size = GENOME_BUFFER;
}

/**
 * \brief    Copy a string of pearls
 * \details  --
 * \param    const string_of_pearls* model
 * \return   \e void
 */
void Genome::copy_string_of_pearls( const string_of_pearls* model )
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
void Genome::delete_string_of_pearls( void )
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
void Genome::load_string_of_pearls( gzFile backup_file )
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
void Genome::save_string_of_pearls( gzFile backup_file )
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
void Genome::load_pearl( gzFile backup_file, pearl& p )
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
  gzread( backup_file, &p.window,         sizeof(p.window) );
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzread( backup_file, &p.TF_tag, sizeof(p.TF_tag) );
  
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
void Genome::save_pearl( gzFile backup_file, pearl& p )
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
  gzwrite( backup_file, &p.window,         sizeof(p.window) );
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  gzwrite( backup_file, &p.TF_tag, sizeof(p.TF_tag) );
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  gzwrite( backup_file, &p.basal_expression_level, sizeof(p.basal_expression_level) );
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  gzwrite( backup_file, &p.functional, sizeof(p.functional) );
}

/**
 * \brief    Draw a random pearl
 * \details  --
 * \param    pearl& p
 * \return   \e void
 */
void Genome::draw_random_pearl( pearl& p )
{
  /*------------------------------------------------------------------ Global attributes */
  
  double probas[5]    = {1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0};
  p.type              = (pearl_type)_parameters->get_prng()->roulette_wheel(probas, 1.0, 5);
  p.identifier        = 0;
  p.parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  p.s    = _parameters->get_prng()->uniform(1, 100);
  p.p    = _parameters->get_prng()->uniform(1, 100);
  p.km   = pow(10.0, _parameters->get_prng()->uniform()*(KM_MAX_LOG-KM_MIN_LOG)+KM_MIN_LOG);
  p.kcat = pow(10.0, _parameters->get_prng()->uniform()*(KCAT_MAX_LOG-KCAT_MIN_LOG)+KCAT_MIN_LOG);
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  p.BS_tag         = _parameters->get_prng()->uniform(-500, 500);
  p.coE_tag        = _parameters->get_prng()->uniform(1, 100);
  p.free_activity  = (_parameters->get_prng()->uniform() < 0.5 ? true : false);
  p.bound_activity = (_parameters->get_prng()->uniform() < 0.5 ? true : false);
  p.window         = _parameters->get_initial_binding_window();
  
  /*------------------------------------------------------------------ Binding site type (BS) attributes */
  
  p.TF_tag = _parameters->get_prng()->uniform(-500, 500);
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  p.basal_expression_level = _parameters->get_prng()->uniform();
}

/**
 * \brief    Increase buffer size relatively to the new genome size
 * \details  --
 * \param    size_t new_size
 * \return   \e void
 */
void Genome::increase_buffer_size( size_t new_size )
{
  assert(new_size <= _string_of_pearls->size*2);
  _string_of_pearls->buffer_size = (new_size/GENOME_BUFFER+1)*GENOME_BUFFER;
  pearl* new_x = new pearl[_string_of_pearls->buffer_size];
  memcpy(new_x, _string_of_pearls->x, sizeof(pearl)*_string_of_pearls->size);
  delete[] _string_of_pearls->x;
  _string_of_pearls->x = new_x;
}

/**
 * \brief    Decrease buffer size relatively to the genome size
 * \details  --
 * \param    void
 * \return   \e void
 */
void Genome::decrease_buffer_size( void )
{
  if (_string_of_pearls->buffer_size > GENOME_BUFFER)
  {
    _string_of_pearls->buffer_size = (_string_of_pearls->size/GENOME_BUFFER+1)*GENOME_BUFFER;
    assert(_string_of_pearls->buffer_size >= _string_of_pearls->size);
    pearl* new_x = new pearl[_string_of_pearls->buffer_size];
    memcpy(new_x, _string_of_pearls->x, sizeof(pearl)*_string_of_pearls->size);
    delete[] _string_of_pearls->x;
    _string_of_pearls->x = new_x;
  }
}

/**
 * \brief    Check if the genome size reaches the maximum genome size
 * \details  --
 * \param    void
 * \return   \e void
 */
void Genome::check_genome_size( void )
{
  if (_string_of_pearls->size > _parameters->get_maximum_genome_size())
  {
    delete[] _string_of_pearls->x;
    _string_of_pearls->x = new pearl[GENOME_BUFFER];
    _string_of_pearls->size = 0;
    _string_of_pearls->buffer_size = GENOME_BUFFER;
  }
}
