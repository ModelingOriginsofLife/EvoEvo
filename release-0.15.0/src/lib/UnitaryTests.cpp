
/**
 * \file      UnitaryTests.cpp
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     UnitaryTests class definition
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

#include "UnitaryTests.h"


/*----------------------------
 * CONSTRUCTORS
 *----------------------------*/

/**
 * \brief    Constructor
 * \details  --
 * \param    void
 * \return   \e void
 */
UnitaryTests::UnitaryTests( Parameters* parameters )
{
  _parameters = parameters;
  initialize_mutation_rates();
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
UnitaryTests::~UnitaryTests( void )
{
  delete[] _mutation_rates;
  _mutation_rates = NULL;
}

/*----------------------------
 * PUBLIC METHODS
 *----------------------------*/

/**
 * \brief    Run unitary tests
 * \details  --
 * \param    std::string parameters_filename
 * \return   \e void
 */
void UnitaryTests::run_unitary_tests( std::string parameters_filename )
{
  printf("> Run unitary tests:\n");
  
  printf("  > Testing Parameters class ...");
  test_Parameters_class(parameters_filename);
  printf(" ok.\n");
  
  printf("  > Testing MutationVector class ...");
  test_MutationVector_class();
  printf(" ok.\n");
  
  printf("  > Testing MutationEvent class ...");
  test_MutationEvent_class();
  printf(" ok.\n");
  
  printf("  > Testing ReplicationReport class ...");
  test_ReplicationReport_class();
  printf(" ok.\n");
  
  printf("  > Testing Genome class ...");
  test_Genome_class();
  printf(" ok.\n");
  
  printf("  > Testing InheritedProteins class ...");
  test_InheritedProteins_class();
  printf(" ok.\n");
  
  printf("  > Testing SpeciesList class ...");
  test_SpeciesList_class();
  printf(" ok.\n");
  
  printf("  > Testing Prng class ...");
  test_Prng_class();
  printf(" ok.\n");
  
  printf("> End of unitary tests.\n\n");
}

/**
 * \brief    Test Parameters class
 * \details  --
 * \param    std::string filename
 * \return   \e void
 */
void UnitaryTests::test_Parameters_class( std::string filename )
{
  Parameters* parameters1 = new Parameters();
  parameters1->load_parameters_from_file(filename);
  
  /*------------------------------------------*/
  /* test Parameters backup constructor       */
  /*------------------------------------------*/
  mkdir("parameters", 0777);
  mkdir("prng", 0777);
  parameters1->save(1000);
  
  Parameters* parameters2 = new Parameters(1000);
  
  Parameters_isEqualTo(parameters2, parameters1);
  
  delete parameters1;
  parameters1 = NULL;
  
  /*------------------------------------------*/
  /* test Parameters copy constructor         */
  /*------------------------------------------*/
  modify_Parameters(parameters2);
  Parameters* parameters3 = new Parameters(*parameters2);
  
  Parameters_isEqualTo(parameters3, parameters2);
  
  delete parameters2;
  parameters2 = NULL;
  
  /*------------------------------------------*/
  /* test Parameters file writing and loading */
  /*------------------------------------------*/
  parameters3->write("./parameters/tmp.txt");
  
  Parameters* parameters4 = new Parameters();
  parameters4->load_parameters_from_file("./parameters/tmp.txt");
  
  Parameters_isEqualTo(parameters4, parameters3);
  
  system("rm -rf parameters");
  system("rm -rf prng");
  
  delete parameters3;
  parameters3 = NULL;
  delete parameters4;
  parameters4 = NULL;
}

/**
 * \brief    Test MutationVector class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_MutationVector_class( void )
{
  _parameters->get_prng()->set_seed((unsigned long int)time(NULL));
  MutationVector* dX = new MutationVector();
  
  for (size_t i = 0; i < 100; i++)
  {
    while (!dX->draw(_parameters->get_prng(), _mutation_rates));
  }
  
  /*----------------------------------------*/
  /* test struct native copy constructor    */
  /*----------------------------------------*/
  const pearl p = *dX->get_dX();
  pearl p1 = p;
  pearl_isEqualTo(&p1, dX->get_dX());
  
  /*----------------------------------------*/
  /* test MutationVector backup constructor */
  /*----------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  dX->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  MutationVector* dX2 = new MutationVector(backup_file);
  gzclose(backup_file);
  
  MutationVector_isEqualTo(dX2, dX);
  
  system("rm -rf tmp");
  
  delete dX;
  dX = NULL;
  
  /*----------------------------------------*/
  /* test MutationVector copy constructor   */
  /*----------------------------------------*/
  MutationVector* dX3 = new MutationVector(*dX2);
  
  MutationVector_isEqualTo(dX3, dX2);
  
  delete dX2;
  dX2 = NULL;
  delete dX3;
  dX3 = NULL;
}

/**
 * \brief    Test MutationEvent class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_MutationEvent_class( void )
{
  _parameters->get_prng()->set_seed((unsigned long int)time(NULL));
  MutationVector* dX = new MutationVector();
  for (size_t i = 0; i < 100; i++)
  {
    while (!dX->draw(_parameters->get_prng(), _mutation_rates));
  }
  MutationEvent* event1 = new MutationEvent(POINT_MUTATION, 100, dX);
  
  /*---------------------------------------*/
  /* test MutationEvent backup constructor */
  /*---------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  event1->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  MutationEvent* event2 = new MutationEvent(backup_file);
  gzclose(backup_file);
  
  MutationEvent_isEqualTo(event2, event1);
  
  system("rm -rf tmp");
  delete event1;
  event1 = NULL;
  
  /*---------------------------------------*/
  /* test MutationEvent copy constructor   */
  /*---------------------------------------*/
  MutationEvent* event3 = new MutationEvent(*event2);
  
  MutationEvent_isEqualTo(event3, event2);
  
  delete event2;
  event2 = NULL;
  delete event3;
  event3 = NULL;
}

/**
 * \brief    Test ReplicationReport class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_ReplicationReport_class( void )
{
  /* 1) create a mutation event ------*/
  _parameters->get_prng()->set_seed((unsigned long int)time(NULL));
  MutationVector* dX = new MutationVector();
  for (size_t i = 0; i < 100; i++)
  {
    while (!dX->draw(_parameters->get_prng(), _mutation_rates));
  }
  MutationEvent* point_mutation_event = new MutationEvent(POINT_MUTATION, 100, dX);
  MutationEvent* hgt_event            = new MutationEvent(HGT, 50, 25, 5, 5, 5, 5, 5);
  MutationEvent* duplication_event    = new MutationEvent(DUPLICATION, 100, 200, 300, 250, 50, 50, 50, 50, 50);
  MutationEvent* deletion_event       = new MutationEvent(DELETION, 100, 200, 300, 250, 50, 50, 50, 50, 50);
  MutationEvent* translocation_event  = new MutationEvent(TRANSLOCATION, 100, 200, 300, 250, 50, 50, 50, 50, 50);
  MutationEvent* inversion_event      = new MutationEvent(INVERSION, 100, 200, 300, 250, 50, 50, 50, 50, 50);
  
  /* 2) create the replication report */
  ReplicationReport* report1 = new ReplicationReport();
  report1->add_mutation_event(point_mutation_event);
  report1->add_mutation_event(hgt_event);
  report1->add_mutation_event(duplication_event);
  report1->add_mutation_event(deletion_event);
  report1->add_mutation_event(translocation_event);
  report1->add_mutation_event(inversion_event);
  report1->compute_mean();
  
  /*-------------------------------------------*/
  /* test ReplicationReport atributes          */
  /*-------------------------------------------*/
  
  /*------------------------------------------------------------------ global data */
  
  assert(report1->get_number_of_events() == 6);
  
  /*------------------------------------------------------------------ point mutations data */
  
  assert(report1->get_nb_point_mutations() == 1);
  
  size_t nb_NC_point_mutations = 0;
  size_t nb_E_point_mutations  = 0;
  size_t nb_TF_point_mutations = 0;
  size_t nb_BS_point_mutations = 0;
  size_t nb_P_point_mutations  = 0;
  
  size_t nb_NC_to_E_transitions  = 0;
  size_t nb_NC_to_TF_transitions = 0;
  size_t nb_NC_to_BS_transitions = 0;
  size_t nb_NC_to_P_transitions  = 0;
  
  size_t nb_E_to_NC_transitions = 0;
  size_t nb_E_to_TF_transitions = 0;
  size_t nb_E_to_BS_transitions = 0;
  size_t nb_E_to_P_transitions  = 0;
  
  size_t nb_TF_to_NC_transitions = 0;
  size_t nb_TF_to_E_transitions  = 0;
  size_t nb_TF_to_BS_transitions = 0;
  size_t nb_TF_to_P_transitions  = 0;
  
  size_t nb_BS_to_NC_transitions = 0;
  size_t nb_BS_to_E_transitions  = 0;
  size_t nb_BS_to_TF_transitions = 0;
  size_t nb_BS_to_P_transitions  = 0;
  
  size_t nb_P_to_NC_transitions = 0;
  size_t nb_P_to_E_transitions  = 0;
  size_t nb_P_to_TF_transitions = 0;
  size_t nb_P_to_BS_transitions = 0;
  
  pearl_type ptype = point_mutation_event->get_mutation_vector()->get_dX()->type;
  switch (ptype)
  {
    case NON_CODING:
      nb_NC_point_mutations = 1;
      break;
    case ENZYME:
      nb_E_point_mutations = 1;
      break;
    case TRANSCRIPTION_FACTOR:
      nb_TF_point_mutations = 1;
      break;
    case BINDING_SITE:
      nb_BS_point_mutations = 1;
      break;
    case PROMOTER:
      nb_P_point_mutations = 1;
      break;
    case NC_TO_E_TRANSITION:
      nb_NC_to_E_transitions = 1;
      break;
    case NC_TO_TF_TRANSITION:
      nb_NC_to_TF_transitions = 1;
      break;
    case NC_TO_BS_TRANSITION:
      nb_NC_to_BS_transitions = 1;
      break;
    case NC_TO_P_TRANSITION:
      nb_NC_to_P_transitions = 1;
      break;
    case E_TO_NC_TRANSITION:
      nb_E_to_NC_transitions = 1;
      break;
    case E_TO_TF_TRANSITION:
      nb_E_to_TF_transitions = 1;
      break;
    case E_TO_BS_TRANSITION:
      nb_E_to_BS_transitions = 1;
      break;
    case E_TO_P_TRANSITION:
      nb_E_to_P_transitions = 1;
      break;
    case TF_TO_NC_TRANSITION:
      nb_TF_to_NC_transitions = 1;
      break;
    case TF_TO_E_TRANSITION:
      nb_TF_to_E_transitions = 1;
      break;
    case TF_TO_BS_TRANSITION:
      nb_TF_to_BS_transitions = 1;
      break;
    case TF_TO_P_TRANSITION:
      nb_TF_to_P_transitions = 1;
      break;
    case BS_TO_NC_TRANSITION:
      nb_BS_to_NC_transitions = 1;
      break;
    case BS_TO_E_TRANSITION:
      nb_BS_to_E_transitions = 1;
      break;
    case BS_TO_TF_TRANSITION:
      nb_BS_to_TF_transitions = 1;
      break;
    case BS_TO_P_TRANSITION:
      nb_BS_to_P_transitions = 1;
      break;
    case P_TO_NC_TRANSITION:
      nb_P_to_NC_transitions = 1;
      break;
    case P_TO_E_TRANSITION:
      nb_P_to_E_transitions = 1;
      break;
    case P_TO_TF_TRANSITION:
      nb_P_to_TF_transitions = 1;
      break;
    case P_TO_BS_TRANSITION:
      nb_P_to_BS_transitions = 1;
      break;
    default:
      break;
  }
  
  assert(report1->get_nb_NC_point_mutations() == nb_NC_point_mutations);
  assert(report1->get_nb_E_point_mutations() == nb_E_point_mutations);
  assert(report1->get_nb_TF_point_mutations() == nb_TF_point_mutations);
  assert(report1->get_nb_BS_point_mutations() == nb_BS_point_mutations);
  assert(report1->get_nb_P_point_mutations() == nb_P_point_mutations);
  
  assert(report1->get_nb_NC_to_E_transitions() == nb_NC_to_E_transitions);
  assert(report1->get_nb_NC_to_TF_transitions() == nb_NC_to_TF_transitions);
  assert(report1->get_nb_NC_to_BS_transitions() == nb_NC_to_BS_transitions);
  assert(report1->get_nb_NC_to_P_transitions() == nb_NC_to_P_transitions);
  
  assert(report1->get_nb_E_to_NC_transitions() == nb_E_to_NC_transitions);
  assert(report1->get_nb_E_to_TF_transitions() == nb_E_to_TF_transitions);
  assert(report1->get_nb_E_to_BS_transitions() == nb_E_to_BS_transitions);
  assert(report1->get_nb_E_to_P_transitions() == nb_E_to_P_transitions);
  
  assert(report1->get_nb_TF_to_NC_transitions() == nb_TF_to_NC_transitions);
  assert(report1->get_nb_TF_to_E_transitions() == nb_TF_to_E_transitions);
  assert(report1->get_nb_TF_to_BS_transitions() == nb_TF_to_BS_transitions);
  assert(report1->get_nb_TF_to_P_transitions() == nb_TF_to_P_transitions);
  
  assert(report1->get_nb_BS_to_NC_transitions() == nb_BS_to_NC_transitions);
  assert(report1->get_nb_BS_to_E_transitions() == nb_BS_to_E_transitions);
  assert(report1->get_nb_BS_to_TF_transitions() == nb_BS_to_TF_transitions);
  assert(report1->get_nb_BS_to_P_transitions() == nb_BS_to_P_transitions);
  
  assert(report1->get_nb_P_to_NC_transitions() == nb_P_to_NC_transitions);
  assert(report1->get_nb_P_to_E_transitions() == nb_P_to_E_transitions);
  assert(report1->get_nb_P_to_TF_transitions() == nb_P_to_TF_transitions);
  assert(report1->get_nb_P_to_BS_transitions() == nb_P_to_BS_transitions);
  
  if (point_mutation_event->get_mutation_vector()->get_dX()->type == ENZYME)
  {
    assert(report1->get_mean_s_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->s);
    assert(report1->get_mean_p_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->p);
    assert(report1->get_mean_km_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->km);
    assert(report1->get_mean_kcat_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->kcat);
    assert(report1->get_mean_BS_tag_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->BS_tag);
    assert(report1->get_mean_coE_tag_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->coE_tag);
    assert(report1->get_mean_TF_tag_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->TF_tag);
    assert(report1->get_mean_basal_expression_level_mutation_size() == point_mutation_event->get_mutation_vector()->get_dX()->basal_expression_level);
  }
  
  /*------------------------------------------------------------------ HGT data */
  
  assert(report1->get_nb_HGT() == 1);
  assert(report1->get_mean_HGT_size() == 25.0);
  assert(report1->get_nb_NC_HGT() == 5);
  assert(report1->get_nb_E_HGT() == 5);
  assert(report1->get_nb_TF_HGT() == 5);
  assert(report1->get_nb_BS_HGT() == 5);
  assert(report1->get_nb_P_HGT() == 5);
  
  /*------------------------------------------------------------------ rearrangements data */
  
  assert(report1->get_nb_rearrangements() == 4);
  assert(report1->get_nb_duplicated_NC() == 50);
  assert(report1->get_nb_duplicated_E() == 50);
  assert(report1->get_nb_duplicated_TF() == 50);
  assert(report1->get_nb_duplicated_BS() == 50);
  assert(report1->get_nb_duplicated_P() == 50);
  assert(report1->get_nb_deleted_NC() == 50);
  assert(report1->get_nb_deleted_E() == 50);
  assert(report1->get_nb_deleted_TF() == 50);
  assert(report1->get_nb_deleted_BS() == 50);
  assert(report1->get_nb_deleted_P() == 50);
  assert(report1->get_nb_duplications() == 1);
  assert(report1->get_nb_deletions() == 1);
  assert(report1->get_nb_translocations() == 1);
  assert(report1->get_nb_inversions() == 1);
  assert(report1->get_mean_rearrangement_size() == 250.0);
  assert(report1->get_mean_duplication_size() == 250.0);
  assert(report1->get_mean_deletion_size() == 250.0);
  assert(report1->get_mean_translocation_size() == 250.0);
  assert(report1->get_mean_inversion_size() == 250.0);
  
  /*-------------------------------------------*/
  /* test ReplicationReport backup constructor */
  /*-------------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  report1->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  ReplicationReport* report2 = new ReplicationReport(backup_file);
  gzclose(backup_file);
  
  ReplicationReport_isEqualTo(report2, report1);
  
  system("rm -rf tmp");
  delete report1;
  report1 = NULL;
  
  /*-------------------------------------------*/
  /* test ReplicationReport copy constructor   */
  /*-------------------------------------------*/
  ReplicationReport* report3 = new ReplicationReport(*report2);
  
  ReplicationReport_isEqualTo(report3, report2);
  
  delete report2;
  report2 = NULL;
  delete report3;
  report3 = NULL;
}

/**
 * \brief    Test Genome class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_Genome_class( void )
{
  /* 1) build the genome */
  ReplicationReport* rep_report1 = new ReplicationReport();
  ReplicationReport* rep_report2 = new ReplicationReport();
  
  Genome* genome1 = new Genome(_parameters, rep_report1);
  
  Parameters* param_copy  = new Parameters(*_parameters);
  Genome*     genome_copy = new Genome(param_copy, rep_report2);
  
  pearl* p        = new pearl;
  initialize_pearl(p);
  bool suitable   = false;
  while (!suitable)
  {
    for (size_t i = 0; i < 1000; i++)
    {
      genome1->add_pearl(p);
      genome_copy->add_pearl(p);
    }
    for (size_t i = 0; i < 10; i++)
    {
      genome1->mutate(_mutation_rates);
      genome_copy->mutate(_mutation_rates);
    }
    if (genome1->get_size() > 0)
    {
      suitable = true;
    }
  }
  delete p;
  p = NULL;
  
  /* 2) test random driven equality */
  ReplicationReport_isEqualTo(rep_report1, rep_report2);
  Genome_isEqualTo(genome_copy, genome1);
  delete genome_copy;
  genome_copy = NULL;
  delete param_copy;
  param_copy = NULL;
  
  /*---------------------------------------*/
  /* test Genome backup constructor        */
  /*---------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  genome1->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  Genome* genome2 = new Genome(_parameters, rep_report1, backup_file);
  gzclose(backup_file);
  
  Genome_isEqualTo(genome2, genome1);
  
  system("rm -rf tmp");
  delete genome1;
  genome1 = NULL;
  
  /*---------------------------------------*/
  /* test Genome copy constructor          */
  /*---------------------------------------*/
  Genome* genome3 = new Genome(*genome2, rep_report1);
  
  Genome_isEqualTo(genome3, genome2);
  
  delete genome2;
  genome2 = NULL;
  delete genome3;
  genome3 = NULL;
}

/**
 * \brief    Test InheritedProteins class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_InheritedProteins_class( void )
{
  /* 1) build the inherited proteins class */
  InheritedProteins* inhprot1 = new InheritedProteins(_parameters);
  
  pearl* p = new pearl;
  initialize_pearl(p);
  for (size_t i = 0; i < 1000; i++)
  {
    if (_parameters->get_prng()->uniform() < 0.5)
    {
      p->type = ENZYME;
    }
    else
    {
      p->type = TRANSCRIPTION_FACTOR;
    }
    inhprot1->add_pearl(p);
  }
  delete p;
  p = NULL;
  inhprot1->initialize_concentration_vector();
  inhprot1->build_index_list();
  
  /*---------------------------------------*/
  /* test Genome backup constructor        */
  /*---------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  inhprot1->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  InheritedProteins* inhprot2 = new InheritedProteins(_parameters, backup_file);
  gzclose(backup_file);
  
  InheritedProteins_isEqualTo(inhprot2, inhprot1);
  
  system("rm -rf tmp");
  delete inhprot1;
  inhprot1 = NULL;
  
  /*---------------------------------------*/
  /* test Genome copy constructor          */
  /*---------------------------------------*/
  InheritedProteins* inhprot3 = new InheritedProteins(*inhprot2);
  
  InheritedProteins_isEqualTo(inhprot3, inhprot2);
  
  delete inhprot2;
  inhprot2 = NULL;
  delete inhprot3;
  inhprot3 = NULL;
}

/**
 * \brief    Test SpeciesList class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_SpeciesList_class( void )
{
  SpeciesList* spl1 = new SpeciesList();
  spl1->set(1, 10.0);
  spl1->set(2, 20.0);
  spl1->set(3, 30.0);
  spl1->set(4, 40.0);
  spl1->set(5, 50.0);
  
  /*---------------------------------------*/
  /* test SpeciesList backup constructor   */
  /*---------------------------------------*/
  gzFile backup_file = gzopen("tmp", "w");
  spl1->save(backup_file);
  gzclose(backup_file);
  
  backup_file = gzopen("tmp", "r");
  SpeciesList* spl2 = new SpeciesList(backup_file);
  gzclose(backup_file);
  
  SpeciesList_isEqualTo(spl1, spl2);
  
  system("rm -rf tmp");
  
  delete spl1;
  spl1 = NULL;
  
  /*---------------------------------------*/
  /* test SpeciesList copy constructor     */
  /*---------------------------------------*/
  SpeciesList* spl3 = new SpeciesList(*spl2);
  SpeciesList_isEqualTo(spl2, spl3);
  delete spl2;
  spl2 = NULL;
  delete spl3;
  spl3 = NULL;
}

/**
 * \brief    Test Prng class
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::test_Prng_class( void )
{
  mkdir("prng", 0777);
  unsigned long int seed = (unsigned long int)time(NULL);
  Prng* prng1 = new Prng();
  prng1->set_seed(seed);
  
  /*---------------------------------------*/
  /* test Prng seed constructor            */
  /*---------------------------------------*/
  Prng* prng2 = new Prng();
  prng2->set_seed(seed);
  Prng_isEqualTo(prng1, prng2);
  delete prng1;
  prng1 = NULL;
  
  /*---------------------------------------*/
  /* test Prng backup constructor          */
  /*---------------------------------------*/
  prng2->save(10000);
  Prng* prng3 = new Prng(10000);
  Prng_isEqualTo(prng2, prng3);
  delete prng2;
  prng2 = NULL;
  
  /*---------------------------------------*/
  /* test Prng copy constructor            */
  /*---------------------------------------*/
  Prng* prng4 = new Prng(*prng3);
  Prng_isEqualTo(prng3, prng4);
  delete prng3;
  prng3 = NULL;
  delete prng4;
  prng4 = NULL;
  system("rm -rf prng");
}

/*----------------------------
 * PROTECTED METHODS
 *----------------------------*/

/**
 * \brief    Initialize the mutation rates vector
 * \details  --
 * \param    void
 * \return   \e void
 */
void UnitaryTests::initialize_mutation_rates( void )
{
  _mutation_rates = new double[NUMBER_OF_MUTATION_RATES];
  _mutation_rates[POINT_MUTATION_RATE]                    = _parameters->get_point_mutation_rate();
  _mutation_rates[DUPLICATION_RATE]                       = _parameters->get_duplication_rate();
  _mutation_rates[DELETION_RATE]                          = _parameters->get_deletion_rate();
  _mutation_rates[TRANSLOCATION_RATE]                     = _parameters->get_translocation_rate();
  _mutation_rates[INVERSION_RATE]                         = _parameters->get_inversion_rate();
  _mutation_rates[TRANSITION_RATE]                        = _parameters->get_transition_rate();
  _mutation_rates[SUBSTRATE_TAG_MUTATION_SIZE]            = _parameters->get_substrate_tag_mutation_size();
  _mutation_rates[PRODUCT_TAG_MUTATION_SIZE]              = _parameters->get_product_tag_mutation_size();
  _mutation_rates[KM_MUTATION_SIZE]                       = _parameters->get_km_mutation_size();
  _mutation_rates[KCAT_MUTATION_SIZE]                     = _parameters->get_kcat_mutation_size();
  _mutation_rates[BINDING_SIZE_TAG_MUTATION_SIZE]         = _parameters->get_binding_size_tag_mutation_size();
  _mutation_rates[CO_ENZYME_TAG_MUTATION_SIZE]            = _parameters->get_co_enzyme_tag_mutation_size();
  _mutation_rates[TRANSCRIPTION_FACTOR_TAG_MUTATION_SIZE] = _parameters->get_transcription_factor_tag_mutation_size();
  _mutation_rates[BASAL_EXPRESSION_LEVEL_MUTATION_SIZE]   = _parameters->get_basal_expression_level_mutation_size();
}

/**
 * \brief    Initialize a default pearl
 * \details  --
 * \param    pearl* gene
 * \return   \e void
 */
void UnitaryTests::initialize_pearl( pearl* obj )
{
  /*------------------------------------------------------------------ Global attributes */
  
  obj->type              = NON_CODING;
  obj->identifier        = 0;
  obj->parent_identifier = 0;
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  obj->s    = 1;
  obj->p    = 1;
  obj->km   = pow(10.0, KM_MIN_LOG);
  obj->kcat = pow(10.0, KCAT_MIN_LOG);
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  obj->BS_tag         = 0;
  obj->coE_tag        = 1;
  obj->free_activity  = false;
  obj->bound_activity = false;
  obj->window         = 0;
  
  /*------------------------------------------------------------------ Bidnign site type (BS) attributes */
  
  obj->TF_tag = 0;
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  obj->basal_expression_level = 0.0;
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  obj->functional = false;
}

/**
 * \brief    Test Parameters struct equality
 * \details  --
 * \param    Parameters* obj1
 * \param    Parameters* obj2
 * \return   \e void
 */
void UnitaryTests::Parameters_isEqualTo( Parameters* obj1, Parameters* obj2 )
{
  /*------------------------------------------------------------------ prng parameters */
  
  assert(obj1->get_seed() == obj2->get_seed());
  
  /*------------------------------------------------------------------ parallel computing */
  
  assert(obj1->get_parallel_computing() == obj2->get_parallel_computing());
  
  /*------------------------------------------------------------------ schemes parameters */
  
  assert(obj1->get_energy_constraints() == obj2->get_energy_constraints());
  assert(obj1->get_membrane_permeability() == obj2->get_membrane_permeability());
  assert(obj1->get_metabolic_inheritance() == obj2->get_metabolic_inheritance());
  assert(obj1->get_enzymatic_inheritance() == obj2->get_enzymatic_inheritance());
  assert(obj1->get_co_enzyme_activity() == obj2->get_co_enzyme_activity());
  assert(obj1->get_score_scheme() == obj2->get_score_scheme());
  assert(obj1->get_selection_threshold() == obj2->get_selection_threshold());
  
  /*------------------------------------------------------------------ space parameters */
  
  assert(obj1->get_width() == obj2->get_width());
  assert(obj1->get_height() == obj2->get_height());
  
  /*------------------------------------------------------------------ output parameters */
  
  assert(obj1->get_experiment_backup_step() == obj2->get_experiment_backup_step());
  assert(obj1->get_tree_backup_step() == obj2->get_tree_backup_step());
  assert(obj1->get_figures_generation_step() == obj2->get_figures_generation_step());
  
  /*------------------------------------------------------------------ genomic parameters */
  
  assert(obj1->get_metabolite_tag_initial_range()->law == obj2->get_metabolite_tag_initial_range()->law);
  assert(obj1->get_metabolite_tag_initial_range()->min == obj2->get_metabolite_tag_initial_range()->min);
  assert(obj1->get_metabolite_tag_initial_range()->max == obj2->get_metabolite_tag_initial_range()->max);
  assert(obj1->get_metabolite_tag_initial_range()->mu == obj2->get_metabolite_tag_initial_range()->mu);
  assert(obj1->get_metabolite_tag_initial_range()->sigma == obj2->get_metabolite_tag_initial_range()->sigma);
  assert(obj1->get_metabolite_tag_initial_range()->lambda == obj2->get_metabolite_tag_initial_range()->lambda);
  
  assert(obj1->get_binding_site_tag_initial_range()->law == obj2->get_binding_site_tag_initial_range()->law);
  assert(obj1->get_binding_site_tag_initial_range()->min == obj2->get_binding_site_tag_initial_range()->min);
  assert(obj1->get_binding_site_tag_initial_range()->max == obj2->get_binding_site_tag_initial_range()->max);
  assert(obj1->get_binding_site_tag_initial_range()->mu == obj2->get_binding_site_tag_initial_range()->mu);
  assert(obj1->get_binding_site_tag_initial_range()->sigma == obj2->get_binding_site_tag_initial_range()->sigma);
  assert(obj1->get_binding_site_tag_initial_range()->lambda == obj2->get_binding_site_tag_initial_range()->lambda);
  
  assert(obj1->get_co_enzyme_tag_initial_range()->law == obj2->get_co_enzyme_tag_initial_range()->law);
  assert(obj1->get_co_enzyme_tag_initial_range()->min == obj2->get_co_enzyme_tag_initial_range()->min);
  assert(obj1->get_co_enzyme_tag_initial_range()->max == obj2->get_co_enzyme_tag_initial_range()->max);
  assert(obj1->get_co_enzyme_tag_initial_range()->mu == obj2->get_co_enzyme_tag_initial_range()->mu);
  assert(obj1->get_co_enzyme_tag_initial_range()->sigma == obj2->get_co_enzyme_tag_initial_range()->sigma);
  assert(obj1->get_co_enzyme_tag_initial_range()->lambda == obj2->get_co_enzyme_tag_initial_range()->lambda);
  
  assert(obj1->get_transcription_factor_tag_initial_range()->law == obj2->get_transcription_factor_tag_initial_range()->law);
  assert(obj1->get_transcription_factor_tag_initial_range()->min == obj2->get_transcription_factor_tag_initial_range()->min);
  assert(obj1->get_transcription_factor_tag_initial_range()->max == obj2->get_transcription_factor_tag_initial_range()->max);
  assert(obj1->get_transcription_factor_tag_initial_range()->mu == obj2->get_transcription_factor_tag_initial_range()->mu);
  assert(obj1->get_transcription_factor_tag_initial_range()->sigma == obj2->get_transcription_factor_tag_initial_range()->sigma);
  assert(obj1->get_transcription_factor_tag_initial_range()->lambda == obj2->get_transcription_factor_tag_initial_range()->lambda);
  
  assert(obj1->get_initial_binding_window() == obj2->get_initial_binding_window());
  assert(obj1->get_initial_number_of_NC_pearls() == obj2->get_initial_number_of_NC_pearls());
  assert(obj1->get_initial_number_of_E_pearls() == obj2->get_initial_number_of_E_pearls());
  assert(obj1->get_initial_number_of_TF_pearls() == obj2->get_initial_number_of_TF_pearls());
  assert(obj1->get_initial_number_of_BS_pearls() == obj2->get_initial_number_of_BS_pearls());
  assert(obj1->get_initial_number_of_P_pearls() == obj2->get_initial_number_of_P_pearls());
  assert(obj1->get_point_mutation_rate() == obj2->get_point_mutation_rate());
  assert(obj1->get_duplication_rate() == obj2->get_duplication_rate());
  assert(obj1->get_deletion_rate() == obj2->get_deletion_rate());
  assert(obj1->get_translocation_rate() == obj2->get_translocation_rate());
  assert(obj1->get_inversion_rate() == obj2->get_inversion_rate());
  assert(obj1->get_transition_rate() == obj2->get_transition_rate());
  assert(obj1->get_substrate_tag_mutation_size() == obj2->get_substrate_tag_mutation_size());
  assert(obj1->get_product_tag_mutation_size() == obj2->get_product_tag_mutation_size());
  assert(obj1->get_km_mutation_size() == obj2->get_km_mutation_size());
  assert(obj1->get_kcat_mutation_size() == obj2->get_kcat_mutation_size());
  assert(obj1->get_binding_size_tag_mutation_size() == obj2->get_binding_size_tag_mutation_size());
  assert(obj1->get_co_enzyme_tag_mutation_size() == obj2->get_co_enzyme_tag_mutation_size());
  assert(obj1->get_basal_expression_level_mutation_size() == obj2->get_basal_expression_level_mutation_size());
  assert(obj1->get_maximum_genome_size() == obj2->get_maximum_genome_size());
  
  /*------------------------------------------------------------------ genetic regulation network level parameters */
  
  assert(obj1->get_genetic_regulation_network_timestep() == obj2->get_genetic_regulation_network_timestep());
  assert(obj1->get_hill_function_theta() == obj2->get_hill_function_theta());
  assert(obj1->get_hill_function_n() == obj2->get_hill_function_n());
  assert(obj1->get_protein_degradation_rate() == obj2->get_protein_degradation_rate());
  
  /*------------------------------------------------------------------ metabolic parameters */
  
  assert(obj1->get_metabolism_timestep() == obj2->get_metabolism_timestep());
  assert(obj1->get_essential_metabolites_toxicity_threshold() == obj2->get_essential_metabolites_toxicity_threshold());
  assert(obj1->get_non_essential_metabolites_toxicity_threshold() == obj2->get_non_essential_metabolites_toxicity_threshold());
  assert(obj1->get_initial_metabolites_amount_in_cells() == obj2->get_initial_metabolites_amount_in_cells());
  assert(obj1->get_energy_reaction_cost_factor() == obj2->get_energy_reaction_cost_factor());
  assert(obj1->get_energy_pumping_cost() == obj2->get_energy_pumping_cost());
  assert(obj1->get_energy_transcription_cost() == obj2->get_energy_transcription_cost());
  assert(obj1->get_energy_toxicity_threshold() == obj2->get_energy_toxicity_threshold());
  
  /*------------------------------------------------------------------ cell parameters */
  
  assert(obj1->get_death_probability() == obj2->get_death_probability());
  assert(obj1->get_variable_lifespan() == obj2->get_variable_lifespan());
  assert(obj1->get_gompert_law()->a == obj2->get_gompert_law()->a);
  assert(obj1->get_gompert_law()->b == obj2->get_gompert_law()->b);
  assert(obj1->get_gompert_law()->c == obj2->get_gompert_law()->c);
  assert(obj1->get_initial_membrane_permeability() == obj2->get_initial_membrane_permeability());
  
  /*------------------------------------------------------------------ population parameters */
  
  assert(obj1->get_migration_rate() == obj2->get_migration_rate());
  assert(obj1->get_hgt_rate() == obj2->get_hgt_rate());
  
  /*------------------------------------------------------------------ environmental parameters */
  
  assert(obj1->get_environment_properties()->number_of_init_cycles == obj2->get_environment_properties()->number_of_init_cycles);
  
  assert(obj1->get_environment_properties()->species_tag_range.law == obj2->get_environment_properties()->species_tag_range.law);
  assert(obj1->get_environment_properties()->species_tag_range.min == obj2->get_environment_properties()->species_tag_range.min);
  assert(obj1->get_environment_properties()->species_tag_range.max == obj2->get_environment_properties()->species_tag_range.max);
  assert(obj1->get_environment_properties()->species_tag_range.mu == obj2->get_environment_properties()->species_tag_range.mu);
  assert(obj1->get_environment_properties()->species_tag_range.sigma == obj2->get_environment_properties()->species_tag_range.sigma);
  assert(obj1->get_environment_properties()->species_tag_range.lambda == obj2->get_environment_properties()->species_tag_range.lambda);
  
  assert(obj1->get_environment_properties()->concentration_range.law == obj2->get_environment_properties()->concentration_range.law);
  assert(obj1->get_environment_properties()->concentration_range.min == obj2->get_environment_properties()->concentration_range.min);
  assert(obj1->get_environment_properties()->concentration_range.max == obj2->get_environment_properties()->concentration_range.max);
  assert(obj1->get_environment_properties()->concentration_range.mu == obj2->get_environment_properties()->concentration_range.mu);
  assert(obj1->get_environment_properties()->concentration_range.sigma == obj2->get_environment_properties()->concentration_range.sigma);
  assert(obj1->get_environment_properties()->concentration_range.lambda == obj2->get_environment_properties()->concentration_range.lambda);
  
  assert(obj1->get_environment_properties()->number_of_species_range.law == obj2->get_environment_properties()->number_of_species_range.law);
  assert(obj1->get_environment_properties()->number_of_species_range.min == obj2->get_environment_properties()->number_of_species_range.min);
  assert(obj1->get_environment_properties()->number_of_species_range.max == obj2->get_environment_properties()->number_of_species_range.max);
  assert(obj1->get_environment_properties()->number_of_species_range.mu == obj2->get_environment_properties()->number_of_species_range.mu);
  assert(obj1->get_environment_properties()->number_of_species_range.sigma == obj2->get_environment_properties()->number_of_species_range.sigma);
  assert(obj1->get_environment_properties()->number_of_species_range.lambda == obj2->get_environment_properties()->number_of_species_range.lambda);
  
  assert(obj1->get_environment_properties()->interaction_scheme == obj2->get_environment_properties()->interaction_scheme);
  assert(obj1->get_environment_properties()->renewal_scheme == obj2->get_environment_properties()->renewal_scheme);
  assert(obj1->get_environment_properties()->variation_scheme == obj2->get_environment_properties()->variation_scheme);
  assert(obj1->get_environment_properties()->variation_localization == obj2->get_environment_properties()->variation_localization);
  
  assert(obj1->get_environment_properties()->introduction_rate == obj2->get_environment_properties()->introduction_rate);
  assert(obj1->get_environment_properties()->diffusion_rate == obj2->get_environment_properties()->diffusion_rate);
  assert(obj1->get_environment_properties()->degradation_rate == obj2->get_environment_properties()->degradation_rate);
  
  /*------------------------------------------------------------------ prime numbers */
  
  for (size_t i = 0; i < 5000; i++)
  {
    assert(obj1->get_prime_numbers()[i] == obj2->get_prime_numbers()[i]);
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test pearl struct equality
 * \details  --
 * \param    pearl* obj1
 * \param    pearl* obj2
 * \return   \e void
 */
void UnitaryTests::pearl_isEqualTo( pearl* obj1, pearl* obj2 )
{
  /*------------------------------------------------------------------ Global attributes */
  
  assert(obj1->type == obj2->type);
  assert(obj1->identifier == obj2->identifier);
  assert(obj1->parent_identifier == obj2->parent_identifier);
  
  /*------------------------------------------------------------------ Enzyme type (E) attributes */
  
  assert(obj1->s == obj2->s);
  assert(obj1->p == obj2->p);
  assert(obj1->km == obj2->km);
  assert(obj1->kcat == obj2->kcat);
  
  /*------------------------------------------------------------------ Transcription factor type (TF) attributes */
  
  assert(obj1->BS_tag == obj2->BS_tag);
  assert(obj1->coE_tag == obj2->coE_tag);
  assert(obj1->free_activity == obj2->free_activity);
  assert(obj1->bound_activity == obj2->bound_activity);
  assert(obj1->window == obj2->window);
  
  /*------------------------------------------------------------------ Bidnign site type (BS) attributes */
  
  assert(obj1->TF_tag == obj2->TF_tag);
  
  /*------------------------------------------------------------------ Promoter type (P) attributes */
  
  assert(obj1->basal_expression_level == obj2->basal_expression_level);
  
  /*------------------------------------------------------------------ Functionality attribute */
  
  assert(obj1->functional == obj2->functional);
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test MutationVector equality
 * \details  --
 * \param    MutationVector* obj1
 * \param    MutationVector* obj2
 * \return   \e void
 */
void UnitaryTests::MutationVector_isEqualTo( MutationVector* obj1, MutationVector* obj2 )
{
  pearl_isEqualTo(obj1->get_dX(), obj2->get_dX());
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test MutationEvent equality
 * \details  --
 * \param    MutationEvent* obj1
 * \param    MutationEvent* obj2
 * \return   \e void
 */
void UnitaryTests::MutationEvent_isEqualTo( MutationEvent* obj1, MutationEvent* obj2 )
{
  assert(obj1->get_mutation_type() == obj2->get_mutation_type());
  assert(obj1->get_point_mutation_location() == obj2->get_point_mutation_location());
  assert(obj1->get_hgt_insert() == obj2->get_hgt_insert());
  assert(obj1->get_nb_NC() == obj2->get_nb_NC());
  assert(obj1->get_nb_E() == obj2->get_nb_E());
  assert(obj1->get_nb_TF() == obj2->get_nb_TF());
  assert(obj1->get_nb_BS() == obj2->get_nb_BS());
  assert(obj1->get_nb_P() == obj2->get_nb_P());
  assert(obj1->get_src_breakpoint1() == obj2->get_src_breakpoint1());
  assert(obj1->get_src_breakpoint2() == obj2->get_src_breakpoint2());
  assert(obj1->get_tgt_breakpoint() == obj2->get_tgt_breakpoint());
  assert(obj1->get_size() == obj2->get_size());
  if (obj1->get_mutation_type() != POINT_MUTATION)
  {
    assert(obj1->get_mutation_vector() == NULL && obj2->get_mutation_vector() == NULL);
  }
  else
  {
    MutationVector_isEqualTo(obj1->get_mutation_vector(), obj2->get_mutation_vector());
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test ReplicationReport equality
 * \details  --
 * \param    ReplicationReport* obj1
 * \param    ReplicationReport* obj2
 * \return   \e void
 */
void UnitaryTests::ReplicationReport_isEqualTo( ReplicationReport* obj1, ReplicationReport* obj2 )
{
  /*------------------------------------------------------------------ Genome size */
  
  assert(obj1->get_old_genome_size() == obj2->get_old_genome_size());
  assert(obj1->get_new_genome_size() == obj2->get_new_genome_size());
  assert(obj1->get_genome_functional_size() == obj2->get_genome_functional_size());
  
  /*------------------------------------------------------------------ GRN data */
  
  assert(obj1->get_nb_functional_regions() == obj2->get_nb_functional_regions());
  assert(obj1->get_nb_enhancers() == obj2->get_nb_enhancers());
  assert(obj1->get_nb_operators() == obj2->get_nb_operators());
  assert(obj1->get_nb_E_regions() == obj2->get_nb_E_regions());
  assert(obj1->get_nb_TF_regions() == obj2->get_nb_TF_regions());
  assert(obj1->get_nb_mixed_regions() == obj2->get_nb_mixed_regions());
  assert(obj1->get_mean_functional_region_size() == obj2->get_mean_functional_region_size());
  assert(obj1->get_mean_E_region_size() == obj2->get_mean_E_region_size());
  assert(obj1->get_mean_TF_region_size() == obj2->get_mean_TF_region_size());
  assert(obj1->get_mean_mixed_region_size() == obj2->get_mean_mixed_region_size());
  assert(obj1->get_mean_enhancer_size() == obj2->get_mean_enhancer_size());
  assert(obj1->get_mean_operator_size() == obj2->get_mean_operator_size());
  assert(obj1->get_mean_operon_size() == obj2->get_mean_operon_size());
  assert(obj1->get_mean_E_operon_size() == obj2->get_mean_E_operon_size());
  assert(obj1->get_mean_TF_operon_size() == obj2->get_mean_TF_operon_size());
  assert(obj1->get_mean_mixed_operon_size() == obj2->get_mean_mixed_operon_size());
  
  /*------------------------------------------------------------------ Genetic redundancy */
  
  assert(obj1->get_mean_regulation_redundancy() == obj2->get_mean_regulation_redundancy());
  assert(obj1->get_mean_metabolic_redundancy() == obj2->get_mean_metabolic_redundancy());
  
  /*------------------------------------------------------------------ Inherited proteins structure */
  
  assert(obj1->get_inherited_size() == obj2->get_inherited_size());
  assert(obj1->get_inherited_nb_E() == obj2->get_inherited_nb_E());
  assert(obj1->get_inherited_nb_TF() == obj2->get_inherited_nb_TF());
  assert(obj1->get_inherited_nb_inner_enzymes() == obj2->get_inherited_nb_inner_enzymes());
  assert(obj1->get_inherited_nb_inflow_pumps() == obj2->get_inherited_nb_inflow_pumps());
  assert(obj1->get_inherited_nb_outflow_pumps() == obj2->get_inherited_nb_outflow_pumps());
  
  /*------------------------------------------------------------------ Phenotype */
  
  assert(obj1->get_id() == obj2->get_id());
  assert(obj1->get_parent_id() == obj2->get_parent_id());
  assert(obj1->get_generation() == obj2->get_generation());
  assert(obj1->get_x() == obj2->get_x());
  assert(obj1->get_y() == obj2->get_y());
  assert(obj1->get_number_of_updates() == obj2->get_number_of_updates());
  assert(obj1->get_number_of_divisions() == obj2->get_number_of_divisions());
  assert(obj1->get_birth_time() == obj2->get_birth_time());
  assert(obj1->get_death_time() == obj2->get_death_time());
  assert(obj1->get_lifespan() == obj2->get_lifespan());
  assert(obj1->get_toxicity() == obj2->get_toxicity());
  assert(obj1->get_inherited_TF_amount() == obj2->get_inherited_TF_amount());
  assert(obj1->get_inherited_E_amount() == obj2->get_inherited_E_amount());
  assert(obj1->get_TF_amount() == obj2->get_TF_amount());
  assert(obj1->get_E_amount() == obj2->get_E_amount());
  assert(obj1->get_inherited_metabolic_amount() == obj2->get_inherited_metabolic_amount());
  assert(obj1->get_min_metabolic_amount() == obj2->get_min_metabolic_amount());
  assert(obj1->get_metabolic_amount() == obj2->get_metabolic_amount());
  assert(obj1->get_max_metabolic_amount() == obj2->get_max_metabolic_amount());
  assert(obj1->get_metabolic_uptake() == obj2->get_metabolic_uptake());
  assert(obj1->get_metabolic_release() == obj2->get_metabolic_release());
  assert(obj1->get_min_energy() == obj2->get_min_energy());
  assert(obj1->get_mean_energy() == obj2->get_mean_energy());
  assert(obj1->get_max_energy() == obj2->get_max_energy());
  assert(obj1->get_min_score() == obj2->get_min_score());
  assert(obj1->get_mean_score() == obj2->get_mean_score());
  assert(obj1->get_max_score() == obj2->get_max_score());
  assert(obj1->get_metabolic_growth_rate() == obj2->get_metabolic_growth_rate());
  assert(obj1->get_Dmetabolic_growth_rate() == obj2->get_Dmetabolic_growth_rate());
  assert(obj1->get_grn_nb_nodes() == obj2->get_grn_nb_nodes());
  assert(obj1->get_grn_nb_edges() == obj2->get_grn_nb_edges());
  assert(obj1->get_metabolic_nb_nodes() == obj2->get_metabolic_nb_nodes());
  assert(obj1->get_metabolic_nb_edges() == obj2->get_metabolic_nb_edges());
  assert(obj1->get_trophic_group() == obj2->get_trophic_group());
  assert(obj1->get_trophic_level() == obj2->get_trophic_level());
  
  /*------------------------------------------------------------------ List of mutation events */
  
  assert(obj1->get_number_of_events() == obj2->get_number_of_events());
  assert(obj1->get_list_of_events()->size() == obj2->get_list_of_events()->size());
  for (size_t i = 0; i < obj1->get_list_of_events()->size(); i++)
  {
    MutationEvent_isEqualTo(obj1->get_list_of_events()->at(i), obj2->get_list_of_events()->at(i));
  }
  
  /*------------------------------------------------------------------ Point mutations data */
  
  assert(obj1->get_nb_point_mutations() == obj2->get_nb_point_mutations());
  assert(obj1->get_nb_NC_point_mutations() == obj2->get_nb_NC_point_mutations());
  assert(obj1->get_nb_E_point_mutations() == obj2->get_nb_E_point_mutations());
  assert(obj1->get_nb_TF_point_mutations() == obj2->get_nb_TF_point_mutations());
  assert(obj1->get_nb_BS_point_mutations() == obj2->get_nb_BS_point_mutations());
  assert(obj1->get_nb_P_point_mutations() == obj2->get_nb_P_point_mutations());
  
  assert(obj1->get_nb_NC_to_E_transitions() == obj2->get_nb_NC_to_E_transitions());
  assert(obj1->get_nb_NC_to_TF_transitions() == obj2->get_nb_NC_to_TF_transitions());
  assert(obj1->get_nb_NC_to_BS_transitions() == obj2->get_nb_NC_to_BS_transitions());
  assert(obj1->get_nb_NC_to_P_transitions() == obj2->get_nb_NC_to_P_transitions());
  
  assert(obj1->get_nb_E_to_NC_transitions() == obj2->get_nb_E_to_NC_transitions());
  assert(obj1->get_nb_E_to_TF_transitions() == obj2->get_nb_E_to_TF_transitions());
  assert(obj1->get_nb_E_to_BS_transitions() == obj2->get_nb_E_to_BS_transitions());
  assert(obj1->get_nb_E_to_P_transitions() == obj2->get_nb_E_to_P_transitions());
  
  assert(obj1->get_nb_TF_to_NC_transitions() == obj2->get_nb_TF_to_NC_transitions());
  assert(obj1->get_nb_TF_to_E_transitions() == obj2->get_nb_TF_to_E_transitions());
  assert(obj1->get_nb_TF_to_BS_transitions() == obj2->get_nb_TF_to_BS_transitions());
  assert(obj1->get_nb_TF_to_P_transitions() == obj2->get_nb_TF_to_P_transitions());
  
  assert(obj1->get_nb_BS_to_NC_transitions() == obj2->get_nb_BS_to_NC_transitions());
  assert(obj1->get_nb_BS_to_E_transitions() == obj2->get_nb_BS_to_E_transitions());
  assert(obj1->get_nb_BS_to_TF_transitions() == obj2->get_nb_BS_to_TF_transitions());
  assert(obj1->get_nb_BS_to_P_transitions() == obj2->get_nb_BS_to_P_transitions());
  
  assert(obj1->get_nb_P_to_NC_transitions() == obj2->get_nb_P_to_NC_transitions());
  assert(obj1->get_nb_P_to_E_transitions() == obj2->get_nb_P_to_E_transitions());
  assert(obj1->get_nb_P_to_TF_transitions() == obj2->get_nb_P_to_TF_transitions());
  assert(obj1->get_nb_P_to_BS_transitions() == obj2->get_nb_P_to_BS_transitions());
  
  assert(obj1->get_mean_s_mutation_size() == obj2->get_mean_s_mutation_size());
  assert(obj1->get_mean_p_mutation_size() == obj2->get_mean_p_mutation_size());
  assert(obj1->get_mean_km_mutation_size() == obj2->get_mean_km_mutation_size());
  assert(obj1->get_mean_kcat_mutation_size() == obj2->get_mean_kcat_mutation_size());
  assert(obj1->get_mean_BS_tag_mutation_size() == obj2->get_mean_BS_tag_mutation_size());
  assert(obj1->get_mean_coE_tag_mutation_size() == obj2->get_mean_coE_tag_mutation_size());
  assert(obj1->get_mean_TF_tag_mutation_size() == obj2->get_mean_TF_tag_mutation_size());
  assert(obj1->get_mean_basal_expression_level_mutation_size() == obj2->get_mean_basal_expression_level_mutation_size());
  
  /*------------------------------------------------------------------ HGT data */
  
  assert(obj1->get_nb_HGT() == obj2->get_nb_HGT());
  assert(obj1->get_mean_HGT_size() == obj2->get_mean_HGT_size());
  assert(obj1->get_nb_NC_HGT() == obj2->get_nb_NC_HGT());
  assert(obj1->get_nb_E_HGT() == obj2->get_nb_E_HGT());
  assert(obj1->get_nb_TF_HGT() == obj2->get_nb_TF_HGT());
  assert(obj1->get_nb_BS_HGT() == obj2->get_nb_BS_HGT());
  assert(obj1->get_nb_P_HGT() == obj2->get_nb_P_HGT());
  
  /*------------------------------------------------------------------ rearrangements data */
  
  assert(obj1->get_nb_rearrangements() == obj2->get_nb_rearrangements());
  assert(obj1->get_nb_duplicated_NC() == obj2->get_nb_duplicated_NC());
  assert(obj1->get_nb_duplicated_E() == obj2->get_nb_duplicated_E());
  assert(obj1->get_nb_duplicated_TF() == obj2->get_nb_duplicated_TF());
  assert(obj1->get_nb_duplicated_BS() == obj2->get_nb_duplicated_BS());
  assert(obj1->get_nb_duplicated_P() == obj2->get_nb_duplicated_P());
  assert(obj1->get_nb_deleted_NC() == obj2->get_nb_deleted_NC());
  assert(obj1->get_nb_deleted_E() == obj2->get_nb_deleted_E());
  assert(obj1->get_nb_deleted_TF() == obj2->get_nb_deleted_TF());
  assert(obj1->get_nb_deleted_BS() == obj2->get_nb_deleted_BS());
  assert(obj1->get_nb_deleted_P() == obj2->get_nb_deleted_P());
  assert(obj1->get_nb_duplications() == obj2->get_nb_duplications());
  assert(obj1->get_nb_deletions() == obj2->get_nb_deletions());
  assert(obj1->get_nb_translocations() == obj2->get_nb_translocations());
  assert(obj1->get_nb_inversions() == obj2->get_nb_inversions());
  assert(obj1->get_mean_rearrangement_size() == obj2->get_mean_rearrangement_size());
  assert(obj1->get_mean_duplication_size() == obj2->get_mean_duplication_size());
  assert(obj1->get_mean_deletion_size() == obj2->get_mean_deletion_size());
  assert(obj1->get_mean_translocation_size() == obj2->get_mean_translocation_size());
  assert(obj1->get_mean_inversion_size() == obj2->get_mean_inversion_size());
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test Genome equality
 * \details  --
 * \param    Genome* obj1
 * \param    Genome* obj2
 * \return   \e void
 */
void UnitaryTests::Genome_isEqualTo( Genome* obj1, Genome* obj2 )
{
  /*------------------------------------------------------------------ string of pearls */
  
  assert(obj1->get_size() == obj2->get_size());
  for (size_t i = 0; i < obj1->get_size(); i++)
  {
    pearl_isEqualTo(obj1->get_pearl(i), obj2->get_pearl(i));
  }
  assert(obj1->get_buffer_size() == obj2->get_buffer_size());
  assert(obj1->get_coding_size() == obj2->get_coding_size());
  assert(obj1->get_non_coding_size() == obj2->get_non_coding_size());
  
  /*------------------------------------------------------------------ concentration vector */
  
  for (size_t i = 0; i < obj1->get_size(); i++)
  {
    assert(obj1->get_concentration_vector()[i] == obj2->get_concentration_vector()[i]);
    assert(obj1->get_concentration_vector()[i] >= 0.0);
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  assert(obj1->get_nb_NC() == obj2->get_nb_NC());
  assert(obj1->get_nb_E() == obj2->get_nb_E());
  assert(obj1->get_nb_TF() == obj2->get_nb_TF());
  assert(obj1->get_nb_BS() == obj2->get_nb_BS());
  assert(obj1->get_nb_P() == obj2->get_nb_P());
  assert(obj1->get_nb_inner_enzymes() == obj2->get_nb_inner_enzymes());
  assert(obj1->get_nb_inflow_pumps() == obj2->get_nb_inflow_pumps());
  assert(obj1->get_nb_outflow_pumps() == obj2->get_nb_outflow_pumps());
  
  /*------------------------------------------------------------------ pearl positions */
  
  for (size_t i = 0; i < obj1->get_nb_TF(); i++)
  {
    assert(obj1->get_TFi()[i] == obj2->get_TFi()[i]);
    assert(obj1->get_TFi()[i] >= 0);
    assert(obj1->get_TFi()[i] < obj1->get_size());
    assert(obj1->get_pearl(obj1->get_TFi()[i])->type == TRANSCRIPTION_FACTOR);
    assert(obj2->get_pearl(obj2->get_TFi()[i])->type == TRANSCRIPTION_FACTOR);
  }
  for (size_t i = 0; i < obj1->get_nb_P(); i++)
  {
    assert(obj1->get_Pi()[i] == obj2->get_Pi()[i]);
    assert(obj1->get_Pi()[i] >= 0);
    assert(obj1->get_Pi()[i] < obj1->get_size());
    assert(obj1->get_pearl(obj1->get_Pi()[i])->type == PROMOTER);
    assert(obj2->get_pearl(obj2->get_Pi()[i])->type == PROMOTER);
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test InheritedProteins equality
 * \details  --
 * \param    InheritedProteins* obj1
 * \param    InheritedProteins* obj2
 * \return   \e void
 */
void UnitaryTests::InheritedProteins_isEqualTo( InheritedProteins* obj1, InheritedProteins* obj2 )
{
  /*------------------------------------------------------------------ string of pearls */
  
  assert(obj1->get_size() == obj2->get_size());
  for (size_t i = 0; i < obj1->get_size(); i++)
  {
    pearl_isEqualTo(obj1->get_pearl(i), obj2->get_pearl(i));
  }
  assert(obj1->get_buffer_size() == obj2->get_buffer_size());
  
  /*------------------------------------------------------------------ concentration vector */
  
  for (size_t i = 0; i < obj1->get_size(); i++)
  {
    assert(obj1->get_concentration_vector()[i] == obj2->get_concentration_vector()[i]);
    assert(obj1->get_concentration_vector()[i] >= 0.0);
  }
  
  /*------------------------------------------------------------------ statistical data */
  
  assert(obj1->get_nb_E() == obj2->get_nb_E());
  assert(obj1->get_nb_TF() == obj2->get_nb_TF());
  assert(obj1->get_nb_inner_enzymes() == obj2->get_nb_inner_enzymes());
  assert(obj1->get_nb_inflow_pumps() == obj2->get_nb_inflow_pumps());
  assert(obj1->get_nb_outflow_pumps() == obj2->get_nb_outflow_pumps());
  
  /*------------------------------------------------------------------ pearl positions */
  
  for (size_t i = 0; i < obj1->get_nb_E(); i++)
  {
    assert(obj1->get_Ei()[i] == obj2->get_Ei()[i]);
    assert(obj1->get_Ei()[i] >= 0);
    assert(obj1->get_Ei()[i] < obj1->get_size());
    assert(obj1->get_pearl(obj1->get_Ei()[i])->type == ENZYME);
    assert(obj2->get_pearl(obj2->get_Ei()[i])->type == ENZYME);
  }
  for (size_t i = 0; i < obj1->get_nb_TF(); i++)
  {
    assert(obj1->get_TFi()[i] == obj2->get_TFi()[i]);
    assert(obj1->get_TFi()[i] >= 0);
    assert(obj1->get_TFi()[i] < obj1->get_size());
    assert(obj1->get_pearl(obj1->get_TFi()[i])->type == TRANSCRIPTION_FACTOR);
    assert(obj2->get_pearl(obj2->get_TFi()[i])->type == TRANSCRIPTION_FACTOR);
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test SpeciesList equality
 * \details  --
 * \param    SpeciesList* obj1
 * \param    SpeciesList* obj2
 * \return   \e void
 */
void UnitaryTests::SpeciesList_isEqualTo( SpeciesList* obj1, SpeciesList* obj2 )
{
  assert(obj1->get_size() == obj2->get_size());
  assert(obj1->get_amount() == obj2->get_amount());
  for (int i = 0; i < (int)obj1->get_size(); i++)
  {
    assert(obj1->get(i+1) == obj2->get(i+1));
    assert(obj1->get_X()[i] == obj2->get_X()[i]);
    assert(obj1->get(i+1) == obj1->get_X()[i]);
    assert(obj2->get(i+1) == obj2->get_X()[i]);
  }
  
  for (int i = 1; i <= 20; i++)
  {
    obj1->add(i, 10.0);
    obj2->add(i, 10.0);
  }
  assert(obj1->get_size() == obj2->get_size());
  assert(obj1->get_amount() == obj2->get_amount());
  for (int i = 0; i < (int)obj1->get_size(); i++)
  {
    assert(obj1->get(i+1) == obj2->get(i+1));
    assert(obj1->get_X()[i] == obj2->get_X()[i]);
    assert(obj1->get(i+1) == obj1->get_X()[i]);
    assert(obj2->get(i+1) == obj2->get_X()[i]);
  }
  
  for (int i = 1; i <= 20; i++)
  {
    obj1->remove(i, 10.0);
    obj2->remove(i, 10.0);
  }
  assert(obj1->get_size() == obj2->get_size());
  assert(obj1->get_amount() == obj2->get_amount());
  for (int i = 0; i < (int)obj1->get_size(); i++)
  {
    assert(obj1->get(i+1) == obj2->get(i+1));
    assert(obj1->get_X()[i] == obj2->get_X()[i]);
    assert(obj1->get(i+1) == obj1->get_X()[i]);
    assert(obj2->get(i+1) == obj2->get_X()[i]);
  }
  
  size_t size = obj1->get_size();
  obj1->increase_size(1000);
  obj2->increase_size(1000);
  obj1->decrease_size(size);
  obj2->decrease_size(size);
  assert(obj1->get_size() == obj2->get_size());
  assert(obj1->get_amount() == obj2->get_amount());
  for (int i = 0; i < (int)obj1->get_size(); i++)
  {
    assert(obj1->get(i+1) == obj2->get(i+1));
    assert(obj1->get_X()[i] == obj2->get_X()[i]);
    assert(obj1->get(i+1) == obj1->get_X()[i]);
    assert(obj2->get(i+1) == obj2->get_X()[i]);
  }
  
  (void)obj1;
  (void)obj2;
}

/**
 * \brief    Test Prng equality
 * \details  --
 * \param    Prng* obj1
 * \param    Prng* obj2
 * \return   \e void
 */
void UnitaryTests::Prng_isEqualTo( Prng* prng1, Prng* prng2 )
{
  size_t test_size = 10000;
  
  /* 1) test uniform law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->uniform() == prng2->uniform());
  }
  /* 2) test integer uniform law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->uniform(-10000, 10000) == prng2->uniform(-10000, 10000));
  }
  /* 3) test bernouilli law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->bernouilli(0.5) == prng2->bernouilli(0.5));
  }
  /* 4) test binomial law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->binomial(1000, 0.5) == prng2->binomial(1000, 0.5));
  }
  /* 5) test multinomial law */
  for (size_t i = 0; i < test_size; i++)
  {
    unsigned int draws1[10];
    unsigned int draws2[10];
    double probas[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    int N = 1000, K = 10;
    prng1->multinomial(draws1, probas, N, K);
    prng2->multinomial(draws2, probas, N, K);
    for (size_t j = 0; j < 10; j++)
    {
      assert(draws1[j] == draws2[j]);
    }
  }
  /* 6) test gaussian law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->gaussian(0, 1) == prng2->gaussian(0, 1));
  }
  /* 7) test exponential law */
  for (size_t i = 0; i < test_size; i++)
  {
    assert(prng1->exponential(10.0) == prng2->exponential(10.0));
  }
  /* 8) test exponential law */
  for (size_t i = 0; i < test_size; i++)
  {
    double probas[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double sum = 5.0;
    int N = 10;
    assert(prng1->roulette_wheel(probas, sum, N) == prng2->roulette_wheel(probas, sum, N));
    
    (void)probas;
    (void)sum;
    (void)N;
  }
}

/**
 * \brief    modify Parameters attributes
 * \details  --
 * \param    Parameters* obj
 * \return   \e void
 */
void UnitaryTests::modify_Parameters( Parameters* obj )
{
  distribution_law range_law;
  range_law.law    = EXPONENTIAL;
  range_law.min    = 0;
  range_law.max    = 0;
  range_law.mu     = 0.0;
  range_law.sigma  = 0.0;
  range_law.lambda = 0.1234;
  
  /*------------------------------------------------------------------ prng */
  
  obj->set_seed(1010101010);
  
  /*------------------------------------------------------------------ parallel computing */
  
  obj->set_parallel_computing(true);
  
  /*------------------------------------------------------------------ simulation schemes */
  
  obj->set_energy_constraints(true);
  obj->set_membrane_permeability(true);
  obj->set_metabolic_inheritance(true);
  obj->set_enzymatic_inheritance(true);
  obj->set_co_enzyme_activity(true);
  obj->set_score_scheme(ESSENTIAL_METABOLITES_COMBINATORIAL_CONTRIBUTION);
  obj->set_selection_threshold(0.1234);
  
  /*------------------------------------------------------------------ space */
  
  obj->set_width(101);
  obj->set_height(101);
  
  /*------------------------------------------------------------------ output */
  
  obj->set_experiment_backup_step(101010);
  obj->set_tree_backup_step(101010);
  obj->set_figures_generation_step(101010);
  
  /*------------------------------------------------------------------ genome level */
  
  obj->set_metabolite_tag_initial_range(&range_law);
  obj->set_binding_site_tag_initial_range(&range_law);
  obj->set_co_enzyme_tag_initial_range(&range_law);
  obj->set_transcription_factor_tag_initial_range(&range_law);
  
  obj->set_initial_binding_window(1234);
  
  obj->set_initial_number_of_NC_pearls(101);
  obj->set_initial_number_of_E_pearls(101);
  obj->set_initial_number_of_TF_pearls(101);
  obj->set_initial_number_of_BS_pearls(101);
  obj->set_initial_number_of_P_pearls(101);
  
  obj->set_point_mutation_rate(0.1234);
  obj->set_duplication_rate(0.1234);
  obj->set_deletion_rate(0.1234);
  obj->set_translocation_rate(0.1234);
  obj->set_inversion_rate(0.1234);
  obj->set_transition_rate(0.1234);
  
  obj->set_substrate_tag_mutation_size(101);
  obj->set_product_tag_mutation_size(101);
  obj->set_km_mutation_size(0.1234);
  obj->set_kcat_mutation_size(0.1234);
  obj->set_binding_size_tag_mutation_size(101);
  obj->set_co_enzyme_tag_mutation_size(101);
  obj->set_transcription_factor_tag_mutation_size(101);
  obj->set_basal_expression_level_mutation_size(0.1234);
  
  obj->set_maximum_genome_size(1010101010);
  
  /*------------------------------------------------------------------ genetic regulation network level */
  
  obj->set_genetic_regulation_network_timestep(0.1234);
  obj->set_hill_function_theta(0.1234);
  obj->set_hill_function_n(0.1234);
  obj->set_protein_degradation_rate(0.1234);
  
  /*------------------------------------------------------------------ metabolic network level */
  
  obj->set_metabolism_timestep(0.1234);
  
  obj->set_essential_metabolites_toxicity_threshold(0.1234);
  obj->set_non_essential_metabolites_toxicity_threshold(0.1234);
  
  obj->set_initial_metabolites_amount_in_cells(0.1234);
  
  obj->set_energy_reaction_cost_factor(0.1234);
  obj->set_energy_pumping_cost(-0.1234);
  obj->set_energy_transcription_cost(-0.1234);
  obj->set_energy_toxicity_threshold(0.1234);
  
  /*------------------------------------------------------------------ cell level */
  
  obj->set_death_probability(0.1234);
  obj->set_variable_lifespan(true);
  obj->set_gompert_law(0.1234, 0.1234, 0.1234);
  obj->set_initial_membrane_permeability(0.1234);
  
  /*------------------------------------------------------------------ population level */
  
  obj->set_migration_rate(0.1234);
  obj->set_hgt_rate(0.1234);
  
  /*------------------------------------------------------------------ environment level */
  
  environment_properties properties;
  properties.number_of_init_cycles   = 12;
  properties.species_tag_range       = range_law;
  properties.concentration_range     = range_law;
  properties.number_of_species_range = range_law;
  properties.interaction_scheme      = NO_INTERACTION;
  properties.renewal_scheme          = CLEAR_MATTER;
  properties.variation_scheme        = RANDOM_SCHEME;
  properties.variation_localization  = RANDOM_LOCALIZATION;
  properties.introduction_rate       = 0.1234;
  properties.diffusion_rate          = 0.1234;
  properties.degradation_rate        = 0.1234;
  obj->set_environment_properties(&properties);
}
