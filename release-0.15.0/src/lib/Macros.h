
/**
 * \file      Macros.h
 * \authors   Charles Rocabert, Carole Knibbe, Guillaume Beslon
 * \date      08-12-2014
 * \copyright Copyright (C) 2014-2016 Charles Rocabert, Carole Knibbe, Guillaume Beslon. All rights reserved
 * \license   This project is released under the GNU General Public License
 * \brief     Definition of macros
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

#ifndef __EVOEVO__Macros__
#define __EVOEVO__Macros__


/*-------------------------------------*/
/* 1) SEED GENERATION MACROS           */
/*-------------------------------------*/

#define MAXIMUM_SEED 10000000 /*!< Maximum drawable seed */

/*-------------------------------------*/
/* 2) MEMORY MANAGEMENT MACROS         */
/*-------------------------------------*/

#define GENOME_BUFFER              500 /*!< Genome buffer size                        */
#define INHERITED_PROTEINS_BUFFER  500 /*!< Inherited proteins buffer size            */
#define SPECIES_LIST_BUFFER        5   /*!< Species list buffer size                  */
#define GRN_REACTIONS_BUFFER       20  /*!< Number of GRN reactions buffer size       */
#define METABOLIC_REACTIONS_BUFFER 10  /*!< Number of metabolic reactions buffer size */

/*-------------------------------------*/
/* 3) MUTATION RATES MACROS            */
/*-------------------------------------*/

#define NUMBER_OF_MUTATION_RATES 14   /*!< Number of mutation rates     */
#define BREAKPOINT_MUTATION_RATE 1e-4 /*!< Mutation rate at breakpoints */

/*-------------------------------------*/
/* 4) HGT MACROS                       */
/*-------------------------------------*/

#define HGT_MIN_SIZE 1  /*!< HGT minimum size */
#define HGT_MAX_SIZE 10 /*!< HGT maximum size */

/*-------------------------------------*/
/* 5) ODE SOLVER MACROS                */
/*-------------------------------------*/

#define ERR_ABS 1e-30 /*!< ODE solver absolute precision */
#define ERR_REL 1e-02 /*!< ODE solver relative precision */
#define H_INIT  1e-01 /*!< ODE solver initial timestep   */
#define H_MAX   5     /*!< ODE solver maximum timestep   */

/*-------------------------------------*/
/* 6) MICHAELIS MENTEN EQUATION MACROS */
/*-------------------------------------*/

#define KM_MIN_LOG    1.0  /*!< Minimum log10 boundarie of Km   */
#define KM_MAX_LOG    3.0  /*!< Maximum log10 boundarie of Km   */
#define KCAT_MIN_LOG  -3.0 /*!< Minimum log10 boundarie of Kcat */
#define KCAT_MAX_LOG  -1.0 /*!< Maximum log10 boundarie of Kcat */

/*-------------------------------------*/
/* 7) CONCENTRATIONS EVOLUTION MACROS  */
/*-------------------------------------*/

#define MINIMUM_CONCENTRATION                  1e-18 /*!< Minimum concentration threshold                    */
#define MINIMUM_SCORE                          1e-03 /*!< Minimum sustainable score to survive               */
#define MINIMUM_HERITABLE_ENZYME_CONCENTRATION 1e-06 /*!< All enzymes below this threshold are not heritable */

/*-------------------------------------*/
/* 8) GRAPHICS MACROS                  */
/*-------------------------------------*/

#define FRAMERATE      0  /*!< Framerate of the graphic display     */
#define CELL_SCALE     5  /*!< Cell's scale in the graphic display  */
#define CELL_SPACE     0  /*!< Space displayed between each cell    */
#define GRADIENT_SCALE 30 /*!< Gradient's scale                     */
#define GRADIENT_SIZE  10 /*!< Gradient's size (nb of color points) */
#define TEXT_SCALE     50 /*!< Text's scale                         */
#define SPAN           5  /*!< Span of the render window            */


#endif /* defined(__EVOEVO__Macros__) */
