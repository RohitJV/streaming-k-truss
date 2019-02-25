/*!
\file
\brief This file contains various constant definitions
\date Started 5/10/17
\author George
\version\verbatim $Id: defs.h 21158 2017-06-10 22:36:08Z karypis $ \endverbatim
*/

#ifndef __DEF_H__
#define __DEF_H__

// #define TEST_CORRECTNESS

#define MAXLINE         1024*128
#define MAX_STRLEN      1024*128


/* command-line options */
#define CMD_KTTYPE              1
#define CMD_IFTYPE              3

#define CMD_IBSIZE              10
#define CMD_JBSIZE              11

#define CMD_SEED                70
#define CMD_DBGLVL              100
#define CMD_HELP                200

#define CMD_STREAM              300

/* kttypes */
#define KTTYPE_BASELINE1        1
#define KTTYPE_BASELINE2        2
#define KTTYPE_BASELINE3        3
#define KTTYPE_BASELINE4        4
#define KTTYPE_BASELINE5        5
#define KTTYPE_BASELINE6        6
#define KTTYPE_BASELINE7        7
#define KTTYPE_BASELINE8        8
#define KTTYPE_BASELINE9        9
#define KTTYPE_BASELINE10       10
#define KTTYPE_BASELINE5b       11
#define KTTYPE_BASELINE10b      12
#define KTTYPE_BASELINE5c       13
#define KTTYPE_TWOPASS          14
#define KTTYPE_SND              15
#define KTTYPE_AND              16
#define KTTYPE_BASELINE10g      17

/* iftype */
#define IFTYPE_TSV              1
#define IFTYPE_METIS            2

/* The text labels for the different tctypes */
static char kttypenames[][20] =
                {"", "baseline1", "baseline2", "baseline3", "baseline4",
                     "baseline5", "baseline6",
                     "baseline7", "baseline8",
                     "baseline9", "baseline10",
                     "baseline5b",
                     "baseline10b",
                     "baseline10g",
                     "baseline5c",
                     "twopass",
                 ""};


/* The text labels for the different iftypes */
static char iftypenames[][10] =
                {"", "tsv", "metis", ""};

#endif
