/* $Header: /data/zender/nco_20150216/nco/src/nco/nco_grp_utl.h,v 1.7 2011-08-01 05:47:46 zender Exp $ */

/* Purpose: Group utilities */

/* Copyright (C) 2011--2011 Charlie Zender
   License: GNU General Public License (GPL) Version 3
   See http://www.gnu.org/copyleft/gpl.html for full license text */

/* Usage:
   #include "nco_var_utl.h" *//* Group utilities */

#ifndef NCO_GRP_UTL_H
#define NCO_GRP_UTL_H

#ifdef HAVE_CONFIG_H
# include <config.h> /* Autotools tokens */
#endif /* !HAVE_CONFIG_H */

/* Standard header files */
#include <stdio.h> /* stderr, FILE, NULL, printf */
#include <stdlib.h> /* strtod, strtol, malloc, getopt, exit */
#include <string.h> /* strcmp. . . */

/* 3rd party vendors */
#include <netcdf.h> /* netCDF definitions and C library */
#include "nco_netcdf.h" /* NCO wrappers for netCDF C library */
#ifdef _OPENMP
#include <omp.h> /* OpenMP pragmas */
#endif /* !_OPENMP */

/* Personal headers */
#include "nco.h" /* netCDF Operator (NCO) definitions */
#include "nco_ctl.h" /* Program flow control functions */
#include "nco_mmr.h" /* Memory management */

/* Dynamic array implementation of group stack */
typedef struct {
  int grp_nbr; /* [nbr] Number of items in stack = number of elements in grp_id array */
  int *grp_id; /* [ID] Group ID */
} grp_stk_sct; 

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

int /* [rcd] Return code */
nco_inq_grps_full /* [fnc] Discover and return IDs of apex and all sub-groups */
(const int grp_id, /* I [ID] Apex group */
 int * const grp_nbr, /* O [nbr] Number of groups */
 int * const grp_ids); /* O [ID] Group IDs of children */
  /* end nco_inq_grps_full() */

void
nco_grp_itr_free /* [fnc] Free group iterator */
(grp_stk_sct * const grp_stk); /* O [sct] Group stack pointer */

int /* [rcd] Return code */
nco_grp_stk_get /* [fnc] Initialize and obtain group iterator */
(const int grp_id, /* I [ID] Apex group */
 grp_stk_sct ** const grp_stk); /* O [sct] Group stack pointer */
  /* end nco_grp_stk_get() */

int /* [rcd] Return code */
nco_grp_stk_nxt /* [fnc] Find and return next group ID */
(grp_stk_sct * const grp_stk, /* O [sct] Group stack pointer */
 int * const grp_id); /* O [ID] Group ID */
  /* end nco_grp_stk_nxt() */

grp_stk_sct * /* O [sct] Group stack pointer */
nco_grp_stk_ntl /* [fnc] Initialize group stack */
(void); 
  /* end nco_grp_stk_ntl() */

void
nco_grp_stk_psh /* [fnc] Push group ID onto stack */
(grp_stk_sct * const grp_stk, /* I/O [sct] Group stack pointer */
 const int grp_id); /* I [ID] Group ID to push */
  /* end nco_grp_stk_psh() */

int /* O [ID] Group ID that was popped */
nco_grp_stk_pop /* [fnc] Remove and return group ID from stack */
(grp_stk_sct * const grp_stk); /* I/O [sct] Group stack pointer */
  /* end nco_grp_stk_pop() */

void
nco_grp_stk_free /* [fnc] Free group stack */
(grp_stk_sct * const grp_stk); /* O [sct] Group stack pointer */
/* end nco_grp_stk_free() */

nm_id_sct * /* O [sct] Variable extraction list */
nco_var4_lst_mk /* [fnc] Create variable extraction list using regular expressions */
(const int nc_id,  /* I [enm] netCDF group ID (root ID of input file) */
 char * const * const var_lst_in, /* I [sng] User-specified list of variable names and rx's */
 const nco_bool EXCLUDE_INPUT_LIST, /* I [flg] Exclude rather than extract */
 const nco_bool EXTRACT_ALL_COORDINATES, /* I [flg] Process all coordinates */
 int * const var_nbr_xtr); /* I/O [nbr] Number of variables in current extraction list */
  /* nco_var4_lst_mk() */

nm_id_sct * /* O [sct] Group extraction list */
nco_grp_lst_mk /* [fnc] Create group extraction list using regular expressions */
(const int nc_id, /* I [enm] netCDF file ID */
 char * const * const grp_lst_in, /* I [sng] User-specified list of group names and rx's */
 const nco_bool EXCLUDE_INPUT_LIST, /* I [flg] Exclude rather than extract */
 int * const grp_nbr); /* I/O [nbr] Number of groups in current extraction list */

int /* [rcd] Return code */
nco_grp_dfn /* [fnc] Define groups in output file */
(const int in_id, /* I [enm] netCDF input-file ID */
 const int out_id, /* I [enm] netCDF output-file ID */
 nm_id_sct *grp_xtr_lst, /* [grp] Number of groups to be defined */
 const int grp_nbr); /* I [nbr] Number of groups to be defined */

int /* [rcd] Return code */
nco_def_grp_rcr
(const int in_id, /* I [enm] netCDF input-file ID */
 const int out_id, /* I [enm] netCDF output-file ID */
 const char * const prn_nm, /* I [sng] Parent group name */
 const int rcr_lvl); /* I [nbr] Recursion level */

#ifdef __cplusplus
} /* end extern "C" */
#endif /* __cplusplus */

#endif /* NCO_GRP_UTL_H */
