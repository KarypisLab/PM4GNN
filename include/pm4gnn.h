/*
 * Copyright 1997-2003, Regents of the University of Minnesota
 *
 * p4gnn.h
 *
 * This file contains function prototypes and constrant definitions for 
 * ParMETIS
 *
 * Started 7/21/03
 * George
 *
 */

#ifndef __p4gnn_h__
#define __p4gnn_h__

#include <mpi.h>
#include <metis.h>

#ifndef _MSC_VER
#define __cdecl
#endif

#if IDXTYPEWIDTH == 32
  /*#define IDX_T         MPI_INT32_T */
  #define IDX_T         MPI_INT
  #define KEEP_BIT      0x40000000L
#elif IDXTYPEWIDTH == 64
  /*#define IDX_T         MPI_INT64_T */
  #define IDX_T         MPI_LONG_LONG_INT
  #define KEEP_BIT      0x4000000000000000LL
#else
  #error "Incorrect user-supplied value fo IDXTYPEWIDTH"
#endif


#if REALTYPEWIDTH == 32
  #define REAL_T        MPI_FLOAT
#elif REALTYPEWIDTH == 64
  #define REAL_T        MPI_DOUBLE
#else
  #error "Incorrect user-supplied value fo REALTYPEWIDTH"
#endif



/*************************************************************************
* Constants 
**************************************************************************/
#define PM4GNN_MAJOR_VERSION        0
#define PM4GNN_MINOR_VERSION        0
#define PM4GNN_SUBMINOR_VERSION     1


/*************************************************************************
* Function prototypes
**************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif

int __cdecl PM4GNN_PartKway(
             idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *adjwgt, 
             idx_t *wgtflag, idx_t *ncon, idx_t *nparts, real_t *tpwgts, real_t *ubvec, 
             idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm);

#ifdef __cplusplus
}
#endif



/*************************************************************************
* Various constants used for the different parameters
**************************************************************************/
/* Debug levels (fields should be ORed) */
#define PM4GNN_DBGLVL_TIME        1      /* Perform timing analysis */
#define PM4GNN_DBGLVL_INFO        2      /* Perform timing analysis */
#define PM4GNN_DBGLVL_PROGRESS    4      /* Show the coarsening progress */
#define PM4GNN_DBGLVL_REFINEINFO  8      /* Show info on communication during folding */
#define PM4GNN_DBGLVL_MATCHINFO   16     /* Show info on matching */
#define PM4GNN_DBGLVL_RMOVEINFO   32     /* Show info on communication during folding */
#define PM4GNN_DBGLVL_REMAP       64     /* Determines if remapping will take place */

#endif 
