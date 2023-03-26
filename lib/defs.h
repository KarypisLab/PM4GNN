/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * defs.h
 *
 * This file contains constant definitions
 *
 * Started 8/27/94
 * George
 *
 * $Id: defs.h 10543 2011-07-11 19:32:24Z karypis $
 *
 */


#define NUM_INIT_MSECTIONS              5

#define RIP_SPLIT_FACTOR                8
#define MAX_NPARTS_MULTIPLIER		20

#define REDIST_WGT              2.0
#define MAXNVWGT_FACTOR         2.0

#define N_MOC_REDO_PASSES       10
#define N_MOC_GR_PASSES         8
#define NREMAP_PASSES           8
#define N_MOC_GD_PASSES         6
#define N_MOC_BAL_PASSES        4
#define NMATCH_PASSES           4

#define MAX_NCON_FOR_DIFFUSION  2
#define SMALLGRAPH              10000

#define LTERM                   (void **) 0     /* List terminator for GKfree() */

#define NGD_PASSES		20

#define UNMATCHED		-1
#define MAYBE_MATCHED		-2
#define TOO_HEAVY		-3

#define HTABLE_EMPTY    	-1

#define NGR_PASSES		4	/* Number of greedy refinement passes */
#define NIPARTS			8	/* Number of random initial partitions */
#define NLGR_PASSES		5	/* Number of GR refinement during IPartition */

#define SMALLFLOAT		0.000001


#define COARSEN_FRACTION	        0.75	/* Node reduction between succesive coarsening levels */
#define COARSEN_FRACTION2	        0.65	/* Node reduction between succesive coarsening levels */
#define UNBALANCE_FRACTION		1.05
#define ORDER_UNBALANCE_FRACTION	1.10

#define MAXVWGT_FACTOR		1.4


/* Debug Levels */
#define DBG_TIME	PM4GNN_DBGLVL_TIME 
#define DBG_INFO	PM4GNN_DBGLVL_INFO
#define DBG_PROGRESS   	PM4GNN_DBGLVL_PROGRESS
#define DBG_REFINEINFO	PM4GNN_DBGLVL_REFINEINFO
#define DBG_MATCHINFO	PM4GNN_DBGLVL_MATCHINFO
#define DBG_RMOVEINFO	PM4GNN_DBGLVL_RMOVEINFO
#define DBG_REMAP	PM4GNN_DBGLVL_REMAP
