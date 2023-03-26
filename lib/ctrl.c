/*!
 * Copyright 1997, Regents of the University of Minnesota
 *
 * \file
 * \brief Functions dealing with manipulating the ctrl_t structure
 *
 *
 * \date Started 10/19/1996
 * \author George Karypis
 * \version\verbatim $Id: ctrl.c 10592 2011-07-16 21:17:53Z karypis $ \endverbatime
 */

#include <pm4gnnlib.h>



/*************************************************************************/
/*! This function sets the ctrl_t structure 
*/
/*************************************************************************/
ctrl_t *SetupCtrl(idx_t *options, idx_t ncon, idx_t nparts, real_t *tpwgts, 
            real_t *ubvec, MPI_Comm comm)
{
  idx_t i, j, defopts;
  ctrl_t *ctrl;

  ctrl = (ctrl_t *)gk_malloc(sizeof(ctrl_t), "SetupCtrl: ctrl");
  memset((void *)ctrl, 0, sizeof(ctrl_t));


  /* communicator-related info */
  MPI_Comm_dup(comm, &(ctrl->gcomm));
  ctrl->comm = ctrl->gcomm;
  ctrl->free_comm = 1;
  gkMPI_Comm_rank(ctrl->gcomm, &ctrl->mype);
  gkMPI_Comm_size(ctrl->gcomm, &ctrl->npes);


  /* options[]-related info */
  ctrl->twohop    = GETOPTION(options, METIS_OPTION_TWOHOP, 0);
  ctrl->ondisk    = GETOPTION(options, METIS_OPTION_ONDISK, 0);
  ctrl->dbglvl    = GETOPTION(options, METIS_OPTION_DBGLVL, 0);
  ctrl->dropedges = GETOPTION(options, METIS_OPTION_DROPEDGES, 0);
  ctrl->fast      = GETOPTION(options, METIS_OPTION_FAST, 0);
  ctrl->seed      = GETOPTION(options, METIS_OPTION_SEED, 7);

  ctrl->seed      = (ctrl->seed == 0 ? ctrl->mype : ctrl->seed*ctrl->mype);
  ctrl->sync      = GlobalSEMax(ctrl, ctrl->seed);
  ctrl->pid       = getpid();

  /* common info */
  ctrl->ncon          = ncon;    
  ctrl->nparts        = nparts;    
  ctrl->redist_factor = 1.0;
  ctrl->redist_base   = 1.0;

  /* setup tpwgts */
  ctrl->tpwgts = rmalloc(nparts*ncon, "SetupCtrl: tpwgts");
  if (tpwgts) {
    rcopy(nparts*ncon, tpwgts, ctrl->tpwgts);
  }
  else {
    for (i=0; i<nparts; i++) {
      for (j=0; j<ncon; j++)
        ctrl->tpwgts[i*ncon+j] = 1.0/nparts;
    }
  }

  /* setup ubvec */
  ctrl->ubvec = rsmalloc(ncon, UNBALANCE_FRACTION, "SetupCtrl: ubvec");
  if (ubvec)
    rcopy(ncon, ubvec, ctrl->ubvec);

  /* initialize the various timers */
  InitTimers(ctrl);

  /* initialize the random number generator */
  srand(ctrl->seed);

  return ctrl;
}


/*************************************************************************/
/*! This function computes the invtvwgts of a graph and stores them in ctrl
*/
/*************************************************************************/
void SetupCtrl_invtvwgts(ctrl_t *ctrl, graph_t *graph)
{
  idx_t j, ncon;

  ncon  = graph->ncon;

  ctrl->invtvwgts = rmalloc(ncon, "SetupCtrl_tvwgts: invtvwgts");

  for (j=0; j<ncon; j++) 
    ctrl->invtvwgts[j] = 1.0/GlobalSESum(ctrl, isum(graph->nvtxs, graph->vwgt+j, ncon));
    
}


/*************************************************************************/
/*! This function de-allocates memory allocated for the control structures 
*/
/*************************************************************************/
void FreeCtrl(ctrl_t **r_ctrl)
{
  ctrl_t *ctrl = *r_ctrl;

  FreeWSpace(ctrl);

  if (ctrl->free_comm)
    gkMPI_Comm_free(&(ctrl->gcomm));

  gk_free((void **)&ctrl->invtvwgts, 
      &ctrl->ubvec, &ctrl->tpwgts, 
      &ctrl->sreq, &ctrl->rreq, &ctrl->statuses,
      &ctrl, 
      LTERM);

  *r_ctrl = NULL;
}

