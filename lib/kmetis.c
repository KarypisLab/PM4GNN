/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * kmetis.c
 *
 * This is the entry point of ParMETIS_PartKway
 *
 * Started 10/19/96
 * George
 *
 * $Id: kmetis.c 10757 2011-09-15 22:07:47Z karypis $
 *
 */

#include <pm4gnnlib.h>

#define RIFN(x) \
  if ((x) == NULL) {\
    printf("PM4GNN ERROR " #x " is NULL.\n");\
    return 0;\
  }

#define RIFNP(x) \
  if ((*x) <= 0) {\
    printf("PM4GNN ERROR " #x " is <= 0.\n");\
    return 0;\
  }


/***********************************************************************************
* This function is the entry point of the parallel k-way multilevel partitionioner. 
* This function assumes nothing about the graph distribution.
* It is the general case.
************************************************************************************/
int PM4GNN_PartKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
        idx_t *adjwgt, idx_t *wgtflag, idx_t *ncon, idx_t *nparts, real_t *tpwgts, 
        real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm)
{
  idx_t h, i, status, nvtxs, npes, mype, seed, dbglvl;
  ctrl_t *ctrl;
  graph_t *graph;
  idx_t moptions[METIS_NOPTIONS];
  size_t curmem;


  /* Check the input parameters and return if an error */
  status = CheckInputsPartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag, 
               ncon, nparts, tpwgts, ubvec, options, edgecut, part, comm);
  if (GlobalSEMinComm(*comm, status) == 0) 
    return METIS_ERROR;

  status = METIS_OK;
  gk_malloc_init();
  curmem = gk_GetCurMemoryUsed();

  /* Set up the control */
  ctrl = SetupCtrl(options, *ncon, *nparts, tpwgts, ubvec, *comm);
  npes = ctrl->npes;
  mype = ctrl->mype;


  /* Take care the nparts == 1 case */
  if (*nparts == 1) {
    iset(vtxdist[mype+1]-vtxdist[mype], 0, part); 
    *edgecut = 0;
    goto DONE;
  }


  /* Take care of npes == 1 case */
  if (npes == 1) {
    nvtxs = vtxdist[1] - vtxdist[0];
    METIS_SetDefaultOptions(moptions);
    moptions[METIS_OPTION_NUMBERING] = 0;

    status = METIS_PartGraphKway(&nvtxs, ncon, xadj, adjncy, vwgt, NULL, 
                 adjwgt, nparts, tpwgts, ubvec, moptions, edgecut, part);
 
    goto DONE;
  }


  /* Setup the graph */
  graph = SetupGraph(ctrl, *ncon, vtxdist, xadj, vwgt, adjncy, adjwgt, *wgtflag);

  /* Setup the workspace */
  AllocateWSpace(ctrl, 10*graph->nvtxs);


  /* Partition the graph */
  STARTTIMER(ctrl, ctrl->TotalTmr);

  ctrl->CoarsenTo = gk_min(vtxdist[npes]+1, 25*(*ncon)*gk_max(npes, *nparts));
  if (vtxdist[npes] < SMALLGRAPH 
      || vtxdist[npes] < npes*20 
      || GlobalSESum(ctrl, graph->nedges) == 0) { /* serially */
    IFSET(ctrl->dbglvl, DBG_INFO, 
        rprintf(ctrl, "Partitioning a graph of size %"PRIDX" serially\n", vtxdist[npes]));
    PartitionSmallGraph(ctrl, graph);
  }
  else { /* in parallel */
    Global_Partition(ctrl, graph);
  }
  ParallelReMapGraph(ctrl, graph);

  icopy(graph->nvtxs, graph->where, part);
  *edgecut = graph->mincut;

  STOPTIMER(ctrl, ctrl->TotalTmr);


  /* Print out stats */
  IFSET(ctrl->dbglvl, DBG_TIME, PrintTimingInfo(ctrl));
  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->gcomm));
  IFSET(ctrl->dbglvl, DBG_INFO, PrintPostPartInfo(ctrl, graph, 0));

  FreeInitialGraphAndRemap(&graph);

DONE:
  FreeCtrl(&ctrl);
  if (gk_GetCurMemoryUsed() - curmem > 0) {
    printf("ParMETIS appears to have a memory leak of %zdbytes. Report this.\n", 
        (ssize_t)(gk_GetCurMemoryUsed() - curmem));
  }
  gk_malloc_cleanup(0);

  return (int)status;
}


/*************************************************************************/
/*! This function checks the validity of the inputs for PartKway
  */
/*************************************************************************/
int CheckInputsPartKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
        idx_t *adjwgt, idx_t *wgtflag, idx_t *ncon, idx_t *nparts, real_t *tpwgts, 
        real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm)
{
  idx_t i, j, mype;
  real_t sum;

  /* Check that the supplied information is actually non-NULL */
  if (comm == NULL) {
    printf("PM4GNN ERROR: comm is NULL. Aborting\n");
    abort();
  }
  gkMPI_Comm_rank(*comm, &mype);

  RIFN(vtxdist);
  RIFN(xadj);
  RIFN(adjncy);
  RIFN(wgtflag);
  RIFN(ncon);
  RIFN(nparts);
  RIFN(tpwgts);
  RIFN(ubvec);
  RIFN(options);
  RIFN(edgecut);
  RIFN(part);

  if (*wgtflag == 2 || *wgtflag == 3) {
    RIFN(vwgt);
    for (j=0; j<*ncon; j++) {
      if (GlobalSESumComm(*comm, isum(vtxdist[mype+1]-vtxdist[mype], vwgt+j, *ncon)) == 0) {
        printf("PM4GNN ERROR: sum weight for constraint %"PRIDX" is zero.\n", j);
        return 0;
      }
    }
  }
  if (*wgtflag == 1 || *wgtflag == 3) 
    RIFN(adjwgt);


  /* Check that the supplied information is actually valid/reasonable */
  if (vtxdist[mype+1]-vtxdist[mype] < 1) {
    printf("PM4GNN ERROR: Poor initial vertex distribution. "
           "Processor %"PRIDX" has no vertices assigned to it!\n", mype);
    return 0;
  }

  RIFNP(ncon);
  RIFNP(nparts);


  for (j=0; j<*ncon; j++) {
    sum = rsum(*nparts, tpwgts+j, *ncon);
    if (sum < 0.999 || sum > 1.001) {
      printf("PM4GNN ERROR: The sum of tpwgts for constraint #%"PRIDX" is not 1.0\n", j);
      return 0;
    }
  }
  for (j=0; j<*ncon; j++) {
    for (i=0; i<*nparts; i++) {
      if (tpwgts[i*(*ncon)+j] < 0.0 || tpwgts[i] > 1.001) {
        printf("PM4GNN ERROR: The tpwgts for constraint #%"PRIDX" and partition #%"PRIDX" is out of bounds.\n", j, i);
        return 0;
      }
    }
  }


  for (j=0; j<*ncon; j++) {
    if (ubvec[j] <= 1.0) {
      printf("PM4GNN ERROR: The ubvec for constraint #%"PRIDX" must be > 1.0\n", j);
      return 0;
    }
  }

  return 1;
}


/*************************************************************************/
/*! This function is the driver to the multi-constraint partitioning 
    algorithm.
*/
/*************************************************************************/
void Global_Partition(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, ncon, nparts;
  real_t ftmp, ubavg, lbavg, *lbvec;

  WCOREPUSH;
 
  ncon   = graph->ncon;
  nparts = ctrl->nparts;
  ubavg  = ravg(graph->ncon, ctrl->ubvec);

  CommSetup(ctrl, graph);

  lbvec = rwspacemalloc(ctrl, ncon);

  if (ctrl->dbglvl&DBG_PROGRESS) {
    idx_t j, nledges=0;
    for (j=0; j<graph->nedges; j++) {
      if (graph->adjncy[j] < graph->nvtxs)
        nledges++;
    }
    rprintf(ctrl, "[%10"PRIDX" %10"PRIDX" %10"PRIDX" %10"PRIDX"] [%.2"PRREAL"] [%"PRIDX"] [", 
        graph->gnvtxs, GlobalSESum(ctrl, graph->nedges),
        GlobalSESumFloat(ctrl, 1.0*nledges/graph->nedges)/ctrl->npes,
	GlobalSEMin(ctrl, graph->nvtxs), GlobalSEMax(ctrl, graph->nvtxs), 
        ctrl->CoarsenTo);
    for (i=0; i<ncon; i++)
      rprintf(ctrl, " %.2"PRREAL"", GlobalSEMinFloat(ctrl,graph->nvwgt[rargmin_strd(graph->nvtxs, graph->nvwgt+i, ncon)*ncon+i]));  
    rprintf(ctrl, "] [");
    for (i=0; i<ncon; i++)
      rprintf(ctrl, " %.2"PRREAL"", GlobalSEMaxFloat(ctrl, graph->nvwgt[rargmax_strd(graph->nvtxs, graph->nvwgt+i, ncon)*ncon+i]));  
    rprintf(ctrl, "]\n");
  }

  if (graph->gnvtxs < 1.3*ctrl->CoarsenTo ||
	(graph->finer != NULL &&
	graph->gnvtxs > graph->finer->gnvtxs*COARSEN_FRACTION)) {

    /* Done with coarsening. Find a partition */
    AllocateRefinementWorkSpace(ctrl, 2*graph->nedges);
    graph->where = imalloc(graph->nvtxs+graph->nrecv, "graph->where");

    InitPartition(ctrl, graph);

    if (ctrl->dbglvl&DBG_PROGRESS) {
      ComputePartitionParams(ctrl, graph);
      ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", 
          graph->gnvtxs, graph->mincut);
      for (i=0; i<graph->ncon; i++) 
        rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
      rprintf(ctrl, "\n");

      /* free memory allocated by ComputePartitionParams */
      gk_free((void **)&graph->ckrinfo, &graph->lnpwgts, &graph->gnpwgts, LTERM);
    }

    /* In case no coarsening took place */
    if (graph->finer == NULL) {
      ComputePartitionParams(ctrl, graph);
      KWayFM(ctrl, graph, NGR_PASSES);
    }
  }
  else {
    Match_Global(ctrl, graph);

    graph_WriteToDisk(ctrl, graph);

    Global_Partition(ctrl, graph->coarser);

    graph_ReadFromDisk(ctrl, graph);

    ProjectPartition(ctrl, graph);

    ComputePartitionParams(ctrl, graph);

    if (graph->ncon > 1 && graph->level < 3) {
      for (i=0; i<ncon; i++) {
        ftmp = rsum(nparts, graph->gnpwgts+i, ncon);
        if (ftmp != 0.0)
          lbvec[i] = (real_t)(nparts) *
          graph->gnpwgts[rargmax_strd(nparts, graph->gnpwgts+i, ncon)*ncon+i]/ftmp;
        else
          lbvec[i] = 1.0;
      }
      lbavg = ravg(graph->ncon, lbvec);

      if (lbavg > ubavg + 0.035) {
        if (ctrl->dbglvl&DBG_PROGRESS) {
          ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
          rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", 
              graph->gnvtxs, graph->mincut);
          for (i=0; i<graph->ncon; i++) 
            rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
          rprintf(ctrl, " [b]\n");
	}

        KWayBalance(ctrl, graph, graph->ncon);
      }
    }

    KWayFM(ctrl, graph, NGR_PASSES);

    if (ctrl->dbglvl&DBG_PROGRESS) {
      ComputeParallelBalance(ctrl, graph, graph->where, lbvec);
      rprintf(ctrl, "nvtxs: %10"PRIDX", cut: %8"PRIDX", balance: ", 
          graph->gnvtxs, graph->mincut);
      for (i=0; i<graph->ncon; i++) 
        rprintf(ctrl, "%.3"PRREAL" ", lbvec[i]);
      rprintf(ctrl, "\n");
    }

    if (graph->level != 0)
      gk_free((void **)&graph->lnpwgts, (void **)&graph->gnpwgts, LTERM);
  }

  WCOREPOP;

}


/*************************************************************************/
/*! This function computes a partitioning of a small graph
  */
/*************************************************************************/
void PartitionSmallGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, h, ncon, nparts, npes, mype;
  idx_t moptions[METIS_NOPTIONS];
  idx_t me;
  idx_t *mypart;
  int lpecut[2], gpecut[2];
  graph_t *agraph;
  idx_t *sendcounts, *displs;
  real_t *gnpwgts, *lnpwgts;

  ncon   = graph->ncon;
  nparts = ctrl->nparts;
  npes   = ctrl->npes;
  mype   = ctrl->mype;

  WCOREPUSH;

  CommSetup(ctrl, graph);
  graph->where = imalloc(graph->nvtxs+graph->nrecv, "PartitionSmallGraph: where");
  agraph       = AssembleGraph(ctrl, graph);
  mypart       = iwspacemalloc(ctrl, agraph->nvtxs);

  METIS_SetDefaultOptions(moptions);
  moptions[METIS_OPTION_SEED] = ctrl->sync + mype;

  METIS_PartGraphKway(&agraph->nvtxs, &ncon, agraph->xadj, agraph->adjncy, 
        agraph->vwgt, NULL, agraph->adjwgt, &nparts, ctrl->tpwgts, NULL,
	moptions, &graph->mincut, mypart);

  lpecut[0] = graph->mincut;
  lpecut[1] = mype;
  gkMPI_Allreduce(lpecut, gpecut, 1, MPI_2INT, MPI_MINLOC, ctrl->comm);
  graph->mincut = gpecut[0];

  if (lpecut[1] == gpecut[1] && gpecut[1] != 0)
    gkMPI_Send((void *)mypart, agraph->nvtxs, IDX_T, 0, 1, ctrl->comm);
  if (lpecut[1] == 0 && gpecut[1] != 0)
    gkMPI_Recv((void *)mypart, agraph->nvtxs, IDX_T, gpecut[1], 1, ctrl->comm, &ctrl->status);

  sendcounts = iwspacemalloc(ctrl, npes);
  displs     = iwspacemalloc(ctrl, npes);

  for (i=0; i<npes; i++) {
    sendcounts[i] = graph->vtxdist[i+1]-graph->vtxdist[i];
    displs[i] = graph->vtxdist[i];
  }

  gkMPI_Scatterv((void *)mypart, sendcounts, displs, IDX_T,
      (void *)graph->where, graph->nvtxs, IDX_T, 0, ctrl->comm);

  lnpwgts = graph->lnpwgts = rmalloc(nparts*ncon, "lnpwgts");
  gnpwgts = graph->gnpwgts = rmalloc(nparts*ncon, "gnpwgts");
  rset(nparts*ncon, 0, lnpwgts);
  for (i=0; i<graph->nvtxs; i++) {
    me = graph->where[i];
    for (h=0; h<ncon; h++)
      lnpwgts[me*ncon+h] += graph->nvwgt[i*ncon+h];
  }
  gkMPI_Allreduce((void *)lnpwgts, (void *)gnpwgts, nparts*ncon, REAL_T, MPI_SUM, ctrl->comm);

  FreeGraph(&agraph);

  WCOREPOP;

  return;
}


/*************************************************************************/
/*! This function assembles the graph into a single processor. */
/*************************************************************************/
graph_t *AssembleGraph(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, j, k, l, gnvtxs, nvtxs, ncon, gnedges, nedges, gsize;
  idx_t *xadj, *vwgt, *adjncy, *adjwgt, *vtxdist, *imap;
  idx_t *axadj, *aadjncy, *aadjwgt, *avwgt, *alabel;
  idx_t *mygraph, *ggraph;
  idx_t *rcounts, *rdispls, mysize;
  real_t *anvwgt;
  graph_t *agraph;

  WCOREPUSH;

  gnvtxs  = graph->gnvtxs;
  nvtxs   = graph->nvtxs;
  ncon    = graph->ncon;
  nedges  = graph->xadj[nvtxs];
  xadj    = graph->xadj;
  vwgt    = graph->vwgt;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  vtxdist = graph->vtxdist;
  imap    = graph->imap;

  /*************************************************************/
  /* Determine the # of idx_t to receive from each processor */
  /*************************************************************/
  rcounts = iwspacemalloc(ctrl, ctrl->npes);
  mysize  = (1+ncon)*nvtxs + 2*nedges;
  gkMPI_Allgather((void *)(&mysize), 1, IDX_T, (void *)rcounts, 1, IDX_T, ctrl->comm);

  rdispls = iwspacemalloc(ctrl, ctrl->npes+1);
  for (rdispls[0]=0, i=1; i<ctrl->npes+1; i++)
    rdispls[i] = rdispls[i-1] + rcounts[i-1];

  /* allocate memory for the recv buffer of the assembled graph */
  gsize   = rdispls[ctrl->npes];
  ggraph  = iwspacemalloc(ctrl, gsize);

  /* Construct the one-array storage format of the assembled graph */
  WCOREPUSH;  /* for freeing mygraph */
  mygraph = iwspacemalloc(ctrl, mysize);

  for (k=i=0; i<nvtxs; i++) {
    mygraph[k++] = xadj[i+1]-xadj[i];
    for (j=0; j<ncon; j++)
      mygraph[k++] = vwgt[i*ncon+j];
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      mygraph[k++] = imap[adjncy[j]];
      mygraph[k++] = adjwgt[j];
    }
  }
  PASSERT(ctrl, mysize == k);

  /**************************************/
  /* Assemble and send the entire graph */
  /**************************************/
  gkMPI_Allgatherv((void *)mygraph, mysize, IDX_T, (void *)ggraph, 
      rcounts, rdispls, IDX_T, ctrl->comm);

  WCOREPOP;  /* free mygraph */


  /*******************************************/
  /* Allocate memory for the assembled graph */
  /*******************************************/
  agraph = CreateGraph();
  agraph->nvtxs  = gnvtxs;
  agraph->ncon   = ncon;
  agraph->nedges = gnedges = (gsize-(1+ncon)*gnvtxs)/2;

  axadj   = agraph->xadj   = imalloc(gnvtxs+1, "AssembleGraph: axadj");
  avwgt   = agraph->vwgt   = imalloc(gnvtxs*ncon, "AssembleGraph: avwgt");
  anvwgt  = agraph->nvwgt  = rmalloc(gnvtxs*ncon, "AssembleGraph: anvwgt");
  aadjncy = agraph->adjncy = imalloc(gnedges, "AssembleGraph: adjncy");
  aadjwgt = agraph->adjwgt = imalloc(gnedges, "AssembleGraph: adjwgt");
  alabel  = agraph->label  = imalloc(gnvtxs, "AssembleGraph: alabel");

  for (k=j=i=0; i<gnvtxs; i++) {
    axadj[i] = ggraph[k++];
    for (l=0; l<ncon; l++)
      avwgt[i*ncon+l] = ggraph[k++];
    for (l=0; l<axadj[i]; l++) {
      aadjncy[j] = ggraph[k++];
      aadjwgt[j] = ggraph[k++];
      j++;
    }
  }

  /*********************************/
  /* Now fix up the received graph */
  /*********************************/
  MAKECSR(i, gnvtxs, axadj);

  for (i=0; i<gnvtxs; i++) {
    for (j=0; j<ncon; j++)
      anvwgt[i*ncon+j] = ctrl->invtvwgts[j]*agraph->vwgt[i*ncon+j];
  }

  iincset(gnvtxs, 0, alabel);

  WCOREPOP;

  return agraph;
}


