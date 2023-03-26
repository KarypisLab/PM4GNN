/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * mmatch.c
 *
 * This file contains code that finds a matching
 *
 * Started 2/22/96
 * George
 *
 * $Id: match.c 10592 2011-07-16 21:17:53Z karypis $
 *
 */

#include <pm4gnnlib.h>

#define LHTSIZE 8192 /* This should be a power of two */
#define MASK    8191 /* This should be equal to LHTSIZE-1 */


/*************************************************************************/
/*! Finds a HEM matching involving both local and remote vertices */
/*************************************************************************/
void Match_Global(ctrl_t *ctrl, graph_t *graph)
{
  idx_t h, i, ii, j, jj, k, v;
  idx_t nnbrs, nvtxs, ncon, cnvtxs, firstvtx, lastvtx, maxi, maxidx, nkept;
  idx_t otherlastvtx, nrequests, nchanged, pass, nmatched, wside;
  idx_t *xadj, *adjncy, *adjwgt, *vtxdist;
  idx_t *match;
  idx_t *peind, *sendptr, *recvptr;
  idx_t *perm, *iperm, *nperm, *changed;
  real_t *nvwgt, maxnvwgt;
  idx_t *nreqs_pe;
  ikv_t *match_requests, *match_granted, *pe_requests;
  idx_t last_unmatched;

  WCOREPUSH;

  maxnvwgt = 0.5/((real_t)(ctrl->CoarsenTo));

  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs   = graph->nvtxs;
  ncon    = graph->ncon;
  xadj    = graph->xadj;
  adjncy  = graph->adjncy;
  adjwgt  = graph->adjwgt;
  nvwgt   = graph->nvwgt;

  vtxdist  = graph->vtxdist;
  firstvtx = vtxdist[ctrl->mype];
  lastvtx  = vtxdist[ctrl->mype+1];

  nnbrs   = graph->nnbrs;
  peind   = graph->peind;
  sendptr = graph->sendptr;
  recvptr = graph->recvptr;

  match  = graph->match = ismalloc(nvtxs+graph->nrecv, UNMATCHED, "GlobalMatch: match");

  /* wspacemalloc'ed arrays */
  nreqs_pe = iset(nnbrs, 0, iwspacemalloc(ctrl, nnbrs));
  perm     = iwspacemalloc(ctrl, nvtxs);
  iperm    = iwspacemalloc(ctrl, nvtxs);
  nperm    = iwspacemalloc(ctrl, nnbrs);
  changed  = iwspacemalloc(ctrl, nvtxs);

  match_requests = ikvwspacemalloc(ctrl, graph->nsend);
  match_granted  = ikvwspacemalloc(ctrl, graph->nrecv);


  /* create the traversal order */
  FastRandomPermute(nvtxs, perm, 1);
  for (i=0; i<nvtxs; i++)
    iperm[perm[i]] = i;
  iincset(nnbrs, 0, nperm);

  /* mark all heavy vertices as TOO_HEAVY so they will not be matched */
  for (nchanged=i=0; i<nvtxs; i++) {
    for (h=0; h<ncon; h++) {
      if (nvwgt[i*ncon+h] > maxnvwgt) {
        match[i] = TOO_HEAVY;
        nchanged++;
        break;
      }
    }
  }

  /* If no heavy vertices, pick one at random and treat it as such so that
     at the end of the matching each partition will still have one vertex.
     This is to eliminate the cases in which once a matching has been 
     computed, a processor ends up having no vertices */
  if (nchanged == 0) 
    match[RandomInRange(nvtxs)] = TOO_HEAVY;

  CommInterfaceData(ctrl, graph, match, match+nvtxs);


  /* set initial value of nkept based on how over/under weight the
     partition is to begin with */
  nkept = graph->gnvtxs/ctrl->npes - nvtxs;


  /* Find a matching by doing multiple iterations */
  for (nmatched=pass=0; pass<NMATCH_PASSES; pass++) {
    wside = (graph->level+pass)%2;
    nchanged = nrequests = 0;
    for (last_unmatched=ii=nmatched; ii<nvtxs; ii++) {
      i = perm[ii];
      if (match[i] == UNMATCHED) {  /* Unmatched */
        maxidx = i;
        maxi   = -1;

        /* Deal with islands. Find another vertex and match it with */
        if (xadj[i] == xadj[i+1]) {
          last_unmatched = gk_max(ii, last_unmatched)+1;
          for (; last_unmatched<nvtxs; last_unmatched++) {
            k = perm[last_unmatched];
            if (match[k] == UNMATCHED) {
              match[i] = firstvtx+k + (i <= k ? KEEP_BIT : 0);
              match[k] = firstvtx+i + (i >  k ? KEEP_BIT : 0);
              changed[nchanged++] = i;
              changed[nchanged++] = k;
              break;
            }
          }
          continue;
        }

        /* Find a heavy-edge matching */
        for (j=xadj[i]; j<xadj[i+1]; j++) {
          k = adjncy[j];
          if (match[k] == UNMATCHED) { 
            if (ncon == 1) {
              if (maxi == -1 || adjwgt[maxi] < adjwgt[j] ||
                  (adjwgt[maxi] == adjwgt[j] && RandomInRange(xadj[i+1]-xadj[i]) == 0)) {
                maxi   = j;
                maxidx = k;
              }
            }
            else {
              if (maxi == -1 || adjwgt[maxi] < adjwgt[j] ||
                  (adjwgt[maxi] == adjwgt[j] && maxidx < nvtxs && k < nvtxs &&
                   BetterVBalance(ncon,nvwgt+i*ncon,nvwgt+maxidx*ncon,nvwgt+k*ncon) >= 0)) {
                maxi   = j;
                maxidx = k;
              }
            }
          }
        }

        if (maxi != -1) {
          k = adjncy[maxi];
          if (k < nvtxs) { /* Take care the local vertices first */
            /* Here we give preference the local matching by granting it right away */
            match[i] = firstvtx+k + (i <= k ? KEEP_BIT : 0);
            match[k] = firstvtx+i + (i >  k ? KEEP_BIT : 0);
            changed[nchanged++] = i;
            changed[nchanged++] = k;
          }
          else { /* Take care any remote boundary vertices */
            match[k] = MAYBE_MATCHED;
            /* Alternate among which vertices will issue the requests */
            if ((wside == 0 && firstvtx+i < graph->imap[k]) || 
                (wside == 1 && firstvtx+i > graph->imap[k])) { 
              match[i] = MAYBE_MATCHED;
              match_requests[nrequests].key = graph->imap[k];
              match_requests[nrequests].val = firstvtx+i;
              nrequests++;
            }
          }
        }
      }
    }


    /***********************************************************
    * Exchange the match_requests, requests for me are stored in
    * match_granted 
    ************************************************************/
    /* Issue the receives first. Note that from each PE can receive a maximum
       of the interface node that it needs to send it in the case of a mat-vec */
    for (i=0; i<nnbrs; i++) {
      gkMPI_Irecv((void *)(match_granted+recvptr[i]), 2*(recvptr[i+1]-recvptr[i]), IDX_T,
                peind[i], 1, ctrl->comm, ctrl->rreq+i);
    }

    /* Issue the sends next. This needs some work */
    ikvsorti(nrequests, match_requests);
    for (j=i=0; i<nnbrs; i++) {
      otherlastvtx = vtxdist[peind[i]+1];
      for (k=j; k<nrequests && match_requests[k].key < otherlastvtx; k++);
      gkMPI_Isend((void *)(match_requests+j), 2*(k-j), IDX_T, peind[i], 1, 
          ctrl->comm, ctrl->sreq+i);
      j = k;
    }

    /* OK, now get into the loop waiting for the operations to finish */
    gkMPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses);
    for (i=0; i<nnbrs; i++) {
      gkMPI_Get_count(ctrl->statuses+i, IDX_T, nreqs_pe+i);
      nreqs_pe[i] = nreqs_pe[i]/2;  /* Adjust for pairs of IDX_T */
    }
    gkMPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);


    /***********************************************************
    * Now, go and service the requests that you received in 
    * match_granted 
    ************************************************************/
    RandomPermute(nnbrs, nperm, 0);
    for (ii=0; ii<nnbrs; ii++) {
      i = nperm[ii];
      pe_requests = match_granted+recvptr[i];
      for (j=0; j<nreqs_pe[i]; j++) {
        k = pe_requests[j].key;
        PASSERTP(ctrl, k >= firstvtx && k < lastvtx, (ctrl, "%"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX"\n", firstvtx, lastvtx, k, j, peind[i]));
        /* myprintf(ctrl, "Requesting a match %"PRIDX" %"PRIDX"\n", pe_requests[j].key, pe_requests[j].val); */
        if (match[k-firstvtx] == UNMATCHED) { /* Bingo, lets grant this request */
          changed[nchanged++] = k-firstvtx;
          if (nkept >= 0) { /* decide who to keep it based on local balance */
            match[k-firstvtx] = pe_requests[j].val + KEEP_BIT;
            nkept--;
          }
          else {
            match[k-firstvtx] = pe_requests[j].val;
            pe_requests[j].key += KEEP_BIT;
            nkept++;
          }
          /* myprintf(ctrl, "Request from pe:%"PRIDX" (%"PRIDX" %"PRIDX") granted!\n", peind[i], pe_requests[j].val, pe_requests[j].key); */ 
        }
        else { /* We are not granting the request */
          /* myprintf(ctrl, "Request from pe:%"PRIDX" (%"PRIDX" %"PRIDX") not granted!\n", peind[i], pe_requests[j].val, pe_requests[j].key); */ 
          pe_requests[j].key = UNMATCHED;
        }
      }
    }


    /***********************************************************
    * Exchange the match_granted information. It is stored in
    * match_requests 
    ************************************************************/
    /* Issue the receives first. Note that from each PE can receive a maximum
       of the interface node that it needs to send during the case of a mat-vec */
    for (i=0; i<nnbrs; i++) {
      gkMPI_Irecv((void *)(match_requests+sendptr[i]), 2*(sendptr[i+1]-sendptr[i]), IDX_T,
                peind[i], 1, ctrl->comm, ctrl->rreq+i);
    }

    /* Issue the sends next. */
    for (i=0; i<nnbrs; i++) {
      gkMPI_Isend((void *)(match_granted+recvptr[i]), 2*nreqs_pe[i], IDX_T, 
                peind[i], 1, ctrl->comm, ctrl->sreq+i);
    }

    /* OK, now get into the loop waiting for the operations to finish */
    gkMPI_Waitall(nnbrs, ctrl->rreq, ctrl->statuses);
    for (i=0; i<nnbrs; i++) {
      gkMPI_Get_count(ctrl->statuses+i, IDX_T, nreqs_pe+i);
      nreqs_pe[i] = nreqs_pe[i]/2;  /* Adjust for pairs of IDX_T */
    }
    gkMPI_Waitall(nnbrs, ctrl->sreq, ctrl->statuses);


    /***********************************************************
    * Now, go and through the match_requests and update local
    * match information for the matchings that were granted.
    ************************************************************/
    for (i=0; i<nnbrs; i++) {
      pe_requests = match_requests+sendptr[i];
      for (j=0; j<nreqs_pe[i]; j++) {
        match[pe_requests[j].val-firstvtx] = pe_requests[j].key;
        if (pe_requests[j].key != UNMATCHED)
          changed[nchanged++] = pe_requests[j].val-firstvtx;
      }
    }

    for (i=0; i<nchanged; i++) {
      ii = iperm[changed[i]];
      perm[ii] = perm[nmatched];
      iperm[perm[nmatched]] = ii;
      nmatched++;
    }

    CommChangedInterfaceData(ctrl, graph, nchanged, changed, match, 
        match_requests, match_granted);
  }


  /* Find a two-hop matching */
  if (ctrl->twohop) {
    WCOREPUSH;

    /* create the "transposed" adjancency list that includes remote vertices */
    idx_t rnvtxs   = nvtxs+graph->nrecv;
    idx_t *rxadj   = iset(rnvtxs+1, 0, iwspacemalloc(ctrl, rnvtxs+1));
    idx_t *radjncy = iwspacemalloc(ctrl, xadj[nvtxs]);
    for (i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        rxadj[adjncy[j]]++;
      }
    }
    MAKECSR(i, rnvtxs, rxadj);
    for (i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++)
        radjncy[rxadj[adjncy[j]]++] = i;
    }
    SHIFTCSR(i, rnvtxs, rxadj);

    /* sort the unmatched nodes in increasing degree and use that as the
       traversal order for matching purposes */
    idx_t ncand = 0;;
    ikv_t *cand = ikvwspacemalloc(ctrl, nvtxs-nmatched); 
    for (ii=nmatched; ii<nvtxs; ii++) {
      i = perm[ii];
      if (match[i] == UNMATCHED) {
        cand[ncand].key = xadj[i+1]-xadj[i];
        cand[ncand].val = i;
        ncand++;
      }
    }
    ikvsorti(ncand, cand);

    for (nchanged=0, ii=0; ii<ncand; ii++) {
      i = cand[ii].val;
      if (match[i] == UNMATCHED) {  
        /* two-hop -- matched pair will be local but intermediate can be remote */ 
        for (maxi=0, j=xadj[i]; j<xadj[i+1]; j++) {
          v = adjncy[j];
          for (jj=rxadj[v]; jj<rxadj[v+1]; jj++) {
            k = radjncy[jj];
            if (k != i && match[k] == UNMATCHED) {
              match[i] = firstvtx+k + (i <= k ? KEEP_BIT : 0);
              match[k] = firstvtx+i + (i >  k ? KEEP_BIT : 0);
              changed[nchanged++] = i;
              changed[nchanged++] = k;
              maxi = 1;
              break;
            }
          }
          if (maxi) 
            break;
        }
      }
    }

    /* fill the gap by matching pairs of low-degree nodes that are not connected */
    if (1) {
      idx_t gap2goal = ((idx_t)(2.0*(1.0-COARSEN_FRACTION2)*(real_t)nvtxs)) - nmatched - nchanged; 

      if (gap2goal > 0) {
        for (k=0, ii=0; ii<ncand; ii++) {
          if (match[cand[ii].val] == UNMATCHED)
            cand[k++] = cand[ii];
        }
        ncand = k;
  
        for (ii=0; ii<ncand; ii++) {
          i = cand[ii].val;
          if (match[i] == UNMATCHED) {  
            for (jj=ii+1; jj<ncand; jj++) {
              k = cand[jj].val;
              if (match[k] == UNMATCHED) {
                match[i] = firstvtx+k + (i <= k ? KEEP_BIT : 0);
                match[k] = firstvtx+i + (i >  k ? KEEP_BIT : 0);
                changed[nchanged++] = i;
                changed[nchanged++] = k;
                gap2goal -= 2;
                break;
              }
            }
            if (gap2goal <= 0)
              break;
          }
        }
      }
    }

    for (i=0; i<nchanged; i++) {
      ii = iperm[changed[i]];
      perm[ii] = perm[nmatched];
      iperm[perm[nmatched]] = ii;
      nmatched++;
    }

    CommChangedInterfaceData(ctrl, graph, nchanged, changed, match, 
        match_requests, match_granted);

    WCOREPOP;
  }

  /* Traverse the vertices and those that were unmatched, match them with themselves */
  cnvtxs = 0;
  for (i=0; i<nvtxs; i++) {
    if (match[i] == UNMATCHED || match[i] == TOO_HEAVY) {
      match[i] = (firstvtx+i) + KEEP_BIT;
      cnvtxs++;
    }
    else if (match[i] >= KEEP_BIT) {  /* A matched vertex which I get to keep */
      cnvtxs++;
    }
  }

  if (ctrl->dbglvl&DBG_MATCHINFO) {
    PrintVector2(ctrl, nvtxs, firstvtx, match, "Match");
    myprintf(ctrl, "Cnvtxs: %"PRIDX"\n", cnvtxs);
    rprintf(ctrl, "Done with matching...\n");
  }

  WCOREPOP;

  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));
  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));

  CreateCoarseGraph_Global(ctrl, graph, cnvtxs);

  if (ctrl->dropedges) 
    DropEdges(ctrl, graph->coarser);

  IFSET(ctrl->dbglvl, DBG_TIME, gkMPI_Barrier(ctrl->comm));
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));
}


/*************************************************************************/
/*! This function creates the coarser graph after a global matching */
/*************************************************************************/
void CreateCoarseGraph_Global(ctrl_t *ctrl, graph_t *graph, idx_t cnvtxs)
{
  idx_t h, i, j, k, l, ii, jj, ll, nnbrs, nvtxs, nedges, ncon;
  idx_t firstvtx, lastvtx, cfirstvtx, clastvtx, otherlastvtx;
  idx_t npes=ctrl->npes, mype=ctrl->mype;
  idx_t cnedges, nsend, nrecv, nkeepsize, nrecvsize, nsendsize, v, u;
  idx_t *xadj, *adjncy, *adjwgt, *vwgt, *vtxdist, *where;
  idx_t *match, *cmap;
  idx_t *cxadj, *cadjncy, *cadjwgt, *cvwgt, *cwhere = NULL, *cvtxdist;
  idx_t *rsizes, *ssizes, *rlens, *slens, *rgraph, *sgraph, *perm;
  idx_t *peind, *recvptr, *recvind;
  real_t *nvwgt, *cnvwgt;
  graph_t *cgraph;
  ikv_t *scand, *rcand;

  WCOREPUSH;

  nvtxs  = graph->nvtxs;
  ncon   = graph->ncon;
  xadj   = graph->xadj;
  vwgt   = graph->vwgt;
  nvwgt  = graph->nvwgt;
  where  = graph->where;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  match  = graph->match;

  vtxdist  = graph->vtxdist;
  firstvtx = vtxdist[mype];
  lastvtx  = vtxdist[mype+1];

  cmap = graph->cmap = ismalloc(nvtxs+graph->nrecv, -1, "Global_CreateCoarseGraph: cmap");

  nnbrs   = graph->nnbrs;
  peind   = graph->peind;
  recvind = graph->recvind;
  recvptr = graph->recvptr;

  /* Initialize the coarser graph */
  cgraph = CreateGraph();
  cgraph->nvtxs  = cnvtxs;
  cgraph->ncon   = ncon;
  cgraph->level  = graph->level+1;
  cgraph->finer  = graph;
  graph->coarser = cgraph;


  /*************************************************************
  * Obtain the vtxdist of the coarser graph 
  **************************************************************/
  cvtxdist = cgraph->vtxdist = imalloc(npes+1, "Global_CreateCoarseGraph: cvtxdist");
  cvtxdist[npes] = cnvtxs;  /* Use last position in the cvtxdist as a temp buffer */

  gkMPI_Allgather((void *)(cvtxdist+npes), 1, IDX_T, (void *)cvtxdist, 1, 
      IDX_T, ctrl->comm);

  MAKECSR(i, npes, cvtxdist);

  cgraph->gnvtxs = cvtxdist[npes];

#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, npes+1, 0, cvtxdist, "cvtxdist");
#endif


  /*************************************************************
  * Construct the cmap vector 
  **************************************************************/
  cfirstvtx = cvtxdist[mype];
  clastvtx  = cvtxdist[mype+1];

  /* Create the cmap of what you know so far locally. */
  for (cnvtxs=0, i=0; i<nvtxs; i++) {
    if (match[i] >= KEEP_BIT) {
      k = match[i] - KEEP_BIT;
      if (k>=firstvtx && k<firstvtx+i)
        continue;  /* Both (i,k) are local and i has been matched via the (k,i) side */

      cmap[i] = cfirstvtx + cnvtxs++;
      if (k != firstvtx+i && (k>=firstvtx && k<lastvtx)) { /* I'm matched locally */
        cmap[k-firstvtx] = cmap[i];
        match[k-firstvtx] += KEEP_BIT;  /* Add the KEEP_BIT to simplify coding */
      }
    }
  }
  PASSERT(ctrl, cnvtxs == clastvtx-cfirstvtx);

  CommInterfaceData(ctrl, graph, cmap, cmap+nvtxs);

  /* Update the cmap of the locally stored vertices that will go away. 
   * The remote processor assigned cmap for them */
  for (i=0; i<nvtxs; i++) {
    if (match[i] < KEEP_BIT) { /* Only vertices that go away satisfy this */
      cmap[i] = cmap[nvtxs+BSearch(graph->nrecv, recvind, match[i])];
    }
  }

  CommInterfaceData(ctrl, graph, cmap, cmap+nvtxs);


#ifndef NDEBUG
  for (i=0; i<nvtxs+graph->nrecv; i++) {
    if (cmap[i] == -1) 
      errexit("cmap[%"PRIDX"] == -1\n", i);
  }
#endif


#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, nvtxs, firstvtx, cmap, "Cmap");
#endif


  /*************************************************************
  * Determine how many adjcency lists you need to send/receive.
  **************************************************************/
  /* first pass: determine sizes */
  for (nsend=0, nrecv=0, i=0; i<nvtxs; i++) {
    if (match[i] < KEEP_BIT) /* This is going away */
      nsend++;
    else {
      k = match[i]-KEEP_BIT;
      if (k<firstvtx || k>=lastvtx) /* This is comming from afar */
        nrecv++;
    }
  }

  scand = ikvwspacemalloc(ctrl, nsend);
  rcand = graph->rcand = ikvmalloc(nrecv, "CreateCoarseGraph: rcand");

  /* second pass: place them in the appropriate arrays */
  nkeepsize = nsend = nrecv = 0;
  for (i=0; i<nvtxs; i++) {
    if (match[i] < KEEP_BIT) { /* This is going away */
      scand[nsend].key = match[i];
      scand[nsend].val = i;
      nsend++;
    }
    else {
      nkeepsize += (xadj[i+1]-xadj[i]);

      k = match[i]-KEEP_BIT;
      if (k<firstvtx || k>=lastvtx) { /* This is comming from afar */
        rcand[nrecv].key = k;
        rcand[nrecv].val = cmap[i] - cfirstvtx;  /* Set it for use during the partition projection */
        PASSERT(ctrl, rcand[nrecv].val>=0 && rcand[nrecv].val<cnvtxs);
        nrecv++;
      }
    }
  }


#ifdef DEBUG_CONTRACT
  PrintPairs(ctrl, nsend, scand, "scand");
  PrintPairs(ctrl, nrecv, rcand, "rcand");
#endif

  /***************************************************************
  * Determine how many lists and their sizes you will send and 
  * receive for each of the neighboring PEs
  ****************************************************************/
  rlens = graph->rlens = imalloc(nnbrs+1, "CreateCoarseGraph: graph->rlens");
  slens = graph->slens = imalloc(nnbrs+1, "CreateCoarseGraph: graph->slens");

  rsizes = iset(nnbrs, 0, iwspacemalloc(ctrl, nnbrs));
  ssizes = iset(nnbrs, 0, iwspacemalloc(ctrl, nnbrs));

  /* Take care the sending data first */
  ikvsortii(nsend, scand);
  slens[0] = 0;
  for (k=i=0; i<nnbrs; i++) {
    otherlastvtx = vtxdist[peind[i]+1];
    for (; k<nsend && scand[k].key < otherlastvtx; k++)
      ssizes[i] += (xadj[scand[k].val+1]-xadj[scand[k].val]);
    slens[i+1] = k;
  }

  /* Take care the receiving data next. You cannot yet determine the rsizes[] */
  ikvsortii(nrecv, rcand);
  rlens[0] = 0;
  for (k=i=0; i<nnbrs; i++) {
    otherlastvtx = vtxdist[peind[i]+1];
    for (; k<nrecv && rcand[k].key < otherlastvtx; k++);
    rlens[i+1] = k;
  }

#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, nnbrs+1, 0, slens, "slens");
  PrintVector(ctrl, nnbrs+1, 0, rlens, "rlens");
#endif

  /***************************************************************
  * Exchange size information
  ****************************************************************/
  /* Issue the receives first. */
  for (i=0; i<nnbrs; i++) {
    if (rlens[i+1]-rlens[i] > 0)  /* Issue a receive only if you are getting something */
      gkMPI_Irecv((void *)(rsizes+i), 1, IDX_T, peind[i], 1, ctrl->comm, ctrl->rreq+i);
  }

  /* Take care the sending data next */
  for (i=0; i<nnbrs; i++) {
    if (slens[i+1]-slens[i] > 0)  /* Issue a send only if you are sending something */
      gkMPI_Isend((void *)(ssizes+i), 1, IDX_T, peind[i], 1, ctrl->comm, ctrl->sreq+i);
  }

  /* OK, now get into the loop waiting for the operations to finish */
  for (i=0; i<nnbrs; i++) {
    if (rlens[i+1]-rlens[i] > 0)  
      gkMPI_Wait(ctrl->rreq+i, &ctrl->status);
  }
  for (i=0; i<nnbrs; i++) {
    if (slens[i+1]-slens[i] > 0)  
      gkMPI_Wait(ctrl->sreq+i, &ctrl->status);
  }


#ifdef DEBUG_CONTRACT
  PrintVector(ctrl, nnbrs, 0, rsizes, "rsizes");
  PrintVector(ctrl, nnbrs, 0, ssizes, "ssizes");
#endif

  /*************************************************************
  * Allocate memory for received/sent graphs and start sending 
  * and receiving data.
  * rgraph and sgraph is a different data structure than CSR
  * to facilitate single message exchange.
  **************************************************************/
  nrecvsize = isum(nnbrs, rsizes, 1);
  nsendsize = isum(nnbrs, ssizes, 1);
  rgraph = iwspacemalloc(ctrl, (2+ncon)*nrecv+2*nrecvsize);

  WCOREPUSH;  /* for freeing sgraph right away */
  sgraph = iwspacemalloc(ctrl, (2+ncon)*nsend+2*nsendsize);

  /* Deal with the received portion first */
  for (l=i=0; i<nnbrs; i++) {
    /* Issue a receive only if you are getting something */
    if (rlens[i+1]-rlens[i] > 0) {
      gkMPI_Irecv((void *)(rgraph+l), (2+ncon)*(rlens[i+1]-rlens[i])+2*rsizes[i], 
          IDX_T, peind[i], 1, ctrl->comm, ctrl->rreq+i);
      l += (2+ncon)*(rlens[i+1]-rlens[i])+2*rsizes[i];
    }
  }


  /* Deal with the sent portion now */
  for (ll=l=i=0; i<nnbrs; i++) {
    if (slens[i+1]-slens[i] > 0) {  /* Issue a send only if you are sending something */
      for (k=slens[i]; k<slens[i+1]; k++) {
        ii = scand[k].val;
        sgraph[ll++] = firstvtx+ii;
        sgraph[ll++] = xadj[ii+1]-xadj[ii];
        for (h=0; h<ncon; h++)
          sgraph[ll++] = vwgt[ii*ncon+h];
        for (jj=xadj[ii]; jj<xadj[ii+1]; jj++) {
          sgraph[ll++] = cmap[adjncy[jj]];
          sgraph[ll++] = adjwgt[jj];
        }
      }

      PASSERT(ctrl, ll-l == (2+ncon)*(slens[i+1]-slens[i])+2*ssizes[i]);

      /*myprintf(ctrl, "Sending to pe:%"PRIDX", %"PRIDX" lists of size %"PRIDX"\n", peind[i], slens[i+1]-slens[i], ssizes[i]); */
      gkMPI_Isend((void *)(sgraph+l), ll-l, IDX_T, peind[i], 1, ctrl->comm, ctrl->sreq+i);
      l = ll;
    }
  }

  /* OK, now get into the loop waiting for the operations to finish */
  for (i=0; i<nnbrs; i++) {
    if (rlens[i+1]-rlens[i] > 0)  
      gkMPI_Wait(ctrl->rreq+i, &ctrl->status);
  }
  for (i=0; i<nnbrs; i++) {
    if (slens[i+1]-slens[i] > 0)  
      gkMPI_Wait(ctrl->sreq+i, &ctrl->status);
  }


#ifdef DEBUG_CONTRACT
  rprintf(ctrl, "Graphs were sent!\n");
  PrintTransferedGraphs(ctrl, nnbrs, peind, slens, rlens, sgraph, rgraph);
#endif

  WCOREPOP;  /* free sgraph */

  /*************************************************************
  * Setup the mapping from indices returned by BSearch to 
  * those that are actually stored 
  **************************************************************/
  perm = iset(graph->nrecv, -1, iwspacemalloc(ctrl, graph->nrecv));
  for (j=i=0; i<nrecv; i++) {
    perm[BSearch(graph->nrecv, recvind, rgraph[j])] = j+1;
    j += (2+ncon)+2*rgraph[j+1];
  }

  /*************************************************************
  * Finally, create the coarser graph
  **************************************************************/
  /* Allocate memory for the coarser graph, and start creating the coarse graph */
  cxadj   = cgraph->xadj  = imalloc(cnvtxs+1, "CreateCoarserGraph: cxadj");
  cvwgt   = cgraph->vwgt  = imalloc(cnvtxs*ncon, "CreateCoarserGraph: cvwgt");
  cnvwgt  = cgraph->nvwgt = rmalloc(cnvtxs*ncon, "CreateCoarserGraph: cnvwgt");
  if (where != NULL)
    cwhere = cgraph->where = imalloc(cnvtxs, "CreateCoarserGraph: cwhere");
  if (ctrl->dropedges)
    cgraph->ndrop = ismalloc(cnvtxs, 0, "CreateCoarserGraph: cgraph->ndrop");

  /* these are just upper bound estimates for now */
  cadjncy = iwspacemalloc(ctrl, nkeepsize+nrecvsize);
  cadjwgt = iwspacemalloc(ctrl, nkeepsize+nrecvsize);


  /* variables used for the htable */
  idx_t kk, m;
  idx_t htsize, maxhtsize=8192, mask, maxclen;
  idx_t *htable=NULL;
  htable = ismalloc(maxhtsize, -1, "htable");

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    if (match[i] >= KEEP_BIT) {
      v = firstvtx+i; 
      u = match[i]-KEEP_BIT;

      if (u>=firstvtx && u<lastvtx && v > u) 
        continue;  /* I have already collapsed it as (u,v) */

      /* determine the maximum length of the combined adjacency list 
         and the size of the required htable */
      maxclen = xadj[i+1]-xadj[i];
      if (v != u) {
        if (u>=firstvtx && u<lastvtx) /* local vertex */
          maxclen += xadj[u-firstvtx+1]-xadj[u-firstvtx];
        else /* remote vertex */
          maxclen += rgraph[perm[BSearch(graph->nrecv, recvind, u)]];

        /* save this, as you will need it later */ 
        if (ctrl->dropedges) 
          cgraph->ndrop[cnvtxs] = maxclen; 
      }

      for (maxclen*=2, htsize=1; htsize<maxclen; htsize*=2);
      mask = htsize-1;
      if (maxhtsize < htsize) {
        maxhtsize = 2*htsize;
        gk_free((void **)&htable, LTERM);
        htable = ismalloc(maxhtsize, -1, "htable");
      }


      /* Collapse the v vertex first, which you know is local */
      for (h=0; h<ncon; h++)
        cvwgt[cnvtxs*ncon+h] = vwgt[i*ncon+h];
      if (where != NULL)
        cwhere[cnvtxs] = where[i];
      nedges = 0;

      /* Collapse the v (i) vertex first */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = cmap[adjncy[j]];
        if (k != cfirstvtx+cnvtxs) {  /* If this is not an internal edge */
          for (kk=k&mask; htable[kk]!=-1 && cadjncy[htable[kk]]!=k; kk=((kk+1)&mask));
          if ((m = htable[kk]) == -1) { /* Seeing this for first time */
            cadjncy[cnedges+nedges] = k; 
            cadjwgt[cnedges+nedges] = adjwgt[j];
            htable[kk] = cnedges+nedges++;
          }
          else {
            cadjwgt[m] += adjwgt[j];
          }
        }
      }

      /* Collapse the u vertex next */
      if (v != u) {
        if (u>=firstvtx && u<lastvtx) { /* Local vertex */
          u -= firstvtx;
          for (h=0; h<ncon; h++)
            cvwgt[cnvtxs*ncon+h] += vwgt[u*ncon+h];

          for (j=xadj[u]; j<xadj[u+1]; j++) {
            k = cmap[adjncy[j]];
            if (k != cfirstvtx+cnvtxs) {  /* If this is not an internal edge */
              for (kk=k&mask; htable[kk]!=-1 && cadjncy[htable[kk]]!=k; kk=((kk+1)&mask));
              if ((m = htable[kk]) == -1) { /* Seeing this for first time */
                cadjncy[cnedges+nedges] = k; 
                cadjwgt[cnedges+nedges] = adjwgt[j];
                htable[kk] = cnedges+nedges++;
              }
              else {
                cadjwgt[m] += adjwgt[j];
              }
            }
          }
        }
        else { /* Remote vertex */
          u = perm[BSearch(graph->nrecv, recvind, u)];
          for (h=0; h<ncon; h++)
            /* Remember that the +1 stores the vertex weight */
            cvwgt[cnvtxs*ncon+h] += rgraph[(u+1)+h];

          for (j=0; j<rgraph[u]; j++) {
            k = rgraph[u+1+ncon+2*j];
            if (k != cfirstvtx+cnvtxs) {  /* If this is not an internal edge */
              for (kk=k&mask; htable[kk]!=-1 && cadjncy[htable[kk]]!=k; kk=((kk+1)&mask));
              if ((m = htable[kk]) == -1) { /* Seeing this for first time */
                cadjncy[cnedges+nedges] = k; 
                cadjwgt[cnedges+nedges] = rgraph[u+1+ncon+2*j+1];
                htable[kk] = cnedges+nedges++;
              }
              else {
                cadjwgt[m] += rgraph[u+1+ncon+2*j+1];
              }
            }
          }
        }

        /* determine the number of edges that you need to drop so that
           nedges = (|adj(u)|+|adj(v)|)/2 */
        if (ctrl->dropedges) 
          cgraph->ndrop[cnvtxs] = gk_max(0, nedges-(cgraph->ndrop[cnvtxs]>>1)); 
      }

      cnedges += nedges;

      /* reset the htable -- reverse order (LIFO) is critical to prevent cadjncy[-1]
       * indexing due to a remove of an earlier entry */
      for (j=cnedges-1; j>=cxadj[cnvtxs]; j--) {
        k = cadjncy[j];
        for (kk=k&mask; cadjncy[htable[kk]]!=k; kk=((kk+1)&mask));
        htable[kk] = -1;
      }

#ifdef XXX
      for (j=cnedges-1; j>=cxadj[cnvtxs]; j--) {
        k = cadjncy[j];
        for (kk=k&mask; ; kk=((kk+1)&mask)) {
          GKASSERT(htable[kk] != -1);
          if (cadjncy[htable[kk]]==k)
            break;
        }
        htable[kk] = -1;
      }
#endif

      cxadj[++cnvtxs] = cnedges;
    }
  }

  gk_free((void **)&htable, LTERM);

  cgraph->nedges = cnedges;

  /* ADD:  In order to keep from having to change this too much */
  /* ADD:  I kept vwgt array and recomputed nvwgt for each coarser graph */
  for (j=0; j<cnvtxs; j++) {
    for (h=0; h<ncon; h++)
      cgraph->nvwgt[j*ncon+h] = ctrl->invtvwgts[h]*cvwgt[j*ncon+h];
  }

  cgraph->adjncy = imalloc(cnedges, "CreateCoarserGraph: cadjncy");
  cgraph->adjwgt = imalloc(cnedges, "CreateCoarserGraph: cadjwgt");
  icopy(cnedges, cadjncy, cgraph->adjncy);
  icopy(cnedges, cadjwgt, cgraph->adjwgt);

  WCOREPOP;

  /* Note that graph->where works fine even if it is NULL */
  gk_free((void **)&graph->where, LTERM);

}


/*************************************************************************/
/*! This function drops some of the edges of the graph to reduce memory
    consumption. */
/*************************************************************************/
void DropEdges(ctrl_t *ctrl, graph_t *graph)
{
  idx_t i, ii, j, k, istart, iend, nvtxs, maxdegree;
  idx_t *xadj, *adjncy, *adjwgt, *imap, *ndrop;
  idx_t *noise, *medianwgts, *keys;

  WCOREPUSH;

  CommSetup(ctrl, graph);

  nvtxs  = graph->nvtxs;
  xadj   = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  ndrop  = graph->ndrop;
  imap   = graph->imap;

  //myprintf(ctrl, ">nedges: %"PRIDX"\n", xadj[nvtxs]);

  STARTTIMER(ctrl, ctrl->AuxTmr2);

  maxdegree = xadj[1];
  for (i=1; i<nvtxs; i++) 
    maxdegree = gk_max(maxdegree, xadj[i+1]-xadj[i]);

  medianwgts = iwspacemalloc(ctrl, nvtxs+graph->nrecv);
  noise      = iwspacemalloc(ctrl, nvtxs+graph->nrecv);
  keys       = iwspacemalloc(ctrl, maxdegree+1);

  for (i=0; i<nvtxs; i++)
    noise[i] = RandomInRange(1021);
  CommInterfaceData(ctrl, graph, noise, noise+nvtxs);

  /* determine the median weight of each adjacency list */
  for (i=0; i<nvtxs; i++) {
    istart = xadj[i];
    iend   = xadj[i+1];

    if (ndrop[i] == 0 || iend-istart < 5) {
      medianwgts[i] = 0;
    }
    else {
      for (k=0, j=istart; j<iend; j++, k++) 
        keys[k] = (adjwgt[j]<<11) + noise[i] + noise[adjncy[j]];

      isorti(k, keys);
      medianwgts[i] = keys[ndrop[i]];
    }
  }
  CommInterfaceData(ctrl, graph, medianwgts, medianwgts+nvtxs);

  /* compact the adjacency structure of the coarser graph to keep only +ve edges */
  for (k=0, i=0; i<nvtxs; i++) {
    istart = xadj[i];
    iend   = xadj[i+1];
    for (j=istart; j<iend; j++) {
      ii = adjncy[j];
      if ((adjwgt[j]<<11) + noise[i] + noise[ii] >= gk_min(medianwgts[i], medianwgts[ii])) {
        adjncy[k]   = imap[adjncy[j]]; /* keep and switch to global ID space */
        adjwgt[k++] = adjwgt[j];
      }
    }
    xadj[i] = k;
  }
  SHIFTCSR(i, nvtxs, xadj);

  graph->nedges = xadj[nvtxs];

  //myprintf(ctrl, "<nedges: %"PRIDX"\n", xadj[nvtxs]);

  gk_free((void **)&graph->ndrop, LTERM);

  FreeCommSetupFields(graph);

  STOPTIMER(ctrl, ctrl->AuxTmr2);
  
  WCOREPOP;
}
