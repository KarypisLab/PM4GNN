/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * main.c
 *
 * This is the entry point of the ILUT
 *
 * Started 10/19/95
 * George
 *
 * $Id: parmetis.c,v 1.5 2003/07/30 21:18:54 karypis Exp $
 *
 */

#include <pm4gnnbin.h>
#define	MAXLINE	64*1024*1024


/*************************************************************************
* This function reads the CSR matrix
**************************************************************************/
void ParallelReadGraph(graph_t *graph, char *filename, MPI_Comm comm)
{
  idx_t i, k, l, pe;
  idx_t npes, mype, ier;
  idx_t gnvtxs, nvtxs, your_nvtxs, your_nedges, gnedges;
  idx_t maxnvtxs = -1, maxnedges = -1;
  idx_t readew = -1, readvw = -1, dummy, edge;
  idx_t *vtxdist, *xadj, *adjncy, *vwgt, *adjwgt;
  idx_t *your_xadj, *your_adjncy, *your_vwgt, *your_adjwgt, graphinfo[4];
  idx_t fmt, ncon, nobj;
  MPI_Status stat;
  char *line = NULL, *oldstr, *newstr;
  FILE *fpin = NULL;

  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  vtxdist = graph->vtxdist = ismalloc(npes+1, 0, "ReadGraph: vtxdist");

  if (mype == npes-1) {
    ier = 0;
    fpin = fopen(filename, "r");

    if (fpin == NULL) {
      printf("COULD NOT OPEN FILE '%s' FOR SOME REASON!\n", filename);
      ier++;
    }

    gkMPI_Bcast(&ier, 1, IDX_T, npes-1, comm);
    if (ier > 0){
      MPI_Finalize();
      exit(0);
    }

    line = gk_cmalloc(MAXLINE+1, "line");

    while (fgets(line, MAXLINE, fpin) && line[0] == '%');

    fmt = ncon = nobj = 0;
    sscanf(line, "%"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX" %"PRIDX"", 
        &gnvtxs, &gnedges, &fmt, &ncon, &nobj);
    gnedges *=2;
    readew = (fmt%10 > 0);
    readvw = ((fmt/10)%10 > 0);
    graph->ncon = ncon = (ncon == 0 ? 1 : ncon);
    graph->nobj = nobj = (nobj == 0 ? 1 : nobj);

    /*printf("Nvtxs: %"PRIDX", Nedges: %"PRIDX", Ncon: %"PRIDX"\n", gnvtxs, gnedges, ncon); */

    graphinfo[0] = ncon;
    graphinfo[1] = nobj;
    graphinfo[2] = readvw;
    graphinfo[3] = readew;
    gkMPI_Bcast((void *)graphinfo, 4, IDX_T, npes-1, comm);

    /* Construct vtxdist and send it to all the processors */
    vtxdist[0] = 0;
    for (i=0,k=gnvtxs; i<npes; i++) {
      l = k/(npes-i);
      vtxdist[i+1] = vtxdist[i]+l;
      k -= l;
    }

    gkMPI_Bcast((void *)vtxdist, npes+1, IDX_T, npes-1, comm);
  }
  else {
    gkMPI_Bcast(&ier, 1, IDX_T, npes-1, comm);
    if (ier > 0){
      MPI_Finalize();
      exit(0);
    }

    gkMPI_Bcast((void *)graphinfo, 4, IDX_T, npes-1, comm);
    graph->ncon = ncon = graphinfo[0];
    graph->nobj = nobj = graphinfo[1];
    readvw = graphinfo[2];
    readew = graphinfo[3];

    gkMPI_Bcast((void *)vtxdist, npes+1, IDX_T, npes-1, comm);
  }

  if ((ncon > 1 && !readvw) || (nobj > 1 && !readew)) {
    printf("fmt and ncon/nobj are inconsistant.  Exiting...\n");
    gkMPI_Finalize();
    exit(-1);
  }


  graph->gnvtxs = vtxdist[npes];
  nvtxs = graph->nvtxs = vtxdist[mype+1]-vtxdist[mype];
  xadj  = graph->xadj  = imalloc(graph->nvtxs+1, "ParallelReadGraph: xadj");
  vwgt  = graph->vwgt  = imalloc(graph->nvtxs*ncon, "ParallelReadGraph: vwgt");

  /*******************************************/
  /* Go through first time and generate xadj */
  /*******************************************/
  if (mype == npes-1) {
    maxnvtxs = vtxdist[1];
    for (i=1; i<npes; i++) 
      maxnvtxs = (maxnvtxs < vtxdist[i+1]-vtxdist[i] ? vtxdist[i+1]-vtxdist[i] : maxnvtxs);

    your_xadj = imalloc(maxnvtxs+1, "your_xadj");
    your_vwgt = ismalloc(maxnvtxs*ncon, 1, "your_vwgt");

    maxnedges = 0;
    for (pe=0; pe<npes; pe++) {
      your_nvtxs = vtxdist[pe+1]-vtxdist[pe];

      for (i=0; i<your_nvtxs; i++) {
        your_nedges = 0;

        while (fgets(line, MAXLINE, fpin) && line[0] == '%'); /* skip lines with '#' */
        oldstr = line;
        newstr = NULL;

        if (readvw) {
          for (l=0; l<ncon; l++) {
            your_vwgt[i*ncon+l] = strtoidx(oldstr, &newstr, 10);
            oldstr = newstr;
          }
        }

        for (;;) {
          edge = strtoidx(oldstr, &newstr, 10) -1;
          oldstr = newstr;

          if (edge < 0)
            break;

          if (readew) {
            for (l=0; l<nobj; l++) {
              dummy  = strtoidx(oldstr, &newstr, 10);
              oldstr = newstr;
            }
          }
          your_nedges++;
        }
        your_xadj[i] = your_nedges;
      }

      MAKECSR(i, your_nvtxs, your_xadj);
      maxnedges = (maxnedges < your_xadj[your_nvtxs] ? your_xadj[your_nvtxs] : maxnedges);

      if (pe < npes-1) {
        gkMPI_Send((void *)your_xadj, your_nvtxs+1, IDX_T, pe, 0, comm);
        gkMPI_Send((void *)your_vwgt, your_nvtxs*ncon, IDX_T, pe, 1, comm);
      }
      else {
        for (i=0; i<your_nvtxs+1; i++)
          xadj[i] = your_xadj[i];
        for (i=0; i<your_nvtxs*ncon; i++)
          vwgt[i] = your_vwgt[i];
      }
    }
    fclose(fpin);
    gk_free((void **)&your_xadj, &your_vwgt, LTERM);
  }
  else {
    gkMPI_Recv((void *)xadj, nvtxs+1, IDX_T, npes-1, 0, comm, &stat);
    gkMPI_Recv((void *)vwgt, nvtxs*ncon, IDX_T, npes-1, 1, comm, &stat);
  }

  graph->nedges = xadj[nvtxs];
  adjncy = graph->adjncy = imalloc(xadj[nvtxs], "ParallelReadGraph: adjncy");
  adjwgt = graph->adjwgt = imalloc(xadj[nvtxs]*nobj, "ParallelReadGraph: adjwgt");

  /***********************************************/
  /* Now go through again and record adjncy data */
  /***********************************************/
  if (mype == npes-1) {
    ier = 0;
    fpin = fopen(filename, "r");

    if (fpin == NULL){
      printf("COULD NOT OPEN FILE '%s' FOR SOME REASON!\n", filename);
      ier++;
    }

    gkMPI_Bcast(&ier, 1, IDX_T, npes-1, comm);
    if (ier > 0){
      gkMPI_Finalize();
      exit(0);
    }

    /* get first line again */
    while (fgets(line, MAXLINE, fpin) && line[0] == '%');

    your_adjncy = imalloc(maxnedges, "your_adjncy");
    your_adjwgt = ismalloc(maxnedges*nobj, 1, "your_adjwgt");

    for (pe=0; pe<npes; pe++) {
      your_nvtxs  = vtxdist[pe+1]-vtxdist[pe];
      your_nedges = 0;

      for (i=0; i<your_nvtxs; i++) {
        while (fgets(line, MAXLINE, fpin) && line[0] == '%');
        oldstr = line;
        newstr = NULL;

        if (readvw) {
          for (l=0; l<ncon; l++) {
            dummy  = strtoidx(oldstr, &newstr, 10);
            oldstr = newstr;
          }
        }

        for (;;) {
          edge   = strtoidx(oldstr, &newstr, 10) -1;
          oldstr = newstr;

          if (edge < 0)
            break;

          your_adjncy[your_nedges] = edge;
          if (readew) {
            for (l=0; l<nobj; l++) {
              your_adjwgt[your_nedges*nobj+l] = strtoidx(oldstr, &newstr, 10);
              oldstr = newstr;
            }
          }
          your_nedges++;
        }
      }
      if (pe < npes-1) {
        gkMPI_Send((void *)your_adjncy, your_nedges, IDX_T, pe, 0, comm);
        gkMPI_Send((void *)your_adjwgt, your_nedges*nobj, IDX_T, pe, 1, comm);
      }
      else {
        for (i=0; i<your_nedges; i++)
          adjncy[i] = your_adjncy[i];
        for (i=0; i<your_nedges*nobj; i++)
          adjwgt[i] = your_adjwgt[i];
      }
    }
    fclose(fpin);
    gk_free((void **)&your_adjncy, &your_adjwgt, &line, LTERM);
  }
  else {
    gkMPI_Bcast(&ier, 1, IDX_T, npes-1, comm);
    if (ier > 0){
      gkMPI_Finalize();
      exit(0);
    }

    gkMPI_Recv((void *)adjncy, xadj[nvtxs], IDX_T, npes-1, 0, comm, &stat);
    gkMPI_Recv((void *)adjwgt, xadj[nvtxs]*nobj, IDX_T, npes-1, 1, comm, &stat);
  }

}


/*************************************************************************
* This function writes out a partition vector
**************************************************************************/
void WritePVector(char *gname, idx_t *vtxdist, idx_t *part, MPI_Comm comm)
{
  idx_t i, j, k, l, rnvtxs, npes, mype, penum;
  FILE *fpin;
  idx_t *rpart;
  char partfile[256];
  MPI_Status status;

  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  if (mype == 0) {
    sprintf(partfile, "%s.part", gname);
    if ((fpin = fopen(partfile, "w")) == NULL) 
      errexit("Failed to open file %s", partfile);

    for (i=0; i<vtxdist[1]; i++)
      fprintf(fpin, "%"PRIDX"\n", part[i]);

    for (penum=1; penum<npes; penum++) {
      rnvtxs = vtxdist[penum+1]-vtxdist[penum];
      rpart = imalloc(rnvtxs, "rpart");
      MPI_Recv((void *)rpart, rnvtxs, IDX_T, penum, 1, comm, &status);

      for (i=0; i<rnvtxs; i++)
        fprintf(fpin, "%"PRIDX"\n", rpart[i]);

      gk_free((void **)&rpart, LTERM);
    }
    fclose(fpin);
  }
  else
    MPI_Send((void *)part, vtxdist[mype+1]-vtxdist[mype], IDX_T, 0, 1, comm); 

}


/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  idx_t i, j, npes, mype, optype, nparts, options[METIS_NOPTIONS];
  idx_t *part=NULL, *sizes=NULL;
  graph_t graph;
  real_t *tpwgts=NULL, ubvec[MAXNCON];
  MPI_Comm comm;
  idx_t wgtflag=0, ndims, edgecut;

  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  gkMPI_Comm_size(comm, &npes);
  gkMPI_Comm_rank(comm, &mype);

  if (argc != 5) {
    if (mype == 0)
      printf("Usage: %s <graph-file> <nparts> <dbglvl> <seed>\n", argv[0]);

    MPI_Finalize();
    exit(0);
  }

  nparts = atoi(argv[2]);

  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_DBGLVL] = atoi(argv[3]);
  options[METIS_OPTION_SEED]   = atoi(argv[4]);

  if (mype == 0) 
    printf("reading file: %s\n", argv[1]);
  ParallelReadGraph(&graph, argv[1], comm);

  rset(graph.ncon, 1.05, ubvec);
  tpwgts = rmalloc(nparts*graph.ncon, "tpwgts");
  rset(nparts*graph.ncon, 1.0/(real_t)nparts, tpwgts);

  if (mype == 0) 
    printf("finished reading file: %s\n", argv[1]);
  
  part  = ismalloc(graph.nvtxs, mype%nparts, "main: part");
  sizes = imalloc(2*npes, "main: sizes");

  wgtflag = 3;
  PM4GNN_PartKway(graph.vtxdist, graph.xadj, graph.adjncy, graph.vwgt, 
      graph.adjwgt, &wgtflag, &graph.ncon, &nparts, tpwgts, ubvec, 
      options, &edgecut, part, &comm);
  WritePVector(argv[1], graph.vtxdist, part, MPI_COMM_WORLD); 

  /* printf("%"PRIDX" %"PRIDX"\n", isum(nvtxs, graph.xadj, 1), isum(nedges, graph.adjncy, 1)); */

  gk_free((void **)&part, &sizes, &tpwgts, &graph.vtxdist, &graph.xadj, &graph.adjncy, 
         &graph.vwgt, &graph.adjwgt, LTERM);

  MPI_Comm_free(&comm);

  MPI_Finalize();

  return 0;
}

