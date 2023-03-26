/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h 10592 2011-07-16 21:17:53Z karypis $
 *
 */

/* ctrl.c */
ctrl_t *SetupCtrl(idx_t *options, idx_t ncon, idx_t nparts, real_t *tpwgts, 
            real_t *ubvec, MPI_Comm comm);
void SetupCtrl_invtvwgts(ctrl_t *ctrl, graph_t *graph);
void FreeCtrl(ctrl_t **r_ctrl);



/* kmetis.c */
int CheckInputsPartKway(idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
        idx_t *adjwgt, idx_t *wgtflag, idx_t *ncon, idx_t *nparts, real_t *tpwgts, 
        real_t *ubvec, idx_t *options, idx_t *edgecut, idx_t *part, MPI_Comm *comm);
void Global_Partition(ctrl_t *, graph_t *);
void PartitionSmallGraph(ctrl_t *, graph_t *);
graph_t *AssembleGraph(ctrl_t *, graph_t *);

/* match.c */
void Match_Global(ctrl_t *, graph_t *);
void CreateCoarseGraph_Global(ctrl_t *, graph_t *, idx_t);
void DropEdges(ctrl_t *ctrl, graph_t *graph);

/* initpart.c */
void InitPartition(ctrl_t *, graph_t *);
void KeepPart(ctrl_t *, graph_t *, idx_t *, idx_t);

/* kwayrefine.c */
void ProjectPartition(ctrl_t *, graph_t *);
void ComputePartitionParams(ctrl_t *, graph_t *);
void KWayFM(ctrl_t *, graph_t *, idx_t);
void KWayBalance(ctrl_t *, graph_t *, idx_t);

/* remap.c */
void ParallelReMapGraph(ctrl_t *, graph_t *);
void ParallelTotalVReMap(ctrl_t *, idx_t *, idx_t *, idx_t, idx_t);
idx_t SimilarTpwgts(real_t *, idx_t, idx_t, idx_t);

/* wspace.c */
void AllocateWSpace(ctrl_t *ctrl, size_t nwords);
void AllocateRefinementWorkSpace(ctrl_t *ctrl, idx_t nbrpoolsize);
void FreeWSpace(ctrl_t *);
void *wspacemalloc(ctrl_t *ctrl, size_t nbytes);
idx_t *iwspacemalloc(ctrl_t *ctrl, size_t n);
real_t *rwspacemalloc(ctrl_t *ctrl, size_t n);
ikv_t *ikvwspacemalloc(ctrl_t *ctrl, size_t n);
rkv_t *rkvwspacemalloc(ctrl_t *ctrl, size_t n);
void cnbrpoolReset(ctrl_t *ctrl);
idx_t cnbrpoolGetNext(ctrl_t *ctrl, idx_t nnbrs);

/* stat.c */
idx_t ComputeSerialEdgeCut(graph_t *graph);
void ComputeSerialBalance(ctrl_t *, graph_t *, idx_t *, real_t *);
void ComputeParallelBalance(ctrl_t *, graph_t *, idx_t *, real_t *);
void Mc_PrintThrottleMatrix(ctrl_t *, graph_t *, real_t *);
void PrintPostPartInfo(ctrl_t *ctrl, graph_t *graph, idx_t movestats);
void ComputeMoveStatistics(ctrl_t *, graph_t *, idx_t *, idx_t *, idx_t *);
void Mc_ComputeMoveStatistics(ctrl_t *, graph_t *, idx_t *, idx_t *, idx_t *);

/* debug.c */
void PrintVector(ctrl_t *, idx_t, idx_t, idx_t *, char *);
void PrintVector2(ctrl_t *, idx_t, idx_t, idx_t *, char *);
void PrintPairs(ctrl_t *, idx_t, ikv_t *, char *);
void PrintGraph(ctrl_t *, graph_t *);
void PrintGraph2(ctrl_t *, graph_t *);
void PrintSetUpInfo(ctrl_t *ctrl, graph_t *graph);
void PrintTransferedGraphs(ctrl_t *, idx_t, idx_t *, idx_t *, idx_t *, idx_t *, idx_t *);
void WriteMetisGraph(idx_t, idx_t *, idx_t *, idx_t *, idx_t *);


/* comm.c */
void CommSetup(ctrl_t *, graph_t *);
void CommUpdateNnbrs(ctrl_t *ctrl, idx_t nnbrs);
void CommInterfaceData(ctrl_t *ctrl, graph_t *graph, idx_t *data, idx_t *recvvector);
void CommChangedInterfaceData(ctrl_t *ctrl, graph_t *graph, idx_t nchanged,
         idx_t *changed, idx_t *data, ikv_t *sendpairs, ikv_t *recvpairs);
idx_t GlobalSEMax(ctrl_t *, idx_t);
idx_t GlobalSEMaxComm(MPI_Comm comm, idx_t value);
idx_t GlobalSEMin(ctrl_t *, idx_t);
idx_t GlobalSEMinComm(MPI_Comm comm, idx_t value);
idx_t GlobalSESum(ctrl_t *, idx_t);
idx_t GlobalSESumComm(MPI_Comm comm, idx_t value);
real_t GlobalSEMaxFloat(ctrl_t *, real_t);
real_t GlobalSEMinFloat(ctrl_t *, real_t);
real_t GlobalSESumFloat(ctrl_t *, real_t);

/* util.c */
void myprintf(ctrl_t *ctrl, char *f_str,...);
void rprintf(ctrl_t *ctrl, char *f_str,...);
void mypridx_tf(ctrl_t *, char *f_str,...);
void rpridx_tf(ctrl_t *, char *f_str,...);
idx_t BSearch(idx_t, idx_t *, idx_t);
void RandomPermute(idx_t, idx_t *, idx_t);
void FastRandomPermute(idx_t, idx_t *, idx_t);
idx_t ispow2(idx_t);
idx_t log2Int(idx_t);
void BucketSortKeysDec(idx_t, idx_t, idx_t *, idx_t *);
real_t BetterVBalance(idx_t, real_t *, real_t *, real_t *);
idx_t IsHBalanceBetterTT(idx_t, real_t *, real_t *, real_t *, real_t *);
idx_t IsHBalanceBetterFT(idx_t, real_t *, real_t *, real_t *, real_t *);
void GetThreeMax(idx_t, real_t *, idx_t *, idx_t *, idx_t *);
size_t rargmax_strd(size_t n, real_t *x, size_t incx);
size_t rargmin_strd(size_t n, real_t *x, size_t incx);
size_t rargmax2(size_t n, real_t *x);
real_t ravg(size_t n, real_t *x);
real_t rfavg(size_t n, real_t *x);

/* graph.c */
graph_t *SetupGraph(ctrl_t *ctrl, idx_t ncon, idx_t *vtxdist, idx_t *xadj,
             idx_t *vwgt, idx_t *adjncy, idx_t *adjwgt, idx_t wgtflag);
void SetupGraph_nvwgts(ctrl_t *ctrl, graph_t *graph);
graph_t *CreateGraph(void);
void InitGraph(graph_t *);
void FreeGraph(graph_t **graph);
void FreeNonGraphFields(graph_t *graph);
void FreeNonGraphNonSetupFields(graph_t *graph);
void FreeCommSetupFields(graph_t *graph);
void FreeInitialGraphAndRemap(graph_t **graph);
void graph_WriteToDisk(ctrl_t *ctrl, graph_t *graph);
void graph_ReadFromDisk(ctrl_t *ctrl, graph_t *graph);

/* timer.c */
void InitTimers(ctrl_t *);
void PrintTimingInfo(ctrl_t *);
void PrintTimer(ctrl_t *, timer, char *);

/* gkmpi.c */
int gkMPI_Comm_size(MPI_Comm comm, idx_t *size);
int gkMPI_Comm_rank(MPI_Comm comm, idx_t *rank);
int gkMPI_Get_count(MPI_Status *status, MPI_Datatype datatype,
        idx_t *count);
int gkMPI_Send(void *buf, idx_t count, MPI_Datatype datatype, idx_t dest,
        idx_t tag, MPI_Comm comm);
int gkMPI_Recv(void *buf, idx_t count, MPI_Datatype datatype,
        idx_t source, idx_t tag, MPI_Comm comm, MPI_Status *status);
int gkMPI_Isend(void *buf, idx_t count, MPI_Datatype datatype, idx_t dest,
        idx_t tag, MPI_Comm comm, MPI_Request *request);
int gkMPI_Irecv(void *buf, idx_t count, MPI_Datatype datatype,
        idx_t source, idx_t tag, MPI_Comm comm, MPI_Request *request);
int gkMPI_Wait(MPI_Request *request, MPI_Status *status);
int gkMPI_Waitall(idx_t count, MPI_Request *array_of_requests, 
        MPI_Status *array_of_statuses);
int gkMPI_Barrier(MPI_Comm comm);
int gkMPI_Bcast(void *buffer, idx_t count, MPI_Datatype datatype,
        idx_t root, MPI_Comm comm);
int gkMPI_Reduce(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, idx_t root, MPI_Comm comm);
int gkMPI_Allreduce(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int gkMPI_Scan(void *sendbuf, void *recvbuf, idx_t count,
        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);
int gkMPI_Allgather(void *sendbuf, idx_t sendcount,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, MPI_Comm comm);
int gkMPI_Alltoall(void *sendbuf, idx_t sendcount,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, MPI_Comm comm);
int gkMPI_Alltoallv(void *sendbuf, idx_t *sendcounts,
        idx_t *sdispls, MPI_Datatype sendtype, void *recvbuf, 
        idx_t *recvcounts, idx_t *rdispls, MPI_Datatype recvtype, 
        MPI_Comm comm);
int gkMPI_Allgatherv(void *sendbuf, idx_t sendcount, MPI_Datatype sendtype, 
        void *recvbuf, idx_t *recvcounts, idx_t *rdispls, 
        MPI_Datatype recvtype, MPI_Comm comm);
int gkMPI_Scatterv(void *sendbuf, idx_t *sendcounts, idx_t *sdispls,
        MPI_Datatype sendtype, void *recvbuf, idx_t recvcount,
        MPI_Datatype recvtype, idx_t root, MPI_Comm comm);
int gkMPI_Gatherv(void *sendbuf, idx_t sendcount, MPI_Datatype sendtype,
        void *recvbuf, idx_t *recvcounts, idx_t *displs, MPI_Datatype recvtype,
        idx_t root, MPI_Comm comm);
int gkMPI_Comm_split(MPI_Comm comm, idx_t color, idx_t key,
        MPI_Comm *newcomm);
int gkMPI_Comm_free(MPI_Comm *comm);
int gkMPI_Finalize();


