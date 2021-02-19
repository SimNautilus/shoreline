#include "ShorelineMPI.h"

// stl includes
#include <numeric>
#include <iostream>
#include <algorithm>

namespace Shoreline {
template <typename T>
int gatherVectors(const std::vector<T>& sendVec,
  std::vector<T>& recvVec, int recvRank, MPI_Comm comm)
{
  MPI_Datatype type = MPI_INT;
  if(std::is_same<int, T>::value) { type = MPI_INT; }
  if(std::is_same<long, T>::value) { type = MPI_LONG; }
  if(std::is_same<double, T>::value) { type = MPI_DOUBLE; }
  int commSize;
  MPI_Comm_size(comm, &commSize);
  // Get the rank of the process
  int commRank;
  MPI_Comm_rank(comm, &commRank);

  // get how many elements are getting sent
  int lengthOnProc = sendVec.size();
  std::vector<int> lengths((commRank==recvRank) ? commSize : 0);
  MPI_Gather(&lengthOnProc, 1, MPI_INT, lengths.data(), 1, MPI_INT, 0, comm);

  // get the total amount of adata being received at recvRank
  // this value will be zero on all other ranks
  int nTotalData = std::accumulate(lengths.begin(), lengths.end(), 0);

  // resize the vector where the data will be received
  recvVec.resize(nTotalData);

  // setup strides for each process
  std::vector<int> recvOffsets(lengths.size());
  std::partial_sum(lengths.begin(), lengths.end(), recvOffsets.begin());
  std::transform(recvOffsets.begin(), recvOffsets.end(), lengths.begin(),
    recvOffsets.begin(), std::minus<int>());

  return MPI_Gatherv(sendVec.data(), lengthOnProc, type,
    recvVec.data(), lengths.data(), recvOffsets.data(), type, recvRank, comm);
}

template int gatherVectors<int>(const std::vector<int>& sendVec,
  std::vector<int>& recvVec, int recvRank, MPI_Comm comm);
template int gatherVectors<double>(const std::vector<double>& sendVec,
  std::vector<double>& recvVec, int recvRank, MPI_Comm comm);
template int gatherVectors<long>(const std::vector<long>& sendVec,
  std::vector<long>& recvVec, int recvRank, MPI_Comm comm);
////////////////////////////////////////////////////////////////////////////////

int scatterField(const std::vector<double>& sendField,
  const std::vector<int>& sendCounts, std::vector<double>& recvField,
  const int sendRank, const MPI_Comm comm)
{
  // resize storage on each process to receive data
  int nDataOnProc = 0;
  MPI_Scatter(sendCounts.data(),1,MPI_INT,&nDataOnProc,1,MPI_INT,sendRank,comm);
  recvField.resize(nDataOnProc);

  // get the offset within the send vector to the data for each process
  std::vector<int> sendOffsets(sendCounts.size());
  std::partial_sum(sendCounts.begin(), sendCounts.end(), sendOffsets.begin());
  std::transform(sendOffsets.begin(), sendOffsets.end(), sendCounts.begin(),
    sendOffsets.begin(), std::minus<int>());

  // send data to processes
  return MPI_Scatterv(sendField.data(), sendCounts.data(), sendOffsets.data(),
    MPI_DOUBLE, recvField.data(), nDataOnProc, MPI_DOUBLE, sendRank, comm);
}
}
