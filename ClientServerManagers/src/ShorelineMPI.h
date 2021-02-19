/* @file ShorelineMPI.h
 * @brief Header file mpi convenience functions
 *
 */

#pragma once
// stl includes
#include <vector>
// mpi include
#include <mpi.h>

namespace Shoreline {

/* @brief gather std::vectors from several processes to one process
 * @param[in] sendVec = on-process vector to send out
 * @param[in] type = the MPI datatype to send. There are a few...
 * @param[in/out] recvVec = storage to receive the vector. Should be passed in
 *                          with length 0. This function will resize the vector
 *                          on the receiving process and not touch it on the
 *                          other processes
 * @param[in] recvRank = the process rank to receive the data
 * @param[in] comm = the MPI communicator to operate on
 */
 template <typename T>
 int gatherVectors(const std::vector<T>& sendVec,
   std::vector<T>& recvVec, int recvRank, MPI_Comm comm);

////////////////////////////////////////////////////////////////////////////////
int scatterField(const std::vector<double>& sendField,
  const std::vector<int>& sendCounts, std::vector<double>& recvField,
  const int sendRank, const MPI_Comm comm);
}
