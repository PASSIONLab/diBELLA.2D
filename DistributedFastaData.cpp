// Created by Saliya Ekanayake on 2019-02-17.

#include "DistributedFastaData.h"
#include "ParallelFastaReader.hpp"


DistributedFastaData::~DistributedFastaData() {
  for (auto &row_seq : row_seqs) {
    delete (row_seq);
  }

  for (auto &col_seq : col_seqs) {
    delete (col_seq);
  }

  for (int i = 0; i < recv_nbrs_count; ++i) {
    delete (recv_fds[i]);
  }
  delete[]recv_fds;
  delete[](to_nbrs_bufss_stat);
  delete[](to_nbrs_buffs_reqs);
  delete[](recv_nbrs_buffs_stats);
  delete[](recv_nbrs_buffs_reqs);
  for (int i = 0; i < recv_nbrs_count; ++i) {
    delete[](recv_nbrs_buffs[i]);
  }
  delete[](recv_nbrs_buffs);
  delete[](recv_nbrs_buff_lengths);
  delete (fd->buffer());
  delete (fd);
  delete[](l_seq_counts);
  delete[](g_seq_offsets);
}

DistributedFastaData::DistributedFastaData(
  const char *file, ushort overlap, ushort k,
  const std::shared_ptr<ParallelOps> &parops)
  : overlap(overlap), k(k), parops(parops) {

  char *buff;
  uint64_t l_start, l_end;
  ParallelFastaReader::read_fasta(file, overlap, parops->world_proc_rank,
                                  parops->world_procs_count, buff, l_start,
                                  l_end);

  fd = new FastaData(buff, k, l_start, l_end);

  l_seq_count = fd->local_count();
  l_seq_counts = new uint64_t[parops->world_procs_count];
  MPI_Allgather(&l_seq_count, 1, MPI_UINT64_T, l_seq_counts,
                1, MPI_UINT64_T, MPI_COMM_WORLD);

  g_seq_offsets = new uint64_t[parops->world_procs_count];
  g_seq_offsets[0] = 0;
  for (int i = 1; i < parops->world_procs_count; ++i) {
    g_seq_offsets[i] = g_seq_offsets[i - 1] + l_seq_counts[i - 1];
  }

  /*! Total sequence count. This could be different from the original
   * input sequence count as we remove sequences less than k-mer length.
   */
  g_seq_count = g_seq_offsets[parops->world_procs_count - 1] +
                l_seq_counts[parops->world_procs_count - 1];

  collect_grid_seqs();
}

void DistributedFastaData::collect_grid_seqs() {
  int pr = parops->grid->GetGridRows();
  int pc = parops->grid->GetGridCols();

  /*! It might seem wrong to use rank_in_col for row_id and rank_in_row for
   * col_id but it is correct. These are ranks when the column and row worlds
   * are viewed separately.
   */
  int grid_row_id = parops->grid->GetRankInProcCol();
  int grid_col_id = parops->grid->GetRankInProcRow();

  uint64_t avg_l_seq_count = g_seq_count / parops->world_procs_count;
  uint64_t avg_grid_seq_count = g_seq_count / pr;

  row_seq_start_idx = avg_grid_seq_count * grid_row_id;
  row_seq_end_idx = (grid_row_id == pr - 1)
                    ? g_seq_count : row_seq_start_idx + avg_grid_seq_count;
  // Subtract 1 to get the zero based inclusive end index;
  --row_seq_end_idx;

  col_seq_start_idx = avg_grid_seq_count * grid_col_id;
  col_seq_end_idx = (grid_col_id == pc - 1)
                    ? g_seq_count : col_seq_start_idx + avg_grid_seq_count;
  // Subtract 1 to get the zero based inclusive end index;
  --col_seq_end_idx;

  /*! Whose sequences will I need */

  // my_nbrs can include me as well. We'll, of course, not send
  // messages to self but I need to know what's needed from locally
  // available squences to satisfy the grid cell I am responsible for.

  find_nbrs(pr, grid_row_id, avg_l_seq_count, avg_grid_seq_count,
            row_seq_start_idx, 1, my_nbrs);
  if (grid_row_id != grid_col_id) {
    find_nbrs(pc, grid_col_id, avg_l_seq_count, avg_grid_seq_count,
              col_seq_start_idx, 0, my_nbrs);
  }

#ifndef NDEBUG
  {
    std::string title = "My neighbors";
    std::string msg;
    for (auto &nbr : my_nbrs) {
      msg += +" (rc_flag: " + std::to_string(nbr.rc_flag)
             + ", nbr_rank: " + std::to_string(nbr.nbr_rank)
             + ", nbr_seq_start_idx: " + std::to_string(nbr.nbr_seq_start_idx)
             + ", nbr_seq_end_idx: " + std::to_string(nbr.nbr_seq_end_idx) +
             ")\n";
    }
    DebugUtils::print_msg(title, msg, parops);
  }
#endif

  // Receive neighbors are my neighbors excluding me, so receive neighbor
  // indexes will keep track of the indexes of my neighbors to receive from
  recv_nbrs_count = 0;
  for (int i = 0; i < my_nbrs.size(); ++i) {
    if (my_nbrs[i].nbr_rank == parops->world_proc_rank) {
      // No need to receive from self
      continue;
    }
    recv_nbrs_idxs.push_back(i);
    ++recv_nbrs_count;
  }

  /*! Issue receive requests to get the buffer lengths from my neighbors
  * for the sequences I need from them.
  */
  auto *recv_nbrs_buff_lengths_reqs = new MPI_Request[recv_nbrs_count];
  recv_nbrs_buff_lengths = new uint64_t[recv_nbrs_count];
  for (int i = 0; i < recv_nbrs_count; ++i) {
    MPI_Irecv(recv_nbrs_buff_lengths + i, 1, MPI_UINT64_T,
              my_nbrs[recv_nbrs_idxs[i]].nbr_rank, 99,
              MPI_COMM_WORLD, recv_nbrs_buff_lengths_reqs + i);
  }

  /*! Who will need my sequences */
  int block_length[5] = {1, 1, 1, 1, 1};
  MPI_Aint displacement[5] = {offsetof(NbrData, rc_flag),
                              offsetof(NbrData, owner_rank),
                              offsetof(NbrData, nbr_rank),
                              offsetof(NbrData, nbr_seq_start_idx),
                              offsetof(NbrData, nbr_seq_end_idx)};
  MPI_Datatype types[] = {MPI_UNSIGNED_SHORT, MPI_INT, MPI_INT, MPI_UINT64_T,
                          MPI_UINT64_T};
  MPI_Datatype MPI_NbrData;
  MPI_Type_create_struct(5, block_length, displacement, types, &MPI_NbrData);
  MPI_Type_commit(&MPI_NbrData);

  /*! Number of neighbors will surely fit in int
   * and moreover we need this to be int to be used with
   * allgatherv, which didn't like when counts and displacements
   * were uint_64 arrays */
  int my_nbrs_count = static_cast<int>(my_nbrs.size());
  auto *all_nbrs_counts = new int[parops->world_procs_count];
  MPI_Allgather(&my_nbrs_count, 1, MPI_INT, all_nbrs_counts, 1, MPI_INT,
                MPI_COMM_WORLD);

  int all_nbrs_count = 0;
  auto *all_nbrs_displas = new int[parops->world_procs_count];
  all_nbrs_displas[0] = 0;
  for (int i = 0; i < parops->world_procs_count; ++i) {
    all_nbrs_count += all_nbrs_counts[i];
    if (i > 0) {
      all_nbrs_displas[i] = all_nbrs_displas[i - 1] + all_nbrs_counts[i - 1];
    }
  }

  auto *all_nbrs = new NbrData[all_nbrs_count];
  MPI_Allgatherv(&my_nbrs[0], my_nbrs_count, MPI_NbrData, all_nbrs,
                 all_nbrs_counts, all_nbrs_displas, MPI_NbrData,
                 MPI_COMM_WORLD);

#ifndef NDEBUG
  {
    std::string title = "All neighbor counts and displacements at the root";
    std::string msg;
    for (int i = 0; i < parops->world_procs_count; ++i) {
      msg += "(rank: " + std::to_string(i)
             + " count: " + std::to_string(all_nbrs_counts[i])
             + " displ:" + std::to_string(all_nbrs_displas[i]) + ")\n";
    }
    DebugUtils::print_msg_on_rank(title, msg, parops, 0);
  }
#endif

#ifndef NDEBUG
  {
    std::string title = "Allgathered neighbors at the root";
    std::string msg;
    int rank = 0;
    for (int i = 0; i < all_nbrs_count; ++i) {
      NbrData *n = all_nbrs + i;
      msg += "rank: " + std::to_string(rank)
             + " (rc_flag: " + std::to_string(n->rc_flag)
             + ", owner_rank: " + std::to_string(n->owner_rank)
             + ", nbr_rank: " + std::to_string(n->nbr_rank)
             + ", nbr_seq_start_idx: " + std::to_string(n->nbr_seq_start_idx)
             + ", nbr_seq_end_idx: " + std::to_string(n->nbr_seq_end_idx) +
             ")\n";
      if (i + 1 == all_nbrs_displas[rank + 1]) {
        ++rank;
      }
    }
    DebugUtils::print_msg_on_rank(title, msg, parops, 0);
  }
#endif

  std::vector<uint64_t> send_lengths;
  std::vector<uint64_t> send_start_offsets;
  std::vector<int> to_nbrs_idxs;
  to_nbrs_count = 0;
  std::sort(all_nbrs, all_nbrs + all_nbrs_count,
            [](const NbrData &a, const NbrData &b) -> bool {
              return a.nbr_rank < b.nbr_rank;
            });
  for (int i = 0; i < all_nbrs_count; ++i) {
    if (all_nbrs[i].nbr_rank != parops->world_proc_rank) {
      // Not requesting my sequences, so continue.
      continue;
    }

    if (all_nbrs[i].owner_rank == parops->world_proc_rank) {
      // Available locally (me requesting my sequences), so continue.
      continue;
    }
    NbrData nbr = all_nbrs[i];
    uint64_t len, start_offset, end_offset;
    fd->buffer_size(nbr.nbr_seq_start_idx, nbr.nbr_seq_end_idx,
                    len, start_offset, end_offset);

    send_lengths.push_back(len);
    send_start_offsets.push_back(start_offset);
    to_nbrs_idxs.push_back(i);
    ++to_nbrs_count;
  }

  auto *to_nbrs_send_reqs = new MPI_Request[to_nbrs_count];
  for (int i = 0; i < to_nbrs_count; ++i) {
    MPI_Isend(&send_lengths[i], 1, MPI_UINT64_T,
              all_nbrs[to_nbrs_idxs[i]].owner_rank, 99, MPI_COMM_WORLD,
              to_nbrs_send_reqs + i);
  }

  /* Wait for length transfers */
  auto recv_stats = new MPI_Status[recv_nbrs_count];
  auto send_stats = new MPI_Status[to_nbrs_count];
  MPI_Waitall(recv_nbrs_count, recv_nbrs_buff_lengths_reqs, recv_stats);
  MPI_Waitall(to_nbrs_count, to_nbrs_send_reqs, send_stats);

#ifndef NDEBUG
  {
    std::string title = "My recv neighbors (excludes me) buff lengths";
    std::string msg;
    for (int i = 0; i < recv_nbrs_count; ++i) {
      NbrData nbr = my_nbrs[recv_nbrs_idxs[i]];
      msg += +" (rc_flag: " + std::to_string(nbr.rc_flag)
             + ", nbr_rank: " + std::to_string(nbr.nbr_rank)
             + ", nbr_seq_start_idx: " + std::to_string(nbr.nbr_seq_start_idx)
             + ", nbr_seq_end_idx: " + std::to_string(nbr.nbr_seq_end_idx) +
             +", nbr buff length: " +
             std::to_string(recv_nbrs_buff_lengths[i]) +
             ")\n";
    }
    DebugUtils::print_msg(title, msg, parops);
  }
#endif

  // TODO - do multiple send/recvs below if buff length > int max

#ifndef NDEBUG
  {
    std::string title = "recv_nbrs and to_nbrs count";
    std::string msg = "(rnc: " + std::to_string(recv_nbrs_count) + " tnc: " +
                      std::to_string(to_nbrs_count) + ")";
    DebugUtils::print_msg(title, msg, parops);
  }
#endif

  recv_nbrs_buffs = new char *[recv_nbrs_count];
  for (int i = 0; i < recv_nbrs_count; ++i) {
    recv_nbrs_buffs[i] = new char[recv_nbrs_buff_lengths[i]];
  }
  recv_nbrs_buffs_reqs = new MPI_Request[recv_nbrs_count];
  recv_nbrs_buffs_stats = new MPI_Status[recv_nbrs_count];
  for (int i = 0; i < recv_nbrs_count; ++i) {
    MPI_Irecv(recv_nbrs_buffs[i], static_cast<int>(recv_nbrs_buff_lengths[i]),
              MPI_CHAR, my_nbrs[recv_nbrs_idxs[i]].nbr_rank, 77, MPI_COMM_WORLD,
              recv_nbrs_buffs_reqs + i);
  }

  to_nbrs_buffs_reqs = new MPI_Request[to_nbrs_count];
  to_nbrs_bufss_stat = new MPI_Status[to_nbrs_count];
  for (int i = 0; i < to_nbrs_count; ++i) {
    MPI_Isend(fd->buffer() + send_start_offsets[i],
              static_cast<int>(send_lengths[i]), MPI_CHAR,
              all_nbrs[to_nbrs_idxs[i]].owner_rank, 77, MPI_COMM_WORLD,
              to_nbrs_buffs_reqs + i);
  }


  delete[](to_nbrs_send_reqs);
  delete[](all_nbrs);
  delete[](all_nbrs_displas);
  delete[](all_nbrs_counts);
  MPI_Type_free(&MPI_NbrData);
  delete[](recv_nbrs_buff_lengths_reqs);

}

void DistributedFastaData::find_nbrs(const int grid_rc_procs_count,
                                     const int grid_rc_id,
                                     const uint64_t avg_l_seq_count,
                                     const uint64_t avg_grid_seq_count,
                                     const uint64_t rc_seq_start_idx,
                                     const ushort rc_flag,
                                     std::vector<NbrData> &my_nbrs) {
  // Note, this rank might not have the sequence, if so we'll search further.
  int start_rank
    = static_cast<int>((rc_seq_start_idx + 1) / avg_l_seq_count);
  while (g_seq_offsets[start_rank] > rc_seq_start_idx) {
    /*! This loop is unlikely to happen unless the ranks above the original
     * start rank were heavily pruned during sequence removal step, which
     * removes sequences in lenght less than k-mer length.
     */
    --start_rank;
  }

  assert(start_rank >= 0 && start_rank < parops->world_procs_count);

  while (g_seq_offsets[start_rank] + l_seq_counts[start_rank] <
         rc_seq_start_idx) {
    /* Another highly unlikely loop. This could happen if the above start_rank
     * contains much less number of sequences than the avg_l_seq_count due to
     * sequence pruning, which means we have to search in the ranks ahead of it
     * to find the sequence we are interested in. */
    ++start_rank;
  }

  assert(start_rank >= 0 && start_rank < parops->world_procs_count);

  uint64_t rc_seq_count_needed = avg_grid_seq_count;
  if (grid_rc_id == grid_rc_procs_count - 1) {
    /*! Last rank in the row or column world, so it might
     * need more than the avg_grid_seq_count */
    rc_seq_count_needed = g_seq_count - (avg_grid_seq_count * grid_rc_id);
  }

  int nbr_rank;
  uint64_t nbr_seq_start_idx, nbr_seq_end_idx;
  uint64_t count = 0;
  uint64_t seq_start_idx = rc_seq_start_idx;
  while (count < rc_seq_count_needed) {
    nbr_rank = start_rank;
    nbr_seq_start_idx = seq_start_idx - g_seq_offsets[start_rank];
    /* See how many sequences to grab from start_rank */
    uint64_t remaining_needed = rc_seq_count_needed - count;
    /* remaining_in_start_rank includes the current
     * seq pointed by rc_seq_start_idx */
    uint64_t remaining_in_start_rank = l_seq_counts[start_rank] -
                                       (seq_start_idx -
                                        g_seq_offsets[start_rank]);
    if (remaining_needed >= remaining_in_start_rank) {
      count += remaining_in_start_rank;
      nbr_seq_end_idx = ((seq_start_idx + remaining_in_start_rank - 1) -
                         g_seq_offsets[start_rank]);
    } else {
      count += remaining_needed;
      nbr_seq_end_idx =
        (seq_start_idx + remaining_needed - 1) - g_seq_offsets[start_rank];
    }
    my_nbrs.emplace_back(rc_flag, parops->world_proc_rank, nbr_rank,
                         nbr_seq_start_idx, nbr_seq_end_idx);
    ++start_rank;
    seq_start_idx = g_seq_offsets[start_rank];
  }

  assert(count == rc_seq_count_needed);
}

uint64_t DistributedFastaData::global_count() {
  return g_seq_count;
}

bool DistributedFastaData::is_ready() {
  return ready;
}

void
DistributedFastaData::push_seqs(int rc_flag, FastaData *fd, uint64_t seqs_count,
                                uint64_t seq_start_idx) {
  ushort len;
  uint64_t start_offset, end_offset_inclusive;
  for (uint64_t i = 0; i < seqs_count; ++i) {
    char *buff = fd->get_sequence(seq_start_idx + i, len, start_offset,
                                  end_offset_inclusive);
    seqan::Peptide *seq = new seqan::Peptide(buff + start_offset, len);
    if (rc_flag == 1) {
      /*! grid row sequence */
      row_seqs.push_back(seq);
    } else {
      /*! grid col sequence */
      col_seqs.push_back(seq);
    }
  }
}

void DistributedFastaData::wait() {
  MPI_Waitall(recv_nbrs_count, recv_nbrs_buffs_reqs, recv_nbrs_buffs_stats);
  MPI_Waitall(to_nbrs_count, to_nbrs_buffs_reqs, to_nbrs_bufss_stat);

#ifndef NDEBUG
  {
    std::string title = "First character of received data from each neighbor";
    std::string msg;
    for (int i = 0; i < recv_nbrs_count; ++i) {
      msg += std::string(recv_nbrs_buffs[i], 1);
    }
    DebugUtils::print_msg(title, msg, parops);
  }
#endif

  int recv_nbr_idx = 0;
  for (auto &nbr : my_nbrs) {
    uint64_t nbr_seqs_count = (nbr.nbr_seq_end_idx - nbr.nbr_seq_start_idx) + 1;
    if (nbr.nbr_rank == parops->world_proc_rank) {
      /*! Local data, so create SeqAn sequences from <tt>fd</tt> */
      push_seqs(nbr.rc_flag, fd, nbr_seqs_count, nbr.nbr_seq_start_idx);
    } else {
      /*! Foreign data in a received buffer, so create a FastaData instance
       * and create SeqAn sequences from it. For foreign data, char buffers
       * only contain the required sequences, so sequence start index is 0.
       */
      uint64_t recv_nbr_l_end = recv_nbrs_buff_lengths[recv_nbr_idx] - 1;
      FastaData *recv_fd = new FastaData(recv_nbrs_buffs[recv_nbr_idx], k, 0,
                                         recv_nbr_l_end);
      push_seqs(nbr.rc_flag, recv_fd, nbr_seqs_count, 0);
      ++recv_nbr_idx;
    }
  }

  if (parops->grid->GetRankInProcRow() == parops->grid->GetRankInProcCol()) {
    /*! Diagonal cell, so col_seqs should point to the same sequences
     * as row_seqs. Also, it should not have received any col sequences
     * at this point */
    assert(col_seqs.size() == 0);
    col_seqs.assign(row_seqs.begin(), row_seqs.end());
  }

  assert(row_seqs.size() == (row_seq_end_idx - row_seq_start_idx) + 1 &&
         col_seqs.size() == (col_seq_end_idx - col_seq_start_idx) + 1);

  ready = true;
}
