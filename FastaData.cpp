// Created by Saliya Ekanayake on 1/7/19.

#include <iostream>
#include <mpi.h>
#include <cassert>
#include "FastaData.hpp"

FastaData::~FastaData() = default;

FastaData::FastaData(std::shared_ptr<char> data, uint64_t l_start,
                     uint64_t l_end, ushort k) {

  id_starts = new uvec_64();
  seq_starts = new uvec_64();

  l_seq_count = 0;
  char c;
  bool in_name = false;
  bool in_seq = false;
  ushort seq_len=0;
  /*! No character count. This includes new line and * characters
   * It also includes entire sequences that are less than k-mer length*/
  ushort nc_count = 0;
  uint64_t idx;
  char *buff = data.get();
  /*! Assume the FASTA content is valid */
  for (uint64_t i = l_start; i <= l_end; ++i) {
    c = buff[i];
    idx = i - nc_count;
    buff[idx] = c;
    if (c == '>') {
      id_starts->push_back(idx);
      seq_len=0;
      in_name = true;
      in_seq = false;
      ++l_seq_count;
    } else if (c == '\n'){
      if (in_name && i + 1 <= l_end){
        seq_starts->push_back(idx+1);
        in_name = false;
        in_seq = true;
      } else if (in_seq && i + 1 <= l_end){
        if (buff[i+1] != '>'){
          ++nc_count;
        } else if (buff[i+1] == '>'){
          if (seq_len < k){
            ++seq_len; // capture the new line character too for removal
            uint64_t seq_id_start = id_starts->back();
            uint64_t seq_start = seq_starts->back();
            /*! + 1 to capture the new line between
             * sequence id and sequence start */
            nc_count += seq_len + (seq_start - seq_id_start) + 1;
            id_starts->pop_back();
            seq_starts->pop_back();
            --l_seq_count;
          }
        }
      }
    } else if (c == '*'){
      if (in_seq){
        ++nc_count;
      }
    } else {
      if (in_seq){
        ++seq_len;
      }
    }
  }

  // Remove the last sequence as well if it's shorter than k
  if (seq_len < k){
    // The last sequence doesn't have a new line at the end
    // as our l_end points to the last real character of the block.
    uint64_t seq_id_start = id_starts->back();
    uint64_t seq_start = seq_starts->back();
    // + 1 to capture the new line between sequence id and sequence start
    nc_count += seq_len + (seq_start - seq_id_start) + 1;
    id_starts->pop_back();
    seq_starts->pop_back();
    --l_seq_count;
  }

  /*! Things can go wrong if you don't end up having at least one sequence,
   * which is unlikely unless the total number of sequences are close to
   * the total number of processes.
   */
  assert(l_seq_count > 0);

  MPI_Allreduce(&l_seq_count, &g_seq_count, 1,
                MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

  this->data = data;
  this->l_start = l_start;
  this->l_end = l_end - nc_count;

  g_seq_offset = 0;
  /*! Use MPI's exclusive scan collective to compute partial prefix sums */
  MPI_Exscan(&l_seq_count, &g_seq_offset, 1,
             MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
}

void FastaData::print() {
  uint64_t idx;
  bool is_global = false;
  ushort len;
  uint64_t start_offset;
  uint64_t end_offset_inclusive;
  std::cout<<"g_seq_offset "<<g_seq_offset<<std::endl;
  for (uint64_t i = 0; i < id_starts->size(); ++i){
    char *beg = data.get() + (*id_starts)[i];
    char *end = data.get() + ((*seq_starts)[i] - 1);
    std::cout.write(beg, end - beg);
    std::cout<<std::endl;
    std::shared_ptr<char> data = get_sequence(i, is_global, len, start_offset, end_offset_inclusive);
    beg = data.get() + start_offset;
    end = data.get() + end_offset_inclusive;
    std::cout.write(data.get()+start_offset, (end - beg)+1);
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
}

std::shared_ptr<char> FastaData::get_sequence(
    uint64_t idx, bool is_global, ushort &len,
    uint64_t &start_offset, uint64_t &end_offset_inclusive) {
  uint64_t l_idx =  is_global ? (idx - g_seq_offset) : idx;
  /*! ((*id_starts)[l_idx+1] - 1) points to the position of the newline
   * character in the l_idx's sequence content */
  len = static_cast<ushort>((l_idx + 1 < id_starts->size()
        ? ((*id_starts)[l_idx+1] - 1)
        : l_end+1) - (*seq_starts)[l_idx]);
  start_offset = (*seq_starts)[l_idx];
  end_offset_inclusive =
      (l_idx+1 < id_starts->size() ? ((*id_starts)[l_idx+1] - 2) : l_end);
  return data;
}

uint64_t FastaData::global_count() {
  return g_seq_count;
}

uint64_t FastaData::count() {
  return l_seq_count;
}

uint64_t FastaData::offset() {
  return g_seq_offset;
}



