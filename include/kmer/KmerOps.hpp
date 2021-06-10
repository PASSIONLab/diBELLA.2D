// Created by Saliya Ekanayake on 10/15/19.

#ifndef DIBELLA_KMEROPS_HPP
#define DIBELLA_KMEROPS_HPP

#include <unordered_set>
#include "MatrixEntry.hpp"
#include "../Types.hpp"
#include "../ParallelOps.hpp"
#include "../Alphabet.hpp"
#include "../ScoreMat.hpp"
#include "../Utils.hpp"
#include "../Defines.hpp"
#include "../Buffer.h"
#include "../DistributedFastaData.hpp"
#include "Kmer.hpp"
#include <stdint.h>

namespace dibella {
  class KmerOps {
  public:
    // GGGG: PSpMat needs CombBLAS 
    static PSpMat<PosInRead>::MPI_DCCols GenerateA(uint64_t seq_count,
        std::shared_ptr<DistributedFastaData>& dfd, ushort k, ushort s,
        Alphabet &alph, const std::shared_ptr<ParallelOps>& parops,
        const std::shared_ptr<TimePod>& tp, int nthreads); //, std::unordered_set<Kmer, Kmer>& local_kmers);

  //   static uint64_t add_kmers(const char *seq, ushort len, uint64_t start_offset,
  //                             uint64_t end_offset_inclusive, ushort k, ushort s,
  //                             Alphabet &alp, uvec_64 &lcol_ids, std::vector<MatrixEntry> &lvals,
  //                             const std::shared_ptr<ParallelOps>& parops,
  //                             std::unordered_set<Kmer, Kmer>& local_kmers)
  //   {

  //     auto num_kmers = static_cast<ushort>((floor((len - k) * 1.0 / s)) + 1);

  //     // TODO: Saliya - this can be improved using bit operation
  //     ushort base = alp.size;
  //     uint64_t kcode = 0;
  //     uint64_t count = 0;
  //     char cap_c;

  //     std::unordered_set<uint64_t> kmers_in_sequence;

  //     for (uint64_t i = start_offset; i <= ((end_offset_inclusive - k) + 1); i += s)
  //     {
  //       kcode = 0;
  //       if (count == num_kmers) break;
  //       std::string kmer_str;

  //       for (uint64_t j = i; j < i + k; ++j)
  //       {
  //         /*! Efficient than using pow() */
  //         cap_c = *(seq + j);
  //         if (cap_c > 96 && cap_c < 123) {
  //           // small case character, so make it uppercase.
  //           cap_c = cap_c - 32;
  //         }
  //         kmer_str += cap_c;
  //         kcode = kcode * base + alp.char_to_code[cap_c];
  //       }

  //       ++count;
  //       lcol_ids.push_back(kcode);
  //       local_kmers.emplace(kmer_str, kcode, alp, false);

  //       /*! Offset is relative to the sequence start, so unsigned short is
  //        * good enough. */
  //       lvals.emplace_back(0, static_cast<ushort &&>(i - start_offset));
  //     }

  //     if (count != num_kmers)
  //     {
  //       fprintf(stderr,
  //               "ERROR: kmerop: add_kmers(): rank: %d, count:%d numk: %d len: %d k: %d s: %d soff: %llu eoffinc: %llu\n",
  //               parops->world_proc_rank, count, num_kmers, len, k, s,
  //               start_offset, end_offset_inclusive);
  //       fflush(stderr);
  //     }
  //     return count;
  //   }
  };
}

#endif // DIBELLA_KMEROPS_HPP
