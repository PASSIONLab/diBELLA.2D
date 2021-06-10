// Created by Saliya Ekanayake on 2019-07-05 and modified by Giulia Guidi on 09/01/2020.

#include "../../include/pw/SeedExtendXdrop.hpp"

uint minOverlapLen = 10000;

void SeedExtendXdrop::PostAlignDecision(const AlignmentInfo& ai, bool& passed, float& ratioScoreOverlap, 
	uint32_t& overhang, uint32_t& overhangT, uint32_t& overlap, const bool noAlign)
{
	auto maxseed = ai.seed;	// returns a seqan:Seed object

	// {begin/end}Position{V/H}: Returns the begin/end position of the seed in the query (vertical/horizonral direction)
	// these four return seqan:Tposition objects
	int begpV = beginPositionV(maxseed);
	int endpV = endPositionV  (maxseed);
	int begpH = beginPositionH(maxseed);
	int endpH = endPositionH  (maxseed);

	unsigned short int overlapLenH = ai.seq_h_seed_length;
	unsigned short int overlapLenV = ai.seq_v_seed_length;

	unsigned short int rlenH = ai.seq_h_length;
	unsigned short int rlenV = ai.seq_v_length;

	unsigned short int minLeft  = min(begpV, begpH);
	unsigned short int minRight = min(rlenV - endpV, rlenH - endpH);

	overlap = minLeft + minRight + (overlapLenV + overlapLenH) / 2;

#ifndef FIXEDTHR
	float myThr = (1 - DELTACHERNOFF) * (ratioScoreOverlap * (float)overlap);

	// Contained overlaps removed for now, reintroduce them later
	// @GGGG-TODO: identify chimeric sequences
	bool contained = false;
	bool chimeric  = false; 

	// seqH is column entry and seqV is row entry, for each column, we iterate over seqHs
	int seqH = ai.seq_v_g_idx, seqV = ai.seq_h_g_idx;
	
	// Reserve length/position if rc [x]
	if(ai.rc)
	{
		uint tmp = begpV;
		begpV = rlenV - endpV;
		endpV = rlenV - tmp;
	}

	if((begpH == 0 & rlenH-endpH == 0) || (begpV == 0 & rlenV-endpV == 0))
		contained = true;
	
	if(!contained)
	{
		// If noAlign is false, set passed to false if the score isn't good enough
		if(!noAlign)
		{
			if((float)ai.xscore < myThr || overlap < minOverlapLen) passed = false;
			else passed = true;
		}

		if(passed)
		{
			uint32_t direction, directionT;
			uint32_t suffix, suffixT;

			// !reverse complement
			if(!ai.rc)
			{
				if(begpH > begpV)
				{
					direction  = 1;
					directionT = 2;

					suffix  = rlenV - endpV;
					suffixT = begpH;
				}	
				else
				{
					direction  = 2;
					directionT = 1;

					suffix  = rlenH - endpH;
					suffixT = begpV;
				} 
			}
			else
			{
				if((begpV > 0) & (begpH > 0) & (rlenV-endpV == 0) & (rlenV-endpV == 0))
				{
					direction  = 0;
					directionT = 0;

					suffix  = begpV; 	// seqV == 2
					suffixT = begpH; 	// seqV == 2			
				}
				else
				{
					direction  = 3;
					directionT = 3;

					suffix  = rlenV - endpV;	// seqV == 1, seqH == 2	
					suffixT = rlenH - endpH;	// seqV == 1, seqH == 2		
				}
			}
			overhang  = suffix  << 2 | direction;
			overhangT = suffixT << 2 | directionT;
		} // if(passed)
	} // if(!contained)
		
#else
	if(ai.xscore >= FIXEDTHR)
		passed = true;
#endif
}

SeedExtendXdrop::SeedExtendXdrop(
    ScoringScheme scoring_scheme,
    ushort seed_length, int xdrop, int seed_count):
    PairwiseFunction(),
    scoring_scheme(scoring_scheme),
    seed_length(seed_length), xdrop(xdrop), seed_count(seed_count){
}

void SeedExtendXdrop::apply(
    uint64_t l_col_idx, uint64_t g_col_idx,
    uint64_t l_row_idx, uint64_t g_row_idx,
    seqan::Dna5String *seqH, seqan::Dna5String *seqV, ushort k,
    dibella::CommonKmers &cks, std::stringstream& ss)
{
  AlignmentInfo ai;

  for (int count = 0; count < seed_count; ++count)
  {
#ifdef TWOSEED
    // row sequence is the same thing as vertical sequence
    ushort LocalSeedVOffset = (count == 0) ? cks.first.first
                                                  : cks.second.first;
	// l_row_seed_start_offset = LocalSeedVOffset											  
    // col sequence is the same thing as horizontal sequence
    ushort LocalSeedHOffset = (count == 0) ? cks.first.second 
                                                  : cks.second.second;
	// l_col_seed_start_offset = LocalSeedHOffset
#else
    // row sequence is the same thing as vertical sequence
    ushort LocalSeedVOffset = cks.pos[0].first;
    // col sequence is the same thing as horizontal sequence
    ushort LocalSeedHOffset = cks.pos[0].second;
#endif

	seqan::Dna5String seedH;
	seqan::Dna5String seedV;

	auto start_time = std::chrono::system_clock::now();
	auto end_time   = std::chrono::system_clock::now();

    // Seed creation params are:
    // horizontal seed start offset, vertical seed start offset, length
    TSeed seed(LocalSeedHOffset, LocalSeedVOffset, seed_length);

	seedH = infix(*seqH, beginPositionH(seed), endPositionH(seed));
	seedV = infix(*seqV, beginPositionV(seed), endPositionV(seed));

	seqan::Dna5StringReverseComplement twin(seedH);

	seqan::Align<seqan::Dna5String> align;
	resize(rows(align), 2);

	if(twin == seedV)
	{
		ai.rc = true;
		seqan::Dna5String twinseqH = *seqH;
		seqan::Dna5StringReverseComplement twinRead(twinseqH);
		LocalSeedHOffset = length(twinseqH) - LocalSeedHOffset - seed_length;

		setBeginPositionH(seed, LocalSeedHOffset);
		setBeginPositionV(seed, LocalSeedVOffset);
		setEndPositionH(seed, LocalSeedHOffset + seed_length);
		setEndPositionV(seed, LocalSeedVOffset + seed_length);

		/* Perform match extension */
		start_time = std::chrono::system_clock::now();
		ai.xscore  = extendSeed(seed, twinseqH, *seqV, seqan::EXTEND_BOTH, scoring_scheme, xdrop, (int)k, seqan::GappedXDrop());
		end_time   = std::chrono::system_clock::now();
		add_time("XA:ExtendSeed", (ms_t(end_time - start_time)).count());

	#ifdef STATS
		assignSource(row(align, 0), infix(twinseqH, beginPositionH(seed),
										endPositionH(seed)));
		assignSource(row(align, 1), infix(*seqV, beginPositionV(seed),
										endPositionV(seed)));
	#endif

	}
	else
	{
		ai.rc = false;
		start_time = std::chrono::system_clock::now();
		ai.xscore = extendSeed(seed, *seqH, *seqV, seqan::EXTEND_BOTH, scoring_scheme, xdrop, (int)k, seqan::GappedXDrop());
		end_time = std::chrono::system_clock::now();
		add_time("XA:ExtendSeed", (ms_t(end_time - start_time)).count());

	#ifdef STATS
		assignSource(row(align, 0), infix(*seqH, beginPositionH(seed),
										endPositionH(seed)));
		assignSource(row(align, 1), infix(*seqV, beginPositionV(seed),
										endPositionV(seed)));
	#endif
	} 

    /*! Note. This aligns the extended seeds globally, NOT the original
     * two sequences.
     *
     * It seems kind of a waste to have to do the alignment
     * again after xdrop seed extension but that's the only
     * way to get the alignment info in SeqAn.
     * See https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/SeedExtension.html
     */
#ifdef STATS
    start_time = std::chrono::system_clock::now();
    globalAlignment(align, scoring_scheme);
    end_time = std::chrono::system_clock::now();
    add_time("XA:global_alignment", (ms_t(end_time - start_time)).count());

    // Compute the statistics of the alignment.
    start_time = std::chrono::system_clock::now();
    computeAlignmentStats(ai[count].stats, align, scoring_scheme);
    end_time = std::chrono::system_clock::now();
    add_time("XA:ComputeStats", (ms_t(end_time - start_time)).count());
#endif

    ai.seq_h_length = length(*seqH);
    ai.seq_v_length = length(*seqV);

	ai.seed = seed;
    ai.seq_h_seed_length = static_cast<ushort>(seed._endPositionH -
                                               seed._beginPositionH);
    ai.seq_v_seed_length = static_cast<ushort>(seed._endPositionV -
                                               seed._beginPositionV);
    ai.seq_h_g_idx = g_col_idx;
    ai.seq_v_g_idx = g_row_idx;
  }

#ifdef STATS
  if (seed_count > 2)
  {
    max_ai = ai[0].stats.alignmentIdentity > ai[1].stats.alignmentIdentity
    ? ai[0] : ai[1];
  }
  double alen_minus_gapopens = (max_ai.stats.alignmentLength - max_ai.stats.numGapOpens) * 1.0;
  ss << g_col_idx << "," << g_row_idx << "," << max_ai.stats.alignmentIdentity
     << "," << max_ai.seq_h_length << "," << max_ai.seq_v_length
     << "," << max_ai.seq_h_seed_length  << "," << max_ai.seq_v_seed_length
     << "," << max_ai.stats.numGapOpens
     << "," << alen_minus_gapopens / max_ai.seq_h_length
     << "," << alen_minus_gapopens / max_ai.seq_v_length
	 << "," << cks.count
	 << std::endl;
#endif
}

// @NOTE This is hard-coded to the number of seeds being <= 2
void
SeedExtendXdrop::apply_batch
(
    seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsh,
	seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsv,
	uint64_t *lids,
	uint64_t col_offset,
	uint64_t row_offset,
    PSpMat<dibella::CommonKmers>::ref_tuples *mattuples,
    std::ofstream &lfs,
	const bool noAlign,
	ushort k,
	uint64_t nreads,
    float ratioScoreOverlap, // GGGG: this is my ratioScoreOverlap variable change name later
    int debugThr
)
{
	seqan::ExecutionPolicy<seqan::Parallel, seqan::Vectorial> exec_policy;

	int numThreads = 1;
	#ifdef THREADED
	#pragma omp parallel
    {
      	numThreads = omp_get_num_threads();
    }
	#endif

	uint64_t npairs = seqan::length(seqsh);
	setNumThreads(exec_policy, numThreads);
	
	lfs << "processing batch of size " << npairs << " with "
		<< numThreads << " threads " << std::endl;

	// for multiple seeds we store the seed with the highest identity
	AlignmentInfo *ai = new AlignmentInfo[npairs];
	std::pair<ushort, ushort> *seedlens = new std::pair<ushort, ushort>[npairs];

	bool *strands = new bool[npairs];
	int  *xscores = new int[npairs];
	TSeed  *seeds = new TSeed[npairs];

	/* GGGG: seed_count is hardcoded here (2) */
	for(int count = 0; count < seed_count; ++count)
	{
		auto start_time = std::chrono::system_clock::now();

		seqan::StringSet<seqan::Gaps<seqan::Dna5String>> seqsh_ex;
		seqan::StringSet<seqan::Gaps<seqan::Dna5String>> seqsv_ex;
		resize(seqsh_ex, npairs, seqan::Exact{});
		resize(seqsv_ex, npairs, seqan::Exact{});
		
	// extend the current seed and form a new gaps object
	#pragma omp parallel for
		for (uint64_t i = 0; i < npairs; ++i)
		{
			dibella::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);

		#ifdef TWOSEED
			ushort LocalSeedVOffset =
				(count == 0) ? cks->first.first : cks->second.first;
			ushort LocalSeedHOffset =
				(count == 0) ? cks->first.second : cks->second.second;
		#else
			// GGGG: TODO check reverse complement
			ushort LocalSeedVOffset = cks.pos[0].first;
			ushort LocalSeedHOffset = cks.pos[0].second;
		#endif

			seqan::Dna5String seedH;
			seqan::Dna5String seedV;

			auto start_time = std::chrono::system_clock::now();
			auto end_time   = std::chrono::system_clock::now();

			// Seed creation params are:
			// horizontal seed start offset, vertical seed start offset, length
			TSeed seed(LocalSeedHOffset, LocalSeedVOffset, seed_length);

			seedH = infix(seqan::source(seqsh[i]), beginPositionH(seed), endPositionH(seed));
			seedV = infix(seqan::source(seqsv[i]), beginPositionV(seed), endPositionV(seed));

			seqan::Dna5StringReverseComplement twin(seedH);

		#ifdef STATS
			seqan::Align<seqan::Dna5String> align;
			resize(rows(align), 2);
		#endif

			if(twin == seedV)
			{
				strands[i] = true;
				seqan::Dna5String twinseqH = seqan::source(seqsh[i]);
				seqan::Dna5StringReverseComplement twinRead(twinseqH);
				LocalSeedHOffset = length(twinseqH) - LocalSeedHOffset - seed_length;

				setBeginPositionH(seed, LocalSeedHOffset);
				setBeginPositionV(seed, LocalSeedVOffset);
				setEndPositionH(seed, LocalSeedHOffset + seed_length);
				setEndPositionV(seed, LocalSeedVOffset + seed_length);

				if(!noAlign) 
				{
					/* Perform match extension */
					start_time = std::chrono::system_clock::now();
					xscores[i] = extendSeed(seed, twinRead, seqan::source(seqsv[i]), seqan::EXTEND_BOTH, scoring_scheme,
							xdrop, (int)k,
							seqan::GappedXDrop());

					end_time = std::chrono::system_clock::now();
					add_time("XA:ExtendSeed", (ms_t(end_time - start_time)).count());
				}
				else
				{
					xscores[i] = 0;
				}

			#ifdef STATS
				assignSource(row(align, 0), infix(twinRead, beginPositionH(seed),
												endPositionH(seed)));
				assignSource(row(align, 1), infix(*seqV, beginPositionV(seed),
												endPositionV(seed)));
			#endif
			}
			else
			{
				strands[i] = false;
				if(!noAlign) 
				{
					start_time = std::chrono::system_clock::now();
					xscores[i] = extendSeed(seed, seqan::source(seqsh[i]), seqan::source(seqsv[i]), seqan::EXTEND_BOTH, scoring_scheme,
							xdrop, (int)k, 
							seqan::GappedXDrop());
					end_time = std::chrono::system_clock::now();
					add_time("XA:ExtendSeed", (ms_t(end_time - start_time)).count());

				#ifdef STATS
					assignSource(row(align, 0), infix(*seqH, beginPositionH(seed),
													endPositionH(seed)));
					assignSource(row(align, 1), infix(*seqV, beginPositionV(seed),
													endPositionV(seed)));
				#endif
				}
				else
				{
					xscores[i] = 0;
				}
			}

		#ifdef STATS
			xscores[i] = extendSeed(seed, seqan::source(seqsh[i]), seqan::source(seqsv[i]),
					   seqan::EXTEND_BOTH, scoring_scheme,
					   xdrop, (int)k, seqan::GappedXDrop());
			assignSource(seqsh_ex[i],
						 infix(seqan::source(seqsh[i]),
							   beginPositionH(seed), endPositionH(seed)));
			assignSource(seqsv_ex[i],
						 infix(seqan::source(seqsv[i]),
							   beginPositionV(seed), endPositionV(seed)));
		#endif
			seeds[i] = seed;
			seedlens[i].first  = static_cast<ushort>(seed._endPositionH -
													 seed._beginPositionH);
			seedlens[i].second = static_cast<ushort>(seed._endPositionV -
													 seed._beginPositionV);
		}

		auto end_time = std::chrono::system_clock::now();
    	add_time("XA:ExtendSeed", (ms_t(end_time - start_time)).count());

	#ifdef STATS
		start_time = std::chrono::system_clock::now();
		// alignment
		globalAlignment(exec_policy, seqsh_ex, seqsv_ex, scoring_scheme);
		
		end_time = std::chrono::system_clock::now();
    	add_time("XA:global_alignment", (ms_t(end_time - start_time)).count());
	#endif
		start_time = std::chrono::system_clock::now();
		
		// Compute stats
		if (count == 0)	// overwrite in the first seed
		{
		#pragma omp parallel for
			for (uint64_t i = 0; i < npairs; ++i)
			{
			#ifdef STATS
				computeAlignmentStats(ai[i].stats, seqsh_ex[i], seqsv_ex[i],
									  scoring_scheme);
			#endif
				ai[i].xscore = xscores[i];
				ai[i].rc     = strands[i];
				ai[i].seed   =   seeds[i];

				ai[i].seq_h_length = seqan::length(seqan::source(seqsh[i]));
				ai[i].seq_v_length = seqan::length(seqan::source(seqsv[i]));

				ai[i].seq_h_seed_length = seedlens[i].first;
				ai[i].seq_v_seed_length = seedlens[i].second;

				ai[i].seq_h_g_idx = col_offset + std::get<1>(mattuples[lids[i]]);
    			ai[i].seq_v_g_idx = row_offset + std::get<0>(mattuples[lids[i]]);
			}
		}
		else
		{
		#pragma omp parallel for
			for (uint64_t i = 0; i < npairs; ++i)
			{

			#ifdef STATS
				seqan::AlignmentStats stats;
				computeAlignmentStats(stats, seqsh_ex[i], seqsv_ex[i],
									  scoring_scheme);
		
				if (stats.alignmentIdentity > ai[i].stats.alignmentIdentity)
				{
					ai[i].stats				= stats;
					ai[i].seq_h_seed_length = seedlens[i].first;
					ai[i].seq_v_seed_length = seedlens[i].second;
				}
			#else
				if (xscores[i] > ai[i].xscore) // GGGG: TODO double check this logic with fresh neurons
				{
					ai[i].xscore = xscores[i];
					ai[i].rc     = strands[i];
					ai[i].seed   =   seeds[i];
					ai[i].seq_h_seed_length = seedlens[i].first;
					ai[i].seq_v_seed_length = seedlens[i].second;
				}
			#endif
			}
		}

		end_time = std::chrono::system_clock::now();
    	add_time("XA:ComputeStats", (ms_t(end_time - start_time)).count());
	}

	delete [] seedlens;
	delete [] xscores;
	delete [] strands;

	auto start_time = std::chrono::system_clock::now();

	std::vector<ushort> perprocessarray(nreads, 0);

	// Dump alignment info
	#pragma omp parallel
	{
	    #pragma omp for
		for (uint64_t i = 0; i < npairs; ++i)
		{
			// Only keep alignments that meet BELLA criteria
			bool passed = false;

			dibella::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);
			PostAlignDecision(ai[i], passed, ratioScoreOverlap, cks->overhang, cks->overhangT, cks->overlap, noAlign);

			if (passed)
			{
				// GGGG: store updated seed start/end position in the CommonKmers pairs (the semantics of these pairs change wrt the original semantics but that's okay)
				cks->first.first   = beginPositionV(ai[i].seed); 	// start on vertical sequence
				cks->first.second  = endPositionV(ai[i].seed); 		// end on vertical sequence
				cks->second.first  = beginPositionH(ai[i].seed);	// start on horizonal sequence
				cks->second.second = endPositionH(ai[i].seed);		// end on horizonal sequence

				cks->lenv 	= ai[i].seq_v_length;
				cks->lenh 	= ai[i].seq_h_length;
				cks->score  = ai[i].xscore;
				cks->passed = passed;	// keep this
			}
		}
	}

	auto end_time = std::chrono::system_clock::now();
  	add_time("XA:StringOp",
			 (ms_t(end_time - start_time)).count());

	delete [] ai;

	return;
}
