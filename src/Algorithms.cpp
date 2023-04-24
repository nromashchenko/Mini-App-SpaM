/**
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * author: Matthias Blanke
 * mail  : matthias.blanke@biologie.uni-goettingen.de
 */

#include <iostream>
#include <fstream>
#include "Algorithms.h"
#include "Scoring.h"
#include "SubstitutionMatrix.h"
#include "Match.h"
#include "MatchManager.h"


/**
 * Calculate fswm distance between reads and genomes considering all spaced words.
 */
bool Algorithms::fswm_complete(const BucketManager &genomeBucketManager,
                               const BucketManager &readBucketManager,
                               Scoring &fswm_distances) {
	SubstitutionMatrix substMat;
	std::ofstream histogramFile;
	if (fswm_params::g_writeHistogram) { histogramFile.open(fswm_params::g_outfoldername + "histogram.txt", std::ios_base::app); }

    // Iterate over all possible matches
    for (size_t matches = 0; matches < genomeBucketManager.buckets.size(); matches++)
    {
        const auto& genomeIndex = genomeBucketManager.buckets[matches];
        const auto& readIndex = readBucketManager.buckets[matches];

        // Loop through genome words with same spaced k-mer
        for (const auto& [genomeSeqID, genomeDontCares] : genomeIndex)
        {
            for (auto genomeDontCare : genomeDontCares)
            {
                // same for read words
                for (const auto& [readSeqID, readDontCares] : readIndex)
                {
                    for (auto readDontCare : readDontCares)
                    {
                        int score = 0;
                        int mismatches = 0;

                        for (int i = 0; i < fswm_params::g_spaces; i++) {
                            score += substMat.chiaromonte[genomeDontCare & 0x03][readDontCare & 0x03];
                            mismatches += substMat.mismatch[genomeDontCare & 0x03][readDontCare & 0x03];
                            genomeDontCare = genomeDontCare >> 2;
                            readDontCare = readDontCare >> 2;
                        }

                        if (fswm_params::g_writeHistogram) {
                            histogramFile << readSeqID << "\t" << genomeSeqID << "\t" << score << std::endl;
                        }

                        if (score > fswm_params::g_filteringThreshold) {
                            fswm_distances.scoringMap[readSeqID][genomeSeqID] += score;
                            fswm_distances.mismatchCount[readSeqID][genomeSeqID] += mismatches;
                            fswm_distances.spacedWordMatchCount[readSeqID][genomeSeqID] += 1;
                        }

                    }
                }
            }
        }
	}

	if (fswm_params::g_writeHistogram) { histogramFile.close(); }

	return true;
}
