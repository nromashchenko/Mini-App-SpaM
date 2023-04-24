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
#ifndef FSWM_BucketManager_H_
#define FSWM_BucketManager_H_

#include <unordered_map>
#include "GlobalParameters.h"


struct BucketManager {
    // TODO: We may need to store this in a more sophisticated data structure
    // (for compression)
    using SeqMismatchIndex = std::vector<mismatch_t>;

    using seqToIndex = std::unordered_map<seq_id_t, SeqMismatchIndex>;

    std::vector<seqToIndex> buckets;

    BucketManager();

    void insert(match_t match, seq_id_t seqID, mismatch_t mismatch);
};
#endif
