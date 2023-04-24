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

#include "BucketManager.h"


BucketManager::BucketManager() :
    buckets(std::numeric_limits<match_t>::max() + 1)
{
}

void BucketManager::insert(match_t match, seq_id_t seqID, mismatch_t mismatch)
{
    buckets[match][seqID].push_back(mismatch);
}
