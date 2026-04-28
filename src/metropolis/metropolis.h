#pragma once

#include "../gauge/GaugeField.h"
#include "../su3/utils.h"

namespace metropolis {
void sweep(GaugeField& field, const GeometryCB& geo, double beta, const std::vector<SU3>& set, size_t& accepted, size_t& proposed, std::mt19937_64& rng);
}
