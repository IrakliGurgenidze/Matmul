//
// Public interface for estimated product size.
//

#pragma once

#include <vector>
#include "Types.h"

#ifndef ESTIMATOR_H
#define ESTIMATOR_H

double estimateProductSize(const std::vector<R1Tuple>& R1Tuples,
                                    const std::vector<R2Tuple>& R2Tuples,
                                    double epsilon = 0.1);

#endif //ESTIMATOR_H
