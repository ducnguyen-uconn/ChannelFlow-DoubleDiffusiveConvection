/**
 * System class of double-diffusive equations for standard channel flows
 *
 * Original author: Duc Nguyen
 */

#ifndef DDCALGO_H
#define DDCALGO_H

#include <memory>
#include "cfbasics/cfvector.h"
#include "cfbasics/mathdefs.h"
#include "channelflow/chebyshev.h"
#include "channelflow/flowfield.h"
#include "channelflow/symmetry.h"
#include "modules/ddc/ddcflags.h"
#include "modules/ddc/macros.h"
#include "modules/ddc/dde.h"
using namespace std;
namespace chflow {

class DDCAlgo : public DNSAlgorithm {
   public:
    void advance(std::vector<FlowField>& fields, int nSteps = 1);

};
}  // namespace chflow
#endif