#ifndef INS_MESH_H
#define INS_MESH_H

#include "state/state.h"
#include "util/measurement.h"

namespace Ins
{
    State::Pva insMech(const State::Pva& prePva,const Measurement::IMU& preImu,const Measurement::IMU& currImu);
};
#endif
