// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GPUdEdx.h
/// \author David Rohr

#ifndef GPUDEDX_H
#define GPUDEDX_H

#include "GPUDef.h"
#include "GPUTPCGeometry.h"
#include "GPUCommonMath.h"

namespace GPUCA_NAMESPACE
{
namespace gpu
{
struct GPUdEdxInfo;
struct GPUParam;

#ifdef GPUCA_ALIROOT_LIB

class GPUdEdx
{
 public:
  GPUd() void clear() {}
  GPUd() void fillCluster(float qtot, float qmax, int padRow, float trackSnp, float trackTgl) {}
  GPUd() void fillSubThreshold(int padRow) {}
  GPUd() void computedEdx(GPUdEdxInfo& output, const GPUParam& param) {}
};

#else

class GPUdEdx
{
 public:
  // The driver must call clear(), fill clusters row by row outside-in, then run computedEdx() to get the result
  GPUd() void clear();
  GPUd() void fillCluster(float qtot, float qmax, int padRow, float trackSnp, float trackTgl);
  GPUd() void fillSubThreshold(int padRow);
  GPUd() void computedEdx(GPUdEdxInfo& output, const GPUParam& param);

 private:
  GPUd() float GetSortTruncMean(float* array, int count, int trunclow, int trunchigh);
  GPUd() void checkSubThresh(int roc);

  static constexpr int MAX_NCL = GPUCA_ROW_COUNT; // Must fit in mNClsROC (unsigned char)!

  float mChargeTot[MAX_NCL]; // No need for default, just some memory
  float mChargeMax[MAX_NCL]; // No need for default, just some memory
  float mSubThreshMinTot = 0.f;
  float mSubThreshMinMax = 0.f;
  unsigned char mNClsROC[4] = { 0 };
  unsigned char mNClsROCSubThresh[4] = { 0 };
  unsigned char mCount = 0;
  unsigned char mLastROC = 255;
  char mNSubThresh = 0;
};

GPUdi() void GPUdEdx::checkSubThresh(int roc)
{
  if (roc != mLastROC) {
    if (mNSubThresh) {
      for (int i = 0; i < mNSubThresh; i++) {
        mChargeTot[mCount] = mSubThreshMinTot;
        mChargeMax[mCount++] = mSubThreshMinMax;
      }
      mNClsROC[mLastROC] += mNSubThresh;
      mNClsROCSubThresh[mLastROC] += mNSubThresh;
    }
    mNSubThresh = 0;
    mSubThreshMinTot = 1e10f;
    mSubThreshMinMax = 1e10f;
  }

  mLastROC = roc;
}

GPUdi() void GPUdEdx::fillCluster(float qtot, float qmax, int padRow, float trackSnp, float trackTgl)
{
  const int roc = GPUTPCGeometry::GetROC(padRow);
  checkSubThresh(roc);
  float factor = CAMath::Sqrt((1 - trackSnp * trackSnp) / (1 + trackTgl * trackTgl));
  factor /= GPUTPCGeometry::PadHeight(padRow);
  qtot *= factor;
  qmax *= factor / GPUTPCGeometry::PadWidth(padRow);

  mChargeTot[mCount] = qtot;
  mChargeMax[mCount++] = qmax;
  mNClsROC[roc]++;
  if (qtot < mSubThreshMinTot) {
    mSubThreshMinTot = qtot;
  }
  if (qmax < mSubThreshMinMax) {
    mSubThreshMinMax = qmax;
  }
}

GPUdi() void GPUdEdx::fillSubThreshold(int padRow)
{
  const int roc = GPUTPCGeometry::GetROC(padRow);
  checkSubThresh(roc);
  mNSubThresh++;
}

#endif // GPUCA_ALIROOT_LIB
} // namespace gpu
} // namespace GPUCA_NAMESPACE

#endif
