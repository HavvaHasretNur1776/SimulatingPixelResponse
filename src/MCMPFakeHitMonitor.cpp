/*****************************************************************************\
* (c) Copyright 2000-2023 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
// from Gaudi
#include "Gaudi/Accumulators/Histogram.h"
#include "LHCbAlgs/Consumer.h"

// from Event
#include "Event/MCHit.h"
#include "Event/MCMPFakeHit.h"

// MPDet
//#include "MPDet/DeMPDetector.h"

/** @class MCMPFakeHitMonitor
 *
 *  @author Manuel Schiller
 *  @date   2023-07-21
 */
class MCMPFakeHitMonitor : public LHCb::Algorithm::Consumer<void( const LHCb::MCMPFakeHits& )> {

public:
  MCMPFakeHitMonitor( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override;
  void       operator()( const LHCb::MCMPFakeHits& deposits ) const override;

private:
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_histX;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_histY;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_histZ;
  mutable std::optional<Gaudi::Accumulators::Histogram<2>> m_histXY;
  mutable std::optional<Gaudi::Accumulators::Histogram<2>> m_histXZ;
  mutable std::optional<Gaudi::Accumulators::Histogram<2>> m_histYZ;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_histErrX;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_histErrY;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_histCorrelXY;
  mutable std::optional<Gaudi::Accumulators::Histogram<2>> m_histResidXY;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_histPullX;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_histPullY;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( MCMPFakeHitMonitor )

MCMPFakeHitMonitor::MCMPFakeHitMonitor( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator, KeyValue{"DepositLocation", LHCb::MCMPFakeHitLocation::Default} ) {}

StatusCode MCMPFakeHitMonitor::initialize() {
  return Consumer::initialize().andThen( [&] {
    using Axis = Gaudi::Accumulators::Axis<double>;
    m_histX.emplace( this, "x", "x;x [mm];entries", Axis{100, -3e3, 3e3} );
    m_histY.emplace( this, "y", "y;y [mm];entries", Axis{100, -2.5e3, 2.5e3} );
    m_histZ.emplace( this, "z", "z;z [mm];entries", Axis{300, 7e3, 10e3} );
    m_histXY.emplace( this, "xy", "xy;x [mm];y [mm];entries", Axis{100, -3e3, 3e3}, Axis{100, -2.5e3, 2.5e3} );
    m_histXZ.emplace( this, "xz", "xz;x [mm];z [mm];entries", Axis{100, -3e3, 3e3}, Axis{300, 7e3, 10e3} );
    m_histYZ.emplace( this, "yz", "yz;y [mm];z [mm];entries", Axis{100, -2.5e3, 2.5e3}, Axis{300, 7e3, 10e3} );
    m_histErrX.emplace( this, "errX", "errX [mm];entries", Axis{100, 0., 1.} );
    m_histErrY.emplace( this, "errY", "errY [mm];entries", Axis{100, 0., 1.} );
    m_histCorrelXY.emplace( this, "correlXY", "correlXY;entries", Axis{100, -1., 1.} );
    m_histResidXY.emplace( this, "residXY", "resid;residual x [mm];residual y [mm];entries", Axis{100, -1., 1.},
                           Axis{100, -1., 1.} );
    m_histPullX.emplace( this, "pullX", "pullX;entries", Axis{100, -5., 5.} );
    m_histPullY.emplace( this, "pullY", "pullY;entries", Axis{100, -5., 5.} );
  } );
}

void MCMPFakeHitMonitor::operator()( const LHCb::MCMPFakeHits& hits ) const {
  for ( const LHCb::MCMPFakeHit* phit : hits ) {
    const LHCb::MCMPFakeHit& hit        = *phit;
    const auto&              mcmidpoint = hit.mchit()->midPoint();
    ++( *m_histX )[hit.x()];
    ++( *m_histY )[hit.y()];
    ++( *m_histZ )[hit.z()];
    ++( *m_histXY )[{hit.x(), hit.y()}];
    ++( *m_histXZ )[{hit.x(), hit.z()}];
    ++( *m_histYZ )[{hit.y(), hit.z()}];
    ++( *m_histErrX )[hit.errx()];
    ++( *m_histErrY )[hit.erry()];
    ++( *m_histCorrelXY )[hit.rhoxy()];
    ++( *m_histResidXY )[{hit.x() - mcmidpoint.x(), hit.y() - mcmidpoint.y()}];
    ++( *m_histPullX )[( hit.x() - mcmidpoint.x() ) / hit.errx()];
    ++( *m_histPullY )[( hit.y() - mcmidpoint.y() ) / hit.erry()];
  }
}
