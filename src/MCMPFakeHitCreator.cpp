/*****************************************************************************\
 * (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      *
 *                                                                             *
 * This software is distributed under the terms of the GNU General Public      *
 * Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
 *                                                                             *
 * In applying this licence, CERN does not waive the privileges and immunities *
 * granted to it by virtue of its status as an Intergovernmental Organization  *
 * or submit itself to any jurisdiction.                                       *
\*****************************************************************************/
//#include "RadDamage.h"
//#include <DetDesc/DetectorElement.h>
//#include <Event/MCVPDigit.h>
//#include <Kernel/VPDataFunctor.h>
//#include <LHCbMath/LHCbMath.h>
//#include <VPDet/DeVP.h>

#include <array>
#include <cassert>
#include <cmath>
#include <optional>
#include <stdexcept>

// form Gaudi
#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiKernel/RndmGenerators.h"
#include "LHCbAlgs/Transformer.h"
#include <Gaudi/Accumulators/Histogram.h>
#include <GaudiKernel/PhysicalConstants.h>

#include "Detector/MP/MPChannelID.h"
#include "Event/MCMPDigit.h"
#include "Event/MCMPFakeHit.h"
#include "Event/MPDigit.h"
#include "LHCbMath/LHCbMath.h"
#include <Event/GenHeader.h>
#include <Event/MCHit.h>

using namespace LHCb;

/**
 *  Algorithm for creating MCMPFakeHits from MCHits. The following processes
 *  are simulated:
 *  - spatial smearing
 */
class MCMPFakeHitCreator : public LHCb::Algorithm::Transformer<LHCb::MCMPFakeHits( const LHCb::MCHits& mcMain )> {

public:
  MCMPFakeHitCreator( const std::string& name, ISvcLocator* pSvcLocator );
  virtual StatusCode initialize() override;
  LHCb::MCMPFakeHits operator()( const LHCb::MCHits& mcMain ) const override;

  static std::optional<std::array<double, 3>> decompose( const std::array<double, 3>& M ) noexcept;

private:
  constexpr static double s_isqrt12 = 0.28867513459481288225;
  Gaudi::Property<double> m_xError{this, "xError", 0.056 * s_isqrt12};
  Gaudi::Property<double> m_yError{this, "yError", 0.165 * s_isqrt12};
  Gaudi::Property<double> m_xyCorrel{this, "xyCorrel", 0.};
  Gaudi::Property<double> m_hitEfficiency{this, "hitEfficiency", 1.0};
  Gaudi::Property<double> m_minEnergy{this, "minEnergy", 0.};

  mutable Rndm::Numbers m_gauss;
  mutable Rndm::Numbers m_uniform;
};

DECLARE_COMPONENT( MCMPFakeHitCreator )

MCMPFakeHitCreator::MCMPFakeHitCreator( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer{name,
                  pSvcLocator,
                  {KeyValue{"InputLocationMain", LHCb::MCHitLocation::MP}},
                  KeyValue{"OutputLocation", LHCb::MCMPFakeHitLocation::Default}} {}

StatusCode MCMPFakeHitCreator::initialize() {
  return Transformer::initialize()
      .andThen( [&] {
        if ( m_hitEfficiency < 0. || m_hitEfficiency > 1. ) {
          error() << "hitEfficiency must be in interval [0., 1.]";
          return StatusCode::FAILURE;
        }
        if ( !std::isfinite( m_minEnergy ) || m_minEnergy < 0. ) {
          error() << "minEnergy must not be negative";
          return StatusCode::FAILURE;
        }
        return StatusCode::SUCCESS;
      } )
      .andThen( [&] { return m_gauss.initialize( randSvc(), Rndm::Gauss( 0.0, 1.0 ) ); } )
      .andThen( [&] { return m_uniform.initialize( randSvc(), Rndm::Flat( 0., 1. ) ); } );
}

std::optional<std::array<double, 3>> MCMPFakeHitCreator::decompose( const std::array<double, 3>& M ) noexcept {
  if ( M[0] < 0 ) return {};
  const auto l11  = std::sqrt( M[0] );
  const auto l21  = M[1] / l11;
  const auto tmp0 = M[2] - l21 * l21;
  if ( tmp0 < 0 ) return {};
  return {{l11, l21, std::sqrt( tmp0 )}};
}

LHCb::MCMPFakeHits MCMPFakeHitCreator::operator()( const LHCb::MCHits& mcMain ) const {
  // decompose the error matrix to smear correctly in x and y
  const std::array<double, 3> cov{{m_xError * m_xError, m_xyCorrel * m_xError * m_yError, m_yError * m_yError}};
  const auto                  maybeL = decompose( cov );
  if ( !maybeL.has_value() ) {
    throw std::runtime_error(
        "MCMPFakeHitCreator::operator()(const LHCb::MCHits&): error matrix not positive definite!" );
  }
  const auto& L = *maybeL;

  // process MCHits
  LHCb::MCMPFakeHits retVal;
  retVal.reserve( mcMain.size() );
  const double hiteff = m_hitEfficiency;
  const double minE   = m_minEnergy;
  for ( const auto& mchit : mcMain ) {
    if ( mchit->energy() < minE ) continue;
    if ( hiteff <= 0. ) continue;
    if ( hiteff < 1. && m_uniform() > hiteff ) continue;
    const auto&                 midpoint = mchit->midPoint();
    const std::array<double, 2> u{{m_gauss(), m_gauss()}};
    const auto                  xsmear = midpoint.x() + L[0] * u[0];
    const auto                  ysmear = midpoint.y() + L[1] * u[0] + L[2] * u[1];
    retVal.push_back( new LHCb::MCMPFakeHit{mchit, xsmear, ysmear, midpoint.z(), m_xError, m_yError, m_xyCorrel} );
  }
  return retVal;
}
