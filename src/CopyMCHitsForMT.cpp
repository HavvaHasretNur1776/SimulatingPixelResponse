/*****************************************************************************\
* (c) Copyright 2000-2018 CERN for the benefit of the LHCb Collaboration      *
*                                                                             *
* This software is distributed under the terms of the GNU General Public      *
* Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   *
*                                                                             *
* In applying this licence, CERN does not waive the privileges and immunities *
* granted to it by virtue of its status as an Intergovernmental Organization  *
* or submit itself to any jurisdiction.                                       *
\*****************************************************************************/

#include <memory>

// from Gaudi
#include "Gaudi/Accumulators/Histogram.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "LHCbAlgs/Transformer.h"

// from Event
#include "Event/MCHit.h"

using namespace Gaudi::Functional;

/**
 *  This small algorithm puts the MCHits from all spills into an array
 *
 *  @author Manuel Schiller
 *  @date   2023-06-28
 */
class CopyMCHitsForMT : public LHCb::Algorithm::Transformer<LHCb::MCHits( const LHCb::MCHits& )> {

public:
  CopyMCHitsForMT( const std::string& name, ISvcLocator* pSvcLocator );
  LHCb::MCHits operator()( const LHCb::MCHits& mcHitsVector ) const override;

private:
  static LHCb::MCHit* clone( const LHCb::MCHit& hit );
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( CopyMCHitsForMT )

CopyMCHitsForMT::CopyMCHitsForMT( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator, {"InputLocation", "/Event/" + LHCb::MCHitLocation::FT},
                   {"OutputLocation", "/Event/MC/MT/Hits"} ) {}

LHCb::MCHits CopyMCHitsForMT::operator()( const LHCb::MCHits& mcHits ) const {
  const auto isNotFT = []( const LHCb::MCHit* hit ) {
    const auto z = hit->midPoint().z();
    return z < 7800. || ( 8050. < z && z < 8475. ) || ( 8750. < z && z < 9150. ) || z > 9450.;
  };
  const std::size_t sz = std::count_if( mcHits.begin(), mcHits.end(), isNotFT );
  LHCb::MCHits      retVal;
  retVal.reserve( sz );

  for ( const LHCb::MCHit* mcHit : mcHits ) {
    if ( isNotFT( mcHit ) ) {
      LHCb::MCHit* newhit = clone( *mcHit );
      retVal.push_back( newhit );
    }
  }
  return retVal;
}

LHCb::MCHit* CopyMCHitsForMT::clone( const LHCb::MCHit& other ) {
  auto retVal = std::make_unique<LHCb::MCHit>();
  retVal->setSensDetID( other.sensDetID() )
      .setEntry( other.entry() )
      .setDisplacement( other.displacement() )
      .setEnergy( other.energy() )
      .setTime( other.time() )
      .setP( other.p() )
      .setMCParticle( other.mcParticle() );
  return retVal.release();
}
