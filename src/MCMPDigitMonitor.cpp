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

#include <boost/container/flat_set.hpp>

// from Gaudi
#include "Gaudi/Accumulators/Histogram.h"
#include "LHCbAlgs/Consumer.h"

// from Event
#include "Event/MCHit.h"
#include "Event/MCMPDeposit.h"
#include "Event/MCMPDigit.h"

// MPDet
//#include "MPDet/DeMPDetector.h"

/** @class MCMPDigitMonitor
 *
 *  @author Lucas Foreman
 *  @date   2023-08-02
 */
class MCMPDigitMonitor : public LHCb::Algorithm::Consumer<void( const LHCb::MCMPDigits& )> {

public:
  MCMPDigitMonitor( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override;
  void       operator()( const LHCb::MCMPDigits& digits ) const override;

private:
  uint nLayers  = 6;
  uint nModules = 28;

  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_nDigits;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_nSignalDigits;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_nNoiseDigits;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_chargePerDigit;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_digitsPerLayer;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_digitsPerModule;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_hitsPerDigit;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( MCMPDigitMonitor )

MCMPDigitMonitor::MCMPDigitMonitor( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator, KeyValue{"MCDigitLocation", LHCb::MCMPDigitLocation::Default} ) {}

StatusCode MCMPDigitMonitor::initialize() {
  return Consumer::initialize().andThen( [&] {
    using Axis = Gaudi::Accumulators::Axis<double>;
    m_nDigits.emplace( this, "nDigits", "Number of digits; Digits/event; Events", Axis{50, 0., 2e6} );
    m_nSignalDigits.emplace( this, "nSignalDigits", "Number of signal digits; Digits/event; Events",
                             Axis{50, 0., 7.5e5} );
    m_nNoiseDigits.emplace( this, "nNoiseDigits", "Number of noise digits; Digits/event; Events", Axis{50, 0., 1.5e6} );
    m_chargePerDigit.emplace( this, "ChargePerDigit", "Charge per digit; Charge/Digit [e]; Digits",
                              Axis{100, 0., 2.5e3} );
    m_hitsPerDigit.emplace( this, "HitsPerDigit", "Hits per digit; Hits/Digit; Digits", Axis{8, -0.5, 7.5} );
    m_digitsPerLayer.emplace( this, "DigitsPerLayer", "Digits per layer; Digits/layer; Layer",
                              Axis{
                                  nLayers,
                                  -0.5,
                                  nLayers - 0.5,
                              } );
    m_digitsPerModule.emplace( this, "DigitsPerModule", "Digits per module; Depoists/Module; Modules",
                               Axis{nModules * nLayers, -0.5, nModules * nLayers - 0.5} );
  } );
}

void MCMPDigitMonitor::operator()( const LHCb::MCMPDigits& digits ) const {
  // Number of (All/Noise/Real/Spillover) Digits per event Histograms
  int                                            nNoise      = 0;
  int                                            nSignal     = 0;
  float                                          digitCharge = 0.;
  boost::container::flat_set<const LHCb::MCHit*> mcHits;
  mcHits.reserve( 16 ); // will normally contain 0 or 1 MCHit

  // Loop over MCMPDigits
  for ( const LHCb::MCMPDigit* pdigit : digits ) {
    const LHCb::MCMPDigit& digit = *pdigit;

    // Build mcHits vector
    mcHits.clear();
    digitCharge = 0.;
    for ( const auto& deposit : digit.deposits() ) {
      for ( const auto triple : *deposit ) {
        const auto& mcHit = triple.mchit();
        if ( mcHit != nullptr ) {
          mcHits.insert( mcHit );
          digitCharge += triple.charge();
        }
      }
    }
    // count empty mcHits vectors
    if ( mcHits.empty() ) {
      ++nNoise;
    } else {
      ++nSignal;
      const auto& digitChannelID = digit.channelID();
      const auto  digitLayer     = digitChannelID.layer();
      const auto  digitModule    = digitChannelID.module();

      ++( *m_digitsPerLayer )[(float)digitLayer];
      ++( *m_digitsPerModule )[(float)digitLayer * nModules + (float)digitModule];
      ++( *m_chargePerDigit )[digitCharge];
      ++( *m_hitsPerDigit )[mcHits.size()];
    }
  }
  ++( *m_nDigits )[digits.size()];
  ++( *m_nNoiseDigits )[nNoise];
  ++( *m_nSignalDigits )[nSignal];
}
