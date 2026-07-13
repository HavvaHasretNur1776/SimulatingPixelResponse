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
#include "Event/MPDigit.h"

// MPDet
//#include "MPDet/DeMPDetector.h"

/** @class MPDigitMonitor
 *
 *  @author Lucas Foreman
 *  @date   2023-08-02
 */
class MPDigitMonitor : public LHCb::Algorithm::Consumer<void( const LHCb::MPDigits& )> {

public:
  MPDigitMonitor( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override;
  void       operator()( const LHCb::MPDigits& digits ) const override;

private:
  static constexpr uint nLayers  = 6;
  static constexpr uint nModules = 28;

  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_nDigits;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_digitsPerLayer;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>> m_digitsPerModule;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( MPDigitMonitor )

MPDigitMonitor::MPDigitMonitor( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator, KeyValue{"DigitLocation", LHCb::MPDigitLocation::Default} ) {}

StatusCode MPDigitMonitor::initialize() {
  return Consumer::initialize().andThen( [&] {
    using Axis = Gaudi::Accumulators::Axis<double>;
    m_nDigits.emplace( this, "nDigits", "Number of digits; Digits/event; Events", Axis{50, 0., 1e4} );
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

void MPDigitMonitor::operator()( const LHCb::MPDigits& digits ) const {

  // Loop over MPDigits
  for ( const LHCb::MPDigit* pdigit : digits ) {
    const LHCb::MPDigit& digit          = *pdigit;
    const auto&          digitChannelID = digit.channelID();
    const auto           digitLayer     = digitChannelID.layer();
    const auto           digitModule    = digitChannelID.module();

    ++( *m_digitsPerLayer )[(float)digitLayer];
    ++( *m_digitsPerModule )[(float)digitLayer * nModules + (float)digitModule];
  }

  ++( *m_nDigits )[digits.size()];
}
