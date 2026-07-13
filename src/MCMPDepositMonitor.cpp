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

#include "NameFormat.h"

// from Gaudi
#include "Gaudi/Accumulators/Histogram.h"
#include "LHCbAlgs/Consumer.h"

// from Event
#include "Event/MCHit.h"
#include "Event/MCMPDeposit.h"

// MPDet
//#include "MPDet/DeMPDetector.h"

/** @class MCMPDepositMonitor
 *
 *  @author Lucas Foreman
 *  @date   2023-08-02
 */

class MCMPDepositMonitor : public LHCb::Algorithm::Consumer<void( const LHCb::MCMPDeposits& )> {

public:
  MCMPDepositMonitor( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override;
  void       operator()( const LHCb::MCMPDeposits& deposits ) const override;

private:
  static constexpr uint nLayers  = 6;
  static constexpr uint nModules = 28;
  static constexpr uint nStaves  = 20;
  static constexpr uint nSensors = 14;

  mutable std::optional<Gaudi::Accumulators::Histogram<1>>                      m_nDeposits;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>>                      m_nSignalDeposits;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>>                      m_nNoiseDeposits;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>>                      m_chargePerDeposit;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>>                      m_chargeArrivalTime;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>>                      m_depositsPerLayer;
  mutable std::optional<Gaudi::Accumulators::Histogram<1>>                      m_depositsPerModule;
  mutable std::array<std::optional<Gaudi::Accumulators::Histogram<1>>, nLayers> m_depositsPerStavePerLayer;
  mutable std::array<std::optional<Gaudi::Accumulators::Histogram<1>>, nLayers * nModules>
      m_depositsPerSensorPerLayerPerModule;
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( MCMPDepositMonitor )

MCMPDepositMonitor::MCMPDepositMonitor( const std::string& name, ISvcLocator* pSvcLocator )
    : Consumer( name, pSvcLocator, KeyValue{"DepositLocation", LHCb::MCMPDepositLocation::Default} ) {}

StatusCode MCMPDepositMonitor::initialize() {
  using MPDigitisation::formatter;
  return Consumer::initialize().andThen( [&] {
    using Axis = Gaudi::Accumulators::Axis<double>;
    m_nDeposits.emplace( this, "nDeposits", "Number of deposits; Deposits/event; Events", Axis{75, 0., 2e6} );
    m_nSignalDeposits.emplace( this, "nSignalDeposits", "Number of signal deposits; Deposits/event; Events",
                               Axis{200, 0., 7.5e5} );
    m_nNoiseDeposits.emplace( this, "nNoiseDeposits", "Number of noise deposits; Deposits/event; Events",
                              Axis{200, 0., 1.5e6} );
    m_chargePerDeposit.emplace( this, "ChargePerDeposit", "Charge per deposit; Charge/Deposit [e]; Deposits",
                                Axis{100, 0., 2.5e3} );
    m_chargeArrivalTime.emplace( this, "ArrivalTimeofCharges", "Charge Arrival time; #it{t} [ns]; Charge [e]",
                                 Axis{100, -25., 150.} );
    m_depositsPerLayer.emplace( this, "DepositsPerLayer", "Deposits per layer; Deposits/layer; Layer",
                                Axis{
                                    nLayers,
                                    -0.5,
                                    nLayers - 0.5,
                                } );
    m_depositsPerModule.emplace( this, "DepositsPerModule", "Deposits per module; Deposits/Module; Modules",
                                 Axis{nModules * nLayers, -0.5, nModules * nLayers - 0.5} );
    for ( uint layer_n = 0; nLayers != layer_n; ++layer_n ) {

      for ( uint module_n = 0; nModules != module_n; ++module_n ) {
        // Deposits Per Sensor for each Module///////////////////////////////////
        std::string name  = formatter()( module_n, "/DepositsPerSensorInModule", layer_n, "Layer" );
        std::string title = formatter()( "; Deposits/Sensor; Sensor", module_n, "Deposits per sensor in Module " );
        m_depositsPerSensorPerLayerPerModule[layer_n * nModules + module_n].emplace(
            this, name.c_str(), title.c_str(), Axis{nStaves * nSensors, -0.5, nStaves * nSensors - 0.5} );
        ////////////////////////////////////////////////////////////////////////
      }
      // Deposits Per Stave for each layer///////////////////////////////////////
      std::string name  = formatter()( layer_n, "DepositsPerStaveInLayer" );
      std::string title = formatter()( "; Deposits/Stave; Stave", layer_n, "Deposits per Stave in Layer " );
      m_depositsPerStavePerLayer[layer_n].emplace( this, name.c_str(), title.c_str(),
                                                   Axis{nModules * nStaves, -0.5, nModules * nStaves - 0.5} );
      //////////////////////////////////////////////////////////////////////////
    }
  } );
}

void MCMPDepositMonitor::operator()( const LHCb::MCMPDeposits& deposits ) const {
  int nNoise  = 0;
  int nSignal = 0;
  // Loop over MCMPDeposits
  for ( const LHCb::MCMPDeposit* pdeposit : deposits ) {
    const LHCb::MCMPDeposit& deposit = *pdeposit;

    for ( const auto triple : deposit ) {
      const auto& mcHit         = triple.mchit();
      const auto  depositCharge = triple.charge();
      const auto  arrivalTime   = triple.time();
      if ( mcHit != nullptr ) {
        ++( *m_chargePerDeposit )[depositCharge];
        ( *m_chargeArrivalTime )[arrivalTime] += depositCharge;
        ++nSignal;
      } else {
        ++nNoise;
      }
    }
    const auto& depositChannelID = deposit.channelID();
    const int   depositLayer     = (int)( depositChannelID.layer() );
    const int   depositModule    = (int)( depositChannelID.module() );
    const int   depositStave     = (int)( depositChannelID.stave() );
    const int   depositSensor    = (int)( depositChannelID.sensor() );

    ++( *m_depositsPerLayer )[depositLayer];
    ++( *m_depositsPerModule )[depositLayer * nModules + depositModule];
    ++( *( m_depositsPerStavePerLayer[depositLayer] ) )[depositModule * nStaves + depositStave];
    ++( *( m_depositsPerSensorPerLayerPerModule[( depositLayer * nModules +
                                                  depositModule )] ) )[depositStave * nSensors + depositSensor];
  }

  ++( *m_nDeposits )[deposits.size()];
  ++( *m_nNoiseDeposits )[nNoise];
  ++( *m_nSignalDeposits )[nSignal];
}
