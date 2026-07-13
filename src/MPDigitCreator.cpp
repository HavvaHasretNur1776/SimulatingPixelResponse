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

// #include "Kernel/VPDataFunctor.h"
// #include "LHCbAlgs/Transformer.h"

#include <cmath>
#include <map>
#include <numeric>
#include <yaml-cpp/yaml.h>

// form Gaudi
#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiKernel/RndmGenerators.h"
#include "LHCbAlgs/Transformer.h"
#include <Gaudi/Accumulators/Histogram.h>

#include "Detector/MP/MPChannelID.h"
#include "Event/MCMPDeposit.h"
#include "Event/MCMPDigit.h"
#include "Event/MPDigit.h"

using namespace LHCb;

// [FIXME] Manuel Schiller, 20230706 - use Gaudi Property for threshold for
// now to get the code to run...
#define THRESHOLD_FROM_PROPERTY

namespace {
  using Threshold = std::array<double, 208>;
}

/** @class MCMPDigitCreator MCMPDigitCreator.h
 *
 *  From the list of MCMPDeposit (deposited energy in each MPChannel),
 *  this algorithm converts the deposited energy in ADC charge according
 *  to the mean number of photoelectron per MIP
 *  + the mean ADC count per photoelectron.
 *  Created digits are put in transient data store.
 *
 *
 *  @author Ellinor Eckstein, Piet Nogga, Hannah Schmitz
 *  @date   2022-10-25
 *
 *
 * TODO: replace electronic noise number,
 * TODO: check if this is the correct convention
 */
class MPDigitCreator
#ifndef THRESHOLD_FROM_PROPERTY // [FIXME] Manuel Schiller, 20230706
    : public LHCb::Algorithm::Transformer<LHCb::MPDigits( const LHCb::MCMPDigits&, const Threshold& )> {
#else  // THRESHOLD_FROM_PROPERTY
    : public LHCb::Algorithm::Transformer<LHCb::MPDigits( const LHCb::MCMPDigits& )> {
#endif // THRESHOLD_FROM_PROPERTY
public:
  /// Standard constructor
  MPDigitCreator( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override; ///< Algorithm initialization

#ifndef THRESHOLD_FROM_PROPERTY // [FIXME] Manuel Schiller, 20230706
  LHCb::MPDigits operator()( LHCb::MCMPDigits const&, Threshold const& ) const override;
#else  // THRESHOLD_FROM_PROPERTY
  LHCb::MPDigits operator()( LHCb::MCMPDigits const& ) const override;
#endif // THRESHOLD_FROM_PROPERTY

private: // TODO: calculate noise!!!!
  /// Noise in number of electrons
  Gaudi::Property<double> m_electronicNoise{this, "ElectronicNoise", 130.0};

  /// Window in which pixel rise time is registered
  Gaudi::Property<double> m_timingWindow{this, "TimingWindow", 25.};

  /// Offset from z/c arrival that timing window starts
  Gaudi::Property<double> m_timingOffset{this, "TimingOffset", 0.};

  /// Option to simulate masked pixels
  Gaudi::Property<bool> m_maskedPixels{this, "MaskedPixels", true};

  /// Option to simulate noisy pixels
  Gaudi::Property<bool> m_noisyPixels{this, "NoisyPixels", false};

  /// Fraction of masked pixels (0.05%)
  Gaudi::Property<double> m_fractionMasked{this, "FractionMasked", 5e-4};

  /// Fraction of noisy pixels (100e, with a 1000 threshold is <3e-16 noise pixels per event)
  Gaudi::Property<double> m_fractionNoisy{this, "FractionNoisy", 0.};

  /// Flag to activate monitoring histograms or not
  Gaudi::Property<bool> m_monitoring{this, "Monitoring", true};

#ifdef THRESHOLD_FROM_PROPERTY // [FIXME] Manuel Schiller, 20230706
  // in electrons, [FIXME] 2.4 ke is not the final value
  Gaudi::Property<double> m_uniformThreshold{this, "UniformThreshold", 2.4e3};
#endif // THRESHOLD_FROM_PROPERTY

  /// Number of signal MPDigits created
  mutable Gaudi::Accumulators::StatCounter<> m_numSignalDigitsCreated{this, "1 Signal MPDigits"};

  /// Number of noise MPDigits created
  mutable Gaudi::Accumulators::StatCounter<> m_numNoiseDigitsCreated{this, "2 Noise MPDigits"};

  /// Number of MPDigits missed due to masked pixels
  mutable Gaudi::Accumulators::StatCounter<> m_numDigitsKilledMasked{this, "3 Lost Random bad pixels"};

  /// Number of MPDigits created due to spillover
  mutable Gaudi::Accumulators::StatCounter<> m_numSpillover{this, "4 Spillover MPDigits"};

  /// Number of MPDigits lost due to timing window
  mutable Gaudi::Accumulators::StatCounter<> m_numLostTiming{this, "5 Out of time MPDigits"};

  mutable Gaudi::Accumulators::Histogram<1> m_hist_mainChargeBefore{
      this, "MainChargeBefore", "Pixel Charge in main BX, e- before threshold", {100, 0., 5e4}};
  mutable Gaudi::Accumulators::Histogram<1> m_hist_mainChargeAfter{
      this, "MainChargeAfter", "Pixel Charge in main BX, e- after threshold", {100, 0., 5e4}};

  mutable Gaudi::Accumulators::Histogram<2> m_hist_sensorVsPixelChargeBefore{
      this,
      "SenorVsPixelChargeBefore",
      "Sensor verse Pixel Charge before threshold (e-)",
      {208, -.5, 207.5},
      {100, 0., 5e4}};
  mutable Gaudi::Accumulators::Histogram<2> m_hist_sensorVsPixelChargeAfter{
      this,
      "SenorVsPixelChargeAfter",
      "Sensor verse Pixel Charge after threshold (e-)",
      {208, -.5, 207.5},
      {100, 0., 5e4}};

  /// Location of the condition with the threshold in e- per sensor //TODO: change this to fix property
  Gaudi::Property<std::string> m_thresholdConditionName{this, "ThresholdCondition",
                                                        "Conditions/Calibration/MP/MPDigitisationParam/Thresholds"};

  /// Random number generators
  mutable Rndm::Numbers m_gauss;
  mutable Rndm::Numbers m_uniform;
  mutable Rndm::Numbers m_poisson;
};

// /**
//  *  Using MCMPDigits as input, this algorithm simulates the response of the
//  *  MightyPix ASIC and produces a MPDigit for each MCMPDigit above threshold.
//  *
//  *  @author Ellinor Eckstein. Piet Nogga, Hannah Schmitz
//  *  @date   2022-10-25
//  */

DECLARE_COMPONENT( MPDigitCreator )

MPDigitCreator::MPDigitCreator( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer{name,
                  pSvcLocator,
                  {KeyValue{"InputLocation", LHCb::MCMPDigitLocation::Default}
#ifndef THRESHOLD_FROM_PROPERTY // [FIXME] Manuel Schiller, 20230706
                   ,
                   KeyValue{"Threshold", name + "_Threshold"}
#endif // THRESHOLD_FROM_PROPERTY
                  },
                  KeyValue{"OutputLocation", LHCb::MPDigitLocation::Default}} {
}

//=============================================================================
// Initialization
//=============================================================================
StatusCode MPDigitCreator::initialize() {
  return Transformer::initialize()
      .andThen( [&] { return m_gauss.initialize( randSvc(), Rndm::Gauss( 0., 1. ) ); } )
      .andThen( [&] { return m_uniform.initialize( randSvc(), Rndm::Flat( 0., 1. ) ); } )
      .andThen( [&] {
        if ( !m_noisyPixels ) return StatusCode( StatusCode::SUCCESS ); // nothing else to do
        const double averageNoisy =
            m_fractionNoisy * LHCb::Detector::MPChannelID::nLayers<>() * LHCb::Detector::MPChannelID::nModules<>() *
            LHCb::Detector::MPChannelID::nStaves<>() * LHCb::Detector::MPChannelID::nSensors<>() *
            LHCb::Detector::MPChannelID::nColumns<>() * LHCb::Detector::MPChannelID::nRows<>();
        if ( averageNoisy < 1 ) {
          m_noisyPixels = false;
          return StatusCode( StatusCode::SUCCESS );
        }
        return m_poisson.initialize( randSvc(), Rndm::Poisson( averageNoisy ) );
      } )
      // .andThen( [&] {
      //   addConditionDerivation( {m_thresholdConditionName}, inputLocation<Threshold>(),
      //                           [&]( YAML::Node const& cond ) -> Threshold {
      //                             return cond["ThresholdPerSensor"].as<std::array<double, 208>>();
      //                           } );
      //   return StatusCode( StatusCode::SUCCESS );
      // } )
      .andThen( [&] { /*setHistoTopDir( "MP/" );*/ } ); // FIXME - how to do
                                                        // this with
                                                        // Gaudi::Accumulators::Histogram?
}

//=============================================================================
// Main execution
//=============================================================================
#ifndef THRESHOLD_FROM_PROPERTY // [FIXME] Manuel Schiller, 20230706
LHCb::MPDigits MPDigitCreator::operator()( LHCb::MCMPDigits const& mcdigits, Threshold const& threshold ) const {
#else  // THRESHOLD_FROM_PROPERTY
LHCb::MPDigits MPDigitCreator::operator()( LHCb::MCMPDigits const& mcdigits ) const {
#endif // THRESHOLD_FROM_PROPERTY
  // Create a container for the digits and transfer ownership to the TES.
  LHCb::MPDigits digits;
  unsigned int   numSpillover           = 0;
  unsigned int   numLostTiming          = 0;
  unsigned int   numDigitsKilledMasked  = 0;
  unsigned int   numSignalDigitsCreated = 0;
  unsigned int   numNoiseDigitsCreated  = 0;

  // Loop over the MC digits.
  for ( const LHCb::MCMPDigit* mcdigit : mcdigits ) {
    // Sum up the collected charge from all deposits, in the time windows
    double chargeMain = 0.;
    for ( const LHCb::MCMPDeposit* dep : mcdigit->deposits() ) {
      for ( const auto triple : *dep ) {
        auto depTime = triple.time() + m_timingOffset;
        if ( depTime > 0. && depTime < m_timingWindow ) {
          // is this in the main bunch
          chargeMain += triple.charge();
          // if no MCHit then assume spillover
          numSpillover += nullptr == triple.mchit();
        } else {
          numLostTiming += nullptr == triple.mchit();
        }
      }
    }
    const LHCb::Detector::MPChannelID channel = mcdigit->channelID();

    if ( m_monitoring ) {
      ++m_hist_mainChargeBefore[chargeMain];
      ++m_hist_sensorVsPixelChargeBefore[{to_unsigned( channel.sensor() ), chargeMain}];
    }
    // Check if the collected charge exceeds the threshold.
#ifndef THRESHOLD_FROM_PROPERTY // [FIXME] Manuel Schiller, 20230706
    if ( chargeMain < threshold[to_unsigned( channel.sensor() )] + m_electronicNoise * m_gauss() ) continue;
#else  // THRESHOLD_FROM_PROPERTY
    if ( chargeMain < m_uniformThreshold + m_electronicNoise * m_gauss() ) continue;
#endif // THRESHOLD_FROM_PROPERTY

    if ( m_monitoring ) {
      ++m_hist_mainChargeAfter[chargeMain];
      ++m_hist_sensorVsPixelChargeAfter[{to_unsigned( channel.sensor() ), chargeMain}];
    }
    if ( m_maskedPixels ) {
      // Draw a random number to choose if this pixel is masked.
      if ( m_uniform() < m_fractionMasked ) {
        ++numDigitsKilledMasked;
        continue;
      }
    }
    // Create a new digit.
    digits.insert( new LHCb::MPDigit{channel} );
    ++numSignalDigitsCreated; // all digits from tracks (in any bunch)
  }
  if ( m_noisyPixels ) {
    const auto numToAdd = m_poisson();
    // Add additional digits without a corresponding MC digit.
    for ( unsigned i = 0; i < numToAdd; ++i ) {
      // Sample the layer, module (mod), ladder, sensor, column and row of the noise pixel.
      const auto layer =
          LHCb::Detector::MPChannelID::Layer( std::floor( m_uniform() * LHCb::Detector::MPChannelID::nLayers<>() ) );
      const auto mod =
          LHCb::Detector::MPChannelID::Module( std::floor( m_uniform() * LHCb::Detector::MPChannelID::nModules<>() ) );
      const auto stave =
          LHCb::Detector::MPChannelID::Stave( std::floor( m_uniform() * LHCb::Detector::MPChannelID::nStaves<>() ) );
      const auto sensor =
          LHCb::Detector::MPChannelID::Sensor( std::floor( m_uniform() * LHCb::Detector::MPChannelID::nSensors<>() ) );
      const auto col =
          LHCb::Detector::MPChannelID::Column( std::floor( m_uniform() * LHCb::Detector::MPChannelID::nColumns<>() ) );
      const auto row =
          LHCb::Detector::MPChannelID::Row( std::floor( m_uniform() * LHCb::Detector::MPChannelID::nRows<>() ) );
      const LHCb::Detector::MPChannelID channel{
          LHCb::Detector::MPChannelID::ForVersion<>{}, layer, mod, stave, sensor, col, row};
      // Make sure the channel has not yet been added to the container.
      if ( digits.object( channel ) ) continue;
      digits.insert( new LHCb::MPDigit{channel} );
      ++numNoiseDigitsCreated;
    }
  }

  // add the values from this event to the running totals
  m_numSpillover += numSpillover;
  m_numLostTiming += numLostTiming;
  m_numDigitsKilledMasked += numDigitsKilledMasked;
  m_numSignalDigitsCreated += ( numSignalDigitsCreated - numSpillover );
  m_numNoiseDigitsCreated += numNoiseDigitsCreated;

  std::sort( digits.begin(), digits.end(),
             []( const LHCb::MPDigit* a, const LHCb::MPDigit* b ) { return a->channelID() < b->channelID(); } );
  return digits;
}
