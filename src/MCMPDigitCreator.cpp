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
/** @class MCFTDigitCreator MCFTDigitCreator.h
 *
 *  From the list of MCFTDeposit (deposited energy in each FTChannel),
 *  this algorithm converts the deposited energy in ADC charge according
 *  to the mean number of photoelectron per MIP
 *  + the mean ADC count per photoelectron.
 *  Created digits are put in transient data store.
 *
 *  TO DO :
 *  - add noise
 *
 *  THINK ABOUT :
 *  - dead channels
 *
 *  @author COGNERAS Eric, Luca Pescatore
 *  @date   2012-04-04
 */

// from Gaudi
#include "GaudiKernel/RndmGenerators.h"
#include "LHCbAlgs/Transformer.h"

#include "Detector/MP/MPChannelID.h"

#include "Event/MCMPDeposit.h"
#include "Event/MCMPDigit.h"

class MCMPDigitCreator : public LHCb::Algorithm::Transformer<LHCb::MCMPDigits( const LHCb::MCMPDeposits& )> {

public:
  /// Standard constructor

  MCMPDigitCreator( const std::string& name, ISvcLocator* pSvcLocator );

  LHCb::MCMPDigits operator()( const LHCb::MCMPDeposits& deposits ) const override;

private:
};

//-----------------------------------------------------------------------------
// Implementation file for class : MCMPDigitCreator
//
// 2012-04-04 : COGNERAS Eric
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( MCMPDigitCreator )

MCMPDigitCreator::MCMPDigitCreator( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer( name, pSvcLocator, KeyValue{"InputLocation", LHCb::MCMPDepositLocation::Default},
                   KeyValue{"OutputLocation", LHCb::MCMPDigitLocation::Default} ) {}

//=============================================================================
// Main execution
//=============================================================================
LHCb::MCMPDigits MCMPDigitCreator::operator()( const LHCb::MCMPDeposits& deposits ) const {
  using LHCb::MCMPDeposit;
  using LHCb::MCMPDeposits;
  using LHCb::MCMPDigit;
  using LHCb::MCMPDigits;
  using LHCb::Detector::MPChannelID;
  // Define digits container and register it in the transient data store
  MCMPDigits digits;
  digits.reserve( deposits.size() );

  auto       it   = deposits.begin();
  const auto last = deposits.end();
  while ( last != it ) {
    const auto first     = it;
    const auto channelID = ( *it )->channelID();
    it = std::find_if( first + 1, last, [channelID]( const auto el ) -> bool { return channelID != el->channelID(); } );
    // make an MCMPDigit
    MCMPDigit* digit = new MCMPDigit( first, it );
    digits.insert( digit );
  }

  const auto lessByMPChannelID = []( const MCMPDigit* d1, const MCMPDigit* d2 ) -> bool {
    return d1->channelID() < d2->channelID();
  };

  // Digits are sorted according to their ChannelID to prepare the clusterisation
  std::stable_sort( digits.begin(), digits.end(), lessByMPChannelID );

  return digits;
}
