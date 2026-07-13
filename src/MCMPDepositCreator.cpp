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

#include <cassert>
#include <cmath>
// #include <cstdio>
#include <boost/container/flat_map.hpp>
#include <limits>
#include <numeric>
#include <string_view>
#include <unordered_map>
#include <vector>
#include <yaml-cpp/yaml.h>

// form Gaudi
#include "GaudiAlg/GaudiTupleAlg.h"
#include "GaudiKernel/RndmGenerators.h"
#include "LHCbAlgs/Transformer.h"
#include <Gaudi/Accumulators/Histogram.h>
#include <GaudiKernel/PhysicalConstants.h>

#include "Detector/MP/MPChannelID.h"
#include "Event/MCMPDeposit.h"
#include "Event/MCMPDigit.h"
#include "Event/MPDigit.h"
#include "LHCbMath/LHCbMath.h"
#include <Event/GenHeader.h>
#include <Event/MCHit.h>

using namespace LHCb;

// [FIXME] fake DeMP class for now
class DeMP {};
// [FIXME] fake DeMPLocation for now
namespace DeMPLocation {
  constexpr std::string_view Default{"DeMPLocation::Default"};
}

namespace MCMPDepositCreatorConditions {
  struct ConditionsCache {
    // [FIXME] getting diffusion coefficient and dead time parameter somewhere
  };
} // namespace MCMPDepositCreatorConditions

/**
 *  Algorithm for creating MCMPDeposits from MCHits. The following processes
 *  are simulated:
 */
class MCMPDepositCreator
    : public LHCb::Algorithm::Transformer<
          LHCb::MCMPDeposits( const LHCb::MCHits& mcPrevPrev, const LHCb::MCHits& mcPrev, const LHCb::MCHits& mcMain,
                              const LHCb::MCHits& mcNext, const LHCb::GenHeader& genHeader
#if 0
                    // [FIXME] get detector to work
                    , const DeMP&
#endif
#if 0
                    // [FIXME] get this to work
                    , const MCMPDepositCreatorConditions::ConditionsCache&
#endif
                              ),
          LHCb::DetDesc::usesBaseAndConditions<GaudiTupleAlg, DeMP, MCMPDepositCreatorConditions::ConditionsCache>> {

public:
  MCMPDepositCreator( const std::string& name, ISvcLocator* pSvcLocator );
  StatusCode         initialize() override; ///< Algorithm initialization
  LHCb::MCMPDeposits operator()( const LHCb::MCHits& mcPrevPrev, const LHCb::MCHits& mcPrev, const LHCb::MCHits& mcMain,
                                 const LHCb::MCHits& mcNext, const LHCb::GenHeader& genHeader
#if 0
                                 // [FIXME] get detector to work
                                 , const DeMP&
#endif
#if 0
                                 // [FIXME] get this to work
                                 , const MCMPDepositCreatorConditions::ConditionsCache&
#endif
                                 ) const override;

private:
  using PixelMap = boost::container::flat_map<LHCb::Detector::MPChannelID, float>;
  using PreDeposits =
      std::unordered_map<LHCb::Detector::MPChannelID,
                         std::tuple<SmartRefVector<LHCb::MCHit>, std::vector<float>, std::vector<float>>>;

  void createPreDeposits( const LHCb::MCHit& hit, const PixelMap& pixels, PreDeposits& predeposits,
                          const float hitTimeOffset ) const;
  void createCharges( const LHCb::MCHit& hit, const double path, const double charge,
                      std::vector<double>& charges ) const;
  void createPixels( const LHCb::MCHit& hit
#if 0
            // [FIXME] get detector working
            , const DeMP& det
#endif
                     ,
                     const std::vector<double>& charges, PixelMap& pixels ) const;

  /// Draw random number from 1/q^2 distribution
  double randomTail( const double qmin, const double qmax ) const;
  /// "fake" an MPChannelID without detector geometry
  static LHCb::Detector::MPChannelID p3d2MPID( double x, double y, double z ) noexcept;

  /// Associated time offset of input containers (ns)
  Gaudi::Property<std::vector<float>> m_hitTimeOffset{
      this,
      "HitTimeOffsets",
      {-50. * Gaudi::Units::ns, -25.0 * Gaudi::Units::ns, 0. * Gaudi::Units::ns, 25. * Gaudi::Units::ns}};
  /// which MCHit containers to include
  Gaudi::Property<std::vector<bool>> m_useMCHitContainer{this, "useMCHitContainer", {true, true, true, true}};

  /// Conversion factor between energy deposit and number of electrons
  Gaudi::Property<double> m_eVPerElectron{this, "eVPerElectron", 3.6};

  /// Max. number of points on hit trajectory
  Gaudi::Property<unsigned> m_nMaxSteps{this, "MaxNumSteps", 150};

  /// Distance between points on hit trajectory
  Gaudi::Property<double> m_stepSize{this, "StepSize", 5 * Gaudi::Units::micrometer};

  /// Number of electrons per micron (MPV)
  Gaudi::Property<double> m_chargeUniform{this, "ChargeUniform", 70.0};

  /// Lower limit of 1/q^2 distribution
  Gaudi::Property<double> m_minChargeTail{this, "MinChargeTail", 10.0};

  /// Temperature of the sensor [K]
  Gaudi::Property<double> m_temperature{this, "Temperature", 253.15 * Gaudi::Units::kelvin};

  /// Discrimination threshold in electrons, for time walk
  Gaudi::Property<double> m_threshold{this, "ChargeThreshold", 1000.0};

  Gaudi::Property<double> m_minElectrons{this, "minElectrons", 100.};
  Gaudi::Property<double> m_minPathLen{this, "minPathLen", 1e-6 * Gaudi::Units::mm};
  Gaudi::Property<double> m_sensorThickness{this, "sensorThickness", .2 * Gaudi::Units::mm};
  Gaudi::Property<double> m_highVoltage{this, "highVoltage", 100. * Gaudi::Units::volt};
  // split charges into nSplit packets for diffusion simulation
  Gaudi::Property<unsigned> m_nSplit{this, "nSplit", 5u};

  mutable Rndm::Numbers m_gauss;
  mutable Rndm::Numbers m_uniform;
};

DECLARE_COMPONENT( MCMPDepositCreator )

MCMPDepositCreator::MCMPDepositCreator( const std::string& name, ISvcLocator* pSvcLocator )
    : Transformer{name,
                  pSvcLocator,
                  {
                      KeyValue{"InputLocationPrevPrev", "/Event/PrevPrev/" + LHCb::MCHitLocation::MP},
                      KeyValue{"InputLocationPrev", "/Event/Prev/" + LHCb::MCHitLocation::MP},
                      KeyValue{"InputLocationMain", LHCb::MCHitLocation::MP},
                      KeyValue{"InputLocationNext", "/Event/Next/" + LHCb::MCHitLocation::MP},
                      KeyValue{"InputLocationGenHeader", LHCb::GenHeaderLocation::Default},
#if 0
                     // [FIXME] Detector element pending...
                     KeyValue{ "DeMPLocation", DeMPLocation::Default },
#endif
#if 0
                     // [FIXME] get this to work
                     KeyValue{ "ConditionsCache", "AlgorithmSpecific-" + name + "-ConditionsCache" }
#endif
                  },
                  KeyValue{"OutputLocation", LHCb::MCMPDepositLocation::Default}} {
}

StatusCode MCMPDepositCreator::initialize() {
  return Transformer::initialize()
      .andThen( [&] {
        bool properties_ok = true;
        if ( m_hitTimeOffset.size() != 4 ) {
          error() << "Must have four HitTimeOffsets";
          properties_ok = false;
        }
        if ( m_useMCHitContainer.size() != 4 ) {
          error() << "Must have four entries in useMCHitContainer";
          properties_ok = false;
        }
        const auto propertyFinitePositive = [&]( const Gaudi::Property<double>& prop ) {
          if ( !std::isfinite( double( prop ) ) || double( prop ) <= 0. ) {
            error() << prop.name() << " must be finite and positive";
            return false;
          }
          return true;
        };
        const std::vector<std::reference_wrapper<const Gaudi::Property<double>>> dbl_properties_to_check{
            m_eVPerElectron, m_stepSize,     m_chargeUniform, m_minChargeTail,   m_temperature,
            m_threshold,     m_minElectrons, m_minPathLen,    m_sensorThickness, m_highVoltage};
        for ( const auto& prop : dbl_properties_to_check ) {
          if ( !propertyFinitePositive( prop ) ) properties_ok = false;
        }
        const auto propertyNonZero = [&]( const Gaudi::Property<unsigned>& prop ) {
          if ( unsigned( prop ) == 0 ) {
            error() << prop.name() << " must be non-zero";
            return false;
          }
          return true;
        };
        const std::vector<std::reference_wrapper<const Gaudi::Property<unsigned>>> unsigned_properties_to_check{
            m_nMaxSteps, m_nSplit};
        for ( const auto& prop : unsigned_properties_to_check ) {
          if ( !propertyNonZero( prop ) ) properties_ok = false;
        }
        if ( !properties_ok ) {
          return StatusCode::FAILURE;
        } else {
          return StatusCode::SUCCESS;
        }
      } )
      .andThen( [&] { return m_gauss.initialize( randSvc(), Rndm::Gauss( 0.0, 1.0 ) ); } )
      .andThen( [&] { return m_uniform.initialize( randSvc(), Rndm::Flat( 0.0, 1.0 ) ); } )
#if 0
    // [FIXME] add per-sensor/... conditions at some later point, just a
    // single number for now
    .andThen( [&] {
        addConditionDerivation(
			       {inputLocation<DeMP>(), "Conditions/Calibration/MP/MPDeadTimeParam"},
			       inputLocation<MPDepositCreatorConditions::ConditionsCache>(),
			       [=]( DeMP const& det, YAML::Node const& deadTimeParam ) -> MPDepositCreatorConditions::ConditionsCache {
				 // Calculate the diffusion coefficient.
				 const double           kt = 2. * m_temperature * Gaudi::Units::k_Boltzmann / Gaudi::Units::eV;
				 std::array<double, 52> tmpDiffusionCoefficient{};
				 det.runOnAllSensors( [&]( const DeMPSensor& sensor ) {
				     tmpDiffusionCoefficient[sensor.module()] =
				       sqrt( kt * sensor.siliconThickness() / ( sensor.sensorHV() ) );
				   } );
				 return {tmpDiffusionCoefficient,
				     deadTimeParam["FunctionParam"].as<std::array<std::array<double, 4>, 52>>()};
				 } );
			       setHistoTopDir( "MP/" );
			       } )
#endif
      ;
}

LHCb::MCMPDeposits MCMPDepositCreator::operator()( const LHCb::MCHits& mcPrevPrev, const LHCb::MCHits& mcPrev,
                                                   const LHCb::MCHits& mcMain, const LHCb::MCHits& mcNext,
                                                   const LHCb::GenHeader& genHeader
#if 0
                                // [FIXME] get detector to work
                                , const DeMP& det
#endif
#if 0
                                // [FIXME] get the conditions cache to work
                                , const MCMPDepositCreatorConditions::ConditionsCache& conditions
#endif
                                                   ) const {
  // figure out over which containers of MCHits we should run
  using MCHitsAndTimeOffset = std::pair<const LHCb::MCHits*, float>;
  std::array<MCHitsAndTimeOffset, 4> mchitcontainers{
      {{&mcPrevPrev, m_useMCHitContainer[0] ? m_hitTimeOffset[0] : std::numeric_limits<float>::quiet_NaN()},
       {&mcPrev, m_useMCHitContainer[1] ? m_hitTimeOffset[1] : std::numeric_limits<float>::quiet_NaN()},
       {&mcMain, m_useMCHitContainer[2] ? m_hitTimeOffset[2] : std::numeric_limits<float>::quiet_NaN()},
       {&mcNext, m_useMCHitContainer[3] ? m_hitTimeOffset[3] : std::numeric_limits<float>::quiet_NaN()}}};
  const auto last_mchitcontainer = std::remove_if( mchitcontainers.begin(), mchitcontainers.end(),
                                                   []( const auto& el ) { return !std::isfinite( el.second ); } );

  PreDeposits predeposits;
  // reserve adequate space in predeposits
  predeposits.reserve(
      std::accumulate( mchitcontainers.begin(), last_mchitcontainer, std::size_t( 0 ),
                       []( std::size_t tmp, const auto& cont ) { return tmp + cont.first->size(); } ) );
  // scratch area for createDeposits
  std::vector<double> charges;
  charges.reserve( m_nMaxSteps );
  PixelMap pixels;
  pixels.reserve( m_nSplit * m_nMaxSteps );
  // create deposits for all containers
  for ( auto it = mchitcontainers.begin(); last_mchitcontainer != it; ++it ) {
    const LHCb::MCHits& hits          = *( it->first );
    const auto          hitTimeOffset = it->second;
    for ( const auto& hit : hits ) {
      // Calculate the total number of electrons based on the G4 energy deposit.
      const double charge = ( hit->energy() / Gaudi::Units::eV ) / m_eVPerElectron;
      // Skip very small deposits.
      if ( charge < m_minElectrons ) continue;
      const double path = hit->pathLength();
      // Skip very small path lengths.
      if ( path < m_minPathLen ) {
        warning() << "Path length in active silicon: " << path << endmsg;
        continue;
      }
      createCharges( *hit, path, charge, charges );

      createPixels( *hit
#if 0
              // [FIXME] get detector working
              , det
#endif
                    ,
                    charges, pixels );

      createPreDeposits( *hit, pixels, predeposits, hitTimeOffset );
    }
  }

  // [FIXME] if you want to simulate dead time, this goes here...

  // convert predeposits to deposits
  LHCb::MCMPDeposits deposits;
  deposits.reserve( predeposits.size() );
  for ( auto&& el : predeposits ) {
    auto dep =
        std::make_unique<MCMPDeposit>( el.first, std::move( std::get<0>( el.second ) ),
                                       std::move( std::get<1>( el.second ) ), std::move( std::get<2>( el.second ) ) );
    deposits.insert( dep.get() );
    // successfully added to deposits, so drop dep's ownership of the
    // MCMPDeposit
    static_cast<void>( dep.release() );
  }
  return deposits;
}

void MCMPDepositCreator::createCharges( const LHCb::MCHit& mchit, const double path, const double charge,
                                        std::vector<double>& charges ) const {
  charges.clear();
  // Calculate the number of points to distribute the deposited charge on.
  unsigned nPoints = std::ceil( path / m_stepSize );
  if ( nPoints > m_nMaxSteps ) nPoints = m_nMaxSteps;
  // Calculate the charge on each point.
  const double mpv   = m_chargeUniform * ( path / Gaudi::Units::micrometer ) / nPoints;
  const double sigma = std::sqrt( std::abs( mpv ) );
  std::generate_n( std::back_inserter( charges ), nPoints, [mpv, sigma, this]() {
    double q;
    do { q = mpv + sigma * m_gauss(); } while ( q <= 0. ); // be paranoid
    return q;
  } );
  auto sum = std::accumulate( charges.begin(), charges.end(), 0. );

  while ( charge > sum + m_minChargeTail ) {
    // Add additional charge sampled from an 1 / n^2 distribution.
    const double   q = randomTail( m_minChargeTail, charge - sum );
    const unsigned i = unsigned( LHCb::Math::round( m_uniform() * ( nPoints - 1 ) ) );
    charges[i] += q;
    sum += q;
  }
  // Calculate scaling factor to match total deposited charge.
  const double scale = charge / sum;
  std::transform( charges.begin(), charges.end(), charges.begin(), [scale]( const auto q ) { return q * scale; } );
}

void MCMPDepositCreator::createPixels( const LHCb::MCHit& hit
#if 0
        // [FIXME] get detector working
        , const DeMP& det
#endif
                                       ,
                                       const std::vector<double>& charges, PixelMap& pixels ) const {
  pixels.clear();
  const auto nPoints = charges.size();
  // [FIXME] get sensor, and its thickness from detector description
  const double thickness = m_sensorThickness;

  const double diffusionCoefficient =
      std::sqrt( thickness / m_highVoltage * 2. * m_temperature * Gaudi::Units::k_Boltzmann / Gaudi::Units::eV );
  // Calculate the distance between two points on the trajectory.
  // [FIXME] This should be in the local frame
  const Gaudi::XYZPoint entry  = ( hit.entry().z() < hit.exit().z() ) ? hit.entry() : hit.exit();
  const Gaudi::XYZPoint exit   = ( hit.entry().z() < hit.exit().z() ) ? hit.exit() : hit.entry();
  const auto            nSplit = m_nSplit;
  for ( unsigned i = 0; nPoints != i; ++i ) {
    const double          f = ( i + 0.5 ) / double( nPoints );
    const Gaudi::XYZPoint point{f * entry.x() + ( 1. - f ) * exit.x(), f * entry.y() + ( 1. - f ) * exit.y(),
                                f * entry.z() + ( 1. - f ) * exit.z()};
    const auto            dzstrip = ( point.z() - entry.z() ) / ( exit.z() - entry.z() );
    const auto            sigmaD  = diffusionCoefficient * std::sqrt( dzstrip );
    const auto            q       = charges[i] / double( nSplit );
    for ( unsigned j = 0; nSplit != j; ++j ) {
      Gaudi::XYZPoint point_smeared( point.x() + sigmaD * m_gauss(), point.y() + sigmaD * m_gauss(), point.z() );
      LHCb::Detector::MPChannelID channelID =
          p3d2MPID( point_smeared.x(), point_smeared.y(), point_smeared.z() ); // [FIXME]
      if ( !true ) { // [FIXME] !isInActiveArea(point_smeared)) {
        if ( msgLevel( MSG::VERBOSE ) )
          verbose() << "Killed MCMPdeposit in channel " << channelID << " due to bad pixel" << endmsg;
        continue;
      }
      auto [it, ok] = pixels.emplace( channelID, q ); // [FIXME]
      if ( !ok ) it->second += q;
    }
  }
}

void MCMPDepositCreator::createPreDeposits( const LHCb::MCHit& hit, const PixelMap& pixels,
                                            MCMPDepositCreator::PreDeposits& predeposits, const float timeOff ) const {
  for ( const auto& [id, charge] : pixels ) {
    // work out time of flight from origin to mid-point of deposit for an
    // infinite momentum particle
    const Gaudi::XYZPoint point{0.5 * ( hit.entry().x() + hit.exit().x() ), 0.5 * ( hit.entry().y() + hit.exit().y() ),
                                0.5 * ( hit.entry().z() + hit.exit().z() )};
    const auto dist  = std::sqrt( float( point.x() * point.x() + point.y() * point.y() + point.z() * point.z() ) );
    const auto t_tof = dist / float( Gaudi::Units::c_light );
    auto       time  = timeOff + static_cast<float>( hit.time() ) - t_tof;

    // FIXME: all Velo
    // add time offset due to timewalk from Larissa's talk 25/5/18
    // using testpulses the function is
    // delta(T) = p0/(Q-p1)**p2 + p3
    // p1 approx threshold level and p3 corrected for in time alignment
    // p0 = 5*25 in ns (plot in BX)
    // p2 = 0.5 mostly (see page 22, some pixels at 0.8) ignore for now
    // scale to e- i.e scale by 1/72 (1000e- -> 13.5 in testpulse a.u.)
    // Then subtract 2ns flat due to pixel to pixel variation
    float timeWalk = 999.; // default to not appearing in next two BX if too small
    if ( charge > m_threshold ) {
      timeWalk = ( 5. * 25. ) / std::sqrt( ( charge - m_threshold ) / 72. ) - 2. * m_uniform();
    }
    time += timeWalk; // apply time walk

    const LHCb::MCHit* mchit =
        ( std::abs( timeOff ) < 12.5 * Gaudi::Units::ns ) ? &hit : static_cast<const LHCb::MCHit*>( nullptr );

    auto& [mchits, charges, times] = predeposits[id];
    mchits.push_back( mchit );
    charges.push_back( charge );
    times.push_back( time );
  }
}

double MCMPDepositCreator::randomTail( const double qmin, const double qmax ) const {
  // Sample charge from 1 / n^2 distribution.
  const double offset = 1. / qmax;
  const double range  = ( 1. / qmin ) - offset;
  const double u      = offset + m_uniform() * range;
  return 1. / u;
}

/** @brief turn a 3d point into a MPChannelID
 *
 * This routines builds periodic tiles of modules around the beam pipe which
 * extend to +/- infinity in x and y. Row and sensor number increase with
 * increasing x, the lowest bit of the stave number and the low two bits
 * increase with increasing x. Column number, the remaining stave number
 * bits, and the remaining module number bits increase with increasing Y.
 * Layer numbers increase with increasing z.
 *
 * @note This routine is a quick and dirty hack to get something to work
 * quickly to enable further development. This should be replaced asap with
 * routines that do things properly, and position sensors, staves, modules
 * using the geometry.
 *
 * @author Manuel Schiller <Manuel.Schiller@glasgow.ac.uk>
 * @date 2022-10-24
 */
LHCb::Detector::MPChannelID MCMPDepositCreator::p3d2MPID( double x, double y, double z ) noexcept {
  using LHCb::Detector::MPChannelID;
  constexpr unsigned nRows     = 320u;
  constexpr unsigned nSensors  = 14u;
  constexpr unsigned nStavesX  = 2u;
  constexpr unsigned nModulesX = 4u;

  constexpr unsigned nColumns  = 116u;
  constexpr unsigned nStavesY  = 10u;
  constexpr unsigned nModulesY = 7u;

  constexpr double x_pitch  = 0.056; // mm
  constexpr double y_pitch  = 0.165; // mm
  constexpr double x_period = 0.5 * x_pitch * nRows * nSensors * nStavesX * nModulesX;
  constexpr double y_period = 0.5 * y_pitch * nColumns * nStavesY * nModulesY;
  // std::printf("[DEBUG] xperiod %+8.2f %+8.2f yperiod %+8.2f %+8.2f\n",
  //         -x_period, +x_period, -y_period, +y_period);
  //  sqash into -x_period <= x < x_period, -y_period <= y < y_period
  // std::printf("[DEBUG] before (%g, %g, %g)\n", x, y, z);
  x = fmod( x, 2. * x_period );
  y = fmod( y, 2. * y_period );
  // std::printf("[DEBUG] during (%g, %g, %g)\n", x, y, z);
  if ( x < x_period ) x += 2. * x_period;
  if ( y < y_period ) y += 2. * y_period;
  if ( x_period <= x ) x -= 2. * x_period;
  if ( y_period <= y ) y -= 2. * y_period;
  // std::printf("[DEBUG]  after (%g, %g, %g)\n", x, y, z);
  //  calculate pixel number in x and y (-x_period, -y_period)
  unsigned ix = ( x + x_period ) / x_pitch;
  unsigned iy = ( y + y_period ) / y_pitch;
  assert( ix < nRows * nSensors * nStavesX * nModulesX );
  assert( iy < nColumns * nStavesY * nModulesY );
  // std::printf("[DEBUG] ix=%9u iy=%9u\n", ix, iy);
  //  convert to MPChannelID
  const unsigned row = ix % nRows;
  ix /= nRows;
  const unsigned column = iy % nColumns;
  iy /= nColumns;
  const unsigned sensor = ix % nSensors;
  ix /= nSensors;
  const unsigned stave = ( ix % nStavesX ) + nStavesX * ( iy % nStavesY );
  ix /= nStavesX;
  iy /= nStavesY;
  const unsigned module = ( ix % nModulesX ) + nModulesX * ( iy % nModulesY );
  ix /= nModulesX;
  iy /= nModulesY;
  assert( 0 == ix );
  assert( 0 == iy );
  // std::printf(
  //         "[DEBUG] row %3u column %3u sensor %2u stave %2u module %2u\n",
  //         row, column, sensor, stave, module);
  z = ( z - 7810. ) / 160.;
  if ( z < 0 ) z = 0;
  if ( z > 10 ) z = 10;
  // assert(0 <= z && z < 10);
  unsigned iz = unsigned( z );
  iz -= iz / 2;
  const unsigned layer = iz;
  // std::printf("[DEBUG] iz=%u\n", iz);
  return MPChannelID( MPChannelID::ForVersion<>{}, MPChannelID::Layer{layer}, MPChannelID::Module{module},
                      MPChannelID::Stave{stave}, MPChannelID::Sensor{sensor}, MPChannelID::Column{column},
                      MPChannelID::Row{row} );
}
