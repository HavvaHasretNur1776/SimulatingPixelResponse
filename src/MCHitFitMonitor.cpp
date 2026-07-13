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
#include "GaudiKernel/RndmGenerators.h"
//tuples
#include "GaudiAlg/GaudiTupleAlg.h"

// from Event
#include "Event/MCHit.h"

//From PolynomialFitter
//#include "PolynomialFitter/PolyFitter.h"
#include "include/PolyFitter.h"

#include <unordered_map>
#include <random>

/** @class MCHitFitMonitor
 *
 *  @author Havva Hasret Nur
 *  @date   2024-02-06
 */
class MCHitFitMonitor : public LHCb::Algorithm::Consumer<void( const LHCb::MCHits& ), Gaudi::Functional::Traits::BaseClass_t<GaudiTupleAlg>> {

public:
  MCHitFitMonitor( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode initialize() override;
  void       operator()( const LHCb::MCHits& deposits ) const override;

private:
  mutable Rndm::Numbers m_gauss;

class TupleFiller{
  private:
    std::vector<float> m_xHit;
    std::vector<float> m_yHit;
    std::vector<float> m_zHit;
    std::vector<float> m_xHitTrue;
    std::vector<float> m_yHitTrue;
    std::vector<float> m_zHitTrue;
    std::vector<float> m_sigmaxHit;
    std::vector<float> m_sigmayHit;
    std::vector<float> m_x_fitted;
    std::vector<float> m_y_fitted;
    std::vector<float> m_x_resi;
    std::vector<float> m_y_resi;
    std::vector<float> m_sigmaX_gaus;
    std::vector<float> m_sigmaY_gaus;
    
    ///additional branches
    std::vector<float> m_paramsX;
    std::vector<float> m_paramsY;
    std::vector<float> m_covarianceX;
    std::vector<float> m_covarianceY;
    std::vector<float> m_matrixX;
    std::vector<float> m_matrixY;
    std::vector<float> m_rhsX;
    std::vector<float> m_rhsY;

    float m_chi2X;
    int m_ndfX;
    float m_chi2Y;
    int m_ndfY;
    float m_probX;
    float m_probY;
    int m_track_number;
    std::size_t m_maxPerTrack;//Max number of hits
    std::size_t m_MatrixSize; //size of matrix/cov: nx(n+1)/2
    std::size_t m_VectorSize; //size of params/rhs : n  

  public:
    TupleFiller(std::size_t maxPerTrack, std::size_t MatrixSize, std::size_t VectorSize) : 
    m_maxPerTrack(maxPerTrack), m_MatrixSize(MatrixSize), m_VectorSize(VectorSize)
    {
        m_xHit.reserve(m_maxPerTrack);
        m_yHit.reserve(m_maxPerTrack);
        m_zHit.reserve(m_maxPerTrack);
        m_xHitTrue.reserve(m_maxPerTrack);
        m_yHitTrue.reserve(m_maxPerTrack);
        m_zHitTrue.reserve(m_maxPerTrack);
        m_sigmaxHit.reserve(m_maxPerTrack);
        m_sigmaxHit.reserve(m_maxPerTrack);
        m_x_fitted.reserve(m_maxPerTrack);
        m_y_fitted.reserve(m_maxPerTrack);
        m_x_resi.reserve(m_maxPerTrack);
        m_y_resi.reserve(m_maxPerTrack);
        // MS20250206: you have three sizes here - per hit things that go
        // with the number of hits, per track things that go with the number
        // of parameters, and per track things that go with the matrix size
        // (n * (n + 1) / 2), and the reserve call should reserve the right
        // amount of space, and the tuple needs a branch for each of these
        // sizes. -> HASRET: why is that an issue ? it is just max number and already within range for nparams and nmatrix
        m_paramsX.reserve(m_VectorSize);
        m_paramsY.reserve(m_VectorSize);
        m_matrixX.reserve(m_MatrixSize);
        m_matrixY.reserve(m_MatrixSize);//it should be nx(n+1)/2
        std::cout<<"m_MatrixSize: "<<m_MatrixSize<<std::endl;
        m_covarianceX.reserve(m_MatrixSize);
        m_covarianceY.reserve(m_MatrixSize);
        m_rhsX.reserve(m_VectorSize);
        m_rhsY.reserve(m_VectorSize);
        m_sigmaX_gaus.reserve(m_maxPerTrack);
        m_sigmaY_gaus.reserve(m_maxPerTrack);
               
    }

    void clear()
    {
        m_xHit.clear();
        m_yHit.clear();
        m_zHit.clear();
        m_xHitTrue.clear();
        m_yHitTrue.clear();
        m_zHitTrue.clear();
        m_sigmaxHit.clear();
        m_sigmayHit.clear();
        m_x_fitted.clear();
        m_y_fitted.clear();
        m_x_resi.clear();
        m_y_resi.clear();
        m_paramsX.clear();
        m_paramsY.clear();
        m_matrixX.clear();
        m_matrixY.clear();
        m_covarianceX.clear();
        m_covarianceY.clear();
        m_rhsX.clear();
        m_rhsY.clear();
        m_sigmaX_gaus.clear();
        m_sigmaY_gaus.clear();
       
    
    }

    void fillPerHit(float x, float y, float z, float x_t, float y_t, float z_t , float sigmax, float sigmay, float x_f, float y_f, float x_r, float y_r, float sigmax_gaus, float sigmay_gaus)
    {
        if (m_xHit.size() >= m_maxPerTrack) {
            std::fprintf(stderr, "too many hits in %s\n", __func__);
            return;
        }
        
        m_xHit.push_back(x);
        m_yHit.push_back(y);
        m_zHit.push_back(z);
        m_xHitTrue.push_back(x_t);
        m_yHitTrue.push_back(y_t);
        m_zHitTrue.push_back(z_t);
        m_sigmaxHit.push_back(sigmax);
        m_sigmayHit.push_back(sigmay);
        m_x_fitted.push_back(x_f);
        m_y_fitted.push_back(y_f);
        m_x_resi.push_back(x_r);
        m_y_resi.push_back(y_r);
        m_sigmaX_gaus.push_back(sigmax_gaus);
        m_sigmaY_gaus.push_back(sigmay_gaus);
    }

    void fillPerTrackX(float chi2, float ndf, float prob, std::vector<float> a_n, std::vector<float> rhs, std::vector<float> matrix, std::vector<float> cov)
    {
        m_chi2X = chi2;
        m_ndfX = ndf;
        m_probX = prob;
        m_paramsX = a_n;
        m_rhsX = rhs;
        m_matrixX = matrix;
        m_covarianceX = cov;
    }
    void fillPerTrackY(float chi2, float ndf, float prob, std::vector<float> a_n, std::vector<float> rhs, std::vector<float> matrix, std::vector<float> cov)
    {
        m_chi2Y = chi2;
        m_ndfY = ndf;
        m_probY = prob;
        m_paramsY = a_n;
        m_rhsY = rhs;
        m_matrixY = matrix;
        m_covarianceY = cov;

    }

    template <typename T>
    void writeToTuple(T& tuple)
    {
        tuple->farray("xHit", m_xHit, "nhits", m_maxPerTrack);
        tuple->farray("yHit", m_yHit,  "nhits", m_maxPerTrack);
        tuple->farray("zHit", m_zHit, "nhits", m_maxPerTrack);
        tuple->farray("xHitTrue", m_xHitTrue, "nhits", m_maxPerTrack);
        tuple->farray("yHitTrue", m_yHitTrue,  "nhits", m_maxPerTrack);
        tuple->farray("zHitTrue", m_zHitTrue, "nhits", m_maxPerTrack);
        tuple->farray( "sigmaxHit", m_sigmaxHit,"nhits", m_maxPerTrack);
        tuple->farray( "sigmayHit", m_sigmayHit,"nhits", m_maxPerTrack);
        tuple->farray( "xHitFitted", m_x_fitted,"nhits", m_maxPerTrack);
        tuple->farray( "yHitFitted", m_y_fitted,"nhits", m_maxPerTrack);
        tuple->farray("X_Residual", m_x_resi, "nhits", m_maxPerTrack);
        tuple->farray("Y_Residual", m_y_resi, "nhits", m_maxPerTrack);
        tuple->farray("sigmaX_gaus", m_sigmaX_gaus, "nhits", m_maxPerTrack);
        tuple->farray("sigmaY_gaus", m_sigmaY_gaus, "nhits", m_maxPerTrack);
       
        //additional branches
        tuple->farray("paramsX", m_paramsX, "nparamsX", m_VectorSize);
        tuple->farray("paramsY", m_paramsY, "nparamsY", m_VectorSize);
        tuple->farray("rhsX", m_rhsX, "nparamsX", m_VectorSize);
        tuple->farray("matrixX", m_matrixX, "nMatX", m_MatrixSize);
        //std::cout<< "HASRET COV SIZE "<< m_covarianceX.size()<<std::endl;
        //std::cout<<"m_MatrixSize: "<<m_MatrixSize<<std::endl;
        tuple->farray("covX", m_covarianceX, "nMatX", m_MatrixSize);
        //std::cout<<"m_covarianceY size: "<<m_MatrixSize<<std::endl;
        tuple->farray("rhsY", m_rhsY, "nparamsY", m_VectorSize);
        tuple->farray("matrixY", m_matrixY, "nMatY", m_MatrixSize);
        tuple->farray("covY", m_covarianceY, "nMatY", m_MatrixSize);

        tuple->column( "chi2X",m_chi2X);
        tuple->column("ndfX",m_ndfX);
        tuple->column( "chi2Y",m_chi2Y);
        tuple->column("ndfY",m_ndfY);
        tuple->column( "Chi2ProbX", m_probX);
        tuple->column( "Chi2ProbY", m_probY);
        //tuple->column("TrackNumber", m_track_number);
        
        tuple->write();
    }
};//tuple subclass


  using Hits = std::vector<const LHCb::MCHit*>;
  using HitsPerMCPart = std::unordered_map<const LHCb::MCParticle*, Hits>;
  struct HitWorkspace {
      std::vector<double> m_x; // 3 pointers, each 8 bytes -> 24 bytes
      std::vector<double> m_y;
      std::vector<double> m_z;
      std::vector<double> m_xNoisy;
      std::vector<double> m_yNoisy;
      std::vector<double> m_zNoisy;
      std::vector<double> m_sigmaX;
      std::vector<double> m_sigmaY;
      std::vector<double> m_sigmaX_gaus;
      std::vector<double> m_sigmaY_gaus;
      unsigned m_degX, m_degY; // 4 bytes each -> 8
      double m_z0; // 8 bytes
      double m_sigmaXvalue, m_sigmaYvalue;
      void clear()
      {
          m_x.clear();
          m_y.clear();
          m_z.clear();
          m_xNoisy.clear();
          m_yNoisy.clear();
          m_zNoisy.clear();
          m_sigmaX.clear();
          m_sigmaY.clear();
          m_sigmaX_gaus.clear();
          m_sigmaY_gaus.clear();
                
    }

    ////Function to eliminate too close z hits and fitting only tracks touching all layers
    unsigned nLayers(const auto& hits) 
    {
    std::array<float, 6> zvals{};
    unsigned n = 0;
    for (const auto& hit: hits) {
        bool layerFound = false;
        for (unsigned i = 0; n != i; ++i) {
            if (std::abs(zvals[i] - hit->midPoint().z()) < 3.f) {
                layerFound = true;
                break;
            }
        }
        if (layerFound) continue;
        if (n >= 6) throw "too many layers";
        zvals[n] = hit->midPoint().z();
        n++;
    }
    return n;
    }

  };//HitWorkspace

  HitsPerMCPart digOutCheatedTracks(const LHCb::MCHits& hits) const;
  std::pair<PolyFitter<double>, PolyFitter<double>> prepareFittersFromHits(unsigned degX, unsigned degY, const Hits& hits,HitWorkspace& wksp) const;
  using Vector = typename PolyFitter<double>::Vector<double>;
  using ConstVector = typename PolyFitter<double>::Vector<const double>;
  using Matrix = typename PolyFitter<double>::Matrix<double>;
  using ConstMatrix = typename PolyFitter<double>::Matrix<const double>;
  static void printFitResult(double chi2, std::size_t ndf, const ConstVector params, const ConstMatrix& cov);
  
};

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( MCHitFitMonitor )

MCHitFitMonitor::MCHitFitMonitor( const std::string& name, ISvcLocator* pSvcLocator )
    //: Consumer( name, pSvcLocator, KeyValue{"DepositLocation", LHCb::MCHitLocation::MP} ),m_tupleFiller(10)  {}
    : Consumer( name, pSvcLocator, KeyValue{"DepositLocation", LHCb::MCHitLocation::MP} )  {}

//    : Consumer( name, pSvcLocator, KeyValue{"DepositLocation", LHCb::MCHitLocation::Default} ) {}

StatusCode MCHitFitMonitor::initialize() {
    std::cout<<"MCHitFitMonitor::initialize()"<<std::endl;
    return Consumer::initialize()
        .andThen( [&] { m_gauss.initialize( randSvc(), Rndm::Gauss( 0.0, 1.0 )); } );

}

void MCHitFitMonitor::operator()( const LHCb::MCHits& hits ) const {
  
  //creating tuple
  std::cout<<"Creating tuple"<<std::endl;
  Tuple tuple = nTuple("TrackFitParameters");
  TupleFiller m_tupleFiller{25, 16, 6};  // Set max hits/properties PerTrack as needed
  std::cout<<"Tuple created"<<std::endl;
  HitWorkspace wksp;
  wksp.m_sigmaXvalue = 1.0;//1.;//0.1;
  wksp.m_sigmaYvalue = 1.0;//1.;//0.1;
  wksp.m_z0 = 8500;//middle of T station
  wksp.m_degX = 5;
  wksp.m_degY = 5;

  //std::cout<<"Preparing workspace"<<std::endl;
  const HitsPerMCPart tracks = digOutCheatedTracks(hits);
    //std::cout<<"Workspace prepared"<<std::endl;
  int track_counter = 0;
  for (const auto& [mcp, hits]: tracks) {
    track_counter += 1;
    if(!mcp) continue;
    if(mcp->goodEndVertex() && mcp->goodEndVertex()->position().z() < 10000.) continue;//excluding early hadronization
    if ( std::abs( mcp->particleID().pid() ) == 11 ) continue;//excluding e-,e+
    //if ( hits.size() < 2 ) continue;//not fitting 2 points
    if ( wksp.nLayers(hits) < 6 ) continue; //eliminate too close z hits and not fitting hits not touching all layers!
    
    auto [fitterX, fitterY] = prepareFittersFromHits(wksp.m_degX, wksp.m_degY, hits, wksp);

    auto fitX = fitterX.fit(std::min(std::max(ptrdiff_t(hits.size())-2,ptrdiff_t(0)),ptrdiff_t(wksp.m_degX)));
      
    auto fitY = fitterY.fit(std::min(std::max(ptrdiff_t(hits.size())-2,ptrdiff_t(0)),ptrdiff_t(wksp.m_degY)));

    if(!fitX){
          std::cout<<"NO RESULT in x plane!"<<std::endl;
    }
    if (!fitY) {
        std::cout<<"NO RESULT in y plane!"<<std::endl;
    }
    if (!fitX || !fitY) continue;

    std::cout<<"Result in X - Plane "<< std::endl;
    auto resultXplane = *fitX ;

    double chi2X = resultXplane.chi2(cppcompat::make_span(wksp.m_zNoisy),cppcompat::make_span(wksp.m_xNoisy),cppcompat::make_span(wksp.m_sigmaX));
    int ndfX = resultXplane.ndf();
    int n_paramsX = resultXplane.nparams();
    //int n_MatX = n_paramsX*(n_paramsX+1.)*0.5;
    ///calculate chi2Prob 
    double probX = TMath::Prob(chi2X, ndfX);

    printFitResult(chi2X, ndfX, resultXplane.params(),resultXplane.cov());

    //std::vector<float> x_fitted = resultXplane.fittedCoordinate(cppcompat::make_span(wksp.m_zNoisy),cppcompat::make_span(wksp.m_xNoisy));

    std::vector<float> a_nX;
    for(auto& a_n : resultXplane.params()){
        a_nX.push_back(a_n);
    }

    std::vector<float> covarianceX ;//as an array of LU of Matrix
    for(unsigned i = 0; i <= n_paramsX ; ++i){
        for(unsigned j = 0; j <= i; ++j){
            covarianceX.push_back(resultXplane.cov()[i][j]);
        }
    }

    std::vector<float> matrixX ;//as an array of LU of Matrix
    for(unsigned i = 0; i < n_paramsX ; ++i){
        for(unsigned j = 0; j <= i; ++j){
            matrixX.push_back(resultXplane.matrix()[i][j]);
        }
    }

    std::vector<float> rhsX ;
    for (auto& r_n : resultXplane.rhs()){
        rhsX.push_back(r_n);
    }
            
    std::cout<<"Result in Y - Plane "<< std::endl;
    auto& resultYplane = *fitY ;

    double chi2Y = resultYplane.chi2(cppcompat::make_span(wksp.m_zNoisy), cppcompat::make_span(wksp.m_yNoisy), cppcompat::make_span(wksp.m_sigmaY));
    int ndfY = resultYplane.ndf();
    int n_paramsY = resultYplane.nparams();
    //int n_MatY = n_paramsY*(n_paramsY+1)*0.5;
    double probY = TMath::Prob(chi2Y, ndfY);

    //std::vector<float> y_fitted = resultYplane.fittedCoordinate(cppcompat::make_span(wksp.m_zNoisy),cppcompat::make_span(wksp.m_yNoisy));

    std::vector<float> a_nY;
    for(auto& a_n : resultYplane.params()){
        a_nY.push_back(a_n);
    }
    std::vector<float> covarianceY ;//as an array of LU of Matrix
    for(unsigned i = 0; i < n_paramsY ; ++i){
        for(unsigned j = 0; j <= i; ++j){
            covarianceY.push_back(resultYplane.cov()[i][j]);
        }
    }
    std::vector<float> matrixY ;//as an array of LU of Matrix
    for(unsigned i = 0; i <= n_paramsY ; ++i){
        for(unsigned j = 0; j <= i; ++j){
            matrixY.push_back(resultYplane.matrix()[i][j]);
        }
    }

    std::vector<float> rhsY ;
    for (auto& r_n : resultYplane.rhs()){
        rhsY.push_back(r_n);
    }

    printFitResult(chi2Y, ndfY, resultYplane.params(),resultYplane.cov());

    // Clear the tuple data
    m_tupleFiller.clear();
     
    for (unsigned i = 0; wksp.m_x.size() != i; ++i){
        
        auto x = wksp.m_xNoisy[i];
        auto y = wksp.m_yNoisy[i];
        auto z = wksp.m_zNoisy[i];

        auto x_true = wksp.m_x[i];
        auto y_true = wksp.m_y[i];
        auto z_true = wksp.m_z[i];

        ////save x/y/z_tracks 
        //auto x_f = x_fitted[i];
        //auto y_f = y_fitted[i];
        auto x_f = resultXplane.fittedCoordinate(z);
        auto y_f = resultYplane.fittedCoordinate(z);
        /*double x_fitted = 0;
        double y_fitted = 0;
        for (unsigned j = 0; n_paramsX != j; ++j){
            x_fitted += a_nX[j]*pow((wksp.m_zNoisy[i]-wksp.m_z0),j);

        }
        for (unsigned j = 0; n_paramsY != j; ++j){
            y_fitted += a_nY[j]*pow((wksp.m_zNoisy[i]-wksp.m_z0),j);

        }*/
        auto sigmaX_gaus = wksp.m_sigmaX_gaus[i];
        auto sigmaY_gaus = wksp.m_sigmaY_gaus[i];
        auto sigmax = wksp.m_sigmaX[i];
        auto sigmay = wksp.m_sigmaY[i];
        auto x_Resi = x_true - x_f;
        auto y_Resi = y_true - y_f;
        std::cout<<" hit z "<< z<<" x "<<x<<" y "<<y<<std::endl;
        m_tupleFiller.fillPerHit(x, y, z, x_true, y_true, z_true, sigmax, sigmay, x_f, y_f, x_Resi, y_Resi, sigmaX_gaus, sigmaY_gaus);  
    }//loop over Hits in each track
 
    m_tupleFiller.fillPerTrackX(chi2X, ndfX, probX, a_nX, rhsX, matrixX, covarianceX);
    m_tupleFiller.fillPerTrackY(chi2Y, ndfY, probY, a_nY, rhsY, matrixY, covarianceY);
    m_tupleFiller.writeToTuple(tuple);
  }//loop over tracks

}//end of MCHits algorithm 

void MCHitFitMonitor::printFitResult(double chi2, std::size_t ndf,
                              const MCHitFitMonitor::ConstVector params,
                              const MCHitFitMonitor::ConstMatrix& cov)
{
    std::printf("FIT RESULT: chi2=%g, ndf=%zu\n\n", chi2, ndf);

    std::printf("params: [");
    for (const auto el : params) {
        std::printf("%s%g\n", (el == params.front()) ? " " : ", ", el);
    };
    std::printf(" ]\n\n");

    std::printf("   cov: [");
    for (std::size_t i = 0, n = params.size(); n != i; ++i) {
        
        for (std::size_t j = 0; n != j; ++j) {
            std::printf("%s%g", !j ? " " : ", ", cov[i][j]);
        }
        std::printf(" ]\n");
    }
}

MCHitFitMonitor::HitsPerMCPart MCHitFitMonitor::digOutCheatedTracks(const LHCb::MCHits& hits) const
{
    HitsPerMCPart tracks;
    for (const LHCb::MCHit* phit : hits) {
        const LHCb::MCHit& hit = *phit;
        const auto* mc_particle = hit.mcParticle();
        auto& hits = tracks[mc_particle];
        hits.push_back(&hit);
    }
    return tracks;
}

std::pair<PolyFitter<double>, PolyFitter<double>> MCHitFitMonitor::prepareFittersFromHits(unsigned degX, unsigned degY,
                                 				const MCHitFitMonitor::Hits& hits,
                                 				MCHitFitMonitor::HitWorkspace& wksp) const
{
    wksp.clear();

    for (const auto* hit: hits) {
        const auto mp = hit->midPoint();
        wksp.m_x.push_back(mp.x());
        wksp.m_y.push_back(mp.y());
        wksp.m_z.push_back(mp.z());

        double noiseX = m_gauss()*wksp.m_sigmaXvalue;
        double noiseY = m_gauss()*wksp.m_sigmaYvalue;
        
        /*wksp.m_xNoisy.push_back(mp.x()+m_gauss()*wksp.m_sigmaXvalue);
        wksp.m_yNoisy.push_back(mp.y()+m_gauss()*wksp.m_sigmaYvalue);*/
        wksp.m_xNoisy.push_back(mp.x()+noiseX);
        wksp.m_yNoisy.push_back(mp.y()+noiseY);
        wksp.m_zNoisy.push_back(mp.z());

        wksp.m_sigmaX_gaus.push_back(noiseX);
        wksp.m_sigmaY_gaus.push_back(noiseY);
    }

    wksp.m_sigmaX.resize(wksp.m_x.size(), wksp.m_sigmaXvalue);
    wksp.m_sigmaY.resize(wksp.m_y.size(), wksp.m_sigmaYvalue);

    std::pair<PolyFitter<double>, PolyFitter<double>> retVal{{degX, wksp.m_z0}, {degY, wksp.m_z0}};
    
    retVal.first.addMeas(cppcompat::make_span(wksp.m_z),
                         cppcompat::make_span(wksp.m_xNoisy),
                         cppcompat::make_span(wksp.m_sigmaX));
    retVal.second.addMeas(cppcompat::make_span(wksp.m_z),
                          cppcompat::make_span(wksp.m_yNoisy),
                          cppcompat::make_span(wksp.m_sigmaY));


    return retVal;
}


