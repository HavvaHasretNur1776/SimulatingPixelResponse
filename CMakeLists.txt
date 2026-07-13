###############################################################################
# (c) Copyright 2000-2021 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
#[=======================================================================[.rst:
MP/MPDigitisation
-----------------
#]=======================================================================]

gaudi_add_module(MPDigitisation
    SOURCES
	src/MPDigitCreator.cpp
        src/MPDigitMonitor.cpp
        src/MCMPDigitMonitor.cpp
        src/MCMPDepositMonitor.cpp
        src/MCMPDigitCreator.cpp
	src/MCMPDepositCreator.cpp
	src/MCMPFakeHitCreator.cpp
	src/CopyMCHitsForMT.cpp
	src/MCMPFakeHitMonitor.cpp
        src/MCHitFitMonitor.cpp
    LINK
        Gaudi::GaudiAlgLib
        Gaudi::GaudiKernel
	Detector::DetectorLib
        LHCb::LHCbKernel
        LHCb::DigiEvent
        LHCb::MCEvent
        LHCb::GenEvent
)
