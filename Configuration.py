###############################################################################
# (c) Copyright 2000-2020 CERN for the benefit of the LHCb Collaboration      #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
"""
High level configuration tools for Boole
"""
__author__ = "Marco Cattaneo <Marco.Cattaneo@cern.ch>"

from functools import cached_property
from Gaudi.Configuration import (
    importOptions,
    getConfigurable,
    log,
    INFO,
    WARNING,
    ERROR,
)

import GaudiKernel.ProcessJobOptions
from Configurables import (
    LHCbConfigurableUser,
    LHCbApp,
    ProcessPhase,
    DigiConf,
    SimConf,
    RichDigiSysConf,
    ApplicationMgr,
    GaudiSequencer,
    EventSelector,
    HistogramPersistencySvc,
    OutputStream,
    RecordStream,
    LHCb__UnpackRawEvent,
)

from DDDB.CheckDD4Hep import UseDD4Hep


class Boole(LHCbConfigurableUser):

    ## Default main sequences
    DefaultSequence = [
        "ProcessPhase/Init", "ProcessPhase/Digi", "ProcessPhase/Link",
        "ProcessPhase/Moni", "ProcessPhase/Filter"
    ]

    __slots__ = {
        "DetectorInit": {
            "DATA": ['Data'],
            "MUON": ['Muon']
        },
        "DetectorDigi": ['VP', 'UT', 'FT', 'Rich', 'Calo', 'Muon', 'MP'],
        "DetectorLink": ['VP', 'UT', 'FT', 'Tr', 'Rich', 'Calo', 'Muon', 'MP'],
        "DetectorMoni": ['VP', 'UT', 'FT', 'Rich', 'Calo', 'Muon', 'MC', 'MP'],
        "EvtMax": -1,
        "SkipEvents": 0,
        "UseSpillover": False,
        "SpilloverPaths": [],
        "TAEPrev": 0,
        "TAENext": 0,
        "TAESubdets": ["Calo", "Muon"],
        "Outputs": ["DIGI"],
        "InputDataType": "SIM",
        "DigiType": "Default",
        "Histograms": "Default",
        "NoWarnings": False,
        "ProductionMode": False,
        "OutputLevel": INFO,
        "DatasetName": "Boole",
        "DataType": "Upgrade",
        "DDDBtag": "",
        "CondDBtag": "",
        "Monitors": [],
        "MainSequence": [],
        "InitSequence": [],
        "DigiSequence": [],
        "LinkSequence": [],
        "MoniSequence": [],
        "FilterSequence": [],
        "EnablePack": True,
        "SiG4EnergyDeposit": True,
        "SplitRawEventFormat": None  #Where the raw event sits, DAQ/RawEvent!
        ,
        "WriteFSR": True,
        "RetinaCluster": True,
        "PreserveSP": True,
        "MergeGenFSR": False,
        "DisableTiming": False
    }

    _propertyDocDct = {
        "DetectorInit":
        """ Dictionary specifying the detectors to take into account in initialisation """,
        "DetectorDigi":
        """ Dictionary specifying the detectors to take into account in digitization """,
        "DetectorLink":
        """ Dictionary specifying the detectors to make linkers """,
        "DetectorMoni":
        """ Dictionary specifying the detectors to monitoring """,
        'EvtMax':
        """ Maximum number of events to process """,
        'SkipEvents':
        """ Number of events to skip """,
        'UseSpillover':
        """ Flag to enable spillover (default False) """,
        'SpilloverPaths':
        """ Paths to fill when spillover is enabled """,
        'TAEPrev':
        """ Number of Prev Time Alignment Events to generate """,
        'TAENext':
        """ Number of Next Time Alignment Events to generate """,
        'TAESubdets':
        """ Subdetectors for which TAE are enabled """,
        'Outputs':
        """ List of outputs: ['MDF','DIGI'] (default 'DIGI') """,
        'InputDataType':
        """ Input Data Type 'XDST' (default 'SIM') """,
        'DigiType':
        """ Defines content of DIGI file: ['Minimal','Default','Extended'] """,
        'Histograms':
        """ Type of histograms: ['None','Default','Expert'] """,
        'NoWarnings':
        """ OBSOLETE, kept for Dirac compatibility. Please use ProductionMode """,
        'ProductionMode':
        """ Enables special settings for running in production """,
        'OutputLevel':
        """ The printout level to use (default INFO) """,
        'DatasetName':
        """ String used to build output file names """,
        'DataType':
        """ Data type. Default '2015' (use 'Upgrade' for LHCb Upgrade simulations)""",
        'DDDBtag':
        """ Tag for DDDB """,
        'CondDBtag':
        """ Tag for CondDB """,
        'Monitors':
        """ List of monitors to execute """,
        'MainSequence':
        """ The default main sequence, see self.DefaultSequence """,
        'InitSequence':
        """ List of initialisation sequences, see KnownInitSubdets """,
        'DigiSequence':
        """ List of subdetectors to digitize, see KnownDigiSubdets """,
        'LinkSequence':
        """ List of MC truth link sequences, see KnownLinkSubdets  """,
        'MoniSequence':
        """ List of subdetectors to monitor, see KnownMoniSubdets """,
        'FilterSequence':
        """ List of Filter sequences, see KnownFilterSubdets  """,
        'EnablePack':
        """ Turn on/off packing of the data (where appropriate/available) """,
        'SiG4EnergyDeposit':
        """ Modelling of energy loss for silicon trackers from Geant4 or in-house.""",
        'WriteFSR':
        """ Flag whether to write out an FSR.  Default : True """,
        'RetinaCluster':
        """ Flag whether to produce Retina Cluster.  Default : False """,
        'PreserveSP':
        """ Flag whether to preserve SuperPixels when producing Retina Cluster.  Default : False """,
        'MergeGenFSR':
        """ Flag whether to merge the generatore level FSRs. Default False """,
        "DisableTiming":
        """Do not run TimingAuditor, useful for QMTests"""
    }

    KnownFilterSubdets = ["ODIN"]
    KnownHistOptions = ["", "None", "Default", "Expert"]
    KnownSpillPaths = ["Prev", "PrevPrev", "Next", "NextNext"]

    __used_configurables__ = [LHCbApp, DigiConf, SimConf, RichDigiSysConf]

    __detLinkListDigiConf = []

    def defineDB(self):
        # Delegate handling to LHCbApp configurable
        self.setOtherProps(LHCbApp(), ["CondDBtag", "DDDBtag", "DataType"])
        LHCbApp().Simulation = True

    def setLHCbAppDetectors(self):
        from Configurables import LHCbApp
        # If detectors set in LHCbApp then use those
        if hasattr(LHCbApp(), "Detectors"):
            if not LHCbApp().Detectors:
                LHCbApp().Detectors = self.getProp("DetectorDigi")
            else:
                log.warning(
                    "Value of 'LHCbApp().Detectors' already set, using that value: %s"
                    % (LHCbApp().Detectors))
        return

    def defineEvents(self):
        # Delegate handling to LHCbApp configurable
        self.setOtherProps(LHCbApp(), ["EvtMax", "SkipEvents"])
        # Setup SIM input
        self.setOtherProp(DigiConf(), "EnablePack")
        # MC 20160921 Next line breaks if EnablePack is false and SIM is packed. I don't understand why it was there
        #        SimConf().setProp("EnableUnpack",self.getProp("EnablePack"))

        detListSim = []
        if 'VP' in self.getProp('DetectorDigi'): detListSim += ['VP']
        if 'UT' in self.getProp('DetectorDigi'): detListSim += ['UT']
        if 'FT' in self.getProp('DetectorDigi'): detListSim += ['FT']
        if 'Rich' in self.getProp('DetectorDigi'): detListSim += ['Rich']
        if 'Calo' in self.getProp('DetectorDigi'):
            detListSim += ['Ecal', 'Hcal']
        if 'Muon' in self.getProp('DetectorDigi'): detListSim += ['Muon']
        if 'HC' in self.getProp('DetectorDigi'): detListSim += ['HC']
        if 'MP' in self.getProp('DetectorDigi'): detListSim += ['MP']
        SimConf().setProp("Detectors", detListSim)

    def configurePhases(self):
        """
        Set up the top level sequence and its phases
        """
        booleSeq = GaudiSequencer("BooleSequencer")
        ApplicationMgr().TopAlg = [booleSeq]
        mainSeq = self.getProp("MainSequence")
        if len(mainSeq) == 0:
            mainSeq = self.DefaultSequence
        booleSeq.Members += mainSeq

        if self.getProp("TAENext") > 0 or self.getProp("TAEPrev") > 0:
            tae = True
            self.enableTAE()
        else:
            tae = False

        detListInit = []
        detListDigi = []
        detListLink = []
        detListMoni = []

        if 'Data' in self.getProp('DetectorInit')['DATA']:
            detListInit += ['Data']

        for det in ['VP', 'UT', 'FT', 'HC']:
            if det in self.getProp('DetectorDigi'):
                detListDigi += [det]
                if det in self.getProp('DetectorLink'): detListLink += [det]
                if det in self.getProp('DetectorMoni'): detListMoni += [det]

        if 'Tr' in self.getProp('DetectorLink'): detListLink += ['Tr']

        for det in ['Rich', 'Calo', 'Muon']:
            if det in self.getProp('DetectorDigi'):
                detListDigi += [det]
                if det in self.getProp('DetectorLink'): detListLink += [det]
                if det in self.getProp('DetectorMoni'): detListMoni += [det]

        if ('Muon' in self.getProp('DetectorInit')['MUON']
                and 'Muon' in self.getProp('DetectorDigi')):
            detListInit += ['Muon']

        if 'MC' in self.getProp('DetectorMoni'): detListMoni += ['MC']

        # add MP from list of FT hits
        if 'MP' in self.getProp('DetectorDigi'):
            detListDigi += ['MP']
            detListLink += ['MP']
            detListMoni += ['MP']
            detListInit += ['MP']

        # Add requested linkers to Digi output
        for det in detListLink:
            self.__detLinkListDigiConf.append(det)

        # Now that requested linkers are saved, add any additional linkers needed for monitoring
        for det in ['VP', 'UT', 'FT', 'Calo']:
            if not (det in detListLink) and (det in detListMoni):
                detListLink += [det]

        initDets = self._setupPhase("Init", detListInit)
        digiDets = self._setupPhase("Digi", detListDigi)
        linkDets = self._setupPhase("Link", detListLink)
        moniDets = self._setupPhase("Moni", detListMoni)

        self.configureInit(tae, initDets)
        self.configureDigi(digiDets)
        self.configureLink(linkDets, moniDets)
        self.configureMoni(moniDets)
        self.configureFilter()

    def configureInit(self, tae, initDets):
        """
        Set up the initialization sequence
        """
        # Start the DataOnDemandSvc ahead of ToolSvc
        ApplicationMgr().ExtSvc += ["DataOnDemandSvc"]
        ApplicationMgr().ExtSvc += ["ToolSvc"]

        ProcessPhase("Init").DetectorList.insert(
            0, "Boole")  # Always run Boole initialisation first!
        initBoole = GaudiSequencer("InitBooleSeq")
        initBoole.Members += ["BooleInit"]

        if UseDD4Hep:
            from Configurables import LHCb__Det__LbDD4hep__IOVProducer as IOVProducer
            initBoole.Members.append(IOVProducer('ReserveIOVDD4hep'))

        # Kept for Dirac backward compatibility
        if self.getProp("NoWarnings"):
            log.warning(
                "Boole().NoWarnings=True property is obsolete and maintained for Dirac compatibility. Please use Boole().ProductionMode=True instead"
            )
            self.setProp("ProductionMode", True)

        # Special settings for production
        if self.getProp("ProductionMode"):
            self.setProp("OutputLevel", ERROR)
            if not LHCbApp().isPropertySet("TimeStamp"):
                LHCbApp().setProp("TimeStamp", True)

        # OutputLevel
        self.setOtherProp(LHCbApp(), "OutputLevel")
        if self.isPropertySet("OutputLevel"):
            level = self.getProp("OutputLevel")
            if level == ERROR or level == WARNING:
                # Additional information to be kept
                getConfigurable("BooleInit").OutputLevel = INFO

        # Do not print event number at every event (done already by BooleInit)
        EventSelector().PrintFreq = -1

        # Load the spillover branches, then kill those not required to prevent further access
        spillPaths = self.getProp("SpilloverPaths")
        killPaths = []
        if len(spillPaths) == 0:
            spillPaths = self.KnownSpillPaths
            self.setProp("SpilloverPaths", spillPaths)

        if self.getProp("UseSpillover"):
            if tae:
                killPaths = self.KnownSpillPaths
            else:
                self.setOtherProp(SimConf(), "SpilloverPaths")
                # Kill any spillover paths not required
                for spill in self.KnownSpillPaths:
                    if spill not in spillPaths:
                        killPaths.append(spill)
        else:
            # Kill all spillover paths
            killPaths = self.KnownSpillPaths

        from Configurables import EventNodeKiller, TESCheck
        spillLoader = TESCheck("SpilloverLoader")
        spillLoader.Inputs = spillPaths
        spillLoader.Stop = False  # In case no spillover on input file
        spillLoader.OutputLevel = ERROR
        spillKiller = EventNodeKiller("SpilloverKiller")
        spillKiller.Nodes = killPaths
        spillHandler = GaudiSequencer("SpilloverHandler")
        spillHandler.Members += [spillLoader, spillKiller]
        spillHandler.IgnoreFilterPassed = True  # In case no spillover on input file
        initBoole.Members += [spillHandler]

        if "Muon" in initDets:
            # Muon Background
            import MuonBackground as MBG
            initMuonSeq = GaudiSequencer("InitMuonSeq")
            initMuonSeq.Members.append(MBG.makeLowEnergyBG("G4"))
            if not tae:
                initMuonSeq.Members.append(MBG.makeFlatSpilloverBG("G4"))

        if 'MP' in initDets:
            from Configurables import CopyMCHitsForMT
            # Move MP hits to the MP container
            initMPSeq = GaudiSequencer("InitMPSeq")
            #GaudiSequencer('DigiMPSeq').Members += [
            #    CopyMCHitsForMT('CopyMCHitsForMTPrevPrev'),
            #    CopyMCHitsForMT('CopyMCHitsForMTPrev'),
            #    CopyMCHitsForMT('CopyMCHitsForMT'),
            #    CopyMCHitsForMT('CopyMCHitsForMTNext'),
            #    CopyMCHitsForMT('CopyMCHitsForMTNextNext')
            #]
            #for i in ['PrevPrev', 'Prev', '', 'Next', 'NextNext']:
            #    CopyMCHitsForMT(
            #        'CopyMCHitsForMT%s' %
            #        i).InputLocation = '/Event%s/MC/FT/Hits' % (i if '' == i
            #                                                    else '/' + i)
            #    CopyMCHitsForMT(
            #        'CopyMCHitsForMT%s' %
            #        i).OutputLocation = '/Event%s/MC/MP/Hits' % (i if '' == i
            #                                                     else '/' + i)

    def configureDigi(self, digiDets):
        """
        Set up the digitization sequence
        """
        importOptions("$STDOPTS/PreloadUnits.opts")  # needed by VELO and ST
        if "VP" in digiDets:
            self.configureDigiVP(GaudiSequencer("DigiVPSeq"), "")
        if "UT" in digiDets:
            self.configureDigiUT(GaudiSequencer("DigiUTSeq"), "")
        if "FT" in digiDets:
            self.configureDigiFT(GaudiSequencer("DigiFTSeq"), "")
        if "Rich" in digiDets:
            self.configureDigiRich(GaudiSequencer("DigiRichSeq"), "")
        if "Calo" in digiDets:
            self.configureDigiCalo(GaudiSequencer("DigiCaloSeq"), "")
        if "Muon" in digiDets:
            self.configureDigiMuon(GaudiSequencer("DigiMuonSeq"), "")
        if "HC" in digiDets:
            self.configureDigiHC(GaudiSequencer("DigiHCSeq"), "")
        if "MP" in digiDets:
            self.configureDigiMP(GaudiSequencer("DigiMPSeq"), "")

    def configureDigiVP(self, seq, tae):
        # VP digitisation and clustering
        if tae == "":
            from Configurables import VPEnsureBanks, VPDepositCreator, VPDigitCreator
            seq.Members += [
                VPEnsureBanks(),
                VPDepositCreator(),
                VPDigitCreator()
            ]
            from Configurables import VPSuperPixelBankEncoder
            superPixelEncoder = VPSuperPixelBankEncoder()
            seq.Members += [superPixelEncoder]
            if self.getProp("RetinaCluster"):
                from Configurables import VPRetinaClusterCreator, RawEventSimpleCombiner
                unpacker = LHCb__UnpackRawEvent(
                    RawEventLocation="VeloSP/RawEvent",
                    RawBankLocations=["VeloSP/RawBanks"],
                    BankTypes=["VP"])
                clusterCreator = VPRetinaClusterCreator(
                    RawBanks=unpacker.RawBankLocations[0],
                    RetinaRawBanks="VPRetina/RawBanks")
                superPixelEncoder.RawEventLocation = unpacker.RawEventLocation.path(
                )

                Combiner = RawEventSimpleCombiner(
                    EnableIncrementalMode=True,
                    RawBanksToCopy=["VPRetinaCluster"],
                    InputRawEventLocations=[
                        clusterCreator.RetinaClusterLocation.path()
                    ],
                )
                if self.getProp("PreserveSP"):
                    # want both the SuperPixel banks and Retina clusters in the output
                    Combiner.RawBanksToCopy += ["VP"]
                    Combiner.InputRawEventLocations += [
                        unpacker.RawEventLocation.path()
                    ]
                seq.Members += [unpacker, clusterCreator, Combiner]
        else:
            raise RuntimeError("TAE not implemented for VP")

    def configureDigiUT(self, seq, tae):
        # Upstream Tracker digitisation
        from Configurables import (MCUTDepositCreator, MCUTDigitCreator,
                                   UTDigitCreator, UTTightDigitCreator,
                                   UTDigitsToRawBankAlg)
        if tae == "":
            importOptions("$UTDIGIALGORITHMSROOT/python/utDigi.py")

            toolName = "SiDepositedCharge"
            if self.getProp("SiG4EnergyDeposit"):
                toolName = "SiGeantDepositedCharge"
            depCreator = MCUTDepositCreator()
            depCreator.DepChargeTool = toolName
            seq.Members += [depCreator]
            seq.Members += [MCUTDigitCreator()]
            digitCreator = UTDigitCreator()
            digitCreator.Saturation = 31
            seq.Members += [digitCreator]
            seq.Members += [UTTightDigitCreator()]
            seq.Members += [UTDigitsToRawBankAlg()]
        else:
            raise RuntimeError("TAE not implemented for UT")

    def configureDigiFT(self, seq, tae):
        # Fibre Tracker digitisation
        if tae == "":
            from Configurables import (
                FTMCHitSpillMerger,
                MCFTDepositCreator,
                MCFTDigitCreator,
                FTClusterCreator,
                FTRawBankEncoder,
                FTRawBankDecoder,
            )

            unpacker = LHCb__UnpackRawEvent(
                'UnpackRawEvent',
                BankTypes=['FTCluster'],
                RawEventLocation='/Event/DAQ/RawEvent',
                RawBankLocations=['/Event/DAQ/RawBanks/FTCluster'])

            seq.Members += [
                FTMCHitSpillMerger(),
                MCFTDepositCreator(),
                MCFTDigitCreator(),
                FTClusterCreator(),
                FTRawBankEncoder(),
                unpacker,
                FTRawBankDecoder(),
            ]
        else:
            raise RuntimeError("TAE not implemented for SCIFI")

    def configureDigiMP(self, seq, tae):
        # MightyPix digitisation
        if tae == "":
            from Configurables import (MCMPFakeHitCreator, MCMPDepositCreator,
                                       MCMPDigitCreator, MPDigitCreator)

            seq.Members += [
                MCMPFakeHitCreator('MCMPFakeHitCreator'),
                MCMPDepositCreator('MCMPDepositCreator'),
                MCMPDigitCreator('MCMPDigitCreator'),
                MPDigitCreator('MPDigitCreator')
            ]
        else:
            raise RuntimeError("TAE not implemented for MP")

    def configureDigiRich(self, seq, tae):
        if tae == "":
            from Configurables import RichDigiSysConf
            self.setOtherProp(RichDigiSysConf(), "UseSpillover")
            self.setOtherProp(RichDigiSysConf(), "DataType")
            RichDigiSysConf().Sequencer = GaudiSequencer("DigiRichSeq")

        else:
            raise RuntimeError("TAE not implemented for RICH")

    def configureDigiCalo(self, seq, tae):
        # Calorimeter digitisation
        from Configurables import CaloSignalAlg, CaloDigitAlg, CaloFutureFillRawBuffer

        seq.Members += [
            CaloSignalAlg("EcalSignal%s" % tae, IsTAE=bool(tae)),
            CaloSignalAlg("HcalSignal%s" % tae, IsTAE=bool(tae))
        ]

        seq.Members += [
            CaloDigitAlg("EcalDigit%s" % tae),
            CaloDigitAlg("HcalDigit%s" % tae)
        ]

        for cal in ["Ecal", "Hcal"]:
            det = '/world/DownstreamRegion/' + cal + ':DetElement-Info-IOV' if UseDD4Hep else '/dd/Structure/LHCb/DownstreamRegion/' + cal

            encoder = CaloFutureFillRawBuffer(
                cal + "FillRawBuffer" + tae,
                InputBank="Raw/" + cal + "/FullAdcs",
                DetectorLocation=det,
                DataCodingType=5)  #4 for BigEndian, 5 for LittleEndian
            seq.Members += [encoder]

    def configureDigiMuon(self, seq, tae):
        from Configurables import MuonDigitization, MuonDigitToTell40RawBuffer
        seq.Members += [MuonDigitization("MuonDigitization%s" % tae)]
        seq.Members += [
            MuonDigitToTell40RawBuffer("MuonDigitToTell40RawBuffer%s" % tae)
        ]

    def configureDigiHC(self, seq, tae):
        # HC digitisation
        if tae == "":
            from Configurables import HCDigitCreator
            seq.Members += [HCDigitCreator()]
        else:
            raise RuntimeError("TAE not implemented for HC")

    def configureFilter(self):
        """
        Set up the filter sequence to selectively write out events
        """
        filterDets = self.getProp("FilterSequence")
        for det in filterDets:
            if det not in self.KnownFilterSubdets:
                log.warning("Unknown subdet '%s' in FilterSequence" % det)

        filterSeq = ProcessPhase("Filter", ModeOR=True)
        filterSeq.DetectorList += filterDets

        if "ODIN" in filterDets:
            from Configurables import OdinTypesFilter
            odinFilter = OdinTypesFilter()
            GaudiSequencer("FilterODINSeq").Members += [odinFilter]

    def configureLink(self, linkDets, moniDets):
        """
        Set up the MC links sequence
        """

        doWriteTruth = ("DIGI" in self.getProp("Outputs")) and (
            self.getProp("DigiType").capitalize() != "Minimal")

        if "VP" in linkDets:
            from Configurables import VPDigitLinker
            seq = GaudiSequencer("LinkVPSeq")
            seq.Members += [VPDigitLinker()]

        if "UT" in linkDets:
            from Configurables import (
                UTDigit2MCHitLinker, UTDigit2MCParticleLinker,
                UTTightDigit2MCHitLinker, UTTightDigit2MCParticleLinker)
            seq = GaudiSequencer("LinkUTSeq")
            seq.Members += [UTDigit2MCHitLinker("UTDigitLinker")]
            seq.Members += [UTTightDigit2MCHitLinker("UTTightDigitLinker")]
            seq.Members += [UTTightDigit2MCParticleLinker("UTTruthLinker")]

        if "Tr" in linkDets and doWriteTruth:
            from Configurables import BuildMCTrackInfo
            seq = GaudiSequencer("LinkTrSeq")
            buildMCTrackInfo = BuildMCTrackInfo()
            seq.Members += [buildMCTrackInfo]

            if "VP" in linkDets:
                buildMCTrackInfo.WithVP = True
            else:
                raise RuntimeError("Need VP to build MCTrackInfo")

            if "UT" in linkDets:
                buildMCTrackInfo.WithUT = True
            else:
                buildMCTrackInfo.WithUT = False
                log.warning("No UT to build MCTrackInfo")

            if "FT" in linkDets:
                buildMCTrackInfo.WithFT = True

        if "Rich" in linkDets and doWriteTruth:
            seq = GaudiSequencer("LinkRichSeq")
            seq.Members += ["Rich::MC::MCRichDigitSummaryAlg"]

        if "Calo" in linkDets:
            from Configurables import CaloFutureRawToDigits
            seq = GaudiSequencer("LinkCaloSeq")

            unpacker_calo = LHCb__UnpackRawEvent(
                'UnpackRawEventCalo',
                BankTypes=['Calo'],
                RawEventLocation="DAQ/RawEvent",
                RawBankLocations=['/Event/DAQ/RawBanks/Calo'])

            unpacker_caloerr = LHCb__UnpackRawEvent(
                'UnpackRawEventCaloError',
                BankTypes=['CaloError'],
                RawEventLocation="DAQ/RawEvent",
                RawBankLocations=['/Event/DAQ/RawBanks/CaloError'])

            from Configurables import CaloDigitLinker

            seq.Members += [unpacker_calo, unpacker_caloerr]
            decoders = []
            for cal in ["Ecal", "Hcal"]:
                det = '/world/DownstreamRegion/' + cal + ':DetElement-Info-IOV' if UseDD4Hep else '/dd/Structure/LHCb/DownstreamRegion/' + cal

                decoder = CaloFutureRawToDigits(
                    name="Future" + cal + "ZSup",
                    RawBanks="DAQ/RawBanks/Calo",
                    ErrorRawBanks="DAQ/RawBanks/CaloError",
                    OutputDigitData="Raw/" + cal + "/CaloDigits",
                    OutputReadoutStatusData="/Event/Transient/DAQ/Status" +
                    cal,
                    DetectorLocation=det)

                decoders += [decoder]
                seq.Members += [decoder]

                digitLinker = CaloDigitLinker("Future" + cal + "DigitLinker")
                digitLinker.DigitLocation = decoder.OutputDigitData
                digitLinker.MCDigitLocation = "MC/" + cal + "/Digits"
                digitLinker.DetectorLocation = det
                digitLinker.HitLinkLocation = "/Event/Link/Raw/" + cal + "/Digits2MCHits"
                digitLinker.ParticleLinkLocation = "/Event/Link/Raw/" + cal + "/Digits2MCParticles"

                seq.Members += [digitLinker]

        if "Muon" in linkDets and doWriteTruth:
            seq = GaudiSequencer("LinkMuonSeq")
            seq.Members += ["MuonDigit2MCParticleAlg"]
            seq.Members += ["MuonTileDigitInfo"]

        if "MP" in linkDets:
            print('MP link not implemented yet')

    def enableTAE(self):
        """
        switch to generate Time Alignment events (only Prev1 for now).
        """
        taeSlots = []
        mainSeq = GaudiSequencer("BooleSequencer").Members

        taePrev = self.getProp("TAEPrev")
        while taePrev > 0:
            digi = mainSeq.index("ProcessPhase/Digi")
            taePhase = ProcessPhase("DigiPrev%s" % taePrev)
            taePhase.RootInTES = "Prev%s/" % taePrev
            mainSeq.insert(digi, taePhase)
            taeSlots.append("Prev%s" % taePrev)
            taePrev -= 1

        taeNext = self.getProp("TAENext")
        while taeNext > 0:
            digi = mainSeq.index("ProcessPhase/Digi")
            taePhase = ProcessPhase("DigiNext%s" % taeNext)
            taePhase.RootInTES = "Next%s/" % taeNext
            mainSeq.insert(digi + 1, taePhase)
            taeSlots.append("Next%s" % taeNext)
            taeNext -= 1
        GaudiSequencer("BooleSequencer").Members = mainSeq

        for taeSlot in taeSlots:
            taePhase = ProcessPhase("Digi%s" % taeSlot)
            taeDets = self.getProp("TAESubdets")
            taePhase.DetectorList = ["Init"] + taeDets
            from Configurables import BooleInit
            slotInit = BooleInit("Init%s" % taeSlot, RootInTES="%s/" % taeSlot)
            GaudiSequencer("Digi%sInitSeq" % taeSlot).Members = [slotInit]

            if "VP" in taeDets:
                self.configureDigiVP(
                    GaudiSequencer("Digi%sVPSeq" % taeSlot), taeSlot)
            if "UT" in taeDets:
                self.configureDigiUT(
                    GaudiSequencer("Digi%sUTSeq" % taeSlot), "UT", taeSlot)
            if "Rich" in taeDets:
                self.configureDigiRich(
                    GaudiSequencer("Digi%sRichSeq" % taeSlot), taeSlot)
            if "Calo" in taeDets:
                self.configureDigiCalo(
                    GaudiSequencer("Digi%sCaloSeq" % taeSlot), taeSlot)
            if "Muon" in taeDets:
                self.configureDigiMuon(
                    GaudiSequencer("Digi%sMuonSeq" % taeSlot), taeSlot)
            if "HC" in taeDets:
                self.configureDigiHC(
                    GaudiSequencer("Digi%sHCSeq" % taeSlot), taeSlot)

    def defineMonitors(self):
        # get all defined monitors
        monitors = self.getProp("Monitors") + LHCbApp().getProp("Monitors")
        # Currently no Boole specific monitors, so pass them all to LHCbApp
        LHCbApp().setProp("Monitors", monitors)

        # Use TimingAuditor for timing, suppress printout from SequencerTimerTool
        from Configurables import ApplicationMgr, AuditorSvc, SequencerTimerTool
        ApplicationMgr().ExtSvc += ['AuditorSvc']
        ApplicationMgr().AuditAlgorithms = True
        if not self.getProp("DisableTiming"):
            AuditorSvc().Auditors += ['TimingAuditor']
        if not SequencerTimerTool().isPropertySet("OutputLevel"):
            SequencerTimerTool().OutputLevel = WARNING

    def saveHistos(self):
        # ROOT persistency for histograms
        ApplicationMgr().HistogramPersistency = "ROOT"
        histOpt = self.getProp("Histograms").capitalize()
        if histOpt == "None" or histOpt == "":
            # HistogramPersistency still needed to read in Muon background.
            # so do not set ApplicationMgr().HistogramPersistency = "NONE"
            return

        histosName = self.getProp("DatasetName")
        if histosName == "": histosName = "Boole"
        if (self.evtMax() > 0):
            histosName += '-' + str(self.evtMax()) + 'ev'
        if histOpt == "Expert": histosName += '-expert'
        histosName += '-histos.root'
        # add RootHistogramSink as output stream (new style histograms)
        from Configurables import Gaudi__Histograming__Sink__Root
        RootHistSink = Gaudi__Histograming__Sink__Root("RootHistSink")
        ApplicationMgr().ExtSvc += [RootHistSink]
        RootHistSink.FileName = histosName

        # Do not output histograms with this service (use RootHistSink above)
        # Exclude all histograms from the to convert list
        HistogramPersistencySvc().ConvertHistos = ['DO NOT CONVERT']
        HistogramPersistencySvc().OutputFile = ''

    def defineInput(self):
        """
          Setup input data type can be SIM or XDST
          """
        if self.getProp("InputDataType") == "XDST":

            importOptions("$BOOLEROOT/options/Boole-RunFromXDST.py")
        if self.getProp("InputDataType") == "XDIGI":
            importOptions("$BOOLEROOT/options/Boole-RunFromXDIGI.py")

    def defineOutput(self):
        """
        Set up output stream according to output data type
        """

        knownOptions = ["MDF", "DIGI"]
        outputs = []
        for option in self.getProp("Outputs"):
            if option not in knownOptions:
                raise RuntimeError(
                    "Unknown Boole().Outputs value '%s'" % option)
            outputs.append(option)

        if "DIGI" in outputs:
            seq = GaudiSequencer("PrepareDIGI")
            ApplicationMgr().TopAlg += [seq]

            # In Minimal case, filter the MCVertices before writing
            if self.getProp("DigiType").capitalize() == "Minimal":
                seq.Members = ["FilterMCPrimaryVtx"]

            # In packed case, run the packing algorithms
            if self.getProp("EnablePack"):
                DigiConf().PackSequencer = seq

            writerName = "DigiWriter"
            digiWriter = OutputStream(writerName, Preload=False)

            digiWriter.RequireAlgs.append("Filter")
            if self.getProp("NoWarnings"
                            ) and not digiWriter.isPropertySet("OutputLevel"):
                digiWriter.OutputLevel = INFO

            # Set up the Digi content
            DigiConf().Writer = writerName
            DigiConf().OutputName = self.outputName()
            DigiConf().setProp("Detectors", self.__detLinkListDigiConf)
            self.setOtherProps(DigiConf(), [
                "DigiType", "TAEPrev", "TAENext", "UseSpillover", "DataType",
                "WriteFSR"
            ])
            if self.getProp("UseSpillover"):
                self.setOtherProps(DigiConf(), ["SpilloverPaths"])

        if "MDF" in outputs:
            # Make sure that file will have no knowledge of other nodes
            from Configurables import EventNodeKiller, bankKiller

            algs = []

            nodeKiller = EventNodeKiller("MDFKiller")
            nodeKiller.Nodes = [
                "Rec", "Trig", "MC", "Raw", "Gen", "Link", "pSim"
            ]
            taePrev = self.getProp("TAEPrev")
            while taePrev > 0:
                nodeKiller.Nodes += ["Prev%s" % taePrev]
                taePrev -= 1
            taeNext = self.getProp("TAENext")
            while taeNext > 0:
                nodeKiller.Nodes += ["Next%s" % taeNext]
                taeNext -= 1
            algs.append(nodeKiller)

            if self.getProp("RetinaCluster") and self.getProp("PreserveSP"):
                # We must only have one Velo bank type, remove SP
                algs.append(bankKiller("KillVPBanks", BankTypes=["VP"]))

            MyWriter = OutputStream(
                "RawWriter", Preload=False, ItemList=["/Event/DAQ/RawEvent#1"])
            if not MyWriter.isPropertySet("Output"):
                MyWriter.Output = "DATAFILE='PFN:" + self.outputName(
                ) + ".mdf' SVC='LHCb::RawDataCnvSvc' OPT='REC'"
            MyWriter.RequireAlgs.append("Filter")
            if self.getProp("NoWarnings"
                            ) and not MyWriter.isPropertySet("OutputLevel"):
                MyWriter.OutputLevel = INFO
            algs.append(MyWriter)
            ApplicationMgr().OutStream += algs

        # Merge and write the genFSRs
        if self.getProp("WriteFSR"):
            if self.getProp("MergeGenFSR"):
                seqGenFSR = GaudiSequencer("GenFSRSeq")
                ApplicationMgr().TopAlg += [seqGenFSR]
                seqGenFSR.Members += ["GenFSRMerge"]

            FSRWriter = RecordStream("FSROutputStreamDigiWriter")
            if not FSRWriter.isPropertySet("OutputLevel"):
                FSRWriter.OutputLevel = INFO

    def outputName(self):
        """
        Build a name for the output file, based in input options
        """
        outputName = self.getProp("DatasetName")
        if (self.evtMax() > 0): outputName += '-' + str(self.evtMax()) + 'ev'
        if len(self.getProp("FilterSequence")) > 0: outputName += '-filtered'
        if self.getProp("DigiType") != "Default":
            outputName += '-%s' % self.getProp("DigiType")
        return outputName

    def evtMax(self):
        return LHCbApp().evtMax()

    def configureMoni(self, moniDets):
        # Set up monitoring
        histOpt = self.getProp("Histograms").capitalize()
        if histOpt not in self.KnownHistOptions:
            raise RuntimeError("Unknown Histograms option '%s'" % histOpt)

        from Configurables import BooleInit, MemoryTool
        booleInit = BooleInit()
        booleInit.addTool(MemoryTool(), name="BooleMemory")
        booleInit.BooleMemory.HistoTopDir = "Boole/"
        booleInit.BooleMemory.HistoDir = "MemoryTool"

        if "UT" in moniDets:
            from Configurables import (MCUTDepositMonitor, MCUTDigitMonitor,
                                       UTDigitMonitor, UT__UTTightDigitMonitor,
                                       UTEffChecker, MCParticle2MCHitAlg,
                                       MCParticleSelector)
            from GaudiKernel.SystemOfUnits import GeV
            mcDepMoni = MCUTDepositMonitor("MCUTDepositMonitor")
            mcDigitMoni = MCUTDigitMonitor("MCUTDigitMonitor")
            digitMoni = UTDigitMonitor("UTDigitMonitor")
            #tightdigitMoni = UT__UTTightDigitMonitor("UTTightDigitMonitor")
            mcp2MCHit = MCParticle2MCHitAlg(
                "MCP2UTMCHitAlg",
                MCHitPath="MC/UT/Hits",
                OutputData="/Event/Link/MC/Particles2MCUTHits")
            effCheck = UTEffChecker("UTEffChecker")
            effCheck.addTool(MCParticleSelector)
            effCheck.MCParticleSelector.zOrigin = 50.0
            effCheck.MCParticleSelector.pMin = 1.0 * GeV
            effCheck.MCParticleSelector.betaGammaMin = 1.0
            GaudiSequencer("MoniUTSeq").Members += [
                mcDepMoni,
                mcDigitMoni,
                digitMoni,
                mcp2MCHit,
                #mcDepMoni, mcDigitMoni, digitMoni, mcp2MCHit, tightdigitMoni,
                effCheck
            ]
            if False:
                from Configurables import UTSpilloverSubtrMonitor
                GaudiSequencer("MoniUTSeq").Members += [
                    UTSpilloverSubtrMonitor("UTSpilloverSubtrMonitor")
                ]
            if histOpt == "Expert":
                mcDepMoni.FullDetail = True
                mcDigitMoni.FullDetail = True
                clusMoni.FullDetail = True
                effCheck.FullDetail = True

        if "FT" in moniDets:
            from Configurables import MCFTDepositMonitor, FTLiteClusterMonitor
            GaudiSequencer("MoniFTSeq").Members += [
                MCFTDepositMonitor(),
                FTLiteClusterMonitor()
            ]

        if "MP" in moniDets:
            from Configurables import MCMPFakeHitMonitor, MCMPDepositMonitor, MCMPDigitMonitor, MPDigitMonitor, MCHitFitMonitor
            GaudiSequencer("MoniMPSeq").Members += [
                MCMPFakeHitMonitor('MCMPFakeHitMonitor'),
                MCMPDepositMonitor('MCMPDepositMonitor'),
                MCMPDigitMonitor('MCMPDigitMonitor'),
                MPDigitMonitor('MPDigitMonitor'),
                MCHitFitMonitor('MCHitFitMonitor')
            ]

        if "Rich" in moniDets and not UseDD4Hep:
            from Configurables import Rich__MC__Digi__DigitQC
            GaudiSequencer("MoniRichSeq").Members += [
                Rich__MC__Digi__DigitQC("RiDigitQC")
            ]

        if "Calo" in moniDets:
            from Configurables import CaloDigitChecker
            importOptions("$CALOMONIDIGIOPTS/CaloDigitChecker.opts")
            GaudiSequencer("MoniCaloSeq").Members += [
                CaloDigitChecker("EcalCheck")
            ]
            GaudiSequencer("MoniCaloSeq").Members += [
                CaloDigitChecker("HcalCheck")
            ]

        if "Muon" in moniDets:
            from Configurables import MuonDigitChecker
            GaudiSequencer("MoniMuonSeq").Members += [MuonDigitChecker()]

        if "VP" in moniDets:
            from Configurables import VPDepositMonitor, VPRetinaFullClusterDecoder, VPClusterEffSimDQ
            GaudiSequencer("MoniVPSeq").Members += [
                VPDepositMonitor(),
                VPRetinaFullClusterDecoder(RawBanks="VPRetina/RawBanks"),
                VPClusterEffSimDQ()
            ]

    # Problem comes in here !!
    def _setupPhase(self, name, knownDets):
        seq = self.getProp("%sSequence" % name)
        if len(seq) == 0:
            seq = knownDets
            self.setProp("%sSequence" % name, seq)
        else:
            for det in seq:
                if det not in knownDets:
                    log.warning(
                        "Unknown subdet '%s' in %sSequence" % (det, seq))
        ProcessPhase(name).DetectorList += seq
        return seq

    def __apply_configuration__(self):
        self.defineInput()
        GaudiKernel.ProcessJobOptions.PrintOff()
        self.defineDB()
        self.setLHCbAppDetectors()
        self.defineEvents()
        self.configurePhases()
        self.defineOutput()
        self.defineMonitors()
        self.saveHistos()
        GaudiKernel.ProcessJobOptions.PrintOn()
        log.info(self)
        GaudiKernel.ProcessJobOptions.PrintOff()
