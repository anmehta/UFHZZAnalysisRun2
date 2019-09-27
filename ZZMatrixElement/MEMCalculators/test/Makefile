ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) -lTreePlayer -lTree -lCintex

ROOTLIBS  += -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH)/ -lZZMatrixElementMELA -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH)/ -lZZMatrixElementMEKD -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH)/ -lZZMatrixElementMEMCalculators
ROOTGLIBS     = $(shell root-config --glibs)

#INCDIR       = -IGoodRunsLists -IPileupReweighting -IMuonEfficiencyCorrections -IMuonMomentumCorrections -IegammaEvent
INCDIR       = -I$(CMSSW_BASE)/src/ZZMatrixElement/MELA/interface 

CXX           = g++
CXXFLAGS      =-I$(ROOTSYS)/include -O -Wall -fPIC
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(GDBFLAGS) $<

# clean
clean:
	rm -f *~ *.o *.o~ core



testMEMCalc: testMEMCalc.C

	g++ -g -o testMEMCalc $(CXXFLAGS) $(LIBS) $(INCDIR) testMEMCalc.C 

