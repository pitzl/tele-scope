
# user:
# set path to eudaq and GBL:

EUDAQ=/opt/eudaq
# export GBL=/home/YOURID/GBL/V01-17-00/cpp

ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)

ROOTLIBS = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS = $(shell $(ROOTSYS)/bin/root-config --glibs)

# -g for gdb
# -pg for gprof
# -std=c++11

CXXFLAGS = -O2 -Wall -Wextra $(ROOTCFLAGS) -I$(EUDAQ)/main/include

tele: tele.cc
	g++ $(CXXFLAGS) -o tele tele.cc \
	$(ROOTLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: tele'

scopes: scopes.cc
	g++ $(CXXFLAGS) -o scopes scopes.cc \
	$(ROOTLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: scopes'

scopem: scopem.cc
	g++ $(CXXFLAGS) -o scopem scopem.cc \
	$(ROOTLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: scopem'

quad: quad.cc
	g++ $(CXXFLAGS) -I$(GBL)/include -o quad quad.cc \
	-L$(GBL)/lib -lGBL $(ROOTLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: quad'

quad3D: quad3D.cc
	g++ $(CXXFLAGS) -I$(GBL)/include -o quad3D quad3D.cc \
	-L$(GBL)/lib -lGBL $(ROOTLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: quad3D'

evds: evds.cc
	g++ $(CXXFLAGS) -o evds evds.cc \
	$(ROOTGLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: evds'

evd: evd.cc
	g++ $(CXXFLAGS) -o evd evd.cc \
	$(ROOTGLIBS) -L$(EUDAQ)/lib -lEUDAQ
	@echo 'done: evd'
