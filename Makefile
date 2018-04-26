
ROOTCFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)

ROOTLIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS   = $(shell $(ROOTSYS)/bin/root-config --glibs)

# -g for gdb
# -pg for gprof
# -std=c++11

CXXFLAGS = -O2 -Wall -Wextra $(ROOTCFLAGS) -I/home/pitzl/eudaq/main/include -I/home/pitzl/GBL/V01-17-00/cpp/include

# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/pitzl/GBL/V01-17-00/cpp/lib/

tele: tele.cc
	g++ tele.cc $(CXXFLAGS) -fopenmp -o tele \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: tele'

scopes: scopes.cc
	g++ $(CXXFLAGS) scopes.cc -o scopes \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: scopes'

scoped: scoped.cc
	g++ $(CXXFLAGS) scoped.cc -o scoped \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: scoped'

scopesh: scopesh.cc
	g++ $(CXXFLAGS) scopesh.cc -o scopesh \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: scopesh'

scopes1: scopes1.cc
	g++ $(CXXFLAGS) scopes1.cc -o scopes1 \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: scopes1'

scoper: scoper.cc
	g++ $(CXXFLAGS) scoper.cc -o scoper \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: scoper'

scopem: scopem.cc
	g++ $(CXXFLAGS) scopem.cc -o scopem \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: scopem'

shallow: shallow.cc
	g++ $(CXXFLAGS) shallow.cc -o shallow \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: shallow'

scopeg: scopeg.cc
	g++ $(CXXFLAGS) scopeg.cc -o scopeg \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ -L/home/pitzl/GBL/V01-17-00/cpp/lib -lGBL
	@echo 'done: scopeg'

telefit: telefit.cc
	g++ $(CXXFLAGS) telefit.cc -o telefit \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ -L/home/pitzl/GBL/V01-17-00/cpp/lib -lGBL
	@echo 'done: telefit'

evd: evd.cc
	g++ $(CXXFLAGS) evd.cc -o evd \
	$(ROOTGLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: evd'

addal6: addal6.cc
	g++ $(ROOTCFLAGS) addal6.cc \
	-Wall -O2 -o addal6
	@echo 'done: addal6'

scopemod: scopemod.cc
	g++ $(CXXFLAGS) scopemod.cc -o scopemod \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: scopemod'

scoper50: scoper50.cc
	g++ $(CXXFLAGS) scoper50.cc -o scoper50 \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: scoper50'

scoper99: scoper99.cc
	g++ $(CXXFLAGS) scoper99.cc -o scoper99 \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: scoper99'

scope: scope.cc
	g++ $(CXXFLAGS) scope.cc -o scope \
	$(ROOTLIBS) -L/home/pitzl/eudaq/lib -lEUDAQ
	@echo 'done: scope'
