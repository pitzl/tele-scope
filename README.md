# tele-scope
test beam pixel telescope analysis based on eudaq only

## prerequisites:

* you need to have sudo priviliges for your system to install some software

* you need ROOT (5 or 6) 
  ```
  export ROOTSYS=/home/YOU/root (or wherever your root-config sits)
  ```
* you need to have git installed:  
  ```
  sudo apt-get install git
  ```
* you need to have cmake installed:  
  ```
  sudo apt-get install cmake
  ```
* (for GeneralBrokenLines you need to have svn installed)

* for proc600 data you need the USB driver from FTDI:  
  ```
  download from http://www.ftdichip.com/Drivers/D2XX.htm  
  (for Linux: x64 (64 bit))  
  gunzip libftd2xx-x86_64-1.3.6.tgz  
  tar -xf libftd2xx-x86_64-1.3.6.tar  
  (you get a directory called release)  
  cd release  
  sudo cp -p WinTypes.h /usr/local/include  
  sudo cp -p ftd2xx.h /usr/local/include  
  (release/build contains libftd2xx.so.1.3.6)  
  cd /usr/local/lib  
  sudo ln -s /home/YOU/release/build/libftd2xx.so.1.3.6  libftd2xx.so  
  ```
* for proc600 data you need to install the libusb-1.0-0-dev package for your system
  ```
  https://sourceforge.net/projects/libusb/files/libusb-1.0/libusb-1.0.20/
  instruction for installation
  ./configure
  make
  make instal
  ```

* for proc600 data you need pXar:  
  ```
  see https://twiki.cern.ch/twiki/bin/viewauth/CMS/Pxar  
  or like this:  
  git clone https://github.com/psi46/pxar.git  
  cd pxar  
  mkdir build  
  cd build  
  cmake -DBUILD_pxarui=OFF ..  
  #or cmake -DBUILD_pxarui=OFF -DINTERFACE_USB=OFF .. #(without usb, without ftdi)  
  make VERBOSE=1 -j4 install  
  ```
* install eudaq as described in  
  https://github.com/eudaq/eudaq/README.md  
  or like this:  
  ```
  git clone https://github.com/eudaq/eudaq.git  
  cd eudaq
  git tag  (list of tags)
  git checkout v1.6.0
  cd eudaq  
  mkdir build  
  cd build  
  export PXARPATH=/home/YOU/pxar  
  cmake -DBUILD_cmspixel=ON ..  
  make install  
  ```

  Problems might occur on OS X when compiling eudaq in pxar.
  Possible solution is to modify: pxar/core/utils/helper.h
  ```
  -#include "api.h"
  +#include "../api/api.h"	
  ```

## tele-scope
* ceckout the tele-scope package:
  ```
  git clone https://github.com/pitzl/tele-scope.git
  cd tele-scope	
  ```

* adjust the makefile according to your setup:
  ```
  change /home/pitzl/eudaq everywhere to your location of eudaq
  ```
  export LD_LIBRARY_PATH=/home/pitzl/ROOT/root/build/lib:/home/pitzl/eudaq-1p6/lib:
  ```

* step 0:  
  prepare a geo.dat file with the telescope and DUT/REF/MOD planes  
  (see one of the examples)  

  get data:
  ```
  mkdir data
  scp -p YOU@desy-cms010:/data/group/pixel/testbeam/native/run025447.raw data/.
  ```

* step 1: telescope triplet alignment:
  ```
  make tele  
  tele -l 99999 -g geo_2016_04f.dat -p 5.6 25447
  (reads data/run025447.raw  
  (writes align_25447.dat and hot_25447.dat)  
  iterate at least three times (simply re-run)  
  creates tele_25447.root  
  ```
* step 2: telescope with DUT and MOD:  
  update runs.dat with run number, geo, GeV
  ```
  make scopem  
  scope 20833  
  (reads runs.dat, which must a link to geo.dat)  
  (reads align_20833.dat and hot_20833.dat)  
  (write alignDUT_20833.dat)  
  iterate 3 times  
  creates scope_20833.root  
  ```

* for quad module data you need GBL:
  ```
  svn co https://svnsrv.desy.de/desy/GeneralBrokenLines/
  cd GeneralBrokenLines
  cd trunk
  cd cpp
  mkdir build
  cd build
  cmake ..
  emacs CMakeCache.txt
  CMAKE_CXX_FLAGS:STRING=-std=c++11
  make
  make install
  export GBL=/home/YOURID/GeneralBrokenLines/trunk/cpp
  ```

## Trouble shooting

* Problem when running ./tele
  ```
  dyld: Library not loaded: @rpath/libEUDAQ.dylib
  Referenced from: /home/YOUR/tele-scope/./tele
  ```

  Solution:
  ```
  sudo ln -s  /home/YOUR/eudaq/lib/libEUDAQ.dylib /usr/local/lib/.
  ```

## Setup the follwoing before running on the NAF 
  ```
   module load gcc/47
   export ROOTSYS=/nfs/dust/cms/user/schuep/software/ilcsoft/v01-17-05/root/5.34.18
   export EUDAQ=/nfs/dust/cms/user/schuep/software/ilcsoft/v01-17-05/Eutelescope/trunk/external/eudaq/v1.5.1
   cd $ROOTSYS
   source bin/thisroot.sh
   cd -
   export LD_LIBRARY_PATH=/nfs/dust/cms/user/schuep/software/ilcsoft/v01-17-05/Eutelescope/trunk/external/eudaq/v1.5.1/lib:$LD_LIBRARY_PATH
  ```

## Running jobs on the NAF batch system

* To submit batch jobs for tele or scope execute (one  mode is required for running)

  ```
  ./scripts/batchJobSubmit.sh

  Other options and features

  ./scripts/batchJobSubmit.sh -f runs_DATE.dat 
  ./scripts/batchJobSubmit.sh -f runs_DATE.dat -r 2450-2460,2465,2470-2475
  ```
