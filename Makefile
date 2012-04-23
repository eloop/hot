CXXFLAGS = $(shell hcustom -c) -I . -I ./3rdparty/include -I ./3rdparty/osx/include 
LDFLAGS := $(shell hcustom -m) -L./3rdparty/osx/lib/ -lblitz -lfftw3f

VEX_Ocean: VEX_Ocean.o
	echo "CC is " $(CXX)
	echo "LD is " $(LD)
	$(CXX) $(LDFLAGS) VEX_Ocean.o -bundle -o VEX_Ocean.dylib