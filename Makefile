CXXFLAGS = `hcustom -c`
LDFLAGS = `hcustom -m`

SOP_Cleave: SOP_Cleave.o
	echo "CC is " $(CXX)
	echo "LD is " $(LD)
	$(CXX) $(LDFLAGS) SOP_Cleave.o -o SOP_Cleave.dylib