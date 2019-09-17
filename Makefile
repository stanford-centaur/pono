all: cosa2

debug: DEBUG=-g
debug: cosa2

SMTSWITCHDIR="./smt-switch/local"

cosa2: cosa2.cpp btor2tools/build/libbtor2parser.a frontends/btor2_encoder.cpp core/rts.cpp
	$(CXX) $(DEBUG) -std=c++17 -I./utils/ -I./core -Ibtor2tools/src/ -I. -I$(SMTSWITCHDIR)/include -L$(SMTSWITCHDIR)/lib -Wl,-rpath=$(SMTSWITCHDIR)/lib -o cosa2 cosa2.cpp bmc.cpp kinduction.cpp unroller.cpp btor2tools/build/libbtor2parser.a core/rts.cpp core/fts.cpp core/prop.cpp frontends/btor2_encoder.cpp -lsmt-switch -lsmt-switch-btor -lgmpxx -lgmp

clean:
	rm -f cosa2
