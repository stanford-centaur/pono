all: test.out

debug: DEBUG=-g
debug: test.out

SMTSWITCHDIR="./smt-switch/local"

test.out: test.cpp btor2tools/build/libbtor2parser.a frontends/btor2_encoder.cpp core/rts.cpp
	$(CXX) $(DEBUG) -std=c++17 -I./utils/ -I./core -Ibtor2tools/src/ -I. -I$(SMTSWITCHDIR)/include -L$(SMTSWITCHDIR)/lib -Wl,-rpath=$(SMTSWITCHDIR)/lib -o test.out test.cpp btor2tools/build/libbtor2parser.a core/rts.cpp core/fts.cpp frontends/btor2_encoder.cpp -lsmt-switch -lsmt-switch-btor -lgmpxx -lgmp

clean:
	rm -rf *.out
