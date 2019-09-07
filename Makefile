all: test.out

debug: DEBUG=-g
debug: test.out

test.out: test.cpp btor2tools/build/libbtor2parser.a frontends/btor2_encoder.cpp core/fts.cpp
	$(CXX) $(DEBUG) -std=c++17 -I./utils/ -I./core -Ibtor2tools/src/ -I. -o test.out test.cpp btor2tools/build/libbtor2parser.a core/fts.cpp frontends/btor2_encoder.cpp -lsmt-switch -lsmt-switch-btor -lgmpxx -lgmp

clean:
	rm -rf *.out
