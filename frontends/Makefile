smvparser.cpp smvparser.h: smvparser.y
	bison --defines=smvparser.h -o smvparser.cpp smvparser.y

lex.yy.c: smvparser.l smvparser.h
	flex smvparser.l

smvparser: lex.yy.c smvparser.cpp smvparser.h smv_encoder.cpp smv_cosa2.cpp
	g++ -std=c++17 -I../utils/ -I../contrib/ -I../core/ -I../ -I../engines/ -I../options/ -I../deps/smt-switch/local/include -L../build/ -L../deps/smt-switch/local/lib smvparser.cpp lex.yy.c smv_encoder.cpp smv_cosa2.cpp -o smvparser -lcosa2 -lsmt-switch -lsmt-switch-btor -ly -lgmp

clean:
	rm -f smvparser smvparser.cpp smvparser.h lex.yy.c
