smvparser.tab.c smvparser.tab.h: smvparser.y
	bison -d smvparser.y

lex.yy.c: smvparser.l smvparser.tab.h
	flex smvparser.l

smvparser: lex.yy.c smvparser.tab.c smvparser.tab.h encoder.cpp smv_cosa2.cpp
	g++ -std=c++17 -I../utils/ -I../contrib/ -I../core/ -I../ -I../engines/ -I../options/ -I../deps/smt-switch/local/include -L../build/ -L../deps/smt-switch/local/lib smvparser.tab.c lex.yy.c encoder.cpp smv_cosa2.cpp -o smvparser -lcosa2 -lsmt-switch -lsmt-switch-btor -ly -lgmp

clean:
	rm smvparser smvparser.tab.c smvparser.tab.h lex.yy.c