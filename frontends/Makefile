smvparser.tab.cc smvparser.tab.hh: smvparser.y
	bison -d -v smvparser.y

lex.yy.c: smvparser.l
	flex smvparser.l

smvparser: smvparser.tab.cc smvparser.tab.hh lex.yy.c smv_encoder.cpp smv_cosa2.cpp
	g++ -std=c++17 -I../utils/ -I../contrib/ -I../core/ -I../ -I../engines/ -I../options/ -I../deps/smt-switch/local/include -L../build/ -L../deps/smt-switch/local/lib smvparser.tab.cc lex.yy.c smv_encoder.cpp smv_cosa2.cpp -o smvparser -lcosa2 -lsmt-switch -lsmt-switch-btor -lsmt-switch-msat -ly -lgmp

clean:
	rm smvparser smvparser.tab.cc smvparser.tab.hh lex.yy.c