OBJ = Interval.o Variables.o settings.o Matrix.o Geometry.o Constraints.o Continuous.o lex.yy.o modelParser.tab.o

.ONESHELL:
all: graph/graph.cpp
	@cd flowstar
	@$(MAKE)
	@$(foreach file, $(OBJ), cp $(file) ../graph;)
	@cd ../graph
	@g++-8 -O3 -w -g -L /usr/local/lib -o ../saw graph.cpp Interval.o Variables.o settings.o Matrix.o Geometry.o Constraints.o Continuous.o lex.yy.o modelParser.tab.o -lmpfr -lgmp -lgsl -lgslcblas -lm -lglpk -lboost_filesystem -lboost_system -lboost_iostreams -I../flowstar
	@cd ..

clean:
	@rm -f saw
	@cd flowstar
	@$(MAKE) clean
	@cd ..
