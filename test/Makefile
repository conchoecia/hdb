CC = g++
CPPFLAGS = --std=c++11 -Wall -I../
TESTFLAGS = --std=c++11 -fprofile-arcs -ftest-coverage
TESTSUFF  = -L /usr/lib -I/usr/include -I../

all: gcov run

run: a.out
	./main
	gcov test_DBG.cpp
	lcov --c --directory . --output-file main_coverage.info
	genhtml main_coverage.info --output-directory out

test: test_DBG.cpp
	$(CC) $(CPPFLAGS) test_DBG.cpp

gcov: test_DBG.cpp
	$(CC) -o main $(TESTFLAGS) test_DBG.cpp $(TESTSUFF)

clean:
	rm -f a.out
	rm -f main
	rm -f *.gcov
	rm -f *.gcno
	rm -f test_case_HKCs
