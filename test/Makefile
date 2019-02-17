CC = g++
CPPFLAGS = --std=c++11 -Wall -I../src/
TESTFLAGS = --std=c++11 -fprofile-arcs -ftest-coverage
TESTSUFF  = -L /usr/lib -I/usr/include -I../src

all: gcov run

run: main
	./main
	gcov test_DBG.cpp
	lcov --c --directory . --output-file main_coverage.info
	lcov --remove main_coverage.info "*/test/*" "/usr/include/*" "/usr/local/*" -o filtered_coverage.info
	rm main_coverage.info
	genhtml filtered_coverage.info --output-directory out

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
	rm -rf out
	rm -f main_coverage.info
	rm -f *.gcda
	rm -f test_dump.tsv

cleanf:
	rm *coverage.info
