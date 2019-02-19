g++ -o main -fprofile-arcs -ftest-coverage test_DBG.cpp -L /usr/lib -I/usr/include
./main
gcov test_catch.cpp
lcov --c --directory . --output-file main_coverage.info
genhtml main_coverage.info --output-directory out
