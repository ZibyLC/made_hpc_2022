GFLAGS = -g
OPTFLAGS = -O3

all: lib

lib_prepare: 
	mkdir build

lib_compile_debug: lib_prepare src/matrix.h src/matrix.cpp
	g++ $(GFLAGS) -c src/matrix.cpp -o build/matrix.o

lib_compile_optimitzation: lib_prepare src/matrix.h src/matrix.cpp
	g++ $(OPTFLAGS) -c src/matrix.cpp -o build/matrix.o

lib_create: build/matrix.o
	ar crf build/libmatrix.a build/matrix.o
	rm -f build/matrix.o

test_debug: build/libmatrix.a test/test.cpp
	g++ $(GFLAGS) -o build/Test -Isrc test/test.cpp build/libmatrix.a

test_optimitzation: build/libmatrix.a test/test.cpp
	g++ $(OPTFLAGS) -o build/Test -Isrc test/test.cpp build/libmatrix.a

debug: lib_prepare lib_compile_debug lib_create test_debug

opt: lib_prepare lib_compile_optimitzation lib_create test_optimitzation

clean:
	rm -rf build
