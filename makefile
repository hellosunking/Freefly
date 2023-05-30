all: krmdup conventional/ksam_rmdup.pipe

krmdup: src/krmdup.multi-thread.cpp
	@g++ src/krmdup.multi-thread.cpp -std=c++11 -fopenmp -O3 -o krmdup.multi-thread

conventional/ksam_rmdup.pipe: src/ksam_rmdup.pipe.cpp
	@g++ src/ksam_rmdup.pipe.cpp -std=c++11 -O3 -o conventional/ksam_rmdup.pipe

clean:
	rm -f krmdup conventional/ksam_rmdup.pipe
