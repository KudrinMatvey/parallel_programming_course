// Copyright 2019 Kudrin Matvey

#include <tbb/tbb.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <utility>
#include <string>
#include <iostream>
#include <omp.h>
#include <vector>
#include <algorithm>
#include "tbb/concurrent_vector.h" 

using std::vector;
using tbb::blocked_range;
using tbb::concurrent_vector;


void Foo(concurrent_vector<int> *a, size_t z,  int p) {
	int k = static_cast<int>(ceil(pow(2, z)));
	for (size_t j = k % p; j + k < a->size(); j += (k + k))
		for (size_t i = 0; i < a->size() - j - k; i++)
			if ((j + i) / (p + p) == (j + i + k) / (p + p))
				if (a->at(j + i) > a->at(j + i + k)) {
					//std::cout << a->at(j + i) << " " << a->at(j + i + k) << std::endl;
					std::iter_swap(a->begin() + j + i,
						a->begin() + j + i + k);
					//std::cout << a->at(j + i) << " " << a->at(j + i + k) << std::endl;

				}
}

void SerialOddEvenMerge(concurrent_vector<int> *a) {
	size_t r = static_cast<int>(ceil(log2(a->size())));
	size_t p = 1;
	for (size_t p = 1; p < a->size(); p += p) {
		for (size_t i = 0; i < r; i ++) {
			Foo(a, r - i, p);
		}
	}
}

void ParallelOddEvenMerge(concurrent_vector<int> *a) {
	size_t r = static_cast<int>(ceil(log2(a->size())));
	size_t p = 1;
	for (size_t p = 1; p < a->size(); p += p) {
		tbb::parallel_for(size_t(0), r, [&](size_t i) {
			Foo(a, r - i, p);
		});
	}
}

char* getCmdOption(char **begin, char **end, const std::string& option) {
	char **itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
		return *itr;
	return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
	return std::find(begin, end, option) != end;
}

concurrent_vector<int> generateRandomArray(const int en, int min, int max) {
	if (min >= max) {
		max = min + max;
	}
	concurrent_vector<int> arr(en);
	srand(static_cast<unsigned int>(time(NULL)));
	for (int j = 0; j < en; j++)
		arr.at(j) = min + (std::rand() % (max - min));
	return arr;
}

bool check(concurrent_vector<int> *arr, int elementsNumber) {
	bool flag = true;
	int min = arr->at(0);
	for (int i = 1; i < elementsNumber; i++) {
		std::cout << arr->at(i) << " ";
		if (arr->at(i) < min) {
			flag = false;
		}
	}
	return flag;
}

size_t calculateStep(int iter) {
    return static_cast<int>(ceil(pow(2, iter)));
}


void batcherLinear(concurrent_vector<int> *a, const int step) {
	concurrent_vector<int> tmp;
	int start;
	for (start = 0; start < step; start++) {
		for (unsigned int i = start, j = 0; i < a->size(); i += step, j++) {
			tmp.push_back(a->at(i));
		}
		SerialOddEvenMerge(&tmp);
		unsigned int i = start, j = 0;
		for (; i < a->size() - start; i += step, j++) {
			a->at(i) = tmp.at(j);
		}
		tmp.clear();
	} 
}

void batcherTbb(concurrent_vector<int> *a, const int step) {
	concurrent_vector<int> tmp;
	int start;
	for (start = 0; start < step; start++) {
		for (unsigned int i = start, j = 0; i < a->size(); i += step, j++) {
			tmp.push_back(a->at(i));
		}
		ParallelOddEvenMerge(&tmp);
		unsigned int i = start, j = 0;
		for (; i < a->size() - start; i += step, j++) {
			a->at(i) = tmp.at(j);
		}
		tmp.clear();
	}  
}


void mergeTwo(concurrent_vector<int> *a) {
	concurrent_vector<int> tmp(a->size());
	size_t e = 0, o = 1, i = 0;
	while (e < a->size() && o < a->size()) {
		if (a->at(e) < a->at(o)) {
			tmp.at(i++) = (a->at(e));
			e = e + 2;
		}
		else {
			tmp.at(i++) = (a->at(o));
			o = o + 2;
		}
	}
	while (e < a->size()) {
		tmp.at(i++) = (a->at(e));
		e = e + 2;
	}
	while (o < a->size()) {
		tmp.at(i++) = (a->at(o));
		o = o + 2;

	}
	for (size_t k = 0; k < a->size(); k++)
	{
		a->at(k) = tmp.at(k);
	}
}

void mergeTwoTbb(concurrent_vector<int> *a) {
	concurrent_vector<int> tmp(a->size());
	size_t e = 0, o = 1, i = 0;
	while (e < a->size() && o < a->size()) {
		if (a->at(e) < a->at(o)) {
			tmp.at(i++) = (a->at(e));
			e = e + 2;
		}
		else {
			tmp.at(i++) = (a->at(o));
			o = o + 2;
		}
	}
	while (e < a->size()) {
		tmp.at(i++) = (a->at(e));
		e = e + 2;
	}
	while (o < a->size()) {
		tmp.at(i++) = (a->at(o));
		o = o + 2;

	}
	tbb::parallel_for(size_t(0), a->size(), [&](size_t k) {
		a->at(k) = tmp.at(k);
	});
	
}



void shellSortLinear(concurrent_vector<int> *a) {
	int step = 0;
	int iter = 0;
	while (calculateStep(iter++) < a->size() / 3) {
		step = calculateStep(iter);
	}
	while (--iter >= 0) {
		step = calculateStep(iter);
		batcherLinear(a, step);
	}
	mergeTwo(a);
}

void shellSortTbb(concurrent_vector<int> *a) {
	int step = 0;
	int iter = 0;
	while (calculateStep(iter++) < a->size() / 3) {
		step = calculateStep(iter);
	}
	while (--iter >= 0) {
		step = calculateStep(iter);
		batcherTbb(a, step);
	}
	mergeTwoTbb(a);

}

int main(int argc, char *argv[]) {
	int elementsNumber = static_cast<int>(ceil(pow(2, 14)));
	int a = 0;
	int b = 10000000;
	concurrent_vector<int> arr_linear;
	
	if (cmdOptionExists(argv, argv + argc, "-n")) {
		char *wcount = getCmdOption(argv, argv + argc, "-n");
		elementsNumber = atoi(wcount);
	}

	if (cmdOptionExists(argv, argv + argc, "-a")) {
		char *wcount = getCmdOption(argv, argv + argc, "-a");
		a = atoi(wcount);
	}

	if (cmdOptionExists(argv, argv + argc, "-b")) {
		char *wcount = getCmdOption(argv, argv + argc, "-b");
		b = atoi(wcount);
	}

	arr_linear = generateRandomArray(elementsNumber, a, b);
	concurrent_vector<int> arr_tbb(arr_linear);

	double start_linear = omp_get_wtime();
	shellSortLinear(&arr_linear);
	double end_linear = omp_get_wtime();

	tbb::task_scheduler_init init;
	double start_parallel = omp_get_wtime();
	shellSortTbb(&arr_tbb);

	double end_parallel = omp_get_wtime();

	if (check(&arr_tbb, elementsNumber)) {
		printf("\nOK: array is lineary sorted");
	}
	else {
		printf("\n ERROR: array is not sorted");
	}
	printf("\nSorted in %f tbb %f, acceleration %f", end_linear - start_linear,\
		end_parallel - start_parallel, end_parallel - start_parallel / end_linear - start_linear);
	
	return 0;
}