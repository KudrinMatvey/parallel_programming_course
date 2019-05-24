//// Copyright 2019 Kudrin Matvey
//
//#include <tbb/tbb.h>
//#include <math.h>
//#include <omp.h>
//#include <time.h>
//#include "iostream"
//#include <algorithm>
//
//
//using namespace tbb;
//
//
//char* getCmdOption(char **begin, char **end, const std::string& option) {
//    char **itr = std::find(begin, end, option);
//    if (itr != end && ++itr != end)
//        return *itr;
//    return 0;
//}
//
//bool cmdOptionExists(char** begin, char** end, const std::string& option) {
//    return std::find(begin, end, option) != end;
//}
//
//int calculateStep(int iter) {
//    return static_cast<int>(ceil(pow(2, iter)));
//}
//
//
//void Foo(concurrent_vector<float> *a, size_t l, int p) {
//	for (int z = a->size(); z != 0; z--) {
//        int k = static_cast<int>(ceil(pow(2, z)));
//        for (size_t j = k % p; j + k < a->size(); j += (k + k))
//            for (size_t i = 0; i < a->size() - j - k; i++)
//                if ((j + i) / (p + p) == (j + i + k) / (p + p))
//                if(a->begin() + l + j + i - 1 > a->begin() + l + j + i + k - 1){
//                   std::iter_swap(a->begin() + l + j + i - 1,
//						a->begin() + l + j + i + k - 1);
//                }
//    }
//}
//
//void ParallelOddEvenMerge( concurrent_vector<float> *a) {
//    int l = 1;
//    for (size_t p = l; p < a->size(); p += p) {
//        tbb::parallel_for(size_t(0), [&](size_t i) {
//            Foo(a, i, p);
//        });
//    }
//}
//
////void SerialOddEvenMerge( float a[], size_t n ) {
////    int l = 1;
////    for (size_t p = l; p < n; p += p) {
////        for(size_t i = 0;i < n; i ++ )
////        Foo(a, i, n, p);
////    }
////        //    Foo(a, i, n);
////}
//
//void batcherLinear(float *a, const int size ,const int step) {
//    //std::cout << "k \n";
//    int tmpsize = size/step + 1;
//    float *tmp = new float[tmpsize];
//    int start;
//    for (start = 0; start < step; start++) {
//    //std::cout << start << " "<< size << " "<< step << std::endl;
//        for (int i = start, j = 0; i < size; i += step, j++) {
//            tmp[j] = a[i];
//        }
//        SerialOddEvenMerge(tmp, tmpsize);
//        int i = start, j = 0;
//        for (; i < size - start; i += step, j++) {
//            a[i] = tmp[j];
//        }
//
//    }  // end of parallel section
//    //delete[] tmp;
//}
//
//void batcherTbb(concurrent_vector<int> *a, const int step) {
//	concurrent_vector<int> *tmp = new concurrent_vector<int>;
//	int start;
//	for (start = 0; start < step; start++) {
//		for (unsigned int i = start, j = 0; i < a->size(); i += step, j++) {
//			tmp->push_back(a->at(i));
//		}
//		ParallelOddEvenMerge(tmp);
//		unsigned int i = start, j = 0;
//		for (; i < a->size() - start; i += step, j++) {
//			a->at(i) = tmp->at(j);
//		}
//		tmp->clear();
//	}  // end of parallel section
//	delete tmp;
//}
//
//
//void shellSortLinear(float *a, int size) {
//    int step = 0;
//    int iter = 0;
//    while (calculateStep(iter++) < size / 3) {
//        step = calculateStep(iter);
//    }
//    while (--iter >= 0) {
//        step = calculateStep(iter);
//        batcherLinear(a, size, step);
//        // SerialOddEvenMerge(a, size);
//
//        //batcherLinear(a, step);
//    }
//}
//
//void shellSortTbb(float *a, int size) {
//    int step = 0;
//    int iter = 0;
//    while (calculateStep(iter++) < size / 3) {
//        step = calculateStep(iter);
//    }
//    for(int t = iter; t >= 0 ;--t) {
//        step = calculateStep(t);
//    printf("hi %d \n", t);
//        try{
//         batcherTbb(a, size, step);
//
//        } catch (...) {
//            printf("c");
//        }
//        // ParallelOddEvenMerge(a, size);
//        //batcherTbb(a, step);
//    }
//    //printf("h");
//    fflush(stdout);
//}
//
//float* generateRandomArray(const int en, int min, int max) {
//    if (min >= max) {
//        max = min + max;
//    }
//    float *arr = new float[en];
//    srand(static_cast<unsigned int>(time(NULL)));
//    float k = 67.5;
//    for (int j = 0; j < en; j++)
//        arr[j] = min + (std::rand() % (max - min)) * k;
//    return arr;
//}
//
//
//int main(int argc, char *argv[]) {
//    int elementsNumber = 512;
//    int a = 0;
//    int b = 2048000;
//    float *arr_linear, *arr_tbb;
//
//    if (cmdOptionExists(argv, argv + argc, "-n")) {
//        char *wcount = getCmdOption(argv, argv + argc, "-n");
//        elementsNumber = atoi(wcount);
//    }
//
//    if (cmdOptionExists(argv, argv + argc, "-a")) {
//        char *wcount = getCmdOption(argv, argv + argc, "-a");
//        a = atoi(wcount);
//    }
//
//    if (cmdOptionExists(argv, argv + argc, "-b")) {
//        char *wcount = getCmdOption(argv, argv + argc, "-b");
//        b = atoi(wcount);
//    }
//
//    arr_linear = generateRandomArray(elementsNumber, a, b);
//    arr_tbb = new float[elementsNumber];
//
//    // TODO: to arraycopy
//    for (int i = 0; i < elementsNumber; ++i) {
//        arr_tbb[i] = arr_linear[i];
//    }
//
//
//
//    double start_linear = omp_get_wtime();
//    // shellSortLinear(arr_linear, elementsNumber);
//    double end_linear = omp_get_wtime();
//
//    tbb::task_scheduler_init init(4);
//    double start_parallel = omp_get_wtime();
//    shellSortTbb(arr_tbb, elementsNumber);
//    double end_parallel = omp_get_wtime();
//    init.terminate();
//
//
//    printf("\nOK: array is lineary sorted in %f tbb %f", end_linear - start_linear, end_parallel - start_parallel);
//    std::cout << "nOK: array is lineary sorted in  tbb" 
//    // <<  end_linear - start_linear << end_parallel - start_parallel 
//    << std::endl;
//    std::cout.flush();
//    fflush(stdout);
//    return 0;
//}






// Copyright 2019 Kudrin Matvey

#include <tbb/tbb.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <utility>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "tbb/concurrent_vector.h"

using std::vector;
using tbb::blocked_range;
using tbb::concurrent_vector;





class OddEvenSorter {
	concurrent_vector<int> *arr;
	const int p, r, d;

public:
	OddEvenSorter(concurrent_vector<int> *arr, const int p, const int r,
		const int  d) : arr(arr), p(p), r(r), d(d) {}
	void operator()(const blocked_range<int>& range) const {
		int begin = range.begin(), end = range.end();
		for (int i = begin; i != end; i++) {
			if ((p & i) == r) {
				if (arr->at(i) > arr->at(i + d)) {
					int tmp = arr->at(i);
					arr->at(i) = arr->at(i + d);
					arr->at(i + d) = tmp;
				}
			}
		}
	}
};

void Foo(concurrent_vector<int> *a, size_t l, int p) {
	for (int z = a->size(); z >= 0; z--) {
        int k = static_cast<int>(ceil(pow(2, z)));
        for (size_t j = k % p; j + k < a->size(); j += (k + k))
            for (size_t i = 0; i < a->size() - j - k; i++)
                if ((j + i) / (p + p) == (j + i + k) / (p + p))
                if(a->begin() + l + j + i - 1> a->begin() + l + j + i + k - 1){
                   std::iter_swap(a->begin() + l + j + i - 1,
						a->begin() + l + j + i + k - 1);
                }
    }
}

void SerialOddEvenMerge(concurrent_vector<int> *a) {
	int l = 1;
	for (size_t p = l; p < a->size(); p += p) {
		for (size_t i = 1; i < a->size(); i++)
			Foo(a, i + 1, p);
	}
}

class BatcherSorter {
	concurrent_vector<int> *a;
	int step;

public:
	BatcherSorter(concurrent_vector<int> *a, const int step) : a(a),
		step(step) {}
	void operator()(const blocked_range<int>& range) const {
		int begin = range.begin(), end = range.end();
		vector<int> *tmp = new vector<int>(a->size() / step);
		for (int start = begin; start != end; start++) {
			for (unsigned int i = start, j = 0; i < a->size(); i += step, j++) {
				tmp->push_back(a->at(i));
			}
			const int length = tmp->size();
			int t = static_cast<int>(ceil(log2(length)));
			int p = static_cast<int>(pow(2, t - 1));
			while (p > 0) {
				int q = static_cast<int>(pow(2, t - 1));
				int r = 0;
				int d = p;
				int i;
				while (d > 0) {
					for (i = 0; i < length - d; ++i) {
						if ((i & p) == r) {
							if (tmp->at(i) > tmp->at(i + d)) {
								std::iter_swap(tmp->begin() + i,
									tmp->begin() + i + d);
							}
						}
					}
					d = q - p;
					q /= 2;
					r = p;
				}
				p /= 2;
			}

			unsigned int i = start, j = 0;
			for (; i < a->size() - start; i += step, j++) {
				a->at(i) = tmp->at(j);
			}
			tmp->clear();
		}
		delete tmp;
	}
};

//void oddEvenMergeSortLinear(std::vector<int> *arr) {
//	const int length = arr->size();
//	int t = static_cast<int>(ceil(log2(length)));
//	int p = static_cast<int>(pow(2, t - 1));
//	while (p > 0) {
//		int q = static_cast<int>(pow(2, t - 1));
//		int r = 0;
//		int d = p;
//		int i;
//		while (d > 0) {
//			for (i = 0; i < length - d; ++i) {
//				if ((i & p) == r) {
//					if (arr->at(i) > arr->at(i + d)) {
//						std::iter_swap(arr->begin() + i, arr->begin() + i + d);
//					}
//				}
//			}
//			d = q - p;
//			q /= 2;
//			r = p;
//		}
//		p /= 2;
//	}
//}
//
//void oddEvenMergeSortTbb(concurrent_vector<int> *arr) {
//	const int length = arr->size();
//	int t = static_cast<int>(ceil(log2(length)));
//	int p = static_cast<int>(pow(2, t - 1));
//	while (p > 0) {
//		int q = static_cast<int>(pow(2, t - 1));
//		int r = 0;
//		int d = p;
//		while (d > 0) {
//			parallel_for(blocked_range<int>(1, length - d - 1, 40000),
//				OddEvenSorter(arr, p, r, d));
//			d = q - p;
//			q /= 2;
//			r = p;
//		}
//		p /= 2;
//	}
//}

char* getCmdOption(char **begin, char **end, const std::string& option) {
	char **itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
		return *itr;
	return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
	return std::find(begin, end, option) != end;
}

vector<int>* generateRandomArray(const int en, int min, int max) {
	if (min >= max) {
		max = min + max;
	}
	vector<int> *arr = new vector<int>(en);
	srand(static_cast<unsigned int>(time(NULL)));
	for (int j = 0; j < en; j++)
		arr->at(j) = min + (std::rand() % (max - min));
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

int calculateStep(int iter) {
    return static_cast<int>(ceil(pow(2, iter)));
}


void batcherLinear(vector<int> *a, const int step) {
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
	}  // end of parallel section
}

void ParallelOddEvenMerge( concurrent_vector<int> *a) {
    int l = 1;
    for (size_t p = l; p < a->size(); p += p) {
        tbb::parallel_for(size_t(0),a->size(), [&](size_t i) {
            Foo(a, i + 1, p);
        });
    }
}


void batcherTbb(concurrent_vector<int> *a, const int step) {
	concurrent_vector<int> *tmp = new concurrent_vector<int>;
	int start;
	for (start = 0; start < step; start++) {
		for (unsigned int i = start, j = 0; i < a->size(); i += step, j++) {
			tmp->push_back(a->at(i));
		}
		
		ParallelOddEvenMerge(a);
		
		unsigned int i = start, j = 0;
		for (; i < a->size() - start; i += step, j++) {
			a->at(i) = tmp->at(j);
		}
		tmp->clear();
	}  // end of parallel section
	delete tmp;
}


void shellSortLinear(vector<int> *a, int size) {
	int step = 0;
	int iter = 0;
	while (calculateStep(iter++) < size / 3) {
		step = calculateStep(iter);
	}
	while (--iter >= 0) {
		step = calculateStep(iter);
		batcherLinear(a, step);
	}
}

void shellSortTbb(concurrent_vector<int> *a, int size) {
	int step = 0;
	int iter = 0;
	while (calculateStep(iter++) < size / 3) {
		step = calculateStep(iter);
	}
	while (--iter >= 0) {
		step = calculateStep(iter);
		batcherTbb(a, step);
	}
}

int main(int argc, char *argv[]) {
	int elementsNumber = 4;
	int a = 0;
	int b = 10000000;
	vector<int> *arr_linear;
	concurrent_vector<int> *arr_tbb;

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
	arr_tbb = new concurrent_vector<int>(elementsNumber);

	for (int i = 0; i < elementsNumber; ++i) {
		arr_tbb->at(i) = arr_linear->at(i);
	}

//	shellSortLinear(arr_linear, elementsNumber);
	ParallelOddEvenMerge(arr_tbb);
	tbb::task_scheduler_init init;
	//shellSortTbb(arr_tbb, elementsNumber);

	if (check(arr_tbb, elementsNumber)) {
		printf("\nOK: array is lineary sorted");
	}
	else {
		printf("\n ERROR: array is not sorted");
	}
	arr_linear->clear();
	arr_tbb->clear();
	return 0;
}