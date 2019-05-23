// Copyright 2019 Kudrin Matvey

#include <tbb/tbb.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "iostream"
#include <algorithm>


using namespace tbb;


char* getCmdOption(char **begin, char **end, const std::string& option) {
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
        return *itr;
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}

int calculateStep(int iter) {
    return static_cast<int>(ceil(pow(2, iter)));
}


void Foo(float *a, size_t l, size_t n, int p) {
    for (int z = n; z != 0; z--) {
        int k = static_cast<int>(ceil(pow(2, z)));
        for (size_t j = k % p; j + k < n; j += (k + k))
            for (size_t i = 0; i < n - j - k; i++)
                if ((j + i) / (p + p) == (j + i + k) / (p + p))
                if(a[l + j + i - 1] > a[l + j + i + k - 1]){
                    std::swap(a[l + j + i - 1], a[l + j + i + k - 1]);
                }
    }
}

void ParallelOddEvenMerge( float a[], size_t n ) {
    int l = 1;
    for (size_t p = l; p < n; p += p) {
        tbb::parallel_for(size_t(0), n, [&](size_t i) {
            Foo(a, i, n, p);
        });
    }
}

void SerialOddEvenMerge( float a[], size_t n ) {
    int l = 1;
    for (size_t p = l; p < n; p += p) {
        for(size_t i = 0;i < n; i ++ )
        Foo(a, i, n, p);
    }
        //    Foo(a, i, n);
}

void batcherLinear(float *a, const int size ,const int step) {
    //std::cout << "k \n";
    int tmpsize = size/step + 1;
    float *tmp = new float[tmpsize];
    int start;
    for (start = 0; start < step; start++) {
    //std::cout << start << " "<< size << " "<< step << std::endl;
        for (int i = start, j = 0; i < size; i += step, j++) {
            tmp[j] = a[i];
        }
        SerialOddEvenMerge(tmp, tmpsize);
        int i = start, j = 0;
        for (; i < size - start; i += step, j++) {
            a[i] = tmp[j];
        }

    }  // end of parallel section
    //delete[] tmp;
}

void batcherTbb(float *a, const int size, const int step) {
    float *tmp = new float[size/step];
    int start;
    int tmpsize;
            std::cout << step << std::endl;
    for (start = 0; start < step; start++) {
        tmpsize = size/step + 1;
        for (int i = start, j = 0; i < size; i += step, j++) {
            tmp[j] = a[i];
        }
        ParallelOddEvenMerge(tmp, tmpsize);
        int i = start, j = 0;
        for (; i < size - start; i += step, j++) {
            a[i] = tmp[j];
        }
    }  // end of parallel section
        std::cout << step << std::endl;
    //delete[] tmp;
}


void shellSortLinear(float *a, int size) {
    int step = 0;
    int iter = 0;
    while (calculateStep(iter++) < size / 3) {
        step = calculateStep(iter);
    }
    while (--iter >= 0) {
        step = calculateStep(iter);
        batcherLinear(a, size, step);
        // SerialOddEvenMerge(a, size);

        //batcherLinear(a, step);
    }
}

void shellSortTbb(float *a, int size) {
    int step = 0;
    int iter = 0;
    while (calculateStep(iter++) < size / 3) {
        step = calculateStep(iter);
    }
    for(int t = iter; t >= 0 ;--t) {
        step = calculateStep(t);
    printf("hi %d \n", t);

        batcherTbb(a, size, step);
        // ParallelOddEvenMerge(a, size);
        //batcherTbb(a, step);
    }
    printf("h");
    fflush(stdout);
}

float* generateRandomArray(const int en, int min, int max) {
    if (min >= max) {
        max = min + max;
    }
    float *arr = new float[en];
    srand(static_cast<unsigned int>(time(NULL)));
    float k = 67.5;
    for (int j = 0; j < en; j++)
        arr[j] = min + (std::rand() % (max - min)) * k;
    return arr;
}


int main(int argc, char *argv[]) {
    int elementsNumber = 512;
    int a = 0;
    int b = 2048000;
    float *arr_linear, *arr_tbb;

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
    arr_tbb = new float[elementsNumber];

    // TODO: to arraycopy
    for (int i = 0; i < elementsNumber; ++i) {
        arr_tbb[i] = arr_linear[i];
    }



    double start_linear = omp_get_wtime();
    shellSortLinear(arr_linear, elementsNumber);
    double end_linear = omp_get_wtime();

    tbb::task_scheduler_init init(4);
    double start_parallel = omp_get_wtime();
    shellSortTbb(arr_tbb, elementsNumber);
    double end_parallel = omp_get_wtime();
    init.terminate();


    printf("\nOK: array is lineary sorted in %f tbb %f", end_linear - start_linear, end_parallel - start_parallel);
    std::cout << "nOK: array is lineary sorted in  tbb" 
    // <<  end_linear - start_linear << end_parallel - start_parallel 
    << std::endl;
    std::cout.flush();
    fflush(stdout);
    return 0;
}
