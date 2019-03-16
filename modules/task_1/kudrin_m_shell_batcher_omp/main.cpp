#include <math.h>
#include <stdio.h>
#include <iostream>
#include <vector>


using namespace std;

void compex(int &a, int &b)
{
    if (a > b) std::swap(a, b);
};


void oddeven_merge_sort(std::vector<int>& arr)
{
    const int length = arr.size();
    int t = (int)ceil(log2(length));
    int p = (int)pow(2, t - 1);

    while (p > 0) {
        int q = (int)pow(2, t - 1);
        int r = 0;
        int d = p;

        while (d > 0) {
            for (int i = 0; i < length - d; ++i) {
                if ((i & p) == r) {
                    compex(arr[i], arr[i + d]);
                }
            }

            d = q - p;
            q /= 2;
            r = p;
        }
        p /= 2;
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

vector<int> generateRandomArray(const int elementsNumber, const int min, const int max) {
    vector<int> arr(elementsNumber);
    for (int j = 0; j < elementsNumber; j++)
        arr[j] = rand() % (max - min);
    return arr;
}

bool check(vector<int>& arr, int elementsNumber)
{
    bool flag = true;
    int min = arr[0];
    for (int i = 1; i < elementsNumber; i++)
    {
        if (arr[i] < min)
        {
            flag = false;
        }
    }
    return flag;
}

int calculateStep(int iter)
{
    int step = 0;
    if (iter % 2)
    {
        step = (int)(8 * pow(2, iter) - 6 * pow(2, (iter + 1) / 2) + 1);
    }
    else
    {
        step = (int)(9 * pow(2, iter) - 9 * pow(2, iter / 2) + 1);
    }
    return step;
}

void batcher(vector<int>& a, const int step) 
{
    vector<int> tmp((a.size() / step) + 2);
    for(int start = 0; start < step; start++)
    {
        
        for (unsigned int i = start, j = 0; i < a.size(); i += step, j++)
        {
            tmp[j] = a[i];
        }
        oddeven_merge_sort(tmp);
        for (unsigned int i = start, j = 0; i < a.size() - start; i += step, j++)
        {
            a[i] = tmp[j];
        }
    }
}

void shellSort(vector<int>& a, int size)
{
    int step = 0;
    int iter = 0;
    while (calculateStep(iter++) < size / 3)
    {
        step = calculateStep(iter);
    }
    while (--iter >= 0)
    {
        step = calculateStep(iter);
        batcher(a, step);
    }
}

int main(int argc, char *argv[])
{
    int elementsNumber = 1000;
    int a = 0;
    int b = 10000;
    vector<int> arr;

    if (cmdOptionExists(argv, argv + argc, "-n"))
    {
        char *wcount = getCmdOption(argv, argv + argc, "-n");
        elementsNumber = atoi(wcount);
    }

    if (cmdOptionExists(argv, argv + argc, "-a"))
    {
        char *wcount = getCmdOption(argv, argv + argc, "-a");
        a = atoi(wcount);
    }

    if (cmdOptionExists(argv, argv + argc, "-b"))
    {
        char *wcount = getCmdOption(argv, argv + argc, "-b");
        b = atoi(wcount);
    }

    arr = generateRandomArray(elementsNumber, a, b);
    shellSort(arr, elementsNumber);

    for (int i = 0; i < elementsNumber; i++)
    {
        printf("%d ", arr[i]);
    }
    if (check(arr, elementsNumber))
    {
        printf("\nOK: array is sorted");
    }
    else
    {
        printf("\n ERROR: array is not sorted");
    }
    return 0;
}