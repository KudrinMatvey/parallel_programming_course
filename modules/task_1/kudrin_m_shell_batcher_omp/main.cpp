char* getCmdOption(char **begin, char **end, const std::string& option) {
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
        return *itr;
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option) {
    return std::find(begin, end, option) != end;
}

int* generateRandomArray(const int elementsNumber, const int min, const int max) {
	int* arr = new int[elementsNumber];
	for (int j = 0; j < elementsNumber; j++)
		arr[j] = rand() % (max - min);
	return arr;
}

bool check(int* arr, int elementsNumber) 
{
	bool flag = true;
	int min = arr[0];
	for(int i = 1; i < elementsNumber; i++)
	{
		if( arr[i] < min ) 
		{
			flag = false;
		}
	}
	return flag;
}

int main(int argc, char *argv[])
{
	int elementsNumber = 1000;
	int a = 0;
	int b = 10000;
	int* arr;
	
	if (cmdOptionExists(argv, argv + argc, "-n"))
	{
		char *wcount = getCmdOption(argv, argv + argc, "-a1");
		elementsNumber = atoi(wcount);
	}

	if (cmdOptionExists(argv, argv + argc, "-a")) 
	{
		char *wcount = getCmdOption(argv, argv + argc, "-b1");
		a = atoi(wcount);
	}

	if (cmdOptionExists(argv, argv + argc, "-b")) 
	{
		char *wcount = getCmdOption(argv, argv + argc, "-a2");
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