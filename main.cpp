#include <string>
#include <fstream>  // работа с файлами
#include <iostream>
#include <cstdlib> // для получения рандомных чисел
#include <ctime>
#include <chrono> // для получения времени
#include <algorithm>
#include <vector>

using namespace std;


/* СОРТИРОВКА ПОДСЧЁТОМ */

//int getMax(int* arr, int n) {
//    int max = arr[0];
//    for (int i = 1; i < n; i++) {
//        if (arr[i] > max)
//            max = arr[i];
//    }
//
//    return max; //максимальный элемент массива
//}
//
//void countSort(int* arr, int n) {
//    int* output = new int[n]; // выделяем память под временный массив
//    int max = getMax(arr, n); // находим максимальный элемент массива
//    int* count = new int[max + 1]; //выделяем память под массив частот
//
//    for (int i = 0; i < max + 1; i++)
//        count[i] = 0; //инициализируем массив частот нулями
//    for (int i = 0; i < n; i++)
//        count[arr[i]]++; // подсчитываем частоту повторения элементов
//    for (int i = 1; i < max + 1; i++)
//        count[i] += count[i - 1]; // вычисляем накопленные частоты
//    for (int i = n - 1; i >= 0; i--) {
//        output[count[arr[i]] - 1] = arr[i]; // устанавливаем элемент на корректное место во временный массив
//        count[arr[i]]--; // уменьшаем частоту, так как уже записали элемент
//    }
//    // реализация алгоритма без использования накопленных частот 
//    //int p = 0;
//    // for (int j = 0; j <= max; ++j) {
//    //    for (unsigned int i = 0; i < count[j]; ++i)
//    //        output[p++] = j;
//    //    }
//    for (int i = 0; i < n; i++) {
//        arr[i] = output[i]; // копируем отсортированный массив во входной 
//    }
//
//    delete[] count; // освобождаем память
//    delete[] output;
//}

/* СОРТИРОВКА ВСТАВКАМИ */

//void insertionSort(int* arr, int n)
//{
//    int i, key, j;
//    for (i = 1; i < n; i++)
//    {
//        key = arr[i];
//        j = i - 1;
//
//        // Move elements of arr[0..i-1], that are greater than key, to one
//        // position ahead of their current position
//        while (j >= 0 && arr[j] > key)
//        {
//            arr[j + 1] = arr[j];
//            j = j - 1;
//        }
//        arr[j + 1] = key;
//    }
//}

/* СОРТИРОВКА СЛИЯНИЕМ */

//// Merges two subarrays of array[].
//// First subarray is arr[begin..mid]
//// Second subarray is arr[mid+1..end]
//void merge(int* array, int const left, int const mid, int const right)
//{
//    auto const subArrayOne = mid - left + 1;
//    auto const subArrayTwo = right - mid;
//
//    // Create temp arrays
//    auto* leftArray = new int[subArrayOne],
//        * rightArray = new int[subArrayTwo];
//
//    // Copy data to temp arrays leftArray[] and rightArray[]
//    for (auto i = 0; i < subArrayOne; i++)
//        leftArray[i] = array[left + i];
//    for (auto j = 0; j < subArrayTwo; j++)
//        rightArray[j] = array[mid + 1 + j];
//
//    auto indexOfSubArrayOne = 0, // Initial index of first sub-array
//        indexOfSubArrayTwo = 0; // Initial index of second sub-array
//    int indexOfMergedArray = left; // Initial index of merged array
//
//    // Merge the temp arrays back into array[left..right]
//    while (indexOfSubArrayOne < subArrayOne && indexOfSubArrayTwo < subArrayTwo) {
//        if (leftArray[indexOfSubArrayOne] <= rightArray[indexOfSubArrayTwo]) {
//            array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
//            indexOfSubArrayOne++;
//        }
//        else {
//            array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
//            indexOfSubArrayTwo++;
//        }
//        indexOfMergedArray++;
//    }
//    // Copy the remaining elements of
//    // left[], if there are any
//    while (indexOfSubArrayOne < subArrayOne) {
//        array[indexOfMergedArray] = leftArray[indexOfSubArrayOne];
//        indexOfSubArrayOne++;
//        indexOfMergedArray++;
//    }
//    // Copy the remaining elements of
//    // right[], if there are any
//    while (indexOfSubArrayTwo < subArrayTwo) {
//        array[indexOfMergedArray] = rightArray[indexOfSubArrayTwo];
//        indexOfSubArrayTwo++;
//        indexOfMergedArray++;
//    }
//    delete[] leftArray;
//    delete[] rightArray;
//}
//
//// begin is for left index and end is
//// right index of the sub-array
//// of arr to be sorted */
//void mergeSort(int* array, int const begin, int const end)
//{
//    if (begin >= end)
//        return; // Returns recursively
//
//    auto mid = begin + (end - begin) / 2;
//    mergeSort(array, begin, mid);
//    mergeSort(array, mid + 1, end);
//    merge(array, begin, mid, end);
//}

/* БЫСТРАЯ СОРТИРОВКА */

//// A utility function to swap two elements
//void swap(int* a, int* b)
//{
//    int t = *a;
//    *a = *b;
//    *b = t;
//}
//
///* This function takes last element as pivot, places
//the pivot element at its correct position in sorted
//array, and places all smaller (smaller than pivot)
//to left of pivot and all greater elements to right
//of pivot */
//int partition(int* arr, int low, int high)
//{
//    int pivot = arr[high]; // pivot
//    int i
//        = (low
//            - 1); // Index of smaller element and indicates
//    // the right position of pivot found so far
//
//    for (int j = low; j <= high - 1; j++) {
//        // If current element is smaller than the pivot
//        if (arr[j] < pivot) {
//            i++; // increment index of smaller element
//            swap(&arr[i], &arr[j]);
//        }
//    }
//    swap(&arr[i + 1], &arr[high]);
//    return (i + 1);
//}
//
///* The main function that implements QuickSort
//arr[] --> Array to be sorted,
//low --> Starting index,
//high --> Ending index */
//void quickSort(int* arr, int low, int high)
//{
//    if (low < high) {
//        /* pi is partitioning index, arr[p] is now
//        at right place */
//        int pi = partition(arr, low, high);
//
//        // Separately sort elements before
//        // partition and after partition
//        quickSort(arr, low, pi - 1);
//        quickSort(arr, pi + 1, high);
//    }
//}

/* ЦИФРОВАЯ СОРТИРОВКА */

//int getMax(int* arr, int n) // функция для нахождения максимального элемента массива
//{
//    int max = arr[0];
//    for (int i = 1; i < n; i++)
//        if (arr[i] > max)
//            max = arr[i];
//    return max;
//}
//
//void countSort(int* arr, int n, int exp)
//{
//    int* output = new int[n];
//    int count[10] = { 0 }; // выделяем память под массив частот размером 10, так как используем 10-чную систему счисления
//
//    // вычисление разряда элемента arr[i] выполняет выражение 
//    // (arr[i] / exp) % 10
//    for (int i = 0; i < n; i++)
//        count[(arr[i] / exp) % 10]++;  // вычисляем значение разряда и считаем его частоту
//
//    for (int i = 1; i < 10; i++)
//        count[i] += count[i - 1];
//
//    for (int i = n - 1; i >= 0; i--) {
//        output[count[(arr[i] / exp) % 10] - 1] = arr[i];
//        count[(arr[i] / exp) % 10]--;
//    }
//
//    for (int i = 0; i < n; i++)
//        arr[i] = output[i];
//}
//
//void radixSort(int* arr, int n)
//{
//    int m = getMax(arr, n); // находим максимальное число
//    // выполняем сортировку массива, пока не дойдем до старшего разряда наибольшего числа 
//
//    int exp = 1;
//    for (int exp = 1; m / exp > 0; exp *= 10)
//        countSort(arr, n, exp);
//}

/* СОРТИРОВКА ВЫБОРОМ */

//void swap(int* xp, int* yp)
//{
//    int temp = *xp;
//    *xp = *yp;
//    *yp = temp;
//}
//
//void selectionSort(int* arr, int n)
//{
//    int i, j, min_idx;
//
//    // One by one move boundary of
//    // unsorted subarray
//    for (i = 0; i < n - 1; i++)
//    {
//
//        // Find the minimum element in
//        // unsorted array
//        min_idx = i;
//        for (j = i + 1; j < n; j++)
//            if (arr[j] < arr[min_idx])
//                min_idx = j;
//
//        // Swap the found minimum element
//        // with the first element
//        if (min_idx != i)
//            swap(&arr[min_idx], &arr[i]);
//    }
//}

/* TIM СОРТИРОВКА */

//const int RUN = 32;
//
//// Merge function merges the sorted runs
//void merge(int* arr, int l, int m, int r)
//{
//    // Original array is broken in two parts left and right array
//    int len1 = m - l + 1;
//    int len2 = r - m;
//    int* left = new int[len1];
//    int* right = new int[len2];
//    for (int i = 0; i < len1; i++)
//        left[i] = arr[l + i];
//    for (int i = 0; i < len2; i++)
//        right[i] = arr[m + 1 + i];
//
//    int i = 0;
//    int j = 0;
//    int k = l;
//
//    // After comparing, we merge those two array in larger sub array
//    while (i < len1 && j < len2)
//    {
//        if (left[i] <= right[j])
//        {
//            arr[k] = left[i];
//            i++;
//        }
//        else
//        {
//            arr[k] = right[j];
//            j++;
//        }
//        k++;
//    }
//
//    // Copy remaining elements of left, if any
//    while (i < len1)
//    {
//        arr[k] = left[i];
//        k++;
//        i++;
//    }
//
//    // Copy remaining element of right, if any
//    while (j < len2)
//    {
//        arr[k] = right[j];
//        k++;
//        j++;
//    }
//
//    delete[] right;
//    delete[] left;
//}
//
//// This function sorts array from left index to to right index which is of size atmost RUN
//void insertionSort(int* arr, int left, int right)
//{
//    for (int i = left + 1; i <= right; i++)
//    {
//        int temp = arr[i];
//        int j = i - 1;
//        while (j >= left && arr[j] > temp)
//        {
//            arr[j + 1] = arr[j];
//            j--;
//        }
//        arr[j + 1] = temp;
//    }
//}
//
//// Iterative Timsort function to sort the array[0...n-1] (similar to merge sort)
//void timSort(int* arr, int n)
//{
//    // Sort individual subarrays of size RUN
//    for (int i = 0; i < n; i += RUN)
//    {
//        int a = i + RUN - 1;
//        int b = n - 1;
//        insertionSort(arr, i, min(a, b));
//    }
//
//    // Start merging from size RUN (or 32). It will merge to form size 64, then 128, 256 and so on ....
//    for (int size = RUN; size < n;
//        size = 2 * size)
//    {
//        // pick starting point of left sub array. We are going to merge arr[left..left+size-1] and arr[left+size, left+2*size-1]
//        // After every merge, we increase left by 2*size
//        for (int left = 0; left < n;
//            left += 2 * size)
//        {
//            // find ending point of left sub array mid+1 is starting point of right sub array
//            int mid = left + size - 1;
//            int right = min((left + 2 * size - 1), (n - 1));
//
//            // merge sub array arr[left.....mid] & arr[mid+1....right]
//            if (mid < right)
//                merge(arr, left, mid, right);
//        }
//    }
//}

/* СОРТИРОВКА ПУЗЫРЬКОМ */

//void bubbleSort(int* arr, int n)
//{
//    for (int i = 0; i < n - 1; i++)
//    {
//        // Last i elements are already in place
//        for (int j = 0; j < n - i - 1; j++) {
//            if (arr[j] > arr[j + 1])
//                swap(arr[j], arr[j + 1]);
//        }
//    }
//}

/* КАРМАННАЯ СОРТИРОВКА */

//struct bucket {
//    int count; // количество элементов в корзине
//    int* data; // массив элементов корзины
//};
//
//int getMax(int* arr, int n) // функция для нахождения максимального элемента массива
//{
//    int max = arr[0];
//    for (int i = 1; i < n; i++)
//        if (arr[i] > max)
//            max = arr[i];
//    return max;
//}
//
//int getExp(int value) {
//    int exp = 1;
//
//    while (value > 10)
//    {
//        value /= 10;
//        exp *= 10;
//    }
//
//    return exp;
//}
//
//void countSort(int* arr, int n)
//{
//    if (!n)
//        return;
//    int* output = new int[n];
//    int max = getMax(arr, n);
//    size_t* count = new size_t[max + 1];
//
//    for (int i = 0; i < max + 1; i++)
//        count[i] = 0;
//
//    for (int i = 0; i < n; i++)
//        count[arr[i]]++;  // подсчитываем  частоту повторения элементов
//
//    for (int i = 1; i <= max; i++)
//        count[i] += count[i - 1];  // вычисляем накопленные суммы
//
//    for (int i = n - 1; i >= 0; i--) {
//        output[count[arr[i]] - 1] = arr[i];  // устанавливаем элемент на корректное место 
//        count[arr[i]]--; // уменьшаем частоту, так как уже записали элемент
//    }
//
//    for (int i = 0; i < n; i++)
//        arr[i] = output[i];
//
//    delete[] output;
//    delete[] count;
//}
//
//void bucketSort(int* arr, int n)
//{
//    struct bucket buckets[10];
//    // вычисляем значение экспоненты
//    int exp = getExp(getMax(arr, n));
//
//    for (int i = 0; i < 10; i++)
//    {
//        buckets[i].count = 0;
//        buckets[i].data = new int[n];
//    }
//    for (int i = 0; i < n; i++) {
//        int bi = arr[i] / exp; // вычисляем индекс корзины
//        buckets[bi].data[buckets[bi].count++] = arr[i]; // добавляем элемент в корзину
//    }
//
//    for (int i = 0; i <= 10; i++) {
//        cout << buckets[i].count << " ";
//        for (size_t j = 0; j < buckets[i].count; j++)
//        {
//            cout << buckets[i].data[j] << " ";
//        }
//        cout << endl;
//    }
//
//    int index = 0;
//    for (int i = 0; i < n; i++)
//    {
//        countSort(buckets[i].data, buckets[i].count);
//        for (int j = 0; j < buckets[i].count; j++) {
//            arr[index] = buckets[i].data[j];
//            index++;
//        }
//        buckets[i].data = 0;
//    }
//}

int getMax(int* A, int n)
{
    if (n == 1)
        return A[0];
    return max(A[n - 1], getMax(A, n - 1));
}
int getExp(int value) {
    int exp = 1;

    while (value > 10)
    {
        value /= 10;
        exp *= 10;
    }

    return exp;
}
void bucketSort(int* arr, int n)
{
    int max = getMax(arr, n);
    int exp = getExp(getMax(arr, n));
    vector<int>* b = new vector<int>[n];
    //vector<int> b;
    // 2) Put array elements in different buckets 
    for (int i = 0; i < n; i++)
    {
        //int bi = n * arr[i] / (max + 1); // Index in bucket
        int bi = arr[i] / exp; // вычисляем индекс корзины
        b[bi].push_back(arr[i]);
    }

    for (int i = 0; i < n; i++)
    {
        cout << b[i].size() << " ";
        for (int j = 0; j < b[i].size(); j++)
        {
            cout << b[i].at(j) << " ";
        }
        cout << endl;
    }

    // 3) Sort individual buckets 
    for (int i = 0; i < n; i++)
        sort(b[i].begin(), b[i].end());
    // 4) Concatenate all buckets into arr[] 
    int index = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < b[i].size(); j++)
            arr[index++] = b[i][j];

    delete[] b;
}

void write_array(const string& filename, const int* arr, const int n)
{
    fstream fs;

    fs.open(filename, fstream::out);
    if (fs.is_open()) // проверяем, что файл успешно открыт
    {
        fs << n << '\n';  // записываем размер массива
        for (int i = 0; i < n; i++)
            fs << arr[i] << ' '; // записываем значения через пробел
        fs << '\n';

        fs.close(); // закрываем файл
    }
}

void read_array(const string& filename, int*& arr, int& n)
{
    fstream fs;

    fs.open(filename, fstream::in);
    if (fs.is_open()) // проверяем, что файл успешно открыт
    {
        fs >> n;  // читаем размер массива
        arr = new int[n];
        for (int i = 0; i < n; i++)
            fs >> arr[i]; // читаем из файла разделительные символы - пробел и перенос строки

        fs.close(); // закрываем файл
    }
}

int main()
{
    string infile = "input_data.txt";
    string outfile = "output_data.txt";
    int size = 1000; // размер файла
    //int size = 200000; // размер файла
    int* rand_arr = new int[size];

    srand(time(nullptr)); // используем текущее время, чтобы сгенерировать рандомные значения
    int lef_border = 10;
    int range_len = 11; // правая граница = range_len + left_border
    //int lef_border = 1;
    //int range_len = 1000; // правая граница = range_len + left_border

    for (int i = 0; i < size; i++)
        rand_arr[i] = lef_border + rand() % range_len; // генерируем число в указанном диапазоне и записываем в массив

    write_array(infile, rand_arr, size); // записываем массив в файл

    int* array_from_file = nullptr;
    int n = 0;
    read_array(infile, array_from_file, n); // читаем массив из файла

    auto start = chrono::high_resolution_clock::now(); // сохраняем время начала работы алгоритма
    //countSort(array_from_file, n); // запускаем сортировку
    //insertionSort(array_from_file, n);
    //mergeSort(array_from_file, 0, n - 1);
    //quickSort(array_from_file, 0, n - 1);
    //radixSort(array_from_file, n);
    //selectionSort(array_from_file, n);
    //timSort(array_from_file, n);
    //bubbleSort(array_from_file, n);
    bucketSort(array_from_file, n);
    auto finish = chrono::high_resolution_clock::now(); // сохраняем время конца работы алгоритма
    chrono::duration<double> elapsed = finish - start;
    cout << "Elapsed time: " << elapsed.count() << " sec" << endl; // вычисляем продолжительность работы в сек. и выводим на экран

    write_array(outfile, array_from_file, n);

    delete[] rand_arr;
    delete[] array_from_file;
    
    return 0;
}