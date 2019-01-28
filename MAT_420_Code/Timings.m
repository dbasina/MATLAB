%% Initialize
x=randn(20000,1);
n=1000;
%%Matlab Sort function
fprintf('Matlab Sort\t\n')
tic
x=sort(x);
toc
temp_1 = x';
temp_2 = x';
temp_3 = x';
temp_13 = x';

%% a) x
%Insertion
fprintf('Insertion Sort: x\t\n')
tic
insertion(temp_1);
toc

%Mergesort
fprintf('MergeSort: x\t\n')

tic
mergesort(temp_2);
toc

%Quicksort
fprintf('QuickSort: x\t\n')
tic
quicksort(temp_3);
toc

%Bubblesort
fprintf('Bubblesort: x\t\n')
tic
bubblesort(temp_13);
toc

%% b) x(n:-1:1)
y= x(n:-1:1);
temp_4 = y';
temp_5 = y';
temp_6 = y';
temp_14= y';
%Insertion
fprintf('Insertion Sort: x(n:-1:1)\t\n')

tic
insertion(temp_4);
toc

%Mergesort
fprintf('MergeSort: x(n:-1:1)\t\n')
tic
mergesort(temp_5);
toc

%Quicksort
fprintf('QuickSort: x(n:-1:1)\t\n')
tic
quicksort(temp_6);
toc

%Bubblesort
fprintf('Bubblesort: x\t\n')
tic
bubblesort(temp_14);
toc

%% c) x([1:2:n,2:2:n])
z= x([1:2:n,2:2:n]);
temp_7 = z';
temp_8 = z';
temp_9 = z';
temp_15 = z';
%Insertion
fprintf('InsertionSort: x([1:2:n,2:2:n])\t\n')
tic
insertion(temp_7);
toc

%Mergesort
fprintf('MergeSort: x([1:2:n,2:2:n])\t\n')
tic
mergesort(temp_8);
toc

%Quicksort
fprintf('QuickSort: x([1:2:n,2:2:n])\t\n')
tic
quicksort(temp_9);
toc

%Bubblesort
fprintf('Bubblesort: x\t\n')
tic
bubblesort(temp_15);
toc
%% d) x(randperm(n))
p=x(randperm(n));
temp_10 = p';
temp_11 = p';
temp_12 = p';
temp_16= p';
%Insertion
fprintf('InsertionSort:x(randperm(n))\t\n')

tic
insertion(temp_10);
toc

%Mergesort
fprintf('MergeSort:x(randperm(n))\t\n')
tic
mergesort(temp_11);
toc

%Quicksort
fprintf('QuickSort:x(randperm(n))\t\n')
tic
quicksort(temp_12);
toc

%Bubblesort
fprintf('Bubblesort: x\t\n')
tic
bubblesort(temp_16);
toc

