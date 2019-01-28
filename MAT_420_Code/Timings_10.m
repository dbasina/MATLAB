A=randn(10,1);
B=A;
C=A;
D=A;
E=A;

fprintf('quicksort: \t\n')
tic
quicksort(A);
toc

fprintf('bubblesort: \t\n')
tic
bubblesort(B);
toc

fprintf('insertionsort: \t\n')
tic
insertion(C');
toc

fprintf('mergesort: \t\n')
tic
mergesort(D);
toc

fprintf('MATLABsort: \t\n')
tic
sort(E);
toc
