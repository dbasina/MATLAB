function f= bubblesort(a)
   for i= 1:1:length(a)          % loop traverses both trackers throught the array n times.
       for j= 1:1:length(a)-1      % loop traverses both trackers throughtout the array 1 time. 
            if (a(j) > a(j+1))
                temp=a(j+1);
                a(j+1)=a(j);
                a(j)=temp;
         
            end
       end
   end
   f=a;
end
            
                
            
       