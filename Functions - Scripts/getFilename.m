if Computer == 1
    parent = 'C:\Users\John Jr\Dropbox (UNC Charlotte)\Projects\Rocket\TVC_Project\Results\2-12';
else
    parent = 'C:\Users\jbore\Dropbox (UNC Charlotte)\Projects\Rocket\TVC_Project\Results\2-12';
end
date = datestr(now, 'mm-dd');
cond = 1;   test = 0;
while cond == 1
    filename = strcat('TVC_',date,'_Results_C',num2str(CTRL.k_s-1),'_',num2str(test),'.mat');
    if exist(filename,'file')
        test = test + 1;
    else
        cond = 0;
    end
end