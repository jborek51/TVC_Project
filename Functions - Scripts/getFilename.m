if Computer == 1
    parent = 'C:\Users\John Jr\Dropbox (UNC Charlotte)\Graduate\Modern Control Theory\Project\Model\Results';
else
    parent = 'C:\Users\jbore\Dropbox (UNC Charlotte)\Graduate\Modern Control Theory\Project\Model\Results';
end
date = datestr(now, 'mm-dd');
cond = 1;   test = 0;
while cond == 1
    filename = strcat('TVC_',date,'_Results_P',num2str(P1),'_C',num2str(CTRL.k_s-1),'_W',num2str(V_W),'_',num2str(test),'.mat');
    if exist(filename,'file')
        test = test + 1;
    else
        cond = 0;
    end
end