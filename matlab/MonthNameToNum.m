function MonthNum = MonthNameToNum(Month)

if strcmp(lower(Month),'jan')
    MonthNum = 1;
elseif strcmp(lower(Month),'feb')
    MonthNum = 2;
elseif strcmp(lower(Month),'mar')
    MonthNum = 3;
elseif strcmp(lower(Month),'apr')
    MonthNum = 4;
elseif strcmp(lower(Month),'may')
    MonthNum = 5;
elseif strcmp(lower(Month),'jun')
    MonthNum = 6;
elseif strcmp(lower(Month),'jul')
    MonthNum = 7;
elseif strcmp(lower(Month),'aug')
    MonthNum = 8;
elseif strcmp(lower(Month),'sep')
    MonthNum = 9;
elseif strcmp(lower(Month),'oct')
    MonthNum = 10;
elseif strcmp(lower(Month),'nov')
    MonthNum = 11;
elseif strcmp(lower(Month),'dec')
    MonthNum = 12;
else
    disp('Invalid Month Name');
    MonthNum = 0;
end