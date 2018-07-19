function Julian = Grg2Juln(date)
    year = date-mod(date,10000);
    month = date-year-mod((date-year),100);
    day = date-year-month;
    
    switch month
    case 100
        Julian = 0;
    case 200
        Julian = 31;
    case 300
        Julian = 59;
    case 400
        Julian = 90;
    case 500
        Julian = 120;
    case 600
        Julian = 151;
    case 700
        Julian = 181;
    case 800
        Julian = 212;
    case 900
        Julian = 243;
    case 1000
        Julian = 273;
    case 1100
        Julian = 304;
    case 1200
        Julian = 334;
    case 1300
        Julian = 365;
    end
    
    Julian = year/10 + Julian + day;
    
    if (((mod(year,40000)==0)&(mod(year,1000000)~=0))|(mod(year,10000000)==0))&(month>200)
        Julian = Julian + 1;
    end
    