function Gregorian = Juln2Grg(date)
    year = date-mod(date,1000);
    date = date-year;
    year = year * 10;
    
    if (((mod(year,40000)==0)&(mod(year,1000000)~=0))|(mod(year,10000000)==0))&(date>59)
        Gregorian = 0;
        if date==60
            Gregorian=1;
        end
        date = date - 1;
    else
        Gregorian = 0;
    end

    if date <= 31
        Gregorian = Gregorian + year + 100 + date;
    elseif date <= 59
        Gregorian = Gregorian + year + 200 + date - 31;
    elseif date <= 90
        Gregorian = Gregorian + year + 300 + date - 59;
    elseif date <= 120
        Gregorian = Gregorian + year + 400 + date - 90;
    elseif date <= 151
        Gregorian = Gregorian + year + 500 + date - 120;
    elseif date <= 181
        Gregorian = Gregorian + year + 600 + date - 151;
    elseif date <= 212
        Gregorian = Gregorian + year + 700 + date - 181;
    elseif date <= 243
        Gregorian = Gregorian + year + 800 + date - 212;
    elseif date <= 273
        Gregorian = Gregorian + year + 900 + date - 243;
    elseif date <= 304
        Gregorian = Gregorian + year + 1000 + date - 273;
    elseif date <= 334
        Gregorian = Gregorian + year + 1100 + date - 304;
    else
        Gregorian = Gregorian + year + 1200 + date - 334;
    end
    