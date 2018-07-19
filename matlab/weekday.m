function DayName = weekday(InDate)

    % check valid date (or at least close)
    if mod(InDate,100) > 31 | mod(floor(InDate/100),100) > 12
        disp('Invalid Date');
        beep; pause(0.5); beep; pause(0.5); beep;
    elseif InDate < 16000000 | InDate > 30000000
        disp('Questionable Date');
        beep; pause(0.5); beep; pause(0.5); beep;
    end
    

	TheDayIWroteThis = 20060101;
	
	day(1).name = 'Sunday';
	day(2).name = 'Monday';
	day(3).name = 'Tuesday';
	day(4).name = 'Wednesday';
	day(5).name = 'Thursday';
	day(6).name = 'Friday';
	day(7).name = 'Saturday';
	
	if InDate == TheDayIWroteThis
        DaysPast = 0;
	else
        TheYear = floor(InDate/10000);
        TheYearIWroteThis = floor(TheDayIWroteThis/10000);
        
        YearsDiff = TheYear - TheYearIWroteThis;
        
        QYears = floor(abs(YearsDiff-1)/4);
        CYears = floor(abs(YearsDiff-1)/100);
        MYears = floor(abs(YearsDiff-1)/400);
        
        if YearsDiff < 0
            if ceil(TheYear/400)*400 < (TheYearIWroteThis-MYears*400) & ceil(TheYear/400)*400 > TheYear
                MYears = MYears + 1;
            end
            if ceil(TheYear/100)*100 < (TheYearIWroteThis-CYears*100) & ceil(TheYear/100)*100 >= TheYear
                CYears = CYears + 1;
            end
            if ceil(TheYear/4)*4 < (TheYearIWroteThis-QYears*4) & ceil(TheYear/4)*4 >= TheYear
                QYears = QYears + 1;
            end
	
            DaysPast = 365*YearsDiff-(MYears-CYears+QYears);
                    
        elseif YearsDiff > 0
            if floor(TheYear/400)*400 > (TheYearIWroteThis+MYears*400) & floor(TheYear/400)*400 <= TheYear
                MYears = MYears + 1;
            end
            if floor(TheYear/100)*100 > (TheYearIWroteThis+CYears*100) & floor(TheYear/100)*100 < TheYear
                CYears = CYears + 1;
            end
            if floor(TheYear/4)*4 > (TheYearIWroteThis+QYears*4) & floor(TheYear/4)*4 < TheYear
                QYears = QYears + 1;
            end
	
            DaysPast = 365*YearsDiff+(MYears-CYears+QYears);
	
        else
            
            DaysPast = 0;
            
        end
        
        % Days difference within year
        DaysPast = DaysPast - mod(Grg2Juln(TheYear*10000+mod(TheDayIWroteThis,10000)),1000) + mod(Grg2Juln(InDate),1000);
        
	end
	
	%disp(sprintf('%d - %s',InDate,day(mod(DaysPast,7)+1).name));
	DayName = day(mod(DaysPast,7)+1).name;
    