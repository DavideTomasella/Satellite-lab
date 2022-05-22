function test_downConverterFilter(t , reader , f_doppler , t_delay, fsampling , ...
                                    chipRate , f_or_t , inital_doppler , ...
                                            initial_delay, phase)
    down = DownconverterFilter();
    down.configDownConverter(fsampling);
    down.configFilter(1,350);
    down.downFilter(reader,0.6*chipRate,chipRate);
    down.downConverter(reader,f_doppler,t_delay,phase);
    down.configFilter(1,350);
    down.downFilter(reader,0.6*chipRate,chipRate);
    
    if f_or_t == 0
        df = f_doppler - inital_doppler;
        if df >= 0
            type = "-";
        else
            type = "--";
        end
        value = num2str(df);
        plot(t,reader.IQsamples(:,1),type,"DisplayName","\Deltaf = "+value);
    else
        dt = t_delay - initial_delay;
        if dt >= 0
            type = "-";
        else
            type = "--";
        end
        value = num2str(dt);
        plot(t,reader.IQsamples(:,1),type,"DisplayName","\Deltat = "+value);
    end    
end