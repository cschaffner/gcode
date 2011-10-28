%% timing and status functions for loops

% runs 10 different timers 
% call with index of timer and fraction of nr of loops that is already done
% in intervals of 4,8,16,32,64,128,128,128 seconds, status updates are
% given

% by Christian Schaffner, c.schaffner@uva.nl
% 22 October 2011

function output=timing(index,fracdone)

    persistent timers;
    persistent interv;
    persistent lasttm;

    if isempty(timers) 
        disp 'Initializing timers';
        timers=uint64(zeros(10,1));
        lasttm=zeros(10,1);
        interv=4*ones(10,1);
    end

    output=0;

    if fracdone==0
        fprintf('starting timer with index %u\n',index);
        timers(index)=floor(tic);
        lasttm(index)=0;
        interv(index)=4;
    else 
        tm=toc(timers(index));
        if (tm>interv(index)+lasttm(index))
             lasttm(index)=tm;
             if(interv(index)<100)
                 % if interv<100, double time interval
                 interv(index)=interv(index)*2;
             end
             fprintf('index %2.0f: after %6.1f seconds: %2.2g%% done. expected time left:%6.1f seconds\n',index,tm,100*fracdone,tm*((1/fracdone)-1));
             output=1;
        end
    end

end
