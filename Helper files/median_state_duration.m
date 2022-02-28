function [median_SL,SL_all, stateLength]= median_state_duration(binbounds, optplot)

        Ntime=length(binbounds);
        states=zeros([Ntime, 1]);
        stateLength=zeros([Ntime, 1]);
        %convert boundaries into states
        for i=2:Ntime
            states(i)=states(i-1)+binbounds(i);
        end
        %for each timepoint, find the length of the state that occurs at
        %that time
        for i=1:Ntime
            stateLength(i)=sum(states==states(i));
        end
        %for each unique event, find its duration
        SL_all=zeros([max(states),1]);
        for i=1:max(states)
            SL_all(i)=sum(states==i);
        end
        median_SL=median(SL_all);
        
        if exist('optplot', 'var')
            figure; hist(SL_all)
            title(['mean = ' num2str(mean(SL_all)) ', med = ' num2str(median(SL_all))])
        end