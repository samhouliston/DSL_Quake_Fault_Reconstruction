function progScreen(i)
    bs = '\b';
    sp = ' ';
    msg = '%d …\n'; % carriage return included in case error interrupts loop
    msglen = length(msg)-3; % compensate for conversion character and carriage return
    %fprintf(1,sp(ones(1,ceil(log10(m+1))+msglen)));
    fprintf(1,[bs(mod(0:2*(ceil(log10(i+1))+msglen)-1,2)+1) msg],i);
end