for A=(1:1:32);

    B=num2str(A);
    fileID=fopen(['\\nask.man.ac.uk\home$\Desktop\Year 3\Vibrations\Vibration Lab\Test_data_vibrations\Group 20 Vibration Tests Data - Edit\' B '.txt']);
    s=textscan(fileID, '%f %f %f' ,'HeaderLines',9, 'CollectOutput',1);
    fclose(fileID);
    C=s{1};
    
    a=C(:,1);
    b=C(:,2);
    c=C(:,3);
    F=a ;%Frequency in first column
    amp=sqrt(b.^2+c.^2);
    amp=log(amp);
    pha=(180/pi)*atan2(c,b);

    figure(A);
    subplot(2,1,1),plot(F,amp);
    title('Amplitude against Frequency')
    xlabel('Frequency (Hz)') % x-axis label
    ylabel('Amplitude (log)') % y-axis label
    subplot(2,1,2),plot(F,pha);
    title('Phase Angle against Frequency')
    xlabel('Frequency (Hz)') % x-axis label
    ylabel('Phase Angle') % y-axis label

end