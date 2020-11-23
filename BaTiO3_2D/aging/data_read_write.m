clear;clc;close all;
%read data
time1=[0,1,2,4,8,16,24,48,96,2700];
time2=[8,16,24,48,96];
n1=length(time1);
n2=length(time2);
folder='D:\BTMO\BTMNO\';
hyloop=cell(n1,4);
for i=1:n1
    for j=1:2
        fileID=fopen([folder,num2str(time1(i)),'h-',num2str(j),'.txt']);
        data=textscan(fileID,'%f %f %f %f','HeaderLines',54,'Delimiter','\t');
        for k=3:4
            hyloop{i,j}=[data{1,3},data{1,4}];
        end
        fclose(fileID);
    end
end
for i=1:n2
    for j=3:4
        fileID=fopen([folder,num2str(time2(i)),'h2-',num2str(j-2),'.txt']);
        data=textscan(fileID,'%f %f %f %f','HeaderLines',54,'Delimiter','\t');
        for k=3:4
            hyloop{i,j}=[data{1,3},data{1,4}];
        end
        fclose(fileID);
    end
end

% data process
% read the sample thickness and wether the unit of field is kV/vm or V
for i=1:n1
    for j=1:2
        fileID=fopen([folder,num2str(time1(i)),'h-',num2str(j),'.txt']);
        filehead=textscan(fileID,'%q %q %q %q',27,'HeaderLines',25,'Delimiter','\t');
        thickness=cell2mat(filehead{1,2}(1));
        sample_thickness=str2double(thickness(1:4));
        unit=cell2mat(filehead{1,3}(22));
        fclose(fileID);
        if unit=="Field (kV/cm)"
            writedata=hyloop{i,j};
            dlmwrite([folder,num2str(time1(i)),'h-',num2str(j),'.csv'],writedata);
        else
            hyloop{i,j}(:,1)=hyloop{i,j}(:,1)*0.1/sample_thickness;
            writedata=hyloop{i,j};
            dlmwrite([folder,num2str(time1(i)),'h-',num2str(j),'.csv'],writedata);
        end
    end
end

for i=1:n2
    for j=3:4
        hyloop{i,j}(:,1)=hyloop{i,j}(:,1);
        writedata=hyloop{i,j};
        dlmwrite([folder,num2str(time2(i)),'h2-',num2str(j-2),'.csv'],writedata);
    end
end
