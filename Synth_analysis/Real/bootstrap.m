clc
clear all
close all

T_lat=readtable("lt_exp.csv");
T_lat=table2array(T_lat);
T_lat_y=T_lat(:,randperm(size(T_lat,2),300));

xx_5=zeros(5,size(T_lat_y,2));
xx_10=zeros(10,size(T_lat_y,2));
xx_25=zeros(25,size(T_lat_y,2));
xx_50=zeros(50,size(T_lat_y,2));
xx_60=zeros(60,size(T_lat_y,2));
xx_100=zeros(100,size(T_lat_y,2));
xx_200=zeros(200,size(T_lat_y,2));

for ii=1:size(T_lat_y,2)
   xx_5(:,ii)=emprand(T_lat_y(:,ii),5,1); 
   xx_10(:,ii)=emprand(T_lat_y(:,ii),10,1);
   xx_25(:,ii)=emprand(T_lat_y(:,ii),25,1);
   xx_50(:,ii)=emprand(T_lat_y(:,ii),50,1); 
   xx_60(:,ii)=emprand(T_lat_y(:,ii),60,1); 
   xx_100(:,ii)=emprand(T_lat_y(:,ii),100,1); 
   xx_200(:,ii)=emprand(T_lat_y(:,ii),200,1); 
end

xx_5(xx_5<0)=0;
xx_10(xx_10<0)=0;
xx_25(xx_25<0)=0;
xx_50(xx_50<0)=0;
xx_60(xx_60<0)=0;
xx_100(xx_100<0)=0;
xx_200(xx_200<0)=0;

csvwrite('Temp/lt_exp_5.csv',xx_5)
csvwrite('Temp/lt_exp_10.csv',xx_10)
csvwrite('Temp/lt_exp_25.csv',xx_25)
csvwrite('Temp/lt_exp_50.csv',xx_50)
csvwrite('Temp/lt_exp_60.csv',xx_60)
csvwrite('Temp/lt_exp_100.csv',xx_100)
csvwrite('Temp/lt_exp_200.csv',xx_200)
