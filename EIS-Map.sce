clear
clc
funcprot(0)
Des=1;

y=messagebox(["Polytechnic Institute - Rio de Janenrio State University." "Labcor -Corrosion Laboratory." "EIS - Map Generator Software."], "EIS - Map Generator","message",["Start Program" "Finish Program" "Finish Scilab"],'modal');
if y==1 then

messagebox(["Be aware that the spreadsheet must be in accordance with the required standard" "Otherwise the program will acknowledge error and will not start"],"EIS - Map Generator","message",["I am aware"],'modal');

u = uigetfile(["*.*";"*.xls"],".", "Choose a file name", %t);

[fd,SST,Sheetnames,Sheetpos] = xls_open( u );
[Value,TextInd] = xls_read(fd,Sheetpos(1))
linha1 = Value(1,:);
Value(1,:) = [];

cb = SST(TextInd(1,:))
id_current = find(cb== 'Current (A)');
id_frequency = find(cb== 'Frequency (Hz)')
id_z = find(cb=='|Z| (ohms)')
id_phase = find(cb=='Phase of Z (deg)')
id_vezes= find(cb=='Point')
id_ddp= find(cb=='Potential (V)')


labels=["area"];
[ok,Area]=getvalue("Sample area(cm²)",labels,list("vec",1),["0.5"])


col_current = Value(:,id_current);
idc=1;   
numfreq=0
i=1
ind=1
vfreq=max(Value(:,id_frequency))
vnv=max(Value(:,id_vezes))+1;
for ind=1:vnv
    if(Value(ind,id_frequency)>0)
        if (min(vfreq)>=Value(ind,id_frequency)) then
          vfreq(idc)=Value(ind,id_frequency);
    
         if (idc>3) then 
             if(vfreq(idc)==vfreq(idc-1)) then
                 vfreq(idc)=[];
                 idc=idc-1;
         end
     end
     idc=idc+1;
   end
   end
end

numfreq=max(size(vfreq));


nexp=0
idc=1
ind=1

for ind=1:vnv
    if (vfreq(1)==Value(ind,id_frequency)) then
        vddp(idc)=Value(ind,id_ddp);
        idc=idc+1;
        nexp=nexp+1;
    end
end

idc=1

for ind=1:vnv
    if (Value(ind,id_frequency)~=0) then
        vZ(idc)=Value(ind,id_z)*Area;
        vPhas(idc)=abs(Value(ind,id_phase));
        idc=idc+1;
    end
end

I=1

for i=1:nexp
    for j=1:numfreq
        lmatZ(i,j)=log(abs(vZ(I)));
            matZ(i,j)=vZ(I);
            if (abs(vPhas(I))>90) then
                matPhas(i,j)=(abs(vPhas(I))-180);
            else
                matPhas(i,j)=abs(vPhas(I));
            end
            
            I=I+1
        end
    end

wg=messagebox(["Chosse a unit for the x axis"],"EIS - Map Generator","message",["Potencial" "Number of experiments"],'modal'); 

if wg==1 then
    a(:)=log10 (vfreq(:))
    b(:)=vddp(:)
    eix='E (V x SCE)'
    p1=[max(vddp)+0.2 4.5]
    else
a(:)=log10 (vfreq(:))
b(:)=1:nexp
p1=[nexp+10 4.5]
eix='Experiment number'

end


vminaf=min(matPhas)
vmaxaf=max(matPhas)
vminmq=min(matZ)
vmaxmq=max(matZ)
vminlmq=min(lmatZ)
vmaxlmq=max(lmatZ)

labels=["min"; "max"];
[ok,vminlmq,vmaxlmq]=getvalue("Plot range for log",labels,list("vec",1,"vec",1),["1.1";"5"]);


gg=messagebox(["What format do you wish to generate the images?"],"EIS - Map Generator","message",["2D image" "3D image", "Both"],'modal');



if gg==1 | gg==3 then
        

clf(0)
scf(0)
xset("colormap",jetcolormap(512))
grayplot(b,a,matPhas)
colorbar(0,90)
xtitle( 'Mapa de angulo de fase', eix, 'log (f /Hz)', 'angulo' , boxed = 1 )
title('Phase / Degree','position',p1)
title('Phase / Degree','fontsize',3)
filename='phase degree'
//xs2pdf(0,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the angle map"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(0,uiputfile(["*.jpg"[,"C:/"]])); 

clf(1)
scf(1)
xset("colormap",jetcolormap(512))
grayplot(b,a,matZ)
colorbar(vminmq,vmaxmq)
xtitle( 'Mapa de impedancia', eix, 'log (f / Hz)', 'angulo' , boxed = 1 )
title('|Z| / Ohm cm²','position',p1)
title('|Z| / Ohm cm²','fontsize',3);
filename='impedance'
//xs2pdf(1,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the impedance map"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(1,uiputfile(["*.jpg"[,"C:/"]])); 

clf(2)
scf(2)
xset("colormap",jetcolormap(512))
grayplot(b,a,lmatZ)
colorbar(vminlmq,vmaxlmq)
//colorbar(0,6)
//colorbar(1.1,5)
xtitle( 'Mapa de impedancia', eix, 'log (f / Hz)', 'angulo' , boxed = 1 )
title('log(|Z|/Ohm cm²)','position',p1)
title('log(|Z|/Ohm cm²)','fontsize',3)
filename='logimpedance'
//xs2pdf(2,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the impedance map in logarithmic"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(2,uiputfile(["*.jpg"[,"C:/"]]));

messagebox(["The generated images are avaiable"],"EIS - Map Generator","message",["OK"],'modal');
end

if gg==2 | gg==3 then
    
clf(3)
scf(3)
xset("colormap",jetcolormap(512))
plot3d1(b,a,matPhas)
colorbar(0,90)
xtitle( 'Mapa de angulo de fase', 'Experiment number', 'log (f) /Hz', 'angulo' , boxed = 1 )
title('Phase / Degree','position',[160.0 4.5])
title('Phase / Degree','fontsize',3)
filename='phase degree'
//xs2pdf(0,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the angle map"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(0,uiputfile(["*.jpg"[,"C:/"]])); 

clf(4)
scf(4)
xset("colormap",jetcolormap(512))
plot3d1(b,a,matZ)
colorbar(vminmq,100000)
xtitle( 'Mapa de impedancia', 'Experiment number', 'log (f) / Hz', '|Z| / Ohm' , boxed = 1 )
title('|Z| / Ohm','position',[160.0 4.5])
title('|Z| / Ohm','fontsize',3);
filename='impedance'
//xs2pdf(1,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the impedance map"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(1,uiputfile(["*.jpg"[,"C:/"]])); 

clf(5)
scf(5)
xset("colormap",jetcolormap(512))
plot3d1(b,a,lmatZ);
//plot3d1(b,a,lmatZ, flag=[-1,1,3], ebox=[min(b) max(b) min(a) max(a) 1.1 5])
//colorbar(1.1,5)
a=gca();
a.tight_limits = 'on';
a.data_bounds(1,3) = vminlmq;
a.data_bounds(2,3) = vmaxlmq;
colorbar(vminlmq,vmaxlmq)
xtitle( 'Mapa de impedancia', 'Experiment number', 'log (f) / Hz', 'log(|Z|/Ohm)' , boxed = 1 )
title('log(|Z|/Ohm)','position',[165.0 4.5])
title('log(|Z|/Ohm)','fontsize',3)
filename='logimpedance'
//xs2pdf(2,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the impedance map in logarithmic"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(2,uiputfile(["*.jpg"[,"C:/"]]));

messagebox(["The generated images are avaiable"],"EIS - Map Generator","message",["OK"],'modal');
end
else
    disp('fim do programa');
    messagebox(["Fim do programa"],"Gerador de Mapa de Impedância","message",["OK"],'modal');
end

