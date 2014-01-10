clear
clc
funcprot(0)
Des=1;

y=messagebox(["Polytechnic Institute - Rio de Janeiro State University" "Labcor - Corrosion Laboratory" "EIS - Map Generator Software"], "EIS - Map Generator","message",["Start" "Finish" "Finish Scilab"],'modal');
if y==1 then

messagebox(["Be aware that the spreadsheet is in accordance with the required standard" "Otherwise the program will not run"],"EIS - Map Generator","message",["I am aware"],'modal');

u = uigetfile(["*.*";"*.xls"],".", "Choose a file name", %t);

[fd,SST,Sheetnames,Sheetpos] = xls_open( u );
[Value,TextInd] = xls_read(fd,Sheetpos(1))
linha1 = Value(1,:);
Value(1,:) = [];

cb = SST(TextInd(1,:))
id_frequency = find(cb== 'Frequency (Hz)')
id_z = find(cb=='|Z| (ohms)')
id_phase = find(cb=='Phase of Z (deg)')
id_vezes= find(cb=='Point')
id_ddp= find(cb=='Potential (V)')

//////////////////////// Pre-test ///////////////////////////

erros = [];
ind=1
vnv=max(Value(:,id_vezes))+1;
for ind=1:vnv
    if(Value(ind,id_frequency)>0) then
        // Verifying if |Z| > 0
        if (Value(ind,id_z)<=0) then
            erros = [erros "The value of |Z| in line "+string(ind+1)+" is invalid."]
        end
    end
end
if length(erros)>0 then
    erros = ["The following errors were found on the sheet:" erros];
    messagebox(erros,"EIS - Map Generator","message",["OK"],'modal');
end

//////////////////////////////////////////////////////////////

labels=["area"];
[ok,Area]=getvalue("Sample area (cm²)",labels,list("vec",1),["1.0"])

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
        lmatZ(i,j)=log10(abs(vZ(I)));
        matZ(i,j)=vZ(I);
        matPhasTeste(i,j) = vPhas(I);
           if (abs(vPhas(I))>90) then
               matPhas(i,j)=(abs(vPhas(I))-180);
           else
               matPhas(i,j)=abs(vPhas(I));
           end
           I=I+1
    end
end


//////////// Smoothing the surface (remove the noise) /////////////
smoothing=messagebox(["Smoothing of spikes:"],"EIS - Map Generator","message",["No" "Yes"],'modal'); 

if smoothing == 2 then
    numeroSuavizacoes = 2
    toleranciaZ = 2
    toleranciaFase = 2
    
    for k=1:numeroSuavizacoes
         
        // Log Módulo
        for i=nexp:-1:1
            diferencaAnterior = 100000000
            zAnterior = lmatZ(nexp,numfreq)
            for j=numfreq:-1:1
                if (abs(lmatZ(i,j)-zAnterior)>abs(toleranciaZ*diferencaAnterior)) then
                    // Vértice sup-esq
                    if(i==1 & j==1) then
                        lmatZ(i,j) = (lmatZ(i,j+1)+lmatZ(i+1,j)+lmatZ(i+1,j+1))/3
                    // Vértice sup-dir
                    elseif(i==1 & j==numfreq) then
                        lmatZ(i,j) = (lmatZ(i,j-1)+lmatZ(i+1,j)+lmatZ(i+1,j-1))/3
                    // Vértice inf-esq
                    elseif(i==nexp & j==1) then
                        lmatZ(i,j) = (lmatZ(i,j+1)+lmatZ(i-1,j)+lmatZ(i-1,j+1))/3
                    // Vértice inf-dir
                    elseif(i==nexp & j==numfreq) then
                        lmatZ(i,j) = (lmatZ(i,j-1)+lmatZ(i-1,j)+lmatZ(i-1,j-1))/3
                    // Borda superior
                    elseif(i==1) then
                        lmatZ(i,j) = (lmatZ(i,j-1)+lmatZ(i,j+1)+lmatZ(i+1,j)+lmatZ(i+1,j+1)+lmatZ(i+1,j-1))/5
                    //Borda inferior
                    elseif(i==nexp) then
                        lmatZ(i,j) = (lmatZ(i,j-1)+lmatZ(i,j+1)+lmatZ(i-1,j)+lmatZ(i-1,j-1)+lmatZ(i-1,j+1))/5
                    //Borda Esquerda
                    elseif(j==1) then
                        lmatZ(i,j) = (lmatZ(i,j+1)+lmatZ(i-1,j)+lmatZ(i+1,j)+lmatZ(i+1,j+1)+lmatZ(i-1,j+1))/5
                    //Borda Direita
                    elseif(j==numfreq) then                
                        lmatZ(i,j) = (lmatZ(i,j-1)+lmatZ(i-1,j)+lmatZ(i+1,j)+lmatZ(i-1,j-1)+lmatZ(i+1,j-1))/5
                    // Massa interna
                    elseif (j<numfreq & i<nexp & j>1 & i>1) then
                        lmatZ(i,j) = (lmatZ(i,j-1)+lmatZ(i,j+1)+lmatZ(i-1,j)+lmatZ(i+1,j)+lmatZ(i+1,j+1)+lmatZ(i-1,j-1)+lmatZ(i-1,j+1)+lmatZ(i+1,j-1))/8
                    end
                end
                diferencaAnterior = lmatZ(i,j)-zAnterior
                zAnterior = lmatZ(i,j)
            end
        end
        
        // Phase
        for i=nexp:-1:1
            diferencaAnterior = 100000000
            faseAnterior = matPhas(nexp,numfreq)
            for j=numfreq:-1:1
                if (abs(matPhas(i,j)-faseAnterior)>abs(toleranciaFase*diferencaAnterior)) then
                    // Vértice sup-esq
                    if(i==1 & j==1) then
                        matPhas(i,j) = (matPhas(i,j+1)+matPhas(i+1,j)+matPhas(i+1,j+1))/3
                    // Vértice sup-dir
                    elseif(i==1 & j==numfreq) then
                        matPhas(i,j) = (matPhas(i,j-1)+matPhas(i+1,j)+matPhas(i+1,j-1))/3
                    // Vértice inf-esq
                    elseif(i==nexp & j==1) then
                        matPhas(i,j) = (matPhas(i,j+1)+matPhas(i-1,j)+matPhas(i-1,j+1))/3
                    // Vértice inf-dir
                    elseif(i==nexp & j==numfreq) then
                        matPhas(i,j) = (matPhas(i,j-1)+matPhas(i-1,j)+matPhas(i-1,j-1))/3
                    // Borda superior
                    elseif(i==1) then
                        matPhas(i,j) = (matPhas(i,j-1)+matPhas(i,j+1)+matPhas(i+1,j)+matPhas(i+1,j+1)+matPhas(i+1,j-1))/5
                    //Inferior boundary
                    elseif(i==nexp) then
                        matPhas(i,j) = (matPhas(i,j-1)+matPhas(i,j+1)+matPhas(i-1,j)+matPhas(i-1,j-1)+matPhas(i-1,j+1))/5
                    //Left boundary
                    elseif(j==1) then
                        matPhas(i,j) = (matPhas(i,j+1)+matPhas(i-1,j)+matPhas(i+1,j)+matPhas(i+1,j+1)+matPhas(i-1,j+1))/5
                    //Rigth boundary
                    elseif(j==numfreq) then                
                        matPhas(i,j) = (matPhas(i,j-1)+matPhas(i-1,j)+matPhas(i+1,j)+matPhas(i-1,j-1)+matPhas(i+1,j-1))/5
                    // Bulk
                    elseif (j<numfreq & i<nexp & j>1 & i>1) then
                        matPhas(i,j) = (matPhas(i,j-1)+matPhas(i,j+1)+matPhas(i-1,j)+matPhas(i+1,j)+matPhas(i+1,j+1)+matPhas(i-1,j-1)+matPhas(i-1,j+1)+matPhas(i+1,j-1))/8
                    end
                end
                diferencaAnterior = matPhas(i,j)-faseAnterior
                faseAnterior = matPhas(i,j)
            end
        end
    
    end
end
/////////////////////////////////////////////////////////////////////////////////////
//////////////////// Calculating real and imaginary components //////////////////////
matFreq = []
for i=1:nexp
    matFreq = [matFreq vfreq];
end
matFreq = matFreq';

matRealZ = [];
matImagZ = [];
matPhasRad = [];
for i=1:nexp
    for j=1:numfreq
        matPhasRad(i,j)=matPhas(i,j)*%pi/180;
        complexo = matZ(i,j)*exp(matPhasRad(i,j)*(%i));
        matRealZ(i,j) = real(complexo);
        matImagZ(i,j) = imag(complexo);
    end
end
matrizNexp = []
matPot = []
for i=1:nexp
    for j=1:numfreq
        matrizNexp(i,j)=i;
        matPot(i,j) = vddp(i);
    end
end

////////////////////////////////////////////////////////////////////////////////////

wg=messagebox(["Choose a unity of x axis"],"EIS - Map Generator","message",["Potential" "Experiment number"],'modal'); 

if wg==1 then
    a(:)=log10 (vfreq(:))
    b(:)=1:nexp
    //b(:)=vddp(:)
    eix='E (V x SCE)'
    pot = vddp(:);
    pot_string=[];
    for i=1:20:length(pot)
        //arredondando para 3 casas decimais
        potencial = round((pot(i))*10^3)/10^3;
        pot_string=[pot_string string(potencial)];
    end
    vet_nexp = 1:20:nexp
else
    a(:)=log10 (vfreq(:))
    b(:)=1:nexp
    eix='Experiment number'
end


vminaf=min(matPhas)
vmaxaf=max(matPhas)
vminmq=min(matZ)
vmaxmq=max(matZ)
vminlmq=min(lmatZ)
vmaxlmq=max(lmatZ)

labels=["phase min (deg)"; "phase max (deg)"; "log min (Ohm.cm²)"; "log max (Ohm.cm²)"; "3D color mode (-1 for no line)"];
[ok,phasemin, phasemax, vminlmq,vmaxlmq,colormode]=...
getvalue("Plot options", labels, list("vec",1,"vec",1,"vec",1,"vec",1,"vec",1),[string(vminaf);string(vmaxaf);string(vminlmq);string(vmaxlmq);"2"]);

// Nyquist diagrams
ng=messagebox(["Do you wish to generate the Nyquist diagrams?"],"EIS - Map Generator","message",["Yes" "No"],'modal');

gg=messagebox(["What format do you wish to generate the images?"],"EIS - Map Generator","message",["2D image" "3D image", "Both"],'modal');


if gg==1 | gg==3 then
        

clf(0)
scf(0)
xset("colormap",jetcolormap(512))
colorbar(phasemin,phasemax)
cbar = gce();
cbar.parent.title.text = "Phase / Deg";
cbar.parent.title.fill_mode = "on"
cbar.parent.title.font_size = 3
grayplot(b,a,matPhas)
xtitle( 'Phase Map', eix, 'log (f /Hz)', 'Phase' , boxed = 1)
filename='phase deg'
eixos=get("current_axes")
eixos.title.font_size = 3;
eixos.x_label.font_size = 3;
eixos.y_label.font_size = 3;
eixos.z_label.font_size = 3;

if wg==1 then
    eixos.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
end

//xs2pdf(0,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the phase map"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(0,uiputfile(["*.jpg"[,"C:/"]])); 

clf(1)
scf(1)
xset("colormap",jetcolormap(512))
colorbar(vminmq,vmaxmq)
cbar = gce();
cbar.parent.title.text = "|Z| / Ohm cm²";
cbar.parent.title.fill_mode = "on"
cbar.parent.title.font_size = 3
grayplot(b,a,matZ)
xtitle( 'Impedance Map', eix, 'log (f / Hz)', 'Phase' , boxed = 1 )
filename='impedance'
eixos=get("current_axes")
eixos.title.font_size = 3;
eixos.x_label.font_size = 3;
eixos.y_label.font_size = 3;
eixos.z_label.font_size = 3;

if wg==1 then
    eixos.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
end

//xs2pdf(1,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the impedance map"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(1,uiputfile(["*.jpg"[,"C:/"]])); 

clf(2)
scf(2)
xset("colormap",jetcolormap(512))
colorbar(vminlmq,vmaxlmq)
cbar = gce();
cbar.parent.title.text = "log(|Z|/Ohm cm²)";
cbar.parent.title.fill_mode = "on"
cbar.parent.title.font_size = 3
grayplot(b,a,lmatZ)
//colorbar(0,6)
//colorbar(1.1,5)
xtitle( 'Impedance Map', eix, 'log (f / Hz)', 'Phase' , boxed = 1 )
filename='logimpedance'
eixos=get("current_axes")
eixos.title.font_size = 3;
eixos.x_label.font_size = 3;
eixos.y_label.font_size = 3;
eixos.z_label.font_size = 3;

if wg==1 then
    eixos.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
end

//xs2pdf(2,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the impedance map in logarithmic"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(2,uiputfile(["*.jpg"[,"C:/"]]));

if ng == 1 then
    clf(6)
    scf(6)
    
    cores = jetcolormap(nexp);
    
    if wg==1 then
        txt='Potential (V)'
        maxVal=max(vddp)
        minVal=min(vddp)
    else 
        txt='Experiment Number'
        maxVal=nexp
        minVal=1
    end
    
    xset("colormap",cores)
    colorbar(minVal,maxVal);
    cbar = gce();
    cbar.parent.title.text = txt;
    cbar.parent.title.fill_mode = "on"
    cbar.parent.title.font_size = 3

    for i=1:nexp
        plot(matRealZ(i,:),matImagZ(i,:),matPot(:,1), "color",cores(i,:));
    end
    
    xtitle( 'Nyquist Diagram', 'Real (Ohm.cm²)', '-Imag (Ohm.cm²)' , boxed = 1)
    filename='nyquist'
    eixos=get("current_axes")
    eixos.title.font_size = 3;
    eixos.x_label.font_size = 3;
    eixos.y_label.font_size = 3;
    eixos.z_label.font_size = 3;
end

messagebox(["Images are available"],"EIS - Map Generator","message",["OK"],'modal');
end

if gg==2 | gg==3 then
    
clf(3)

scf(3)
xset("colormap",jetcolormap(512))
colorbar(vminaf,vmaxaf)
cbar = gce();
cbar.parent.title.text = "Phase / deg";
cbar.parent.title.fill_mode = "on"
cbar.parent.title.font_size = 3

plot3d1(b,a,matPhas)

e=gce();
e.hiddencolor=-1;
e.color_mode=colormode;
xtitle( 'Phase Map', eix, 'log (f) /Hz', 'Phase' , boxed = 1 )
filename='phase degree'
eixos=get("current_axes")
eixos.title.font_size = 3;
eixos.x_label.font_size = 3;
eixos.y_label.font_size = 3;
eixos.z_label.font_size = 3;
eixos.tight_limits = 'on';
eixos.data_bounds(1,3) = phasemin;
eixos.data_bounds(2,3) = phasemax;

if wg==1 then
    eixos.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
end

//xs2pdf(0,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the angle map"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(0,uiputfile(["*.jpg"[,"C:/"]])); 

clf(4)
scf(4)
xset("colormap",jetcolormap(512))
colorbar(vminmq,vmaxmq)
cbar = gce();
cbar.parent.title.text = "|Z| / Ohm.cm²";
cbar.parent.title.fill_mode = "on"
cbar.parent.title.font_size = 3
plot3d1(b,a,matZ)
e=gce();
e.hiddencolor=-1;
e.color_mode=colormode;
xtitle( 'Impedance Map', eix, 'log (f) / Hz', '|Z| / Ohm.cm²' , boxed = 1 )
filename='impedance'
eixos=get("current_axes")
eixos.title.font_size = 3;
eixos.x_label.font_size = 3;
eixos.y_label.font_size = 3;
eixos.z_label.font_size = 3;

if wg==1 then
    eixos.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
end

//xs2pdf(1,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the impedance map"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(1,uiputfile(["*.jpg"[,"C:/"]])); 

clf(5)
scf(5)
xset("colormap",jetcolormap(512))

colorbar(min(lmatZ),max(lmatZ))
cbar = gce();
cbar.parent.title.text = "log(|Z|/Ohm.cm²)";
cbar.parent.title.fill_mode = "on"
cbar.parent.title.font_size = 3
plot3d1(b,a,lmatZ);
e=gce();
e.hiddencolor=-1;
e.color_mode=colormode;
//plot3d1(b,a,lmatZ, flag=[-1,1,3], ebox=[min(b) max(b) min(a) max(a) 1.1 5])
//colorbar(1.1,5)
ax=gca();
ax.tight_limits = 'on';
ax.data_bounds(1,3) = vminlmq;
ax.data_bounds(2,3) = vmaxlmq;
xtitle( 'Impedance Map', eix, 'log (f) / Hz', 'log(|Z|/Ohm.cm²)' , boxed = 1 )
filename='logimpedance'
ax.title.font_size = 3;
ax.x_label.font_size = 3;
ax.y_label.font_size = 3;
ax.z_label.font_size = 3;

if wg==1 then
    ax.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
end

//xs2pdf(2,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the impedance map in logarithmic"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(2,uiputfile(["*.jpg"[,"C:/"]]));

/////////////////// 3D Nyquist Diagrams ////////////////////
if ng == 1 then
    clf(7)
    scf(7)
    xset("colormap",jetcolormap(64)); 
    colorbar(0,max(matImagZ));
    cbar = gce();
    cbar.parent.title.text = "Imaginary Component";
    cbar.parent.title.fill_mode = "on"
    cbar.parent.title.font_size = 3
    if wg==1 then
        eix='E (V x SCE)'
        surf(matPot,matRealZ,matImagZ);
    else 
        eix='Experiment Number'
        surf(matrizNexp,matRealZ,matImagZ);
    end
    e=gce();
    e.hiddencolor=-1;
    e.color_mode=colormode;
    xtitle( 'Nyquist Diagram', eix, 'Real', 'Imaginary' , boxed = 1)
    filename='nyquist'
    eixos=get("current_axes")
    eixos.title.font_size = 3;
    eixos.x_label.font_size = 3;
    eixos.y_label.font_size = 3;
    eixos.z_label.font_size = 3;
end

messagebox(["Images are available"],"EIS - Map Generator","message",["OK"],'modal');
end
else
    disp('End');
    messagebox(["End"],"EIS - Map Generator","message",["OK"],'modal');
end

// Close all files
file('close',file() )
