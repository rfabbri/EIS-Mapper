clear
clc
funcprot(0)
Des=1;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Funções para Kramers-Kronig
function f = imag_kk(f, vreal)
    calc_imag = [];
    for i = 1:length(f) 
        vfreqaux = f;
        vrealaux = vreal;
        vfreqaux(:,i) = [];
        vrealaux(:,i) = [];
        integrando = ((vrealaux - vreal(i))./(vfreqaux.^2 - f(i)^2));  
        integralKK = intsplin(vfreqaux,integrando);
        calc_imag(i) = integralKK;
    end
    f = ((2*f./%pi)'.*calc_imag)';
endfunction
function f = real_kk(f, vimag, idPot)
    calc_real = [];
    for i = 1:length(f) 
        vfreqaux = f;
        vimagaux = vimag;
        vfreqaux(:,i) = [];
        vimagaux(:,i) = [];
        integrando = (((vfreqaux.*vimagaux - f(i).*vimag(i))./(vfreqaux.^2 - f(i)^2)));
        integralKK = intsplin(vfreqaux,integrando);
        calc_real(i) = integralKK;
    end
    Zreal_inf = matRealZ(idPot,1);
    f = (((2 ./ %pi)'.*calc_real)+Zreal_inf);
endfunction
function y = inverte(x)
    n = length(x);
    for i=1:n
        y(i) = x(n+1-i);
    end
    y=y';
endfunction
function [p]=polyfit(x, y, grauPolinomio)
    // return coefficient vector or poly if fourth string argument given
    x = x(:); 
    y = y(:);
    n = length(x);
    if length(y) <> n, error('x and y must have same length'), end
    v = ones(n,grauPolinomio+1);
    for i=2:grauPolinomio+1
        v(:,i) = x.*v(:,i-1);
    end
    p = (v\y)';
endfunction
function f=extrapolacao(f, vReal, vImag)
// min
    
    grauPol = 10;
    x = f(1:10);
    y = vReal(1:10);
    [coefs] = polyfit(x, y, grauPol);
    polinomio = poly(coefs, 'x', 'coeff');
    
    x_extra = [f(1:10).*0.01 f(1:10).*0.1];
    curva_extra = horner(polinomio,x_extra);
    
    f_extra = [x_extra f];
    vReal_extra = [curva_extra vReal]

    x = f(1:10);
    y = vImag(1:10);
    [coefs] = polyfit(x, y, grauPol);
    polinomio = poly(coefs, 'x', 'coeff');
    
    curva_extra = horner(polinomio,x_extra);
    vImag_extra = [curva_extra vImag];
    
    f =[f_extra; vReal_extra; vImag_extra];
endfunction
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
id_current= find(cb=='Current (A)')

//////////////////////// Pre-test ///////////////////////////

erros = [];
ind=1

s = size(Value);
vnv = s(1);          // Número de pontos
//vnv=max(Value(:,id_vezes))+1;
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
        vCurrent(idc)=Value(ind,id_current);
        lvCurrent(idc)=log10(abs(Value(ind,id_current)*1000));
        idc=idc+1;
        nexp=nexp+1;
    end
end

idc=1

for ind=1:vnv
    if (Value(ind,id_frequency)~=0) then
        vZ(idc)=Value(ind,id_z)*Area;
        vPhas(idc)=-(Value(ind,id_phase));
        idc=idc+1;
    end
end

I=1

for i=1:nexp
    for j=1:numfreq
        lmatZ(i,j)=log10(abs(vZ(I)));
        matZ(i,j)=vZ(I);
        matPhasTeste(i,j) = vPhas(I);
        matPhas(i,j)=(vPhas(I));
//           if (abs(vPhas(I))>90) then
//               matPhas(i,j)=(abs(vPhas(I))-180);
//           else
//               matPhas(i,j)=abs(vPhas(I));
//           end
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

labels=["phase min (deg)"; "phase max (deg)"; "log min (Ohm.cm²)"; "log max (Ohm.cm²)"; "3D color mode (-1 for no line, 0 for surface mesh, 1 for surface with lines)"];
[ok,phasemin, phasemax, vminlmq,vmaxlmq,colormode]=...
getvalue("Plot options", labels, list("vec",1,"vec",1,"vec",1,"vec",1,"vec",1),[string(vminaf);string(vmaxaf);string(vminlmq);string(vmaxlmq);"1"]);

// Validating Values
if vmaxaf < (phasemax-1) then
    messagebox([string(phasemax)+" degrees is higher than the maximum phase."],"EIS - Map Generator","message",["OK"],'modal');
end
if vminaf > (phasemin+1) then
    messagebox([string(phasemin)+" degrees is lower than the minimum phase."],"EIS - Map Generator","message",["OK"],'modal');
end

// Nyquist diagrams
ng=messagebox(["Do you wish to generate the Nyquist diagrams?"],"EIS - Map Generator","message",["Yes" "No"],'modal');

// Polarization superposition
pg=messagebox(["Do you wish to plot the DC values over the maps?"],"EIS - Map Generator","message",["Yes" "No"],'modal');

// Kramers-Kronig transform
kk_validate=messagebox(["Do you wish to validate the data using Kramers-Kronig transform?"],"EIS - Map Generator","message",["Yes" "No"],'modal');

// Verificando se tem reversão de potencial
indiceReversao = [];
for i=2:nexp
    if vddp(i) < vddp(i-1) then
        // tem reversão!
        indiceReversao = i-1;
        break;
    end
end

el = 2;
if indiceReversao ~= [] then
    // Linhas de potencial
    el=messagebox(["Potencial reversion detected. Do you wish to plot equi-potential lines?"],"EIS - Map Generator","message",["Yes" "No"],'modal');
end    


gg=messagebox(["What format do you wish to generate the images?"],"EIS - Map Generator","message",["2D image" "3D image", "Both"],'modal');


if gg==1 | gg==3 then
        

clf(0)
scf(0)
xset("colormap",jetcolormap(512))
colorbar(vminaf,vmaxaf)
cbar = gce();
cbar.parent.title.text = "Phase / Deg";
cbar.parent.title.fill_mode = "off"
cbar.parent.title.font_size = 3
grayplot(b,a,matPhas)
xtitle( 'Phase Map', eix, 'log (f /Hz)', 'Phase' , boxed = 0)
filename='phase deg'
eixos=get("current_axes")
eixos.title.font_size = 3;
eixos.x_label.font_size = 3;
eixos.y_label.font_size = 3;
eixos.z_label.font_size = 3;
if wg==1 then
    eixos.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
end
///////////////////////////////////////////////////////////////////////////
// Polarization Curve Superposition
if pg == 1 then
    maxFreq = max(a);
    maxCurrent = max(vCurrent);
    for i = 1:length(vCurrent);
        auxCurrent(i) = vCurrent(i)*maxFreq/maxCurrent;
    end
    a2=newaxes();
    a2.y_location = 'right'; 
    plot2d(b, abs(vCurrent), logflag =  "nl", leg="DC values", style=[color("black")]);
    curva= a2.children(1).children(1)
    curva.thickness=3;
    
    if  el == 1 then
        paresReversao = [];
        if indiceReversao ~= [] then
            for i=indiceReversao:nexp
                paresReversao = [paresReversao [i;(2*indiceReversao-i)] ]
            end
            ss=size(paresReversao);
            s = ss(2);
            for i = 2:2:s
                x = [paresReversao(1,i), paresReversao(2,i)];
                y = [vCurrent(paresReversao(1,i)), abs(vCurrent(paresReversao(2,i)))];
                plot(x, y,'r');
                eixos=get("current_axes");
                curva= eixos.children(1).children(1)
                curva.thickness=1;
                eixos.axes_visible = ["off","off","on"];
            end
        end
    end
    
    xtitle( '', '', '$\bold{\mathrm{j(A.cm^{-2})}}$', '' , boxed = 0);
    a2.axes_visible = ["off","on","on"];
    a2.filled = "off";
    a2.y_label.font_size = 3;
    a2.thickness=3
    colorbar(vminaf,vmaxaf)
    cbar = gce();
    cbar.parent.title.text = "Phase / Deg";
    cbar.parent.title.fill_mode = "off"
    cbar.parent.title.font_size = 3

end
///////////////////////////////////////////////////////////////////////////


//xs2pdf(0,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the phase map"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(0,uiputfile(["*.jpg"[,"C:/"]])); 

clf(1)
scf(1)
xset("colormap",jetcolormap(512))
colorbar(vminmq,vmaxmq)
cbar = gce();
cbar.parent.title.text = "|Z| / Ohm cm²";
cbar.parent.title.fill_mode = "off"
cbar.parent.title.font_size = 3
grayplot(b,a,matZ)
xtitle( 'Bode Plot', eix, 'log (f / Hz)', 'Phase' , boxed = 0 )
filename='impedance'
eixos=get("current_axes")
eixos.title.font_size = 3;
eixos.x_label.font_size = 3;
eixos.y_label.font_size = 3;
eixos.z_label.font_size = 3;
if wg==1 then
    eixos.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
end
///////////////////////////////////////////////////////////////////////////
// Polarization Curve Superposition
if pg == 1 then
    maxFreq = max(a);
    maxCurrent = max(vCurrent);
    for i = 1:length(vCurrent);
        auxCurrent(i) = vCurrent(i)*maxFreq/maxCurrent;
    end
    a2=newaxes();
    a2.y_location = 'right'; 
    plot2d(b, abs(vCurrent), logflag =  "nl", leg="DC values", style=[color("black")]);
    curva= a2.children(1).children(1)
    curva.thickness=3;
    if  el == 1 then
        paresReversao = [];
        if indiceReversao ~= [] then
        for i=indiceReversao:nexp
            paresReversao = [paresReversao [i;(2*indiceReversao-i)] ]
        end
        ss=size(paresReversao);
        s = ss(2);
        for i = 2:2:s
            x = [paresReversao(1,i), paresReversao(2,i)];
            y = [vCurrent(paresReversao(1,i)), abs(vCurrent(paresReversao(2,i)))];
            plot(x, y,'r');
            eixos=get("current_axes");
            curva= eixos.children(1).children(1)
            curva.thickness=1;
            eixos.axes_visible = ["off","off","off"];
        end
        end
    end
    
    xtitle( '', '', '$\bold{\mathrm{j(A.cm^{-2})}}$', '' , boxed = 0);
    a2.axes_visible = ["off","on","on"];
    a2.filled = "off";
    a2.y_label.font_size = 3;
    a2.thickness=3
    colorbar(vminlmq,vmaxlmq)
    cbar = gce();
    cbar.parent.title.text = "|Z| / Ohm cm²";
    cbar.parent.title.fill_mode = "off"
    cbar.parent.title.font_size = 3
end
///////////////////////////////////////////////////////////////////////////

//xs2pdf(1,filename)

//messagebox(["Choose the directory and name in which you want to save the image of the impedance map"],"EIS - Map Generator","message",["I agree"],'modal');
//xs2jpg(1,uiputfile(["*.jpg"[,"C:/"]])); 

clf(2)
scf(2)
xset("colormap",jetcolormap(512))
colorbar(vminlmq,vmaxlmq)
cbar = gce();
cbar.parent.title.text = "log(|Z|/Ohm cm²)";
cbar.parent.title.fill_mode = "off"
cbar.parent.title.font_size = 3
grayplot(b,a,lmatZ)
//colorbar(0,6)
//colorbar(1.1,5)
xtitle( 'Bode Plot', eix, 'log (f / Hz)', 'Phase' , boxed = 0 )
filename='logimpedance'
eixos=get("current_axes")
eixos.title.font_size = 3;
eixos.x_label.font_size = 3;
eixos.y_label.font_size = 3;
eixos.z_label.font_size = 3;
if wg==1 then
    eixos.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
end
///////////////////////////////////////////////////////////////////////////
// Polarization Curve Superposition
if pg == 1 then
    maxFreq = max(a);
    maxCurrent = max(vCurrent);
    for i = 1:length(vCurrent);
        auxCurrent(i) = vCurrent(i)*maxFreq/maxCurrent;
    end
    a2=newaxes();
    a2.y_location = 'right'; 
    plot2d(b, abs(vCurrent), logflag =  "nl", leg="DC values", style=[color("black")]);
    curva= a2.children(1).children(1)
    curva.thickness=3;
    if  el == 1 then
        paresReversao = [];
        if indiceReversao ~= [] then
        for i=indiceReversao:nexp
            paresReversao = [paresReversao [i;(2*indiceReversao-i)] ]
        end
        ss=size(paresReversao);
        s = ss(2);
        for i = 2:2:s
            x = [paresReversao(1,i), paresReversao(2,i)];
            y = [vCurrent(paresReversao(1,i)), abs(vCurrent(paresReversao(2,i)))];
            plot(x, y,'r');
            eixos=get("current_axes");
            curva= eixos.children(1).children(1)
            curva.thickness=1;
            eixos.axes_visible = ["off","off","off"];
        end
        end
    end
    
    xtitle( '', '', '$\bold{\mathrm{j(A.cm^{-2})}}$', '' , boxed = 0);
    a2.axes_visible = ["off","on","on"];
    a2.filled = "off";
    a2.y_label.font_size = 3;
    a2.thickness=3
    colorbar(vminlmq,vmaxlmq)
    cbar = gce();
    cbar.parent.title.text = "log(|Z|/Ohm cm²)";
    cbar.parent.title.fill_mode = "off"
    cbar.parent.title.font_size = 3
    
end
///////////////////////////////////////////////////////////////////////////

clf(8)
scf(8)

if wg==1 then
    plot2d(vddp, abs(vCurrent), logflag =  "nl");
    xtitle( 'DC values', eix, '$\mathrm{j(A.cm^{-2})}$', '' , boxed = 0)
    filename='dc values'
    eixos=get("current_axes")
    eixos.title.font_size = 3;
    eixos.x_label.font_size = 3;
    eixos.y_label.font_size = 3;
    eixos.z_label.font_size = 3;
//    eixos.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
else
    
    plot2d(b, abs(vCurrent), logflag =  "nl");
    xtitle( 'DC values', eix, '$\mathrm{j(A.cm^{-2})}$', '' , boxed = 0)
    filename='dc values'
    eixos=get("current_axes")
 
   
    if  el == 1 then
        paresReversao = [];
        if indiceReversao ~= [] then
        for i=indiceReversao:nexp
            paresReversao = [paresReversao [i;(2*indiceReversao-i)] ]
        end
        ss=size(paresReversao);
        s = ss(2);
        for i = 2:2:s
            x = [paresReversao(1,i), paresReversao(2,i)];
            y = [vCurrent(paresReversao(1,i)), abs(vCurrent(paresReversao(2,i)))];
            plot(x, y,'r');
            eixos=get("current_axes");
            curva= eixos.children(1).children(1)
            curva.thickness=1;
//            eixos.axes_visible = ["off","off","off"];
        end
        end
    end
    
    eixos.title.font_size = 3;
    eixos.x_label.font_size = 3;
    eixos.y_label.font_size = 3;
    eixos.z_label.font_size = 3;
//    eixos.thickness=3
    
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
    cbar.parent.title.fill_mode = "off"
    cbar.parent.title.font_size = 3

    for i=1:nexp
        plot(matRealZ(i,:),matImagZ(i,:),matPot(:,1), "color",cores(i,:));
    end
    
    xtitle( 'Nyquist Diagram', 'Real (Ohm.cm²)', '-Imag (Ohm.cm²)' , boxed = 0)
    filename='nyquist'
    eixos=get("current_axes")
    eixos.title.font_size = 3;
    eixos.x_label.font_size = 3;
    eixos.y_label.font_size = 3;
    eixos.z_label.font_size = 3;
    
end

/////////////////// fim plots 2d //////////////////////////////////
end

if gg==2 | gg==3 then
    
clf(3)

scf(3)
xset("colormap",jetcolormap(512))
colorbar(vminaf,vmaxaf)
cbar = gce();
cbar.parent.title.text = "Phase / deg";
cbar.parent.title.fill_mode = "off"
cbar.parent.title.font_size = 3

plot3d1(b,a,matPhas);

e=gce();
e.hiddencolor=-1;
e.color_mode=colormode;
xtitle( 'Phase Map', eix, 'log (f) /Hz', 'Phase' , boxed = 0 )
filename='phase degree'
eixos=get("current_axes")
eixos.title.font_size = 3;
eixos.x_label.font_size = 3;
eixos.y_label.font_size = 3;
eixos.z_label.font_size = 3;
eixos.tight_limits = 'on';
eixos.data_bounds(1,3) = phasemin;
eixos.data_bounds(2,3) = phasemax;
eixos.zoom_box =  [min(b), min(a), max(b), max(a), phasemin, phasemax];

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
cbar.parent.title.fill_mode = "off"
cbar.parent.title.font_size = 3
plot3d1(b,a,matZ)
e=gce();
e.hiddencolor=-1;
e.color_mode=colormode;
xtitle( 'Bode Plot', eix, 'log (f) / Hz', '|Z| / Ohm.cm²' , boxed = 0 )
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
cbar.parent.title.fill_mode = "off"
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
ax.zoom_box =  [min(b), min(a), max(b), max(a), vminlmq, vmaxlmq];
xtitle( 'Bode Plot', eix, 'log (f) / Hz', 'log(|Z|/Ohm.cm²)' , boxed = 0 )
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
    cbar.parent.title.text = "-Im";
    cbar.parent.title.fill_mode = "off"
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
    xtitle( 'Nyquist Diagram', eix, 'Real', '-Im' , boxed = 0)
    filename='nyquist'
    eixos=get("current_axes")
    eixos.title.font_size = 3;
    eixos.x_label.font_size = 3;
    eixos.y_label.font_size = 3;
    eixos.z_label.font_size = 3;
end

///////////////////////// fim plots 3d ///////////////////////////////////
end

if kk_validate == 1 then
        matImagKK = [];
    matRealKK = [];
    vfreq = inverte(matFreq(1,:));
    for i = 1:nexp
        Zreal = inverte(matRealZ(i,:));
        Zimag = -1 * abs(inverte(matImagZ(i,:)));
//        imagCalcKK = imag_kk(vfreq, Zreal);
//        realCalcKK = real_kk(vfreq, -Zimag, i);

        aux = extrapolacao(vfreq, Zreal, Zimag);
        vfreq_extra=aux(1,:);
        Zreal_extra=aux(2,:);
        Zimag_extra=aux(3,:);
        imagCalcKK_extra = imag_kk(vfreq_extra, Zreal_extra);
        realCalcKK_extra = real_kk(vfreq_extra, -Zimag_extra,i);
        imagCalcKK = imagCalcKK_extra(21:length(imagCalcKK_extra));
        realCalcKK = realCalcKK_extra(21:length(realCalcKK_extra));
        

        matRealKK(:,i) = inverte(realCalcKK)';
        matImagKK(:,i) = inverte(imagCalcKK)';
        
//        clf(900);
//        scf(900);
//        plot(log10(vfreq), (Zreal),'b*');
//        plot(log10(vfreq), -(Zimag),'r*');
//        plot(log10(vfreq), (realCalcKK),'b');
//        plot(log10(vfreq), -(imagCalcKK),'r');
//        xtitle('E = ' + string(vddp(i)) + " V", 'log (f /Hz)', 'Real e Imag (Exp e KK)', '');
//        eixos=get("current_axes")
//        eixos.title.font_size = 3;
//        eixos.x_label.font_size = 3;
//        eixos.y_label.font_size = 3;
//        eixos.z_label.font_size = 3;
//        pause;
    end
    matRealKK=matRealKK';
    matImagKK=-matImagKK';
//    erKK_real = abs(matRealKK-matRealZ)/max(abs(matRealZ));
//    erKK_imag = abs(matImagKK - matImagZ)/max(abs(matImagZ));
    erKK_real = (abs((matRealKK-matRealZ) ./((matRealZ)))*100);
    erKK_imag = (abs((matImagKK - matImagZ) ./((matImagZ)))*100);
    

    for i=1:nexp
        for j=1:numfreq
            if erKK_real(i,j) > 100
                 erKK_real(i,j)=100;
            end
        end
    end
    
    for i=1:nexp
        for j=1:numfreq
            if erKK_imag(i,j) > 100
                 erKK_imag(i,j)=100;
            end
        end
    end

    clf(9)
    scf(9)
    xset("colormap",jetcolormap(512))
    colorbar(min(erKK_imag),max(erKK_imag))
    cbar = gce();
    cbar.parent.title.text = "% Deviation";
    cbar.parent.title.fill_mode = "off";
    cbar.parent.title.font_size = 3;
//    plot3d1(b,a,erKK_imag);
    grayplot(b,a,(erKK_imag));
    xtitle( 'Kramers-Kronig - Imaginary', eix, 'log (f /Hz)', 'Percentual Error' , boxed = 1);
    eixos=get("current_axes");
    eixos.title.font_size = 3;
    eixos.x_label.font_size = 3;
    eixos.y_label.font_size = 3;
    eixos.z_label.font_size = 3;
    if wg==1 then
        eixos.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
    end
    
    clf(10)
    scf(10)
    xset("colormap",jetcolormap(512));
    colorbar(min(erKK_real),max(erKK_real));
    cbar = gce();
    cbar.parent.title.text = "% Deviation";
    cbar.parent.title.fill_mode = "off";
    cbar.parent.title.font_size = 3;
    grayplot(b,a,(erKK_real));
    xtitle( 'Kramers-Kronig - Real', eix, 'log (f /Hz)', '% Error' , boxed = 1);
    eixos=get("current_axes");
    eixos.title.font_size = 3;
    eixos.x_label.font_size = 3;
    eixos.y_label.font_size = 3;
    eixos.z_label.font_size = 3;
    if wg==1 then
        eixos.x_ticks = tlist(["ticks","locations","labels"], vet_nexp, pot_string);
    end
end

messagebox(["Images are available"],"EIS - Map Generator","message",["OK"],'modal');
else
    disp('End');
    messagebox(["End"],"EIS - Map Generator","message",["OK"],'modal');
end
// Close all files
file('close',file() )
