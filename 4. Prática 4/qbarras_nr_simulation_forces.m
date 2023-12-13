%% Inputs: 
% Angles: t2v, t3v, t4v
% Mechanism lengths: a, b, c, d
% End-effector position: AP, tAP
% Support forces: F12, F14

%% Evaluate positions of points O2,A,B,O4,P
rO2=zeros(length(t2v),2);
rA=a*[cos(t2v) sin(t2v)];
rB=rA+b*[cos(t3v) sin(t3v)];
rO4=[rO2(:,1)+d rO2(:,2)];
rP=rA+AP*[cos(t3v+tAP) sin(t3v+tAP)];

%% Show simulation
outvid=0; % if =1 outputs GIF video animation
figure(51), clf, set(51,'position',[0 0 690 650])
if outvid==1 
    filename='qbarras_nr.gif';
    frame=getframe(gcf); im=frame2im(frame); [A,map]=rgb2ind(im,256);
    imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1e-3);
end
rT=[rO2; rA; rB; rO4; rP]; 
maxfvl=a;
minx=min(rT(:,1))-maxfvl; maxx=max(rT(:,1))+maxfvl; miny=min(rT(:,2))-maxfvl; maxy=max(rT(:,2))+maxfvl; 
axis equal, axis off, axis([minx maxx miny maxy]),

hAP=line([rA(1,1) rP(1,1)],[rA(1,2) rP(1,2)]); set(hAP,'Color','g','LineStyle','-','Marker','o'),
hPB=line([rP(1,1) rB(1,1)],[rP(1,2) rB(1,2)]); set(hPB,'Color','g','LineStyle','-','Marker','o'),
hO2A=line([rO2(1,1) rA(1,1)],[rO2(1,2) rA(1,2)]); set(hO2A,'Color','k','LineStyle','-','Marker','o'),
hAB=line([rA(1,1) rB(1,1)],[rA(1,2) rB(1,2)]); set(hAB,'Color','k','LineStyle','-','Marker','o'),
hBO4=line([rB(1,1) rO4(1,1)],[rB(1,2) rO4(1,2)]); set(hBO4,'Color','k','LineStyle','-','Marker','o'),
hO4O2=line([rO4(1,1) rO2(1,1)],[rO4(1,2) rO2(1,2)]); set(hO4O2,'Color','k','LineStyle','-','Marker','o'),

hO2Ao=line([rO2(1,1) rA(1,1)],[rO2(1,2) rA(1,2)]); set(hO2Ao,'Color',.7*[1 1 1],'LineStyle',':'),
hABo=line([rA(1,1) rB(1,1)],[rA(1,2) rB(1,2)]); set(hABo,'Color',.7*[1 1 1],'LineStyle',':'),
hBO4o=line([rB(1,1) rO4(1,1)],[rB(1,2) rO4(1,2)]); set(hBO4o,'Color',.7*[1 1 1],'LineStyle',':'),
hO4O2o=line([rO4(1,1) rO2(1,1)],[rO4(1,2) rO2(1,2)]); set(hO4O2o,'Color',.7*[1 1 1],'LineStyle',':'),
mx=maxx-minx; my=maxy-miny;
text(rO2(1,1)-mx/100,rO2(1,2)-my/20,'$O_2$')
text(rA(1,1)-mx/100,rA(1,2)+my/20,'$A$')
text(rB(1,1)-mx/100,rB(1,2)+my/20,'$B$')
text(rO4(1,1)-mx/100,rO4(1,2)-my/20,'$O_4$')
text(rP(1,1)-mx/100,rP(1,2)+my/20,'$P$')

% Support forces and torques
maxF=max(max(abs([F12 F14])))/maxfvl; Fs=F12+F14;
hq12=line([rO2(1,1) rO2(1,1)-F12(1,1)/maxF],[rO2(1,2) rO2(1,2)-F12(1,2)/maxF]);
hq14=line([rO4(1,1) rO4(1,1)-F14(1,1)/maxF],[rO4(1,2) rO4(1,2)-F14(1,2)/maxF]);
hqs=line([(rO2(1,1)+rO4(1,1))/2 (rO2(1,1)+rO4(1,1))/2-Fs(1,1)/maxF],[(rO2(1,2)+rO4(1,2))/2 (rO2(1,2)+rO4(1,2))/2-Fs(1,2)/maxF]);
set([hq12 hq14 hqs],'linewidth',3)
 
if outvid==1 dts=5; else dts=1; end
for n=2:dts:length(t)
    set(hO2A,'xdata',[rO2(n,1) rA(n,1)],'ydata',[rO2(n,2) rA(n,2)]);
    set(hAB,'xdata',[rA(n,1) rB(n,1)],'ydata',[rA(n,2) rB(n,2)]);
    set(hBO4,'xdata',[rB(n,1) rO4(n,1)],'ydata',[rB(n,2) rO4(n,2)]);
    set(hAP,'xdata',[rA(n,1) rP(n,1)],'ydata',[rA(n,2) rP(n,2)]);
    set(hPB,'xdata',[rP(n,1) rB(n,1)],'ydata',[rP(n,2) rB(n,2)]);
    hA=line([rA(n-1,1) rA(n,1)],[rA(n-1,2) rA(n,2)]); set(hA,'Color','r','LineStyle',':');
    hB=line([rB(n-1,1) rB(n,1)],[rB(n-1,2) rB(n,2)]); set(hB,'Color','r','LineStyle',':');
    hP=line([rP(n-1,1) rP(n,1)],[rP(n-1,2) rP(n,2)]); set(hP,'Color','m','LineStyle',':');
    set(hq12,'xdata',[rO2(n,1) rO2(n,1)-F12(n,1)/maxF],'ydata',[rO2(n,2) rO2(n,2)-F12(n,2)/maxF]);
    set(hq14,'xdata',[rO4(n,1) rO4(n,1)-F14(n,1)/maxF],'ydata',[rO4(n,2) rO4(n,2)-F14(n,2)/maxF]);
    set(hqs,'xdata',[(rO2(n,1)+rO4(n,1))/2 (rO2(n,1)+rO4(n,1))/2-Fs(n,1)/maxF],'ydata',[(rO2(n,2)+rO4(n,2))/2 (rO2(n,2)+rO4(n,2))/2-Fs(n,2)/maxF]);

    if outvid==1 
        frame=getframe(gcf); im=frame2im(frame); [A,map]=rgb2ind(im,256);
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1e-3);
    else pause(1e-12);
    end
end,
