%--------------------------------------------------------------------------
% auther:Binglong Zhu
% Email: blzhu@buaa.edu.cn; sdszbl@163.com
% Date: 2022年5月14日20:20:51
% Function: Successive Parabolic Interpolation + Axesion transformation
%--------------------------------------------------------------------------
clc;
close all;
clear all;

format long;
figure(1)
%% function 1
%  f=@(x)x.^3-13*x.^2+20*x+100;% best_position= 7.813437226365976,best_fitness= -60.369882183024487
% Vmin=-5;
% Vmax=20;
% fplot(f,[Vmin Vmax],'k');hold on;
% % axis equal;
% xbest=12;
% [x,fo]=successiveparabolicinterpolation(f,Vmin,Vmax,xbest);




%% function 2
% https://github.com/gopherclass/successive-parabolic-interpolation
% best_position= 1.427551777634208 ,best_fitness=-1.775725653147415
f=@(x)x.^2/10-2*sin(x);
Vmin=-20;
Vmax=20;
fplot(f,[Vmin Vmax],'k');hold on;
% axis equal;
xbest=-14;
[x,fo]=successiveparabolicinterpolation(f,Vmin,Vmax,xbest);


saveas( 1, 'Result.jpg')

figure(2)
plot(1:1:numel(fo),fo,'-*')
saveas( 2, 'ft.jpg')

disp(['best_position:' num2str(x(end),'% 20.16g') ',best_fitness:' ...
    num2str(fo(end),'% 20.16g')] );


%% Subfunction
function [xiter,fiter]=successiveparabolicinterpolation(f,Vmin,Vmax,xbest)
if Vmin>xbest || xbest>Vmax
    disp('ERROR! the xbest should be given in range of [Vmin,Vmax]!');
    xiter=xbest;
    fiter=f(xbest);
    return;
end
xtol=1e-4;
ytol=1e-8;
kmax=10;
pmax=3;
for k=1:kmax
    
    beta=(Vmax-Vmin)/5*(kmax-k+1)/kmax % Maximum axial transformation radius，越小越慢的走出下坡，越小越难跳出局部最优，大了利于全局搜索
    betaArr(k)=beta;
    delt=beta;% trust domain radius,越小越利于局部搜索，越大越利于全局搜索
    
    % 生成初始三个点，这三个点要确保在区间内
    xc=[xbest-delt*rand xbest  xbest+delt*rand];% trust domain
    while 1
        if xc(1)>Vmin && xc(3)<Vmax
            break;
        else
            xc=[xbest-delt*rand xbest  xbest+delt*rand];% trust domain
        end
    end
    fc=f(xc);% solve Error value
    %三点按误差函数升序排列
    [fc,fcI]=sort(fc);% ascending order
    xc=xc(fcI);
    
    
    xbest=xc(1);
    ybest=fc(1);
    
    % 绘制出这三个点
    plot(xc,fc,'o','markerfacecolor','w');% plot three points
    
    %求极点位置对称轴xe
    xe=(fc(1)*(xc(2)^2-xc(3)^2)-fc(2)*(xc(1)^2-xc(3)^2)+fc(3)*(xc(1)^2-xc(2)^2))/(2*(fc(1)*(xc(2)-xc(3))-fc(2)*(xc(1)-xc(3))+fc(3)*(xc(1)-xc(2))));
    if isnan(xe)||isinf(xe)
        xe=xc(1);
    end
    
    %求解xe处误差函数fxe
    if xe<Vmax && xe>Vmin
        fxe=f(xe);
        [a,b,c]=parabola(xc,fc);%当前抛物线绘图
        plot(xe,fxe,'o','markerfacecolor','w','color','b');hold on;%plot the vertex of the parabola
        
    else
        fxe= inf ;
    end
    
    if fxe<fc(1)
        xbest=xe; % 每个分支都给xbest和ybest赋值
        ybest=fxe;
        %抛物线搜索
        p=0;
        while 1
            p=p+1
            if(p>pmax)
                break;
            end
            if any(abs(xe-xc)<xtol) || all(abs(fxe-fc)<ytol)
                break;
            else
                if fxe<fc(1)
                    xc(3)=xe;% 用新找的最好的解xe代替之前最差的xc(3)
                    fc(3)=fxe;
                    [xe,fxe]=Parabolic(f,Vmax,Vmin,xc,fc)%有可能fxe变大了
                    if Vmin>xe || xe>Vmax %如果xe不在参数范围内，则给定xbest
                        xe=xbest;
                        fxe=ybest;
                    else %求出新的fxe，新的fxe可能变大了
                        if fxe<ybest
                        xbest=xe;
                        ybest=fxe;
                        end
                    end
                else
                    break;
                end
            end
        end
    else
        %状态转移搜索
        xk=xc(1);%已有的误差最小的x
        
        xkm1=xc(2);%x(k-1) %中间误差xc(2)
        x2new=xk+beta*rand*(xk-xkm1)/norm(xk-xkm1);%x(k+1) TT转移算子
        fx2new=f(x2new);
        
        xkm1=xc(3);%x(k-1) %最大误差的处xc(2)
        x3new=xk+beta*rand*(xk-xkm1)/norm(xk-xkm1);%x(k+1) TT转移算子
        fx3new=f(x3new);
        
        if fx2new<=fx3new && fx2new<ybest
            xbest=x2new;
            ybest=fx2new;
            drawarrow([xc(2), f(xc(2))], [xbest, ybest]);hold on;% 有箭头就表明有抛物线开口朝下了
            plot(xbest,ybest,'p');hold on;% the best point with minimum f
        elseif fx3new<fx2new && fx3new<ybest
            xbest=x3new;
            ybest=fx3new;
            drawarrow([xc(3), f(xc(3))], [xbest, ybest]);hold on;% 有箭头就表明有抛物线开口朝下了
            plot(xbest,ybest,'p');hold on;% the best point with minimum f
        end
        
    end
    %     %如果是最后一次迭代，这特意求出最优解xbest的y值
    %     ybest=f(xbest);
    %
    xiter(k)=xbest;
    fiter(k)=ybest;
    if k>1 && ybest>fiter(k-1)
        disp('error!');
    end
end
plot(xbest,ybest,'*');hold on;% the best point with minimum f
end

%% 抛物线插值求极值
function [xe,fxe]=Parabolic(f,Vmax,Vmin,x,fx)
% 按误差升序排列
[fx,fxI]=sort(fx);% ascending order
x=x(fxI);

% 找到过(xi,fxi)这三点的抛物线的对称轴x，或者说求极点位置对称轴xe
xe=(fx(1)*(x(2)^2-x(3)^2)-fx(2)*(x(1)^2-x(3)^2)+fx(3)*(x(1)^2-x(2)^2))/(2*(fx(1)*(x(2)-x(3))-fx(2)*(x(1)-x(3))+fx(3)*(x(1)-x(2))));
if isnan(xe)||isinf(xe)
    xe=x(1);
end

%求解xe处误差函数fxe
if xe<Vmax && xe>Vmin
    fxe=f(xe);
    plot(xe,fxe,'o','markerfacecolor','w','color','b');hold on;%plot the vertex of the parabola
else
    fxe= inf ;
end

end



% 绘制抛物线，并给出系数a b c
function [a,b,c]=parabola(xc,fc)
x1=xc(1);x2=xc(2);x3=xc(3);
y1=fc(1);y2=fc(2);y3=fc(3);
% x1=p1(1);y1=p1(2);
% x2=p2(1);y2=p2(2);
% x3=p3(1);y3=p3(2);
if x1~=x2 && x2~=x3 && x1~=x3
    a=-((-x2* y1+x3 *y1+x1* y2-x3* y2-x1* y3+x2* y3)/((-x1+x2) *(x2-x3) *(-x1+x3)));
    b=-((x2^2 *y1 - x3^2 *y1 - x1^2 *y2 + x3^2 *y2 + x1^2 *y3 - x2^2 *y3)/((x1 - x2) *(x1 - x3) *(x2 - x3)));
    c=-((-x2^2 *x3 *y1+x2 *x3^2 *y1+x1^2 *x3 *y2-x1 *x3^2 *y2-x1^2 *x2 *y3+x1 *x2^2 *y3)/((x1-x2)* (x1-x3) *(x2-x3)));
else
    a=0;b=0;c=0;
end

% P=[a,b,c];
% x=min(px):0.1:max(px);
px=xc;
x=linspace(min(px)-(max(px)-min(px))/0.5,max(px)+(max(px)-min(px))/0.5);
y=a*x.^2+b*x+c;
plot(x,y,'--');
hold on;
end

function drawarrow(x,y,lineType,ax)
switch nargin
    case 2
        lineType='arrow';
        ax=gca;
    case 3
        ax=gca;
end
% 调整坐标大小以适应箭头长度
xlim=ax.XLim;
ylim=ax.YLim;
xlimmin=xlim(1);xlimmax=xlim(2);
ylimmin=ylim(1);ylimmax=ylim(2);
if xlimmin>min(x(1),y(1)), xlimmin=min(x(1),y(1));end
if xlimmax<max(x(1),y(1)), xlimmax=max(x(1),y(1));end
if ylimmin>min(x(2),y(2)), ylimmin=min(x(2),y(2));end
if ylimmax<max(x(2),y(2)), ylimmax=max(x(2),y(2));end
ax.XLim = [xlimmin,xlimmax];
ax.YLim = [ylimmin,ylimmax];
xlim=ax.XLim;
ylim=ax.YLim;
pos=ax.Position;
x_ratio = pos(3)/(xlim(2)-xlim(1));
y_ratio = pos(4)/(ylim(2)-ylim(1)); % 缩放比例
orig_pos=[-xlim(1)*x_ratio+pos(1),-ylim(1)*y_ratio+pos(2)]; % figure坐标系中的原点坐标
x=x.*[x_ratio,y_ratio];y=y.*[x_ratio,y_ratio];
x=x+orig_pos;y=y+orig_pos;
% annotation(lineType,[x(1),y(1)],[x(2),y(2)])
ar=annotation(lineType,[x(1),y(1)],[x(2),y(2)]);
ar.Color='Red';
ar.LineStyle='--';
ar.LineWidth=1;
end
