% �Զ���MATLAB����ErrorEllipse����ӡ��ά���ݵ�������Բ

% �޸���ԭ����http://���ﲻ���visiondummy��com/2014/04/draw-error-ellipse-representing-covariance-matrix/

% ����ֵΪһ��n��2����ֵ��������Ÿ���p

function [a,b,R,X0,Y0,r_ellipse] = ErrorEllipse(datamatrix,p,colo,style,sign,fill,pl)

data = datamatrix;

% ����Э���������������������ֵ
covariance = cov(data);
[eigenvec, eigenval ] = eig(covariance);

% ��ȡ�����������
[largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);

% ��ȡ�������ֵ
largest_eigenval = max(max(eigenval));

% ������С������������С����ֵ
if(largest_eigenvec_ind_c == 1)
    smallest_eigenval = max(eigenval(:,2));
    smallest_eigenvec = eigenvec(:,2);
else
    smallest_eigenval = max(eigenval(:,1));
    smallest_eigenvec = eigenvec(1,:);
end

% ����X��������������ֱ�ӵļнǣ�ֵ��Ϊ[-pi,pi]
angle = atan2(largest_eigenvec(2), largest_eigenvec(1));

% ���н�Ϊ��ʱ����2pi����ֵ
if(angle < 0)
    angle = angle + 2*pi;
end

% �������ݵ����о�ֵ����ʽΪ2��1�ľ���
avg = mean(data);

% ����������Բ�Ĳ�������������ֵ����ת�Ƕȡ���ֵ���������
chisquare_val = sqrt(chi2inv(p,2));
theta_grid = linspace(0,2*pi);
phi = angle;

X0=avg(1);
Y0=avg(2);
a=chisquare_val*sqrt(largest_eigenval);
b=chisquare_val*sqrt(smallest_eigenval);

% ����ԲͶ�䵽ֱ���������� 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

% ��ת����
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

% ��ˣ���ת��Բ
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

if pl == 1
    % ��ӡ������Բ
    plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'Color',colo,'LineStyle',style,'LineWidth',1)
    
    % ��ӡԭʼ����
    % plot(data(:,1), data(:,2), '.');
    hold on
    
    if fill == 1
%         scatter(data(:,1),data(:,2),25,'filled',colo,'MarkerFaceAlpha',.6)
        hold on
        e = scatter(mean(data(:,1)),mean(data(:,2)),40,'d','filled');
        e.MarkerFaceColor = colo;
        e.MarkerEdgeColor = colo;
        e.LineWidth = 2;
%         scatter(mean(data(:,1)),mean(data(:,2)),80,colo,'d','MarkerEdgeColor',colo,'LineWidth',2)
    else
        %     scatter(data(:,1),data(:,2),10,colo,'MarkerEdgeAlpha',.8)
        %     scatter(mean(data(:,1)),mean(data(:,2)),100,colo,sign)
    end
    
end
% if fill == 1
% %     scatter(data(:,1),data(:,2),10,colo,'filled','MarkerFaceAlpha',1)
%     scatter(mean(data(:,1)),mean(data(:,2)),100,colo,sign)
% else
% %     scatter(data(:,1),data(:,2),10,colo,'MarkerEdgeAlpha',.4)
%     scatter(mean(data(:,1)),mean(data(:,2)),100,colo,sign)
% end


end