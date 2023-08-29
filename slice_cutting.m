function [cell_pit]=slice_cutting(get_p1,if_s,cell_pit)
global x11 y11

for i=1:length(get_p1)
    % 获得修剪后的坐标
    x_p1=get_p1{i}(1,:); % 获取x坐标
    y_p1=get_p1{i}(2,:); % 获取y坐标
    z_p1=get_p1{i}(3,:); % 获取z坐标
    % 找到切片起始点/结束点
    x_p1min=min(x_p1);
    x_p1max=max(x_p1);
    [~,Index1] = min(abs(x11-x_p1min));
    [~,Index2] = min(abs(x11-x_p1max));
    % 修剪
    cell_s={};
    num=1;
    for j=Index1:Index2-1
        order_xx_1=find(x_p1 >= x11(j));
        order_xx_2=find(x_p1 < x11(j+1));
        order_xx=intersect(order_xx_1,order_xx_2);
        slice_matrix_y=y_p1(order_xx);
        slice_matrix_z=z_p1(order_xx);
        slice_matrix=[slice_matrix_y;slice_matrix_z];
        order_clean=find(abs(slice_matrix_z) > y11(j));
        slice_matrix(:,order_clean)=[];
        if ~isempty(slice_matrix)
            order_cell_s(1,num)=j;      % 修剪后切片对应齿廓索引起始
            order_cell_s(2,num)=j+1;           
            cell_s{1,num}=slice_matrix; % 修剪后切片坐标
            num=num+1;
        end
    end    
    

% 对整体切片
    std=1e-05; % delta_now要求
    std1=1e-05; % delta_now1要求
    std2=1e-08; % delta_now2要求
    slice_select1={}; % 筛选后的切片坐标
    order_fig_select=[];  % 筛选后的切片索引
    delta_now1=[];
    delta_now=[];
    delta_now2=[];
    Index_start=1;
    Index_end=Index_start;
    while Index_end <= length(cell_s)
        get_p1_yz=[];
        state_num=1;
        for m=Index_start:Index_end
            get_p1_yz=[get_p1_yz , cell_s{m}];
        end
        %对坐标排序整合
        get_p1_y=get_p1_yz(1,:);
        get_p1_z=get_p1_yz(2,:);
        [p1_y,I]=sort(get_p1_y);
        delta_now=mean(diff(p1_y,1));    % 切片密度指标   
        [p2_z,I2]=sort(get_p1_z,'descend');
        p2_y=get_p1_y(I2);
        delta_now1=abs(p2_z(1)-y11(order_cell_s(1,Index_start))); % 切片边界指标
        if length(p2_y)>=7
            delta_now2=abs(max(p2_y(1:5))-min(p2_y(1:5))); % 切片边界厚度指标 
        else
            delta_now2=abs(max(p2_y(:))-min(p2_y(:)));
            state_num=0;
            disp("点数少于7");
        end
        if delta_now <= std && delta_now1 <= std1 && delta_now2 >= std2 && state_num
            p1_z=get_p1_z(I);       
            slice_matrix1=[];
            if Index_end ~= length(cell_s)
                slice_matrix1(1,:)=p1_y;
                slice_matrix1(2,:)=p1_z;
                slice_matrix1_pre=slice_matrix1;
            else
                slice_matrix1= slice_matrix1_pre;
            end
            order_fig_select1(1,1)=order_cell_s(1,Index_start);
            order_fig_select1(2,1)=order_cell_s(2,Index_end);
            slice_select1{end+1}=slice_matrix1;
            order_fig_select=[order_fig_select,order_fig_select1];       
            Index_start=Index_end+1;
            Index_end=Index_start;   
        else
            Index_end=Index_end+1;  
       end   
    end            


% 切片坐标归类到相应轮廓x坐标
    for n=1:length(slice_select1)
        slice_select_now=slice_select1{n};
        Index_start0=order_fig_select(1,n);
        Index_end0=order_fig_select(2,n);        
        for k=(Index_start0):Index_end0-1
            cell_pit{1,k}{1,end+1}=slice_select_now;
            cell_pit{1,k}{2,length(cell_pit{1,k}(1,:))}=i;
            cell_pit{1,k}{3,length(cell_pit{1,k}(1,:))}=if_s(i);
        end

    end
    % 处理最后一个切片
        nn=length(slice_select1);
        slice_select_now=slice_select1{nn};
        Index_start0=order_fig_select(1,nn);
        Index_end0=order_fig_select(2,nn);
        kk=Index_end0;
        cell_pit{1,kk}{1,end+1}=slice_select_now;
        cell_pit{1,kk}{2,end}=i;
        cell_pit{1,kk}{3,end}=if_s(i);

end

end



