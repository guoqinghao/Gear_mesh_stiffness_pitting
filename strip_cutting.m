function [cell_pit,get_strip_point]=strip_cutting(cell_pit)
global y11  
get_strip_point=[];                                                % 得到切条后坐标
Index_slice=find( ~cellfun('isempty',cell_pit(1,:)));              % 得到涉及点蚀的轮廓索引

for ii=1:length(Index_slice)                                       %%%%% 所有切片
    i=Index_slice(ii);                                             % 对应的轮廓点序号
    slice0=cell_pit{1,i};                                          % 得到当前切片
    
    %% 获得当前切片统一的切条厚度向量
    get_strip_end2=[];                                             % 得到切条z坐标矩阵
    state_strip=ones(1,length(slice0(1,:)));                       % 记录每个点蚀的状态（是否切完）
    strip_start=y11(i);                                            % 切条初始值
    strip_end=strip_start;                                         % 切条终止值

    strip_start1=[];                                               % 判断是否有点超过y11
    while any(state_strip==1)                                      %%%%% 如果切片上还有数据，切一条
        strip_z=[];                                                % 记录该层每个切片“切一条”的切条厚度
        I_state=find(state_strip==1);                              % 找到切片上不为空的点蚀的索引  
        for j=1:length(I_state)                                    %%%%% 遍历切片上每个点蚀
            jj=I_state(j);
            %确定当前点蚀切条分辨率和切条横向要求
            slice_strip1=slice0{1,jj};                             % 获得单个切片坐标
            slice_strip1_y=slice_strip1(1,:);
            slice_strip1_z=slice_strip1(2,:);
            [slice_strip1zz,~]=sort(slice_strip1_z,'descend');     % 按z排序,降序
            delta_strip_z=diff(slice_strip1zz);                    % 求相邻点z坐标差值
            delta_strip_zmean=mean(abs(delta_strip_z));            % 确定切条分辨率
            
            [slice_strip1yy,~]=sort(slice_strip1_y,'descend');     % 按y排序,降序
            delta_strip_y=diff(slice_strip1yy);                    % 求相邻点y坐标差值
            delta_strip_ymean=mean(abs(delta_strip_y));            % 确定切条要求
               
            % 判断是否有点超过y11
            order_strip_test10=find(slice_strip1_z > y11(i));
            order_strip_test20=find(slice_strip1_z <= y11(i)+10);
            order_strip_test0=intersect(order_strip_test10,order_strip_test20);
            strip_matrix_y_test0=slice_strip1_z(order_strip_test0);  % 切条初始数据
            if ~isempty(strip_matrix_y_test0)
                strip_start1=[strip_start1,max(strip_matrix_y_test0)];
%                 disp("有点超过y11");
            end
            
            num_wh=0;                                              % 迭代次数
            while 1                                                %%%%% 获取当前点蚀当前切条厚度
                num_wh=num_wh+1;                                   % 更新迭代次数
                strip_end=strip_end-delta_strip_zmean;             % 切条终止值迭代
                % 单步切条
                order_strip_1=find(slice_strip1_z > strip_end);
                order_strip_2=find(slice_strip1_z <= strip_start);
                order_strip=intersect(order_strip_1,order_strip_2);
                strip_matrix_y=slice_strip1_y(order_strip);        % 切条初始数据
                delta_strip_matrix_y=abs(max(strip_matrix_y)-min(strip_matrix_y)); % 切条长度
                % 处理切条范围内没点的情况
                if isempty(delta_strip_matrix_y)                   
                    delta_strip_matrix_y=0;  
                end
                % 判断当前点蚀切条厚度是否满足要求
                if delta_strip_matrix_y > 5*delta_strip_ymean     
                    strip_z(end+1)=strip_end;                      % 记录当前点蚀切条z坐标
                    break;                    
                end
                % 判断迭代次数,检查是否切完
                if num_wh > 5
                    % 检查之后切条是否有点
                    strip_start_test=strip_end;
                    strip_end_test=0;  
                    order_strip_1_test=find(slice_strip1_z > strip_end_test);
                    order_strip_2_test=find(slice_strip1_z <= strip_start_test);
                    order_strip_test=intersect(order_strip_1_test,order_strip_2_test);
                    if isempty(order_strip_test)                   % 若无点，则该切片遍历完
                        state_strip(jj)=0;                         % 标记该点蚀切完
                        break;
                    end                    
                end
            end
            
        end
        if isempty(strip_z)
            break;
        end
        get_strip_end1=min(strip_z);
        get_strip_end2(end+1)=get_strip_end1;                      % 记录公共切条的z坐标
        strip_start=get_strip_end1;                                % 更新切条初始值
        strip_end=strip_start;

    end
    
    
    %% 用厚度向量进行切条
    delta_strip_matrix1=[];                                        % 切条后的长度,厚度，z坐标
    % 根据是否有点超过y11，决定总体切条的起始点
    if isempty(strip_start1)
        strip_start=y11(i);                                            % 切条初始值
    else
        strip_start=max(strip_start1);
    end

    for m=1:length(get_strip_end2)                                 % 遍历每个切条厚度
        strip_end=get_strip_end2(m);                               % 切条终止值
        delta_strip_matrix0=[];                                    % 切条后的长度分量
        for n=1:length(slice0(1,:))                                %%%%% 遍历切片上每个点蚀
            % 获得单个切片坐标
            slice_strip1=slice0{1,n};                             
            slice_strip1_y=slice_strip1(1,:);
            slice_strip1_z=slice_strip1(2,:);
            % 单步切条
            order_strip_1=find(slice_strip1_z >= strip_end);
            order_strip_2=find(slice_strip1_z <= strip_start);
            order_strip=intersect(order_strip_1,order_strip_2);
            strip_matrix=slice_strip1_y(order_strip);               % 切条数据
            if isempty(strip_matrix)                   
                delta_strip_matrix0=0;
                if strip_start >= y11(i)
                    disp("第一个切片就为空");
                end
            else
                delta_strip_matrix0(end+1)=abs(max(strip_matrix)-min(strip_matrix));             % 切条长度
                get_strip_point(:,end+1)=[max(strip_matrix);(strip_end+strip_start)/2;y11(i)]; % 切条后坐标
                get_strip_point(:,end+1)=[min(strip_matrix);(strip_end+strip_start)/2;y11(i)];
            end
            % 处理切条范围内没点的情况
            if isempty(delta_strip_matrix0)                   
                delta_strip_matrix0=0;
                disp("切条厚度下无对应切条");
            end
        end
        if strip_start>y11(i)
            delta_strip_matrix1(:,end+1)=[sum(delta_strip_matrix0);y11(i)-strip_end;strip_end];%切条后的长度,厚度，z坐标
        else
            delta_strip_matrix1(:,end+1)=[sum(delta_strip_matrix0);strip_start-strip_end;strip_end];%切条后的长度,厚度，z坐标
        end
        strip_start=strip_end;                                      % 更新切条起始点

    end
    cell_pit{2,i}=delta_strip_matrix1;
end
      
