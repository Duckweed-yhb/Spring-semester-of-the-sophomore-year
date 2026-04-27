% =========================================================
% 二维静电场电位分布求解：高斯-赛德尔 & 超松弛迭代法
% 边界条件：
%   上边界 φ = 30V
%   左边界 φ = 0V
%   下边界（右侧）∂φ/∂n = 0（Neumann边界）
%   右边界 φ = 100V
%   左下角凹区域 φ = 0V
% =========================================================
clear; clc; close all;

%% 1. 网格与参数设置
nx = 20;   % x方向节点数
ny = 20;   % y方向节点数
tol = 1e-5; % 收敛误差
max_iter = 10000; % 最大迭代次数
omega = 1.5; % 超松弛因子（可根据需要调整）

phi_GS = zeros(ny, nx);   % 高斯-赛德尔迭代结果
phi_SOR = zeros(ny, nx);   % 超松弛迭代结果

%% 2. 边界条件初始化
% 上边界：φ = 30V
phi_GS(1, :) = 30;
phi_SOR(1, :) = 30;

% 左边界：φ = 0V
phi_GS(:, 1) = 0;
phi_SOR(:, 1) = 0;

% 右边界：φ = 100V
phi_GS(:, end) = 100;
phi_SOR(:, end) = 100;

% 左下角凹区域：φ = 0V
phi_GS(ny/2+1:end, 1:ny/2) = 0;
phi_SOR(ny/2+1:end, 1:ny/2) = 0;

%% 3. 高斯-赛德尔迭代法
fprintf('=== 高斯-赛德尔迭代法 ===\n');
iter_GS = 0;
while iter_GS < max_iter
    max_diff = 0;
    for j = 2:ny-1
        for i = 2:nx-1
            % 跳过左下角凹区域
            if j > ny/2 && i <= ny/2
                continue;
            end
            
            old_val = phi_GS(j, i);
            phi_GS(j, i) = 0.25 * (phi_GS(j-1, i) + phi_GS(j+1, i) + ...
                                    phi_GS(j, i-1) + phi_GS(j, i+1));
            
            % 处理 Neumann 边界（下边界 ∂φ/∂n = 0）
            if j == ny
                phi_GS(j, i) = 0.25 * (2*phi_GS(j-1, i) + phi_GS(j, i-1) + phi_GS(j, i+1));
            end
            
            diff = abs(phi_GS(j, i) - old_val);
            if diff > max_diff
                max_diff = diff;
            end
        end
    end
    iter_GS = iter_GS + 1;
    if max_diff < tol
        break;
    end
end
fprintf('迭代次数：%d\n', iter_GS);

%% 4. 超松弛迭代法（SOR）
fprintf('\n=== 超松弛迭代法（SOR） ===\n');
iter_SOR = 0;
while iter_SOR < max_iter
    max_diff = 0;
    for j = 2:ny-1
        for i = 2:nx-1
            % 跳过左下角凹区域
            if j > ny/2 && i <= ny/2
                continue;
            end
            
            old_val = phi_SOR(j, i);
            new_val = 0.25 * (phi_SOR(j-1, i) + phi_SOR(j+1, i) + ...
                              phi_SOR(j, i-1) + phi_SOR(j, i+1));
            
            % 处理 Neumann 边界（下边界 ∂φ/∂n = 0）
            if j == ny
                new_val = 0.25 * (2*phi_SOR(j-1, i) + phi_SOR(j, i-1) + phi_SOR(j, i+1));
            end
            
            % SOR 更新
            phi_SOR(j, i) = (1 - omega) * old_val + omega * new_val;
            
            diff = abs(phi_SOR(j, i) - old_val);
            if diff > max_diff
                max_diff = diff;
            end
        end
    end
    iter_SOR = iter_SOR + 1;
    if max_diff < tol
        break;
    end
end
fprintf('迭代次数：%d\n', iter_SOR);

%% 5. 结果可视化
figure('Name','电位分布对比');

subplot(1,2,1);
contourf(phi_GS, 20);
colorbar;
title(sprintf('高斯-赛德尔迭代（迭代次数：%d）', iter_GS));
axis equal;

subplot(1,2,2);
contourf(phi_SOR, 20);
colorbar;
title(sprintf('超松弛迭代（ω=%.2f，迭代次数：%d）', omega, iter_SOR));
axis equal;

sgtitle('二维静电场电位分布');