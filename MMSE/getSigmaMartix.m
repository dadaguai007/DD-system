function [Rxxd,distanceEigenDiag,Vd]=getSigmaMartix(dividerArray,M)
% M 为符号数
% 获取协方差矩阵,数据的特征值，数据的特征向量（降序排列）



Rxxd = dividerArray * dividerArray' / M;%得到估计的协方差矩阵
% 计算协方差矩阵Rxxd的特征分解（特征向量矩阵，特征值对角矩阵）
[Vd, distanceEigen] = eig(Rxxd);
% 提取特征值并转换为行向量
distanceEigenDiag = diag(distanceEigen)';
% 对特征值向量进行升序排序
[distanceEigenDiag, eigenMark] = sort(distanceEigenDiag);
% 将升序结果翻转为降序，（最大特征值在前）
distanceEigenDiag = fliplr(distanceEigenDiag);
% 使特征向量与特征值降序对应
Vd = fliplr(Vd(:, eigenMark));


end