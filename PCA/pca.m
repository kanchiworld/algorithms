# GNU Octave and MATLAB script for PCA demo
# runas as
# octave
# octave_prompt> pca
#
#

X = [ 2.5; 0.5; 2.2; 1.9; 3.1; 2.3; 2; 1; 1.5; 1.1 ]
Y = [2.4; 0.7; 2.9; 2.2; 3; 2.7; 1.6; 1.1; 1.6; 0.9 ]

figure(1)
plotmatrix(X, Y, "*")

# cov(X,Y) calculates covariance of X vector and Y vecotr, so the result
# is a single number.
#cov(X,Y)                                                                                                                                                                                
#help eig                                                                                                                                                                                

# cov ([X(:), Y(:)]) always results in a 2x2 matrix
covXY = cov ([X(:), Y(:)])                                                                                                                                                                  
lmda = eig(covXY)                                                                                                                                                                           
[V, lmda] = eig(covXY)                   

#If you need to rearrange V from highest eigenvalue to lowest
#W = [V(:,2),V(:,1)]

Xm = X - mean(X)
Ym = Y - mean(Y)

XmYm = [ Xm, Ym ]

ComponentData = V' * XmYm'
#PC1 = ComponentData(2, :)
PC1 = ComponentData(2:end, :) # All row except row 1, all columns

#V1 = V(:, 2)
V1 = V(:, 2:end) # All rows, all columns except column 1

ReducedMAData = V1 * PC1
ReducedX = ReducedMAData(1, :) + mean(X)
ReducedY = ReducedMAData(1, :) + mean(Y)

figure(2)
plotmatrix(ReducedX', ReducedY', "*")

#To get the original Mean adjusted data
OrigMAData = V * ComponentData



