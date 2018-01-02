clear
close all
clc

%% load the data

% The data will be in a 'sparse matrix' format
tmp = load('ml_1m.mat'); 
Xsp = tmp.X;

%visualize the non-zero elements in Xsp (if there is dot in the plot, there is a
%non-zero element)
spy(Xsp);

%we can obtain the mask matrix from Xsp
M = double(Xsp==1);

clear tmp %release the variable from the memory

[I,J] = size(Xsp);
N = nnz(Xsp); %number of non-zeros

fprintf('Only %2.4f percent of the matrix is full!\n',100*N/(I*J));

%create a list of the elements, which will be useful later
%Xi, Xj : coordonnee
%Xs : value
[Xi,Xj,Xs] = find(Xsp);
Xlist = [Xi,Xj,Xs];

clear Xi
clear Xj
clear Xs


%% initialize the factor matrices

K = 10; %set the rank

W = 1 * randn(I,K); 
H = 1 * randn(K,J);

%% stochastic gradient descent


batchSize = 500; %the number of elements that we will use at each iteration

eta = 0.03; %step-size
numIter = 100;

rmse_sgd = zeros(numIter,1);

for t = 1:numIter
    
    %get a random batch from the data
    % N nb of non zero
    data_index = randperm(N,batchSize); %todo
    %this vector will contain some random numbers between 1 and N. 
    %Its size should be batchSize. You can use the function randperm
    
    %for each element in the data batch, update the corresponding elements
    %in W and H
    for i = 1:batchSize
       % for each element in the batch, find its corresponding 'i', 'j',
       % and value by using the Xlist array
       % data_index[i] : range 
       % 1,2,3 : colone 
       cur_i = Xlist(data_index(i),1);%todo
       cur_j = Xlist(data_index(i),2);%todo
       cur_x = Xlist(data_index(i),3);%todo
       
       %compute the current xhat, for the current i and j
       cur_xhat = W(cur_i,:) * H(:,cur_j);%todo
       
       %compute the gradients for the 'corresponding elements' of W and H
       %not all the elements of W and H will be updated
       %grad_w = (X - WH) * H'
       %grad_h = W' * (X - WH)
       grad_w = - (cur_x - cur_xhat) * H(:, cur_j)'; %todo
       grad_h = - W(cur_i, :)' * (cur_x - cur_xhat); %todo
       
       %take a gradient step
       %cur_i : actual coordonee of x-axis
       %cur_y : actual coordonee of y-axis
       %W = W - learningrate * grad_w
       %H = H - learningrate * grad_h
       W(cur_i,:) = W(cur_i,:) - eta * grad_w ;%todo
       H(:,cur_j) = H(:,cur_j) - eta * grad_h ;%todo
    end
    
    %compute the root-mean-squared error
    %rmse_sgd(t) = sqrt((norm(M.( - W * H)).^2)/N)%todo
    rmse_sgd(t) = sqrt(sum(sum(((M.*(W*H))-Xsp).^2))/N);
    %norm(M.*(W*H)-Xsp,'fro')/sqrt(N)
    fprintf('Iteration %d\n',t);
end

plot(rmse_sgd);

%% Make recommendations

user_index = 11;

%use W, H, and M to compute a 'movie_index' for the user

user_xhat = (1-M(:,user_index)).*(W*H(:,user_index));
[xx,movie_index] = max(user_xhat);%todo

fprintf('Recommend movie %d to user %d\n',movie_index,user_index);



