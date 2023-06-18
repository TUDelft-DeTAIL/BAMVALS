function PhiW = compute_PhiW(u,TN,Vm,Vp,r)
    % PhiW = compute_PhiW(u,TN)
    % --------------------------
    % Computes projection matrix Phi*W from a truncated MIMO Voltera Series in Tensor Network (TN) format
    % given inputs u(:,1),u(:,2),....
    %
    % PhiW      =   matrix, y(:,:) contains contraction over all inputs and first d-1  cores
    %
    % u         =   matrix, u(:,k) contains the kth input,
    %
    % TN        =   Tensor Network.
    %
    % Reference
    % ---------
    %
    % 10/10-2016, Eva Memmel

    % generate N x (pM+1) matrix U with input samples
    [N,p]=size(u);
    D=length(TN.core);
    J=TN.n(1,3); %J = p*M+1
    L=TN.n(1,2); %number of outputs
    M=(J-1)/p; %memory
    U=zeros(N,J);
    u=[zeros(M-1,p);u];
    for i=M:N+M-1            
        temp=ones(1,J);
        for j=1:M
            temp(2+(j-1)*p:2+j*p-1)=u(i-j+1,:);                
        end   
        U(i-M+1,:)=temp;
    end
    if L==1
        PhiW=dotkron(Vm,U,Vp);
    elseif sweepindex == 1
        PhiW=kron(dotkron(U,Vp),Vm);
    else
        PhiW=dotkron(Vm,U,Vp);
        PhiW=reshape(PhiW,[N,L,r(1)*n*r(2)]);
        PhiW=permute(PhiW,[2 1 3]);
        PhiW=reshape(PhiW,[N*L,r(1)*n*r(2)]);
    end

end

