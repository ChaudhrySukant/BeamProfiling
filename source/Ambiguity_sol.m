function [D_solved_shortest, D_solved_avrg, D_solved, N_solved, N_length] = Ambiguity_sol(Phi,Lambda,Range_Max)
%Calculates estimated absolute distance (D_solved_shortest, D_solved_avrg) on a multi-frequency PoA measurement
%using the estimated phases vector (Phi) at modulation wavelengths vector (Lambda)
%via exaustive search (number of full cycles (N_solved) that provide higher agreement between estimated distance pairs (D_1 VS D_i)
%within measuremet range (Range_Max)
%Solving is always by comparison with first wavelengtht: D_solved & N_solved provide pairs solutions [1 i]
%D_solved_shortest is the distance estimation from the shortest wavelength and D_solved_avrg is calculated as average of estimated absolute distances D_solved (both only computed when there is consistency in first wavelength full cycles between pairs, -1 otherwise)

N=cell(size(Lambda));
N_length=ones(size(Lambda));
D=cell(size(Lambda));

% Compute all possible distances per wavelength within Range_Max
for i=1:length(Lambda)
    N{i}=0:1:floor(Range_Max/(Lambda(i)/2));
    N_length(i)=length(N{i});    
    D{i}=Lambda(i)/2*Phi(i)/(2*pi)+N{i}*Lambda(i)/2;
end

% Compute number of cycles per wavelenth that provides closer agreement
% between pairs
D_difs=cell(1,length(Lambda)-1);
Min_dif=cell(1,length(Lambda)-1);
N_solved=ones(length(Lambda)-1,2);
D_solved=ones(length(Lambda)-1,2);
for i=1:length(Lambda)-1
    D_difs{i}=ones(length(N{1}),length(N{2}));
    for k=1:length(N{1})
        for m=1:length(N{i+1})
            D_difs{i}(k,m)=abs(D{1}(k)-D{i+1}(m));
        end
    end
    Min_dif{i}=min(min(D_difs{i}));
    [index_1, index_2, ~]=find(D_difs{i}==Min_dif{i}(1));
    N_solved(i,:)=[index_1(1)-1, index_2(1)-1];
    D_solved(i,:)=[Lambda(1)/2*Phi(1)/(2*pi)+N_solved(i,1)*Lambda(1)/2, Lambda(i+1)/2*Phi(i+1)/(2*pi)+N_solved(i,2)*Lambda(i+1)/2];
end

% Check for consistency between solutions from different wavelength pairs
if nnz(N_solved(:,1)-N_solved(1,1))>0
    warning('Inconsistent solution'); %Different number of cycles for first wavelength
    D_solved_avrg=-1;
    D_solved_shortest=-1;
else
    D_solved_avrg=mean([D_solved(1,1); D_solved(:,2)]);
    D_solved_shortest=D_solved(Lambda==min(Lambda),Lambda==min(Lambda));
end

end