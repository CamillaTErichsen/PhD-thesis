%{
calulate the plsr to predient gradient from genetic expression
Author: Morgan
Ref:    https://www.pnas.org/doi/full/10.1073/pnas.1820754116
Link:   https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md
%}

%%
% addpath(genpath('/n02dat01/users/dyli/Grad_data/support_data/Gene_results'));

% this is the expression for each gene in BNA
gene_regional_expression = load('bna_expression_15633gene_new.txt');
disp(size(gene_regional_expression));

% Note that here we use the left hemisphere only
nregs_lh=246;

X=gene_regional_expression(1:nregs_lh,:); % Predictors
% Y=horzcat(mytstat_Maast_lh,mytstat_Dublin_lh,mytstat_Cobre_lh); % Response variable
Y = load('tmap_A946_L.csv'); % the t map which size is (246,1)      change to L/R to calculate for both hemispheres

% z-score:
X=zscore(X);
Y=zscore(Y);

%%
genes=textread('gene_name.txt','%s');

geneindex=1:15633;

% number of bootstrap iterations:
bootnum=10000;   % test with 1000 iterations first, then move to 10k iterations for the final

% Do PLS in 2 dimensions (with 2 components):
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
 
% store regions' IDs and weights in descending order of weight for both components:
[R1,p1]=corr([XS(:,1),XS(:,2)],Y);
 
% align PLS components with desired direction for interpretability 
if R1(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end
 
[PLS1w,x1] = sort(stats.W(:,1),'descend');
% PLS1ids=genes(x1);
geneindex1=geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
% PLS2ids=genes(x2);
geneindex2=geneindex(x2);
 
% print out results
csvwrite('PLS1_ROIscores_minKL.csv',XS(:,1));
csvwrite('PLS2_ROIscores_minKL.csv',XS(:,2));
 
% define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];
 
%%
% start bootstrap
for i=1:bootnum
    i
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
    
    temp=stats.W(:,2);%extract PLS2 weights
    newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
end
 
% get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');
 
% get bootstrap weights
temp1=PLS1w./PLS1sw';
temp2=PLS2w./PLS2sw';
 
% order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
% PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
[Z2 ind2]=sort(temp2,'descend');
% PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);

% change the name!
save('Z1vals_10000_L.txt', 'Z1', '-ascii');
save('Z2vals_10000_L.txt', 'Z2', '-ascii');
save('ind1vals_10000_L.txt', 'geneindex1', '-ascii');
save('ind2vals_10000_L.txt', 'geneindex2', '-ascii');

% PCA
p1 = (1-normcdf(abs(Z1))) * 2;
p2 = (1-normcdf(abs(Z2))) * 2;
