function [pic50_sigma, hill_spread] = fit_spread_parameters(pic50s, varargin)
% Code to examine the spread of values that 
% a given ion channel assay produces. Input should be
% repeated measurements of the same positive control.
% Example here is hERG IonWorks Quattro Cisapride screen at AZ.
% Data were published with Elkins et al. (2013) J. Pharmacological and
% Toxicological Methods.
%
% First argument is a vector of pIC50s.
% Second argument is a matching vector of Hill coefficients (optional).
%
pic50_sigma = 0;
hill_spread= 0;
close all

if nargin == 1
    dataset.values = {pic50s};
elseif nargin == 2
    dataset.values = {pic50s, varargin{1}};
else
    error('Only two args allowed') 
end
            
num_drugs = 1;

x = linspace(0,9,1000)';

for i=1:length(dataset.values)
    figure
    h = histogram(dataset.values{i});  
    ylabel('Frequency','FontSize',16)
    if i <= num_drugs
        xlim([3 9])
        xlabel('pIC50 (-log_{10} IC50(M))','FontSize',16)
        set(gca,'FontSize',14)
        fit_type = 'logistic';    
    else
        xlim([0 5.05])
        xlabel('Hill coefficient','FontSize',16)
        set(gca,'FontSize',14)
        fit_type = 'loglogistic';
    end
    hold on    
    
    temporary = mle(dataset.values{i},'distribution',fit_type);
    mu(i) = temporary(1);
    sigma(i) = temporary(2);   
    
    bin_width = (max(dataset.values{i})-min(dataset.values{i}))/h.NumBins;
    prob_func(i,:) = pdf(fit_type,x,mu(i),sigma(i));   
    y(i,:) = prob_func(i,:)*bin_width*length(dataset.values{i});
       
    figure(i)
    hold on
    plot(x,y(i,:),'r-','LineWidth',2)
   
    % Print out in latex format...

    if i <= num_drugs
        %fprintf('%s & %s & %s & %s & %i & %s & %4.4g & %4.4g \\\\ \n', temp{1}, temp{2}, temp{3},dataset.source{i},length(dataset.values{i}),fit_type,mu(i), sigma(i))
%         fprintf('Parameter &  N & $\\mu$ & $\\sigma$ & $\\alpha$ & $1/\\beta$ \\\\ \n')
%         fprintf('pIC50 Spread & %i & %4.4g & %4.4g \\\\ \n', length(dataset.values{i}),mu(i),sigma(i))
        pic50_sigma = sigma(i);
    else
%         fprintf('Parameter &  N & $\\alpha$ & $1/\\beta$ \\\\ \n')
%         fprintf('Hill Spread & %i & %4.4g & %4.4g \\\\ \n', length(dataset.values{i}),exp(mu(i)), sigma(i))
        hill_spread = sigma(i);
    end
end

