clear all
%rng(4)
models = {'sms_33', 'ms_33', 'BIOMD0000000201', 'BIOMD0000000202', 'BIOMD0000000212',...
    'BIOMD0000000016', 'BIOMD0000000207', 'Calvin_Ray', 'BIOMD0000000012',...
    'BIOMD0000000219', 'BIOMD0000000239'};

cputimes = zeros(2, numel(models));
cnflenght = [];
for nn=1:numel(models)
    model = models{nn};
    
    InputName = strcat('C:\Users\sedaghat\Dropbox (Personal)\MONET_MetabolicNetworks\Data\RealMetabolicNetworks\',model, '_binEFM.mat');
    
    load(InputName)
    mat=full(sparseEFM);
    mat=mat'; %
    
    mat = Preprocessing( mat );
    mat = Irredundant(mat);
    
    [FK_cnf, FK_cpu, cnf_len]= FK_Dualization( mat );
    %Berge input is mat = [EFM x Reactions]
    [berge_cnf, totaliter, maxiter, berge_cpu]= berge(mat);
    
    if (size(berge_cnf,1) ~= size(FK_cnf))
        display('Something is wrong!');
        stop
    elseif size(intersect(berge_cnf, FK_cnf, 'rows'),1)==size(berge_cnf, 1)
        display('Hooray!! FK algorithm works well!!');
    else
        display('Something is worng!!');
    end
    
    disp(['Berge computation time: ',num2str(berge_cpu),' sec']);
    disp(['FK computation time: ',num2str(FK_cpu),' sec']);
    
    cputimes(1, nn) = FK_cpu;
    cnflenght(nn,:) = cnf_len;
    cputimes(2, nn) = berge_cpu;
end
%
% Coli_AC_comp
% Coli_AC_full
% Coli_glc_comp
% Coli_glc_full

% save(OutputName, 'res', 'totaliter', 'maxiter', 'processtime', '-v7.3')
for nn=1:(numel(models)-1)
%     figure
    plot(cnflenght(nn, 1:(find(cnflenght(nn,:)==0, 1, 'first')-1)), 'b--o');
    name = strrep(models(nn), '_', '-');
    title(name)
    xlabel('Iterations')
    ylabel('Length of CNF')
    saveas(gcf, char(strcat(name, '.pdf')))
end