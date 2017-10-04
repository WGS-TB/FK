clear all
%rng(4)
models = {'BIOMD0000000048','BIOMD0000000228', 'sms_33', 'ms_33', ...
    'BIOMD0000000034', 'BIOMD0000000042', 'BIOMD0000000201', ...
'BIOMD0000000406', 'BIOMD0000000030', 'BIOMD0000000028', 'BIOMD0000000029',...
'BIOMD0000000202', 'BIOMD0000000212', 'BIOMD0000000016',...
'BIOMD0000000207', 'Calvin_Ray', 'BIOMD0000000044','BIOMD0000000045', ...
'BIOMD0000000227', 'BIOMD0000000041','BIOMD0000000036','BIOMD0000000039',... 
'BIOMD0000000043','BIOMD0000000235', 'BIOMD0000000012'};%'ms_44', 'BIOMD0000000404',

cputimes = zeros(numel(models),3);
cnflenght = [];
mcnflenght = [];
Imcnflenght = [];
nn=1;
for nn=1:numel(models)
    model = models{nn};
    InputName = strcat('C:\Users\sedaghat\Dropbox (Personal)\MONET_MetabolicNetworks\Data\RealMetabolicNetworks\',model, '_binEFM.mat');
    
    load(InputName)
    mat=full(sparseEFM);
    mat=mat'; %
    
    mat = Preprocessing( mat );
    mat = Irredundant(mat);
    
    %Berge input is mat = [EFM x Reactions]
    [berge_cnf, totaliter, maxiter, berge_cpu]= berge(mat);
    
    %     rng(4)
    [FK_cnf, FK_cpu, cnf_len]= FK_Dualization( mat );
    
    [MFK_cnf, MFK_cpu, mcnf_len]= MFK_Dualization( mat );
    
    [IMFK_cnf, IMFK_cpu, Imcnf_len]= Improved_MFK_Dualization( mat );
    
    
    if (size(berge_cnf,1) ~= size(FK_cnf))
        display('Something is wrong!');
        stop
    elseif size(intersect(berge_cnf, FK_cnf, 'rows'),1)==size(berge_cnf, 1)
        display('Hooray!! FK algorithm works well!!');
    else
        display('Something is worng!!');
        stop
    end
    
    
    if (size(berge_cnf,1) ~= size(MFK_cnf))
        display('Something is wrong!');
        stop
    elseif size(intersect(berge_cnf, MFK_cnf, 'rows'),1)==size(berge_cnf, 1)
        display('Hooray!! MFK algorithm works well!!');
    else
        display('Something is worng!!');
        stop
    end
    
    if (size(berge_cnf,1) ~= size(IMFK_cnf))
        display('Something is wrong!');
        stop
    elseif size(intersect(berge_cnf, IMFK_cnf, 'rows'),1)==size(berge_cnf, 1)
        display('Hooray!! IMFK algorithm works well!!');
    else
        display('Something is worng!!');
        stop
    end
    
    disp(['Berge computation time: ', num2str(berge_cpu),' sec']);
    disp(['FK computation time: ', num2str(FK_cpu),' sec']);
    disp(['MFK computation time: ', num2str(MFK_cpu),' sec']);
    disp(['IMFK computation time: ', num2str(IMFK_cpu),' sec']);
    
    cputimes(nn,1) = berge_cpu;
    cputimes(nn,2) = FK_cpu;
    cputimes(nn,3) = MFK_cpu;
    cputimes(nn,4) = IMFK_cpu;
    
    cnflenght(nn,:) = cnf_len;
    mcnflenght(nn,:) = mcnf_len;
    Imcnflenght(nn,:) = Imcnf_len;
    
end
%
% Coli_AC_comp
% Coli_AC_full
% Coli_glc_comp
% Coli_glc_full

% save(OutputName, 'res', 'totaliter', 'maxiter', 'processtime', '-v7.3')
for nn=1:(numel(models))
    if(max(cnflenght(nn,:)) > 1)
        figure
        plot(cnflenght(nn, 1:(find(cnflenght(nn,:)==0, 1, 'first')-1)), 'k--o',...
            'MarkerSize',4,...
            'MarkerFaceColor','k');
        hold on;
        plot(mcnflenght(nn, 1:(find(mcnflenght(nn,:)==0, 1, 'first')-1)), 'b--x',...
            'MarkerSize',4,...
            'MarkerFaceColor','b');
        hold on;
        plot(Imcnflenght(nn, 1:(find(Imcnflenght(nn,:)==0, 1, 'first')-1)), 'r--*',...
            'MarkerSize',4,...
            'MarkerFaceColor','r');
        legend('FK','MFK','IMFK','Location','northwest')
        
        name = strrep(models(nn), '_', '-');
        title({char(name);
            ['(FK, MFK, IMFK) time = (' num2str(cputimes(nn,2)), ', ',num2str(cputimes(nn,3)), ', ', num2str(cputimes(nn,4)), ')'];
            })
        
        xlabel('Iterations')
        ylabel('Length of CNF')
        xlim([1 inf])
        ylim([1 inf])
        hold off
        saveas(gcf, char(strcat('C:\Users\sedaghat\Dropbox (Personal)\MONET_MetabolicNetworks\Results\', name, '.pdf')))
    end
end

close all
