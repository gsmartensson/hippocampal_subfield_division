%% Subfield division tool
% Prior to running this script, data must have been processed through
% FreeSurfer 6.0 with -hippocampal-subfields-T1 flag (see https://surfer.nmr.mgh.harvard.edu/fswiki/HippocampalSubfields)
% All subject Freesurfer output folders need to be located in the same root folder extpath:
% extpath/subj1
% extpath/subj2
% ....
%
% This script is described in "Atrophy of the Posterior Subiculum Is
% Associated with Memory Impairment, Tau-and Aβ Pathology in Non-demented
% Individuals" by Lindberg et al.
%
% Author: Gustav Mårtensson 2017
% 
%%
clc, clear all, close all

%% Input Parameters
% Path to folder with freesurfer processed data to be analyzed. Data organized as:
% extpath/subj1/mri/, extpath/subj2/mri/, etc. Script will perform subfield
% division on all subjects in this folder
extpath='/path/to/folders/';%'/home/gustav/Desktop/images/mri/subfield';

% Number of subfield divisions
N=10;

% Filenames for .csv output files.
% text_filename contains raw volumes for each segment of the hippocampal
% subfields, text_filename_prc are the contains the ratios of each segment
% compared to whole subfield.
text_filename = '/path/to/output/volume_output.csv';%'/home/gustav/Desktop/images/mri/data_volumes.csv';
text_filename_prc = '/path/to/output/volume_ratios.csv';%'/home/gustav/Desktop/images/mri/data_ratios.csv';

% Plotting variable, plot_bol = true will illustratively plot the segmented subfields.
plot_bol = false;

% Indices of subfields in the FreeSurfer LUT. Don't need to be changed.
sub_idx=[203,204,205,206,208,209,210,211,212,214,215,226];

%% Declaring LUT names, to map results to specific subfield.
homepath = pwd;

% LUT names from freesurfer
lut_name = {
    'alveus'
    'perforant_pathway'
    'parasubiculum'
    'presubiculum'
    'subiculum'
    'CA1'
    'CA2'
    'CA3'
    'CA4'
    'GC-ML-DG'
    'HATA'
    'fimbria'
    'lateral_ventricle'
    'molecular_layer_HP'
    'hippocampal_fissure'
    'entorhinal_cortex'
    'molecular_layer_subiculum'
    'Amygdala'
    'Cerebral_White_Matter'
    'Cerebral_Cortex'
    'Inf_Lat_Vent'
    'Perirhinal'
    'Cerebral_White_Matter_Edge'
    'Background'
    'Ectorhinal'
    'HP_tail'};

lut_idx = 201:226;

%% Read list of subjects in extpath, remove subjects without ../mri/lh.hippoSfLabels-T1.v10.mgz
% list all subjects in
cd(extpath)
subjlist = dir;
for i=length(subjlist):-1:1
    if subjlist(i).name(1) == '.' || subjlist(i).isdir == false
        subjlist(i) = [];
    end
end

% remove subject if hippocampal subfield files don't exist
for i=length(subjlist):-1:1
    subjname = subjlist(i).name;
    filepath=fullfile(extpath,subjname,'mri','lh.hippoSfLabels-T1.v10.mgz');
    %filepath=[extpath,'/',subjname,'/mri/','lh.hippoSfLabels-T1.v10.mgz'];
    if ~exist(filepath)
        disp([filepath, ' does not exist'])
        subjlist(i) = [];
    end
end
for i=length(subjlist):-1:1
    subjname = subjlist(i).name;
    filepath=fullfile(extpath,subjname,'mri','lh.hippoSfLabels-T1.v10.mgz');
    if ~exist(filepath)
        disp([filepath, ' does not exist'])
        subjlist(i) = [];
    end
end

cd(homepath)

%% Loop over all subjects to process
for is = 1:length(subjlist) % Loop over subjects
    
    %% Read MRI files and extract subfields
    subjname = subjlist(is).name;
    disp(['Processing images ',num2str(is),'/',num2str(length(subjlist)),': ',subjname])
    hs = {'lh','rh'};
    data = [];
    for hsi=1:length(hs)
        filepath=fullfile(extpath,subjname,'/mri/',[hs{hsi},'.hippoSfLabels-T1.v10.mgz']);
        
        % load MRI file and extract subfields specified in sub_idx
        temp = MRIread(filepath);
        resolution = prod(temp.volres); % voxel volume
        
        temp2 = zeros(size(temp.vol(:,:,round(end/2))));
        iplus = length(data);
        for i=1:length(sub_idx)
            data(i+iplus).vol = temp.vol==sub_idx(i);
            data(i+iplus).vox2mm3 = resolution;
            data(i+iplus).roi_name = [hs{hsi},'_',lut_name{lut_idx==sub_idx(i)}];
            data(i+iplus).tot_vol = sum(data(i).vol(:))*data(i).vox2mm3;
            temp2 = temp2 + data(i+iplus).vol(:,:,round(end/2))*i;
        end
    end
    
    cd(homepath)
    
    %% Create output .csv files with relevant labels
    in_vec = 1:length(data);
    
    if exist(text_filename) == 0
        fileID = fopen(text_filename,'a');
        fileID_prc = fopen(text_filename_prc,'a');
        fprintf(fileID,'%s','Subject name');
        fprintf(fileID_prc,'%s','Subject name');
        for in = in_vec
            for i = 1:N
                fprintf(fileID,'\t%s',[data(in).roi_name,'_',num2str(i)]);
                fprintf(fileID_prc,'\t%s',[data(in).roi_name,'_ratio_',num2str(i)]);
            end
            fprintf(fileID,'\t%s',[data(in).roi_name,'_total_volume']);
            fprintf(fileID_prc,'\t%s',[data(in).roi_name,'_sum_ratios']);
        end
        
        fprintf(fileID,'\n');
        fprintf(fileID_prc,'\n');
        fclose(fileID);
        fclose(fileID_prc);
    end
    
    %% loop over all subfields/FS LUT labels of a subject
    for in = in_vec
        clearvars x* y*  z* ix* cla* temp perc volumes leg
        
        fileID = fopen(text_filename,'a');
        fileID_prc = fopen(text_filename_prc,'a');
        if in ==in_vec(1)
            strformat = ['%s'];
        else
            strformat = [''];
        end
        for i = 1:(N+1)
            strformat = [strformat, '\t%6.5f'];
        end
        if in ==in_vec(end)
            strformat = [strformat, '\n'];
        end
        roi_n = in;
        
        % Extract image data of subfield to temporary variable A
        A = (data(roi_n).vol)>0;
        % get indices of subfield voxels
        [x, y, z] = ind2sub(size(A), find(A==1));
        
        % find maximum distance between two edges
        maxl=0;
        ix_vec = [0,0];
        if isempty(x) % if no subfield was segmented, stop processing
            fprintf(fileID,strformat,subjname,zeros(2*N+1,1));
            fprintf(fileID_prc,strformat,subjname,zeros(2*N+1,1));
            fclose(fileID);
            fclose(fileID_prc);
        else
            % loop over all voxel indices, find maximum distance between
            % 2 voxels (assumes cubical voxels)
            for i=1:length(x)
                tt = ((x(i)-x).^2 + (y(i)-y).^2 + (z(i)-z).^2);
                [m,ix] = max(tt);
                if m>maxl
                    ix_vec(1) = i;
                    ix_vec(2) = ix;
                    maxl = m;
                end
            end
            ix_vec = sort(ix_vec);
            xmin = x(ix_vec(1)); ymin = y(ix_vec(1)); zmin = z(ix_vec(1));
            
            tot_num = sum(A(:));
            
            % set edge subfield voxel in origo
            x = x - xmin; y = y - ymin; z = z - zmin;
            
            col_vec = {'k*','g*','c*','m*','y*'}; cl=length(col_vec);
            col_vec = {'k.','g.','c.','m.','b.','y.'}; cl=length(col_vec);
            
            % generate axis line between voxels distanced furthest from
            % each other
            xx = linspace(x(ix_vec(1)),x(ix_vec(2)),N+1);
            yy = linspace(y(ix_vec(1)),y(ix_vec(2)),N+1);
            zz = linspace(z(ix_vec(1)),z(ix_vec(2)),N+1);
            
            % check which of the N segments each voxel belongs to
            dist = [];
            linecoor = [xx;yy;zz];
            for i =length(x):-1:1
                coor = [x(i),y(i),z(i)];
                %     dist(i,j) = Inf;
                for j =1:length(xx)
                    tv = coor-linecoor(:,j)';
                    dist(i,j) = tv*tv';
                end
                [mdist,mix] = min(dist(i,:));
                if mix==1
                    class(i) = 2;
                elseif mix==N+1
                    class(i) = N+1;
                elseif dist(i,mix-1)<=dist(i,mix+1)
                    class(i) = mix;
                elseif dist(i,mix-1)>dist(i,mix+1)
                    class(i) = mix+1;
                end
            end
            
            %% Plot and calculate volumes and ratios
            if plot_bol
                hold off
                figure(77)
            end
            
            % compute volumes and ratios for each segment
            for i =1:N
                ixlist2 = (class == (i+1));
                
                perc(i) = sum(ixlist2)/tot_num;
                volumes(i) = sum(ixlist2)*data(roi_n).vox2mm3;
                leg1{i} = ['R: ', num2str(i),' - ',num2str(perc(i),2)];
                if plot_bol
                    plot3(x(ixlist2), y(ixlist2), z(ixlist2), col_vec{mod(i,cl)+1});
                    hold on
                end
                
            end
            
            titstr = [upper(data(in).roi_name(1)),data(in).roi_name(2:end)]; titstr = strrep(titstr,'_',' ');
            if plot_bol
                disp(titstr)
                title(titstr,'Interpreter','latex')
                grid on
                plot3(xx,yy,zz,'r-o','LineWidth',2)
                plot3(xx(1),yy(1),zz(1),'k-o','LineWidth',2)
                xlabel('x','Interpreter','latex'),ylabel('y','Interpreter','latex'),zlabel('z','Interpreter','latex')
                ll = [-1,1];
                axis equal
                xlim(xlim + ll), ylim(ylim + ll),zlim(zlim + ll)
                drawnow
            end
            
        end
        %% write ratios and volumes to files
        if true
            if in ==in_vec(1)
                fprintf(fileID,strformat,subjname,volumes,tot_num*data(roi_n).vox2mm3);
                fprintf(fileID_prc,strformat,subjname,perc,sum(perc));
            else
                fprintf(fileID,strformat,volumes,tot_num*data(roi_n).vox2mm3);
                fprintf(fileID_prc,strformat,perc,sum(perc));
            end
        end
        fclose(fileID);
        fclose(fileID_prc);
    end
end
disp('Done!')