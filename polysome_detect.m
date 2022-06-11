%% Polysome_detect

% This script is to detect polysomes in cellular tomograms based on ribosome positions and orientations determined during sub-tomogram refinement.
% define polysomes based on the distance between the mRNA exit site of one ribosome to the entry site of another ribosome.

% Written by Liang Xue
% Last update March 2021 

% The script was written in MATLAB 2016b
% Prerequisite: put the two associated fromEuler_RELION.m and rotmat2eulang.m files in MATLAB search Path.





%% Preparing a motl table contain all relevant information 
% The polysome detect script is run with a motivelist(motl) table that combined information(particle coordinates, shifts and rotations, classification, etc.) from previous refinement and classification runs.
% Usually, this motl table is converted from RELION refinement data star or equivalent data.
% In the motl table, each column is to store all information of a particle/ribosome. Each row is to store one item. 
% Here, by default, it is a table of 20 rows and each is defined as below:
%        1         : whether this ribosome belongs to an expressome. Leave empty if no applicable
%        2         : 
%        3         : 
%        4         : particle sequential number to track each particle in the tomogram. Usually adopted from particle's filename 
%        5         : tomoNum. Unique number identifier for each tomogram
%        6         : polysome identifier. In each tomogram, ribosomes belong to the same polysome will have a same non-zero number. 
%        7         : The relative ranking of this ribosome within the polysome it belongs. Note: the number may not be continuous, but always the sequence is from small to large
%        8         : x-coordinate in full unbin tomogram, pixels
%        9         : y-coordinate in full unbin tomogram, pixels
%        10        : z-coordinate in full unbin tomogram, pixels
%        11        : calculated mRNA entry site:  x in full unbin tomogram, pixels
%        12        : calculated mRNA entry site:  y in full unbin tomogram, pixels
%        13        : calculated mRNA entry site:  z in full unbin tomogram, pixels
%        14        : calculated mRNA exit site:  x in full unbin tomogram, pixels
%        15        : calculated mRNA exit site:  y in full unbin tomogram, pixels
%        16        : calculated mRNA exit site:  z in full unbin tomogram, pixels
%        17        : rot in RELION data star. aound z axis
%        18        : tilt in RELION data star.  aound new Y axis
%        19        : psi in RELION data star. aound new z axis
%        20        : ribosome class identifiers to store previous classification results on translation states. Leave empty if not available

% For a dataset of N ribosome particles, tne table will be 20 x N.





%% Prepare a motl table from RELION
% Here it requires the refinement(average) of the ribosome sub-tomograms has been performed, with both data star file (RELION-3 format) and average map available.
% Map back the ribosomes into tomogram and visualize them in 3D if possible (e.g. TOM toolbox/tom_paste2.m or https://github.com/builab/subtomo2Chimera)

% A template motl of 70S ribosomes in Mycoplasma cellular tomograms is provided alongside the script.
% skip this part if using the provided motl_template

% read data star from RELION refinement
fid1 = fopen('run_data.star','r'); 
dataArray1 = textscan(fid1,'%s %f %f %f %f %f %f %s %s %d %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',26); %format might differ
fclose(fid1);

% read ribosome classification results if possible
cid1 = fopen('jobxxx_class001.star','r'); 
riboclass1 = textscan(cid1,'%s %f %f %f %f %f %f %s %[^\n]','headerlines',26); %format might differ
fclose(cid1);

cid2 = fopen('jobxxx_class002.star','r'); 
riboclass2 = textscan(cid2,'%s %f %f %f %f %f %f %s %[^\n]','headerlines',26); %format might differ
fclose(cid2);

cid4 = fopen('jobxxx_class003.star','r'); 
riboclass3 = textscan(cid3,'%s %f %f %f %f %f %f %s %[^\n]','headerlines',26); %format might differ
fclose(cid3);

% ... more class data star if there are


% localize the relative positions of the mRNA entry and exit sites in the ribosome average map
% These sites needs to be carefully and somehow empirically defined, based on known structures and the average map(change the density threshold in Chimera and check where the mRNA density orignates). 
% One possible way: open average in Chimera; open volume eraser and move to the right positions; command line type getcrd; open replylog to get coordinates xyz(in Angstrom) 

% below are based on the Mycoplasma pneumoniae ribosome average
boxsize=[256; 256; 256]; % the local coordinates depends on box size; update
pixelsize=1.7005; % A/pixel; Update!
% mRNA entry check
localxyz_mRNAentry=([128.721; 241.009; 236.059])/pixelsize; %note if these number were directly copied from Chimera log,they are in A. Here should convert to pixels
offset_mRNAentry = localxyz_mRNAentry-(boxsize+1)/2;  % calculate offset relative to the box center
% mRNA exit check
localxyz_mRNAexit=([186.995; 248.341; 300])/pixelsize; %note if these number were directly copied from Chimera log,they are in A. Here should convert to pixels
offset_mRNAexit = localxyz_mRNAexit-(boxsize+1)/2;  % calculate offset relative to the box center



% combine and convert RELION results to motl  
motl_all=zeros(20, size(dataArray1{1,1},1));  

for j = 1:size(dataArray1{1,8},1)
    
    % relative local particle seq from RELION data star
	rpn_dirty=char(dataArray1{1,8}(j,1));   %column number may change
    motl_all(4,j)=str2double(rpn_dirty(58:63)); %extract position may change. 
    
    % get tomonum from RELION data star 
	TN_dirty=regexp(dataArray1{1,1}(j,1), '\d*','Match'); %note the format might change
	motl_all(5,j)= str2double(TN_dirty{1,1});
        
	% rotations still in RELION/xmipp system
    motl_all(17:19,j) = [dataArray1{1,11}(j,1), dataArray1{1,12}(j,1), dataArray1{1,13}(j,1)];
        
	% re-centered xyz. center of the ribosome
	motl_all(8:10,j) = [(dataArray1{1,2}(j,1) - dataArray1{1,14}(j,1)), (dataArray1{1,3}(j,1) - dataArray1{1,15}(j,1)), (dataArray1{1,4}(j,1) - dataArray1{1,16}(j,1))]; %column number may change
    
    % rlnAngleRot, rlnAngleTilt, rlnAnglePsi to calculate the Rotation matrix
    R=fromEuler_RELION(dataArray1{1,11}(j,1), dataArray1{1,12}(j,1), dataArray1{1,13}(j,1)); %column number may change
    
    % calculate the mRNA entry site for this ribosome. xyz in full tomgram
    motl_all(11:13,j) = [(dataArray1{1,2}(j,1) - dataArray1{1,14}(j,1)); (dataArray1{1,3}(j,1) - dataArray1{1,15}(j,1)); (dataArray1{1,4}(j,1) - dataArray1{1,16}(j,1))] + R * offset_mRNAentry;
    
    % calculate the mRNA exit site for this ribosome. xyz in full tomgram
    motl_all(14:16,j) = [(dataArray1{1,2}(j,1) - dataArray1{1,14}(j,1)); (dataArray1{1,3}(j,1) - dataArray1{1,15}(j,1)); (dataArray1{1,4}(j,1) - dataArray1{1,16}(j,1))] + R * offset_mRNAexit;
    
    % annote row 20 (or other rows) based on classification results(if applicable)
    %1
    if ismember(char(dataArray1{1,8}(j,1)),riboclass1{1,8}(:,1))
        motl_all(20,j)=1;    
    end
    %2
    if ismember(char(dataArray1{1,8}(j,1)),riboclass2{1,8}(:,1))
        motl_all(20,j)=2;    
    end
    %3
    if ismember(char(dataArray1{1,8}(j,1)),riboclass3{1,8}(:,1))
        motl_all(20,j)=3;    
    end
    
    % other class annotation if there
    
end


dlmwrite('motl_xxx.txt', motl_all); 






% Either use the provided motl_template.txt file to run below script 
% or use the the newly generated motl_xxx.txt (double check if its format is same to the motl_tenplate)


clear;
%% processing_add_polysome_information

motl=dlmread('motl_template.txt'); % The motl_template.txt file provided along with the script

% split motl by tomoNum
% The polysome detection will be run one cellular tomogram after one tomogram 
% each tomoNum represents one tomogram

tomoNum=unique(motl(5,:));
mkdir 'motl_TomoNum';
for i=1:size(tomoNum,2)
    
    motl_tn=motl(:,motl(5,:)==tomoNum(1,i)); 
    dlmwrite(['./motl_TomoNum/motl_' num2str(tomoNum(1,i)) '.txt'], motl_tn)

end



% The polysome detect will be run one tomogram after one tomogram 
mkdir 'motl_annoted_addpolysome_TomoNum'
cd 'motl_annoted_addpolysome_TomoNum'

infotrack=zeros(6,size(tomoNum,2)); % to store some basic information for quick checking

for t=1:size(tomoNum,2)
    
    disp(['Read motl_annoted of tomogram ' num2str(tomoNum(1,t))]);
    motl=dlmread(['../motl_TomoNum/motl_' num2str(tomoNum(1,t)) '.txt']); 
    polyid=0; % make sure this polyid is initialized. 0:mono-ribosome; >0: belongs to polysome.  
 
    PN=size(motl,2);
    for p=1:PN
        for p2=1:PN
            if p~=p2
                distance=norm(motl(14:16,p)-motl(11:13,p2));
                % Note!! this distance is the most important criterion to determine polysomes. 
                % distance <41 pixels(here pixel size is 0.17nm, so ~7nm) 
                % Note! this distance is related to how the mRNA entry and exit sites are defined, and needs to be tested according to different data.
                % Carefully choose the distance criterion by 3D ribosome mapping, distance distribution plotting, etc.
                if distance < 41  % in pixels
                    %case 1
                    if (motl(6,p)>0) && (motl(6,p2)>0)
                        motl(7,motl(6,:)==motl(6,p2)) = motl(7,motl(6,:)==motl(6,p2)) + motl(7,p);
                        motl(6,motl(6,:)==motl(6,p2)) = motl(6,p);  
                    %case 2
                    elseif (motl(6,p)==0) && (motl(6,p2)>0)
                            motl(6,p)=motl(6,p2);
                            motl(7,motl(6,:)==motl(6,p2)) = motl(7,motl(6,:)==motl(6,p2)) +1;                         
                    %case 3
                    elseif (motl(6,p)>0) && (motl(6,p2)==0) 
                            motl(6,p2)=motl(6,p);
                            motl(7,p2)=motl(7,p)+1;
                    %case 4
                    elseif (motl(6,p)==0) && (motl(6,p2)==0) 
                            polyid=polyid+1;
                            motl(6,p)=polyid;
                            motl(6,p2)=polyid;
                            motl(7,p)=1;
                            motl(7,p2)=2; 
                    end
                  
                end 
            end 
        end 
    end 
    
%   % some quick numbers to see how the script works
    infotrack(1,t)=tomoNum(1,t);  % which tomogram
    infotrack(2,t)=PN; % how many ribosomes in this tomogram
    infotrack(3,t)=sum(motl(6,:)>0); % how many ribosomes annoated as polysomes
    infotrack(4,t)=sum(motl(6,:)>0)/PN; % the percentage of polysomes against all ribosomes
    infotrack(5,t)=size(unique(motl(6,:)),2)-1; % how many polysomes
    infotrack(6,t)=max(motl(7,:)); % the longest polysomes
    
    dlmwrite(['motl_annoted_addpolysome_' num2str(tomoNum(1,t)) '.txt'], motl);
    
end

% Important notes:
% Map the detected polysomes in 3D to how the detection goes
% Tune the "distance criterion" if possible. 
% Smaller distances -> limited to polysomes or patrs of the polysomes that assemble very tightly .
% Large distances -> higher possibility to have wrong annotations(the ribosomes just randomly come close)
% The distance-based polysome annotation is highly related to how the mRNA entry site and exit site are defined. So carefully select these positions.
% Be extremely careful when comparing statistics between different datasets. 
% This script has not been tested for membrane-associated ribosomes, but in principal should also work. 


% combine each tomogram's annotated motl into one motl_annoated_addpolysome if needed
motl_all=[];
for t=1:size(tomoNum,2)
    motl_TN=dlmread(['./motl_annoted_addpolysome_' num2str(tomoNum(1,t)) '.txt']); 
    motl_all=[motl_all motl_TN];
end
    
dlmwrite('motl_annoted_addpolysome.txt', motl_all);















clear;
%% Some basic plotting in MATLAB


motl=dlmread('motl_annoted_addpolysome.txt'); 
tomoNum=unique(motl(5,:));


%% This is to visualize the overall distribution of neighbouring ribosomes.

nn_xyz=[];  % to store the normalized neighbour coordinates xyz matrix; rotated
nn_entry2exit=[]; 

for t=1:size(tomoNum,2)
    
    motl=dlmread(['./motl_annoted_addpolysome_' num2str(tomoNum(1,t)) '.txt']); 
    
%     motl=motl(:,motl(6,:)>0);   % If not comment out, then only plot the annoated polysomes

    for p=1:size(motl,2)
	
		motl_nop=motl; motl_nop(:,p)=[];  %make a motl_nop that only lacks this particle-p
		motl_rest2p=motl_nop(8:10,:) - motl(8:10,p);  %distance of all of rest particles to this particle-p
        motl_restentry2pexit=motl_nop(11:13,:) - motl(14:16,p); 
		% d_all=vecnorm(motl_rest2p); d_mink=mink(d_all, 8); % select the 8 smallest distance
		    
        % rlnAngleRot, rlnAngleTilt, rlnAnglePsi to calculate Rotationmatrix
        R=fromEuler_RELION(motl(17,p), motl(18,p), motl(19,p));
    
		for r=1:size(motl_rest2p,2)
			if norm(motl_rest2p(1:3, r)) < 300   %pixels; for the motl template, roughly a shell constraint of 1.5x 70S ribosome; Update!
                neighbour_rotated=R\motl_rest2p(1:3,r);
				nn_xyz=[nn_xyz neighbour_rotated]; 
                neighbour_rotated2=R\motl_restentry2pexit(1:3,r);
				nn_entry2exit=[nn_entry2exit neighbour_rotated2]; 
                
            end
        end
    end
end

ppc=size(nn_xyz,2);

% before plotting, make sure handedness is correct; also convert to nm 

% plotting coordinates, center to center
scatter3(nn_xyz(1,1:ppc)*0.17005,nn_xyz(2,1:ppc)*0.17005,nn_xyz(3,1:ppc)*0.17005,'.') ;  %0.17005nm is the pixel of template data; update
% rotate the scatter3 plotting, and should be able to see two clusters of points(each point = the center of one neighbouring ribosome).


% plotting entry to exit distance distribution
d=[];
for m=1:ppc
    d=[d (norm(nn_entry2exit(1:3,m)))*0.17005]; %0.17005nm is the pixel of template data; update
end
d=d(d<20); % only show those smaller than 20nm; change if needed
histogram(d,40, 'BinWidth', 0.5, 'FaceColor', [0.5 0.5 0.5], 'Normalization','count')





clear;
%% This is to plot the neighbours based on polysome sequence information

motl=dlmread('motl_annoted_addpolysome.txt'); 
tomoNum=unique(motl(5,:));

nn_xyz=[];  %normalized neighbour coordinates xyz matrix used to store found neighbours for all ribosome. rotated
nn_followecenter2exit_xyz = [];
nn_exitvector=[];
precursor_rot=[];
follower_rot=[];


ddata_center2center=[];  %distance between (mass) center of two neighbouring ribosomes
ddata_exit2entry=[];  %distance between exit and entry of two neighbouring ribosomes
offset_entry = ([128.721; 241.009; 236.059]/1.7005) - ([256; 256; 256]+1)/2;  % Update entry site xyz, pixel size and box size
offset_exit = ([186.995; 248.341; 300]/1.7005) - ([256; 256; 256]+1)/2;  % Update exit site xyz, pixel size and box size


for t=1:size(tomoNum,2)
    
    motl=dlmread(['./motl_annoted_addpolysome_' num2str(tomoNum(1,t)) '.txt']); 
 
    motl=motl(:,motl(6,:)>0);   %if not comment out, then only consider polysomes

    for p=1:size(motl,2) 
        
        % rlnAngleRot, rlnAngleTilt, rlnAnglePsi to calculate Rotationmatrix
        R=fromEuler_RELION(motl(17,p), motl(18,p), motl(19,p));
        
        motl_pl=motl(:,motl(6,:)==motl(6,p));
        precursor=motl_pl(:, motl_pl(7,:)==(motl(7,p)-1));
        follower=motl_pl(:, motl_pl(7,:)==(motl(7,p)+1));
        
        if ~isempty(follower)
            ddata_center2center=[ddata_center2center norm(follower(8:10,1)-motl(8:10,p))];
            v_exit2entry=follower(8:10,1)+(fromEuler_RELION(follower(17,1),follower(18,1),follower(19,1))*offset_entry)-motl(8:10,p)-(R*offset_exit);
            ddata_exit2entry=[ddata_exit2entry norm(v_exit2entry)]; 
            
            neighbour_rotated=R\(follower(8:10,1)-motl(8:10,p)); % rotated to normalize
            nn_xyz=[nn_xyz neighbour_rotated]; 
            
            newcalculated_exit = motl(8:10,p)+ R * offset_exit;
            neighbour2exit_rotated=R\(follower(8:10,1)-newcalculated_exit);
            nn_followecenter2exit_xyz=[nn_followecenter2exit_xyz neighbour2exit_rotated]; 
            %calculate an extra parameter of this follower_ vector from mass center to mRNA exit site
            neighbour_exit_rotated=R\(follower(14:16,1)-follower(8:10,1));
            nn_exitvector=[nn_exitvector neighbour_exit_rotated];
            %calculate an extra parameter of this follower_its rotation vector
            %The order of rotation angles is x-axis, y-axis, z-axis. Also called extrinsic rotations
            follower_rot_rotated=rotmat2eulang(R\(fromEuler_RELION(follower(17,1),follower(18,1), follower(19,1))), 'XYZ')'; %extrinsic, psi, Theta, phi
            follower_rot=[follower_rot follower_rot_rotated];
                
        end
    end
end

ppc=size(nn_xyz,2);

%3D scatter plot
scatter3(nn_xyz(1,1:ppc),nn_xyz(2,1:ppc),nn_xyz(3,1:ppc),'.') ;

% plot the relative position nn_xyz in x-y, x-z, and y-z
subplot(1,3,1);
plot(nn_xyz(1,1:ppc),nn_xyz(2,1:ppc), '.');
subplot(1,3,2);
plot(nn_xyz(1,1:ppc),nn_xyz(3,1:ppc), '.');
subplot(1,3,3);
plot(nn_xyz(2,1:ppc),nn_xyz(3,1:ppc), '.');


% plot rotation (euler angle pairs) of the following ribosomes 
subplot(1,3,1);
plot(follower_rot(1,1:ppc),follower_rot(2,1:ppc), '.');
subplot(1,3,2);
plot(follower_rot(2,1:ppc),follower_rot(3,1:ppc), '.');
subplot(1,3,3);
plot(follower_rot(1,1:ppc),follower_rot(3,1:ppc), '.');

% do K mean cluster of the rotation of following ribosomes based on euler angles
data_forK=follower_rot;
[idx,C] = kmeans(data_forK,2);
% color the euler angle plot to see if make sense
c1=[0.9 0.7 1];
c2=[0 1 0];
subplot(1,3,1);
plot(follower_rot(1,idx==1),follower_rot(2,idx==1), '.', 'MarkerFaceColor', c1);
hold on
plot(follower_rot(1,idx==2),follower_rot(2,idx==2), '.', 'MarkerFaceColor', c2);
subplot(1,3,2);
plot(follower_rot(2,idx==1),follower_rot(3,idx==1), '.', 'MarkerFaceColor', c1);
hold on
plot(follower_rot(2,idx==2),follower_rot(3,idx==2), '.', 'MarkerFaceColor', c2);
subplot(1,3,3);
plot(follower_rot(1,idx==1),follower_rot(3,idx==1), '.', 'MarkerEdgeColor', c2);
hold on
plot(follower_rot(1,idx==2),follower_rot(3,idx==2), '.', 'MarkerEdgeColor', c1);
hold off



% Name the two clusters as t-t and t-b. Usually the larger cluster is t-t. Check in 3D if needed.
% different runs may switch the order
t_b=sum(idx==2);
t_t=sum(idx==1);


% use clustering based on rotations to color relative positions (nn_xyz) of the following ribosomes
data_toshow=nn_xyz*1.7005;
figure;
sz=20;
c1=[0.9 0.7 1];
scatter3(data_toshow(1,idx==1), data_toshow(2,idx==1), -data_toshow(3,idx==1),sz, c1, 'filled');
hold on
c2=[0 1 0];
scatter3(data_toshow(1,idx==2), data_toshow(2,idx==2), -data_toshow(3,idx==2), sz, c2, 'filled');

% also plot the mRNA eentry and exit site of the preceding ribosome if needed
scatter3(offset_entry(1)*1.7005, offset_entry(2)*1.7005,-offset_entry(3)*1.7005, 90, [0 0 1], 'filled'); % emtry site, big blue dot
scatter3(offset_exit(1)*1.7005, offset_exit(2)*1.7005,-offset_exit(3)*1.7005, 90, [0.8 0.8 0.8], 'filled'); % exit site, big grey dot


figure
subplot(1,3,1)
plot(data_toshow(1,idx==1), data_toshow(2,idx==1),'c.');
hold on 
plot(data_toshow(1,idx==2), data_toshow(2,idx==2),'b.');
subplot(1,3,2)
plot(data_toshow(1,idx==1), -data_toshow(3,idx==1),'c.');
hold on 
plot(data_toshow(1,idx==2), -data_toshow(3,idx==2),'b.');
subplot(1,3,3)
plot(data_toshow(2,idx==1), -data_toshow(3,idx==1),'c.');
hold on 
plot(data_toshow(2,idx==2), -data_toshow(3,idx==2),'b.');
hold off


% % manually check in raw tomograms to see if t-t vs t-b seperation based on K clustering is working
% motl_follower_tb=follower_fullinfo(:,idx==1);
% % specify a tomoNum and check each after each
% tomoNum=667;
% mft=motl_follower_tb(:,motl_follower_tb(5,:)==tomoNum);  
% % put a sphere at detected polysomes' centers in bin8 tomograms
% vol=zeros(464,464,225);
% sphere=tom_spheremask(ones(25, 25, 25), 7, 2);
% for j = 1:size(mft,2)
%     xyz = floor([mft(8,j) mft(9,j) mft(10,j)]/8 - size(sphere)./2 + 1);
%     xyz = round(xyz);
%     vol= tom_paste2(vol, sphere, xyz, 'max');
% end
% tom_mrcwrite(vol, 'name', ['check_tb_' num2str(tomoNum) '_bin8.mrc']);
% 
% c=uisetcolor;



% %% calculate some distance parameters 
% mean_center2center=mean(ddata_center2center) * 1.7005;  %in A  =207
% std_center2center=std(ddata_center2center) * 1.7005;  %in A  =18
% mean_exit2entry=mean(ddata_exit2entry) * 1.7005;  %in A  =45
% std_exit2entry=std(ddata_exit2entry) * 1.7005;  %in A  =15
% 
% % for only t-t
% mean_center2center_tt=mean(ddata_center2center(1,idx==1)) * 1.7005;  %in A  =204
% std_center2center_tt=std(ddata_center2center(1,idx==1)) * 1.7005;  %in A  =17
% mean_exit2entry_tt=mean(ddata_exit2entry(1,idx==1)) * 1.7005;  %in A  =42
% std_exit2entry_tt=std(ddata_exit2entry(1,idx==1)) * 1.7005;  %in A  =14
% 
% % for only t-b
% mean_center2center_tb=mean(ddata_center2center(1,idx==2)) * 1.7005;  %in A  =219
% std_center2center_tb=std(ddata_center2center(1,idx==2)) * 1.7005;  %in A  =17
% mean_exit2entry_tb=mean(ddata_exit2entry(1,idx==2)) * 1.7005;  %in A  =54
% std_exit2entry_tb=std(ddata_exit2entry(1,idx==2)) * 1.7005;  %in A  =15





% Written by Liang Xue, EMBL 2020
