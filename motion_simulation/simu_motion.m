ccc
%% simulation parameters
ratio = 1;
flip_angle = 15; % degree
TR = 2; % [ms]
M0 = 1;
heart_rate = 60:10:120; % simulated heart rate [min^-1]
max_disp = 0.5:0.1:2; % simulated maximum myocardium displacement, normalizde by number of slice thickness
T1 = 10:10:3000; % simulated T1 range [ms]

flip_angle = flip_angle*ratio;

%% load slice profile
load('./SliceFast.mat')
Slice_profile = Slice.Profile;

%% simulation with motion
for iheart_rate = heart_rate
    for imax_disp = max_disp 
        
        flip_angle_slice = flip_angle.*Slice_profile;
        
        Mz = ones(3800, 1, 'single');
        Mz_PD = ones(3800, 1, 'single');
        fprintf(sprintf('start simulation for:\nheart rate = %g/min\nmaxium displacement = %g slice thickness\n',iheart_rate, imax_disp))
        tic
        SI_with_motion = zeros(9, 1000, 300, 'single');
        parfor i = 1:length(T1)
%             fprintf([num2str(i),'\n'])
            simu_T1 = T1(i);
            slab_T1 = ones(3800,1)*simu_T1;
            SI_with_motion(:,:,i) = simu_steady_state_with_motion(flip_angle_slice, Mz, M0, TR, slab_T1, imax_disp, iheart_rate);
        end
        toc

        save_name = sprintf('HR_%g_max_disp_%g_TR_%g_FA_%g.mat', iheart_rate, imax_disp, TR, flip_angle);
        save(save_name,'SI_with_motion','T1');
        fprintf(sprintf('done\n'))
    end
end

%% without motion
fprintf('start simulation without motion:')
tic
SI_no_motion = zeros(9, 1000, 300, 'single');
parfor i = 1:length(T1)
    simu_T1 = T1(i);
    slab_T1 = ones(3800,1)*simu_T1;
    SI_no_motion(:,:,i) = simu_steady_state(flip_angle_slice, Mz, M0, TR, slab_T1);
end
toc
save('SI_without_motion','SI_no_motion','T1');
fprintf('all finished')