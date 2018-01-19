function[data,cov_images,params, meta] = generate_data(meta,X,SNR,TR,sim_hemodynamics,outfile,varargin)
%GENERATE_DATA  generate synthetic brain data
%
%USAGE:
%   [data, ground_truth, params, meta] = generate_data(meta,X,SNR,[TR],[sim_hemodynamics],[outfile],[{niiargs}]);
%
%INPUTS:
%
%     meta: (1) a struct with the following fields:
%        nvoxels: total number of voxels containing brain
%     coordToCol: dimx by dimy by dimz matrix of voxel numbers (zeros
%                 indicate no voxel at the corresponding location)
%     colToCoord: nvoxels by D matrix of voxel locations, where D is the
%                 number of dimensions (usually D = 3)
%
%           (2) alternatively, the user may pass in a vector of dimension
%               lengths (e.g. [10 20 30]).  if so, a meta matrix will be
%               automatically constructed for the corresponding
%               (hyper)cube.
%
%        X: design matrix (ntrials by ncovariates) containing the
%           activation of each covariate on each trial.
%
%      SNR: desired signal-to-noise ratio of the synthetic data.
%
%       TR: temporal resolution of images, in seconds.  only used if
%           sim_hemodynamics is set to TRUE.  this should be a positive
%           scalar.  (OPTIONAL.)
%
% sim_hemo: a boolean flag specifying whether the design matrix should be
%           convolved with the hemodynamic response function to more
%           accurately mimic fMRI data across successive images.
%           (OPTIONAL; default is false.)
%
%  outfile: path to store the resulting images in NIFTI format.  this can
%           be useful for testing out a data processing pipeline,
%           particularly with sim_hemo = true.
%
%  niiargs: arguments to pass to make_nii (function for making NIFTI
%           images).  see MAKE_NII for details.
%
%
%OUTPUTS:
%
%         data: a 1 by ntrials cell array.  each cell contains an nvoxels by 1
%               vector of voxel activations
%
% ground_truth: an nvoxels by ncovariates matrix containing the true images
%               for each covariate (i.e. radial basis functions).
%
%       params: a struct with the following fields specifying source
%               parameters for the radial basis functions:
%               centers: a size(X,2) by D matrix of source center locations
%                widths: a size(X,2)-dimensional vector of source widths
%               weights: a copy of X (specifies source weights), convolved
%                        with the HRF if applicable.
%
%         meta: the meta struct.  if the meta struct has been constructed
%               automatically, it can be useful to return it for use in
%               subsequent computations outside of this function.
%
%EXAMPLE USAGE:
%   
%   %TOY EXAMPLES:
%   %1 (unique) source active per trial; 5 trials total; high SNR
%   X1 = eye(5);
%   [data1,ground_truth1] = generate_data(meta,X1,1e5);
%
%   %5 sources, with each source randomly active during each trial; 10
%   %trials total; medium SNR
%   X2 = rand(10,5);
%   [data2,ground_truth2] = generate_data(meta,X2,1e3);
%
%   %2 sources active on alternating trials; 10 trials total; low SNR
%   X3 = repmat(eye(2),5,1);
%   [data3,ground_truth3] = generate_data(meta,X3,1);
%
%   %2 sources, active on alternating trials; 10 trials total; low SNR; 2
%   second TR; convolve with hemodynamic response function and generate
%   NIFTI images.
%   generate_data(meta,X3,1,2,true,'synthetic_images.nii');
%
%   %"REAL" EXAMPLE:
%   %2500 sources, with a unique set of 100 active on each trial.  high
%   %SNR; 4 second TR; convolve with hemodynamic response function and
%   generate NIFTI images.
%   X4 = repmat(eye(25),1,100);
%   generate_data(meta,X4,1e4,4,true,'synthetic_images.nii');
%
%SEE ALSO: MAKE_NII, CONV


% 4/10/12   JRM     wrote it.
% 5/15/12   JRM     support SNR = 0
% 2/26/13   JRM     added optional ability to convolve design matrix with
%                   hemodynamic response function.  also added ability to
%                   output images in NIFTI format.
% 3/7/13    JRM     output source parameters, minor tweaks & fixes.
% 6/30/13   JRM     fixed bug in hemodynamic convolution.
% 7/10/13   JRM     add ability to construct meta matrix automatically

if exist('sim_hemodynamics','var') && sim_hemodynamics
    assert(exist('TR', 'var'), 'Must specify TR duration.')
    assert(isscalar(TR) && TR > 0,'TR must be a positive scalar.');
    X = hrf_conv(X,TR);
end

if isnumeric(meta), meta = construct_meta(meta); end

%create images (RBF's) for each covariate
[cov_images,params.centers,params.widths] = get_cov_images(meta,size(X,2));
params.weights = X;
data = repmat({zeros(1,meta.nvoxels)},1,size(X,1));

%create noisy image for each trial
for i = 1:size(X,1)    
    data{i} = X(i,:)*cov_images'; %each image is a weighted sum of the covariate images, where the weightings are given by X(i,:)
    if SNR <= eps(0)
        data{i} = randn(size(data{i}));
    else
        data{i} = data{i} + (mean(data{i})/SNR)*randn(size(data{i})); %SNR = data_mean/noise_std; so noise_std = data_mean/SNR
    end    
end

if exist('outfile','var')
    assert(ischar(outfile),'output file must be a string');
    
    images = zeros([size(meta.coordToCol) length(data)]);
    for i = 1:length(data)
        images(:,:,:,i) = get_brain_image(data{i},meta);
    end
    
    nii = make_nii(images, varargin{:});
    save_nii(nii, outfile);
end


function[img,centers,widths] = get_cov_images(meta,n)
centers = zeros(n,size(meta.colToCoord,2));
widths = zeros(1,n);
img = zeros(meta.nvoxels,n);
for i = 1:n
    centers(i,:) = meta.colToCoord(randi(meta.nvoxels),:); %pick random voxel as center
    widths(i) = log(min(size(meta.coordToCol))*rand/4);
    dists = pdist2(meta.colToCoord,centers(i,:));
    img(:,i) = exp(-dists./exp(widths(i)));
end


function[C] = hrf_conv(X,TR)
%define hrf 
%(from http://kendrickkay.net/GLMdenoise/doc/GLMdenoise/utilities/getcanonicalhrf.html)
hrf = [0 0.0314738742235483 0.132892311247317 0.312329209862644 0.441154423620173 0.506326320948033 0.465005683404153 0.339291735120426 0.189653785392583 0.0887497190889423 0.0269546540274463 -0.00399259325523179 -0.024627314416849 -0.0476309054781231 -0.0550487226952204 -0.0533213710918957 -0.0543354934559645 -0.053251015547776 -0.0504861257190311 -0.0523878049128595 -0.0480250705100501 -0.0413692129609857 -0.0386230204112975 -0.0309582779400724 -0.0293100898508089 -0.0267610584328128 -0.0231531738458546 -0.0248940860170463 -0.0256090744971939 -0.0245258893783331 -0.0221593630969677 -0.0188920336851537 -0.0205456587473883 -0.0230804062250214 -0.0255724832493459 -0.0200646133809936 -0.0101145804661655 -0.014559191655812];

%resample HRF according to TR
if TR < 1
    hrf = interp(hrf,1/TR);
elseif TR > 1
    %round TR to nearest tenth-second        
    hrf = resample(hrf,10,round(TR*10));
end

C = zeros(size(X));
for i = 1:size(X,2)
    next = conv(X(:,i)', hrf, 'full');
    C(:,i) = next(1:end-(length(hrf)-1));
end

function[img] = get_brain_image(x,meta)
img = meta.coordToCol;
img(img == 0) = nan;
for i = 1:meta.nvoxels
    img(meta.colToCoord(i,1),meta.colToCoord(i,2),meta.colToCoord(i,3)) = x(i);
end

function[meta] = construct_meta(dims)
if size(dims, 1) > 1
    colToCoord = dims;
    dims = cellfun(@(x)(length(unique(x))), slices(dims, 2));    
else
    colToCoord = fullfact(dims);
end
meta = struct('nvoxels', prod(dims), 'coordToCol', zeros(dims), 'colToCoord', colToCoord);




