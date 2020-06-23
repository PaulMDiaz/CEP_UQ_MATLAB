function F1 =  eigenvectors(W, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Plot the estimated eigenvectors for the active subspace analysis with
%   optional bootstrap ranges.
%
%   Inputs:
%          W: m-by-k array that contains k of the estimated eigenvectors
%             from the active subspace analysis
%          W_br: (optional) m-by-(2*k) array that contains the lower and 
%                upper bounds for the components of the eigenvectors
%          opts: (optional) structure array which contain plotting options
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin)
    W_br = [];
    opts = [];
elseif length(varargin) == 1
    if isnumeric(varargin{1})
        W_br = varargin{1};
        opts = [];
    elseif isstruct(varargin{1})
        W_br = [];
        opts = varargin{1};
    else
        error('ERROR: Inappropriate inputs passed')
    end
elseif length(varargin) == 2
    if isnumeric(varargin{1}) && isstruct(varargin{2})
        W_br = varargin{1};
        opts = varargin{2};
    elseif isstruct(varargin{1}) && isnumeric(varargin{2})
        opts = varargin{1};
        W_br = varargin{2};
    else
        error('ERROR: Inappropriate inputs passed')
    end
else
    error('ERROR: Too many inputs.')
end

[m, n] = size(W);

% Determine layout of subplots based on n.
subplot_vert  = [1, 1, 1, 2, 2, 2];
subplot_horiz = [1, 2, 3, 2, 3, 3];
if n > 6
    n = 6;
    warning('WARNING: Only plotting first 6 eigenvectors.')
end
subplot_vert  = subplot_vert(n);
subplot_horiz = subplot_horiz(n);

% Get plotting options.
opts = plot_opts(opts);

%FF  = figure()

F1 = [];
for i = 1:1
    %subplot(subplot_vert, subplot_horiz, i)
    
    % Plot bootstrap errors if provided.
    if ~isempty(W_br) && (max(abs(W_br(:, 2*(i-1)+2) - W_br(:, 2*(i-1)+1))) > 1e-7)
        fill([1:1:m, m:-1:1], [W_br(:, 2*(i-1)+1)', fliplr(W_br(:, 2*(i-1)+2)')], opts.err_color)
        
    end
    
    % Plot eigenvector.
    F1 =  plot(W, ...
         'markeredgecolor', 'k', ...
         'marker', opts.marker, ...
         'markersize', opts.markersize, ...
         'linewidth', opts.linewidth)
    %set(F1, 'markerfacecolor', get(F1, 'color'));
    %hold on
    if isempty(opts.xlabel)
    %xlabel('Index', 'fontsize', opts.fontsize)
    xlabel('i', 'fontsize', opts.fontsize,'interpreter','latex')
    else
        xlabel(opts.xlabel, 'fontsize', opts.fontsize)
    end
    if isempty(opts.ylabel)
        ylabel('$w_{j,i}$', 'fontsize', opts.fontsize,'interpreter','latex')
    else
        ylabel(opts.ylabel, 'fontsize', opts.fontsize)
    end

  

    grid('on')
    grid on;
    set(gca, ...
        'XLim', [1, m], ...
        'XTick', 1:m, ...
        'XScale', 'Linear', ...
        'YScale', 'linear', ...
        'TickLabelInterpreter','latex',...
        'fontsize', opts.fontsize)
    if ~isempty(opts.xticklabel)
        set(gca, 'XTickLabel', opts.xticklabel)
    end
end
L = legend('$\mathbf{w}_1$','$\mathbf{w}_2$','location','best')
set(L,'interpreter','latex','fontsize',opts.fontsize);
end