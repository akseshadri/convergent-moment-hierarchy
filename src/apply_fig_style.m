function apply_fig_style(ax)
% APPLY_FIG_STYLE  Apply the repository's consistent axes style.
%
%   APPLY_FIG_STYLE(ax) sets Arial fonts (12 pt ticks, 14 pt labels, 12 pt
%   legend with no box), a single light major grid (no minor grid), and box
%   on, every figure shares the same format. 
%
%   Input:
%     ax  — axes handle (defaults to gca).
    if nargin < 1, ax = gca; end
    set(ax, 'FontName','Arial', 'FontSize',12, 'Box','on', 'Layer','top', ...
            'XGrid','on', 'YGrid','on', 'XMinorGrid','off', 'YMinorGrid','off', ...
            'GridAlpha',0.15, 'LineWidth',0.75);
    set(get(ax,'XLabel'), 'FontName','Arial', 'FontSize',14);
    set(get(ax,'YLabel'), 'FontName','Arial', 'FontSize',14);
    lg = get(ax,'Legend');
    if ~isempty(lg), set(lg, 'FontName','Arial', 'FontSize',12, 'Box','off'); end
end
