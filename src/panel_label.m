function panel_label(varargin)
% PANEL_LABEL  Place a bold panel letter at the top-left of an axes.
%
%   Inputs:
%     ax   — axes handle (optional; defaults to gca).
%     str  — panel-label text, e.g. '(b)'.
    if nargin == 1, ax = gca; str = varargin{1};
    else,           ax = varargin{1}; str = varargin{2}; end
    text(ax, 0.03, 0.97, str, 'Units','normalized', 'FontName','Arial', ...
         'FontSize',14, 'FontWeight','bold', 'VerticalAlignment','top');
end
