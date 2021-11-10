% dots reach lifetime
l_out = find(rand(ndots, 1) < f_kill);

if any(l_out)
    
    nout_l = length(l_out);
    
    % choose new coordinates
    xy(l_out, :) = [rect(3) * rand(nout_l, 1) rect(4) * rand(nout_l, 1)];
    
    t = 2 * pi * rand(nout_l, 1);
    %                 index = randperm(size(t, 1));
    %                 index = index(:,1: round(size(t, 1) * Coherence));
    %
    %                 if length(index) > 0
    %                     t(index) = angle;
    %                 end
    cs = [cos(t), sin(t)];
    dxdy(l_out, :) = [pfs(l_out) pfs(l_out)] .* cs;
    
end

% check to see which dots have gone beyond the borders of the annuli
x_out = find(xy(:, 1) < 0 | xy(:, 1) > rect(3));	% dots to reposition
y_out = find(xy(:, 2) < 0 | xy(:, 2) > rect(4));

if any(x_out)
    
    nout_x = length(x_out);
    
    % choose new coordinates
    xy(x_out, :) = [rect(3) * rand(nout_x, 1) rect(4) * rand(nout_x, 1)];
    
    t = 2 * pi * rand(nout_x, 1);
    %                 index = randperm(size(t, 1));
    %                 index = index(:,1: round(size(t, 1) * Coherence));
    %                 t(index) = angle;
    cs = [cos(t), sin(t)];
    dxdy(x_out, :) = [pfs(x_out) pfs(x_out)] .* cs;
    
end

if any(y_out)
    
    nout_y = length(y_out);
    
    % choose new coordinates
    xy(y_out, :) = [rect(3) * rand(nout_y, 1) rect(4) * rand(nout_y, 1)];
    
    t = 2 * pi * rand(nout_y, 1);
    %                 index = randperm(size(t, 1));
    %                 index = index(:,1: round(size(t, 1) * Coherence));
    %                 t(index) = angle;
    cs = [cos(t), sin(t)];
    dxdy(y_out, :) = [pfs(y_out) pfs(y_out)] .* cs;
    
end