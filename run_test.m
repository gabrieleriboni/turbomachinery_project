m_max = 1;
main_def;
save('main_workspace.mat');

clear;
m_max = 1;
sensitivity_analysis;
save('sens_workspace.mat');

clear;
M = load('main_workspace.mat');
S = load('sens_workspace.mat');

fields = fieldnames(M);
for i = 1:length(fields)
    f = fields{i};
    if isfield(S, f)
        val_m = M.(f);
        val_s = S.(f);
        if isnumeric(val_m) && isnumeric(val_s) && numel(val_m) == 1 && numel(val_s) == 1
            if abs(val_m - val_s) > 1e-6
                fprintf('DIFFERENCE IN %s: main = %.8f, sens = %.8f\n', f, val_m, val_s);
            end
        end
    end
end
