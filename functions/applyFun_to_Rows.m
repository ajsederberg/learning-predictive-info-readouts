function y_out = applyFun_to_Rows(data, fn, which_argout)

y_out = cell(size(data, 1), 1);

for i_row = 1:size(data, 1)

    if which_argout == 2
        [~, y_out{i_row}] = fn(data(i_row, :));
    else
        y_out{i_row} = fn(data(i_row, :));
    end
end

% try y_out = cell2mat(y_out);
% catch
%     display('returning cell')
% end  