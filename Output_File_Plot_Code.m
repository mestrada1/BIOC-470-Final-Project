uiimport() % Loaded output_file.txt as table
output_file = table2struct(outputfile); % Converted table to structure array


% Plotting p-values for first 15 genes
Genes = {};
p_values = [];
for i = 1:15
    Genes{i} = output_file(i).gene;
    p_values(i) = output_file(i).p; 
end

x = 1:15;
plot(x, p_values, 'bo')
set(gca, 'XTick', 1:15, 'XTickLabel', Genes)
xlabel('Gene');
ylabel('p-value');
title('Mutation Significance for Given Genes');

% Plotting q-values for first 15 genes
q_values = [];
for i = 1:15
    q_values(i) = output_file(i).q;
end

x = 1:15;
plot(x, q_values, 'bo')
set(gca, 'XTick', 1:15, 'XTickLabel', Genes)
xlabel('Gene');
ylabel('q_value');
title('False Discovery Rate for Given Genes');