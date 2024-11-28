QTC_Coeff_file = fopen("A2_Q2_Zout_QTC_Coeff.txt", 'r');


lines = readlines('A2_Q1_Zout_QTC_Coeff.txt');
lineCount = length(lines);
disp(['Total number of lines: ', num2str(lineCount)]);


lines = readlines('A2_Q1_Zout_MDiff.txt');
lineCount = length(lines);
disp(['Total number of lines: ', num2str(lineCount)]);
