function save_data_to_csv(data, genotype, ID, varName, savecsvpath)

T = array2table(data, 'variableNames', varName);
T.Genotype = genotype;
T.AnimalID = ID;
T = movevars(T, 'Genotype', 'Before', 'Session1');
T = movevars(T, 'AnimalID', 'Before', 'Genotype');

writetable(T, savecsvpath);

end