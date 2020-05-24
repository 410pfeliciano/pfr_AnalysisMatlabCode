function savingdata(LFP, tetChannel, savePath)
        prompt = 'Assign a new variable name to your data ? ';
        xx = input(prompt, 's');
        assignin('base',xx,LFP) % assignin the new variable
        fprintf(1, 'Saving Files!\n');
        file_name1 = sprintf('LFPch#%d.mat',tetChannel);
        save ([savePath filesep file_name1], xx, 'time') % saving spike data
        % save ([savePath filesep file_name1], LFP, 'time','-v7.3') % saving spike data
        fprintf(1, 'Done!\n');
    end