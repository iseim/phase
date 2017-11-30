% polymer analysis (on position arrays saved from .gsd files using Tom's script whingdingdilly and
% 'savemat' in scipy.io)

%parameters
num_RNA = 100;
len_RNA = 33;
num_protein = floor(num_RNA * 21.7);
len_protein = 8;
box_size = 22.6 * num_RNA^(1/3);
vect = attract2;

%RNA data
for i = 1:num_RNA
    for j = 1:length(vect)
        for k = 1:len_RNA
            ind = (i-1)*len_RNA + k;
            RNA.(['RNA', num2str(i)]).(['Particle', num2str(k)])(j,1) = vect{1,j}(ind,1);
            RNA.(['RNA', num2str(i)]).(['Particle', num2str(k)])(j,2) = vect{1,j}(ind,2);
            RNA.(['RNA', num2str(i)]).(['Particle', num2str(k)])(j,3) = vect{1,j}(ind,3);
        end
    end
end

%calculate msd for each particle in each polymer
for i = 1:num_RNA
    for j = 1:len_RNA
        [msdRNA.(['RNA', num2str(i)])(j,:), lags] = polymer_msd(RNA.(['RNA', num2str(i)]).(['Particle', num2str(j)]) ...
            , box_size, 1);
    end
end


%protein data
for i = 1:num_protein
    for j = 1:length(vect)
        for k = 1:len_protein
            ind = (i-1)*len_protein + k + num_RNA*len_RNA;
            protein.(['protein', num2str(i)]).(['Particle', num2str(k)])(j,1) = vect{1,j}(ind,1);
            protein.(['protein', num2str(i)]).(['Particle', num2str(k)])(j,2) = vect{1,j}(ind,2);
            protein.(['protein', num2str(i)]).(['Particle', num2str(k)])(j,3) = vect{1,j}(ind,3);
        end
    end
end

%calculate msd for each particle in each polymer
for i = 1:num_protein
    for j = 1:len_protein
        [msdPRO.(['protein', num2str(i)])(j,:), lags] = polymer_msd(protein.(['protein', num2str(i)]).(['Particle', num2str(j)]) ...
            , box_size, 1);
    end
end


            
            
            
            
            
            
            
            
            
            
            