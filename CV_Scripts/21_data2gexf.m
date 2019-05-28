function [frequency_interp] = data2gexf(ninterp, threshold, scale, fout_name)
%data2gexf.m
%This function parses Jia's data to produce a .gexf file, which we use in
%Gephi to make an image and a movie.
%Sample call: data2gexf(11, 0, 0, 'test.gexf');
%
%Assumptions:
%       (1) edge2dom has been compiled to e2d and is in the same directory 
%       as this function.
%       (2) the input data are in 2019cv_genotype_map_v3/
%
% input: 
%        ninterp: number of interpolation points between each generation.
%
%        threshold: minumum frequency required for a genotype to be
%                   included in the genotype network
%
%        scale: Boolean. If 1, then freq = 100.*log10(freq+1), else freq
%               unchanged
%
%        fout_name: name of .gexf output file (e.g., out.gexf)
%
%Author: Joshua L. Payne%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine the set of unique genotypes that appear in any of the replicates 
%of either the V0 or Vc populations. use a hash for this
dict = java.util.Hashtable;
numseqs = 1;
for i = 1:8
    
    if i == 1 %V0, rep 1
        fid = fopen('2019cv_genotype_map_v3/pc1_genotypeMap_cv2019_v3.txt');
    elseif i == 2 %V0, rep 2
        fid = fopen('2019cv_genotype_map_v3/pc2_genotypeMap_cv2019_v3.txt');
    elseif i == 3 %V0, rep 3
        fid = fopen('2019cv_genotype_map_v3/pc3_genotypeMap_cv2019_v3.txt');
    elseif i == 4 %V0, rep 4
        fid = fopen('2019cv_genotype_map_v3/pc4_genotypeMap_cv2019_v3.txt');
    elseif i == 5 %Vc, rep 1
        fid = fopen('2019cv_genotype_map_v3/pm1_genotypeMap_cv2019_v3.txt');
    elseif i == 6 %Vc, rep 2
        fid = fopen('2019cv_genotype_map_v3/pm2_genotypeMap_cv2019_v3.txt');
    elseif i == 7 %Vc, rep 3
        fid = fopen('2019cv_genotype_map_v3/pm3_genotypeMap_cv2019_v3.txt');
    elseif i == 8 %Vc, rep 4
        fid = fopen('2019cv_genotype_map_v3/pm4_genotypeMap_cv2019_v3.txt');
    end
        
    %strip off the first two lines
    line = fgets(fid);
    line = fgets(fid);
    
    %loop over all lines in .txt file
    while ~feof(fid)
    
        %get the next line
        line = fgets(fid);
        
        %get the 2nd token, the amino acid sequence
        [~, remain] = strtok(line);
        [seq, remain] = strtok(remain);
    
        %have we seen this sequence before? If not, add it to the hash
        if isempty(dict.get(seq))
            dict.put(seq, numseqs);
            seqlist{numseqs} = seq;
            numseqs = numseqs + 1;
        end
        
    end
    fclose(fid);
    
end
numseqs = numseqs - 1; %because matlab

%build an adjacency matrix and write an edge list to a .txt file
adj = zeros(numseqs);
fid = fopen('edges.txt','w');
numedges = 1;
for i = 1:numseqs
    
    %get sequence i
    seq1 = seqlist{i};
    
    %loop over the other sequences
    for j = (i+1):numseqs
        
        %get sequence i
        seq2 = seqlist{j};
        
        %calculate the edit distance between the two seqiences
        edit_dist = sum(seq1 ~= seq2);
        
        %if the edit distance is 1, add an edge
        if edit_dist == 1
           adj(i,j) = 1;
           edges(numedges, 1) = i;
           edges(numedges, 2) = j;
           fprintf(fid,'%d\n%d\n', i, j);
           numedges = numedges + 1;
        end
        
    end
    
end
numedges = numedges - 1; %because Matlab
fclose(fid);

%determine which nodes are in the dominant genotype network
%(task outsourced to igraph in C, see edge2dom.c)
eval(['!./e2d ' num2str(numedges) ' edges.txt edges_dom.txt nodes_dom.txt']);
nodes_dom = load('nodes_dom.txt');
nodes_dom_mask = zeros(numseqs,1);
nodes_dom_mask(nodes_dom) = 1;
edges_dom = load('edges_dom.txt');
edges_dom = edges_dom + 1; %because Matlab

fprintf('Total number of nodes %d, number of nodes in dominant network %d\n', numseqs, length(nodes_dom));

%get the frequency of each sequence in each replicate of the V0 and Vc
%populations
frequency = zeros(numseqs,9,8); %numseqs x numgens x numreps
for i = 1:8
    
    if i == 1 %V0, rep 1
        fid = fopen('2019cv_genotype_map_v3/pc1_genotypeMap_cv2019_v3.txt');
    elseif i == 2 %V0, rep 2
        fid = fopen('2019cv_genotype_map_v3/pc2_genotypeMap_cv2019_v3.txt');
    elseif i == 3 %V0, rep 3
        fid = fopen('2019cv_genotype_map_v3/pc3_genotypeMap_cv2019_v3.txt');
    elseif i == 4 %V0, rep 4
        fid = fopen('2019cv_genotype_map_v3/pc4_genotypeMap_cv2019_v3.txt');
    elseif i == 5 %Vc, rep 1
        fid = fopen('2019cv_genotype_map_v3/pm1_genotypeMap_cv2019_v3.txt');
    elseif i == 6 %Vc, rep 2
        fid = fopen('2019cv_genotype_map_v3/pm2_genotypeMap_cv2019_v3.txt');
    elseif i == 7 %Vc, rep 3
        fid = fopen('2019cv_genotype_map_v3/pm3_genotypeMap_cv2019_v3.txt');
    elseif i == 8 %Vc, rep 4
        fid = fopen('2019cv_genotype_map_v3/pm4_genotypeMap_cv2019_v3.txt');
    end
        
    %strip off the first two lines
    line = fgets(fid);
    line = fgets(fid);
    
    %loop over all lines in the .txt file
    while ~feof(fid)
    
        %get the line
        line = fgets(fid);
        
        %get the generation number
        [gen, remain] = strtok(line);
        gen = str2num(gen);
       
        %get the sequence
        [seq, remain] = strtok(remain);
        index = dict.get(seq);
        
        %get the frequency
        freq = strtok(remain);
        freq = str2num(freq);
        
        if scale
            frequency(index,gen+1,i) = 100*log10(freq+1);
        else
            frequency(index,gen+1,i) = freq; %numseqs x numgens x numreps
        end
        
    end
    
    fclose(fid); 
    
end

%Jia didn't explicitly write the frequency of the ancestral sequence (100%)
%in generations 0-4 of V0, so I write it here.
frequency(1, 1:5, 1:4) = 100; %numseqs x numgens x numreps

%interpolate frequencies between generations
if ninterp > 0

    %initialize the new frequency matrix (numseqs x (numgens - 1)*ninterp x numreps)
    frequency_interp = zeros(numseqs, 8*ninterp - 8, 8) + nan;

    %for replicate i
    for i = 1:8
       
        %for sequence j
        for j = 1:numseqs
            
            %for gen k to gen k+1
            for k = 1:8
            
                if k == 1
                    low = 1;
                else
                    temp = frequency_interp(j, :, i);
                    low = min(find(isnan(temp))) - 1;
                end
                
                high = low + ninterp - 1;
                interp = linspace(frequency(j,k,i), frequency(j,k+1,i), ninterp);

                frequency_interp(j, low:high, i) = interp;
                
            end
            
        end

    end
    
else
    
    frequency_interp = frequency;
    
end

%open the gexf file
fid = fopen(fout_name,'w');

%write the preamble of the gexf
fprintf(fid, '%s\n','<gexf xmlns="http://www.gexf.net/1.1draft"');
fprintf(fid, '%s\n','xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"');
fprintf(fid, '%s\n','xsi:schemaLocation="http://www.gexf.net/1.1draft');
fprintf(fid, '%s\n','http://www.gexf.net/1.1draft/gexf.xsd"');
fprintf(fid, '%s\n','version="1.3">');
fprintf(fid, '%s\n', '<graph mode="dynamic" defaultedgetype="undirected">');
fprintf(fid, '%s\n', '<attributes class="node" mode="dynamic">');

fprintf(fid, '%s\n', '<attribute id="V0rep1_freq" title="V0rep1_freq" type="float"/>');
fprintf(fid, '%s\n', '<attribute id="V0rep2_freq" title="V0rep2_freq" type="float"/>');
fprintf(fid, '%s\n', '<attribute id="V0rep3_freq" title="V0rep3_freq" type="float"/>');
fprintf(fid, '%s\n', '<attribute id="V0rep4_freq" title="V0rep4_freq" type="float"/>');
fprintf(fid, '%s\n', '<attribute id="Vcrep1_freq" title="Vcrep1_freq" type="float"/>');
fprintf(fid, '%s\n', '<attribute id="Vcrep2_freq" title="Vcrep2_freq" type="float"/>');
fprintf(fid, '%s\n', '<attribute id="Vcrep3_freq" title="Vcrep3_freq" type="float"/>');
fprintf(fid, '%s\n', '<attribute id="Vcrep4_freq" title="Vcrep4_freq" type="float"/>');

fprintf(fid, '%s\n', '<attribute id="V0gen8_class" title="V0gen8_class" type="float"/>');
fprintf(fid, '%s\n', '<attribute id="Vcgen8_class" title="Vcgen8_class" type="float"/>');

fprintf(fid, '%s\n', '<attribute id="V0gen8_freq" title="V0gen8_freq" type="float"/>');
fprintf(fid, '%s\n', '<attribute id="Vcgen8_freq" title="Vcgen8_freq" type="float"/>');

fprintf(fid, '%s\n', '</attributes>');
fprintf(fid, '%s\n', '<nodes>');

if ninterp > 0
    numgens = size(frequency_interp,2);
else
    numgens = 9;
end

%write the nodes
for i = 1:numseqs
    
    %only include nodes that are in the dominant network
    if nodes_dom_mask(i)
        
        %e.g., <node id="1" label="LFFGFKIKNVKIIY" start="0" end="9" >'
        fprintf(fid, '<node id="%d" label="%s" start="0" end="%d" >\n' , i, seqlist{i}, numgens);
        
        fprintf(fid, '<attvalues>\n');
        
        %V0rep1 frequency and present/absent over all generations 
        for j = 1:numgens %loop over generations
            fprintf(fid, '<attvalue for="V0rep1_freq" value="%f" start="%d" end="%d"/>\n', frequency_interp(i, j, 1), j-1, j); %numseqs x numgens x numreps
        end
        %V0rep2 frequency and present/absent over all generations
        for j = 1:numgens %loop over generations
            fprintf(fid, '<attvalue for="V0rep2_freq" value="%f" start="%d" end="%d"/>\n', frequency_interp(i, j, 2), j-1, j);
        end
        %V0rep3 frequency and present/absent over all generations
        for j = 1:numgens %loop over generations
            fprintf(fid, '<attvalue for="V0rep3_freq" value="%f" start="%d" end="%d"/>\n', frequency_interp(i, j, 3), j-1, j);
        end
        %V0rep4 frequency and present/absent over all generations
        for j = 1:numgens %loop over generations
            fprintf(fid, '<attvalue for="V0rep4_freq" value="%f" start="%d" end="%d"/>\n', frequency_interp(i, j, 4), j-1, j);
        end
        %Vcrep1 frequency and present/absent over all generations
        for j = 1:numgens %loop over generations
            fprintf(fid, '<attvalue for="Vcrep1_freq" value="%f" start="%d" end="%d"/>\n', frequency_interp(i, j, 5), j-1, j);
        end
        %Vcrep2 frequency and present/absent over all generations
        for j = 1:numgens %loop over generations
            fprintf(fid, '<attvalue for="Vcrep2_freq" value="%f" start="%d" end="%d"/>\n', frequency_interp(i, j, 6), j-1, j);
        end
        %Vcrep3 frequency and present/absent over all generations
        for j = 1:numgens %loop over generations
            fprintf(fid, '<attvalue for="Vcrep3_freq" value="%f" start="%d" end="%d"/>\n', frequency_interp(i, j, 7), j-1, j);
        end
        %Vcrep4 frequency and present/absent over all generations
        for j = 1:numgens %loop over generations
            fprintf(fid, '<attvalue for="Vcrep4_freq" value="%f" start="%d" end="%d"/>\n', frequency_interp(i, j, 8), j-1, j);
        end
        
        %determine which sequences are only observed in a single V0
        %replicate in the final generation
        if all(frequency(i, 9, 1:4) == 0) %numseqs x numgens x numreps
           
            fprintf(fid, '<attvalue for="V0gen8_class" value="%d" start="%d" end="%d"/>\n', 0, 0, numgens);
            fprintf(fid, '<attvalue for="V0gen8_freq" value="%f" start="%d" end="%d"/>\n', 0, 0, numgens);
            
        elseif sum(frequency(i, 9, 1:4) > threshold) == 1
            
            temp = frequency(i, 9, 1:4);
            index = find(temp > threshold);
            fprintf(fid, '<attvalue for="V0gen8_class" value="%d" start="%d" end="%d"/>\n', index, 0, numgens);
            fprintf(fid, '<attvalue for="V0gen8_freq" value="%f" start="%d" end="%d"/>\n', frequency(i, 9, index), 0, numgens);
            
        else
        
            temp = frequency(i, 9, 1:4);
            value = max(temp);
            fprintf(fid, '<attvalue for="V0gen8_class" value="%d" start="%d" end="%d"/>\n', 5, 0, numgens);
            fprintf(fid, '<attvalue for="V0gen8_freq" value="%f" start="%d" end="%d"/>\n', value, 0, numgens);
            
        end
        
        %determine which sequences are only observed in a single Vc
        %replicate in the final generation
        if all(frequency(i, 9, 5:8) == 0)
           
            fprintf(fid, '<attvalue for="Vcgen8_class" value="%d" start="%d" end="%d"/>\n', 0, 0, numgens);
            fprintf(fid, '<attvalue for="Vcgen8_freq" value="%f" start="%d" end="%d"/>\n', 0, 0, numgens);
            
        elseif sum(frequency(i, 9, 5:8) > threshold) == 1
            
            temp = frequency(i, 9, 5:8);
            index = find(temp > threshold);
            fprintf(fid, '<attvalue for="Vcgen8_class" value="%d" start="%d" end="%d"/>\n', index, 0, numgens);
            fprintf(fid, '<attvalue for="Vcgen8_freq" value="%f" start="%d" end="%d"/>\n', frequency(i, 9, index+4), 0, numgens);
            
        else
        
            temp = frequency(i, 9, 5:8);
            value = max(temp);
            fprintf(fid, '<attvalue for="Vcgen8_class" value="%d" start="%d" end="%d"/>\n', 5, 0, numgens);
            fprintf(fid, '<attvalue for="Vcgen8_freq" value="%f" start="%d" end="%d"/>\n', value, 0, numgens);
            
        end
        
        fprintf(fid, '%s\n', '</attvalues>');
        fprintf(fid, '%s\n', '</node>');

    end
end
fprintf(fid, '%s\n', '</nodes>');
    
%write the edges
fprintf(fid, '%s\n', '<edges>');
for i = 1:size(edges_dom,1)
    fprintf(fid, '<edge source="%d" target="%d" />\n', nodes_dom(edges_dom(i, 1)), nodes_dom(edges_dom(i, 2)));
end

fprintf(fid, '%s\n', '</edges>');
fprintf(fid, '%s\n', '</graph>');
fprintf(fid, '%s\n', '</gexf>');

fclose(fid);