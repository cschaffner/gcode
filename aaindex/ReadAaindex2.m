clear all;

% reading in aaindex2
count=0;
fid = fopen('aaindex2');
% initialize array
aaindex2=single(zeros(20,20,67));

while ~feof(fid)
    C = textscan(fid, 'H %s\n D %[^\n]');
    if isempty(C{1}) 
           break ;
    end
    D = textscan(fid,'%c %*[^\n]',1);
    while ~(D{1} == 'M' && strcmp(D{2},'rows'))
        D = textscan(fid, '%c %s %s %s %*[^\n]',1);
    end
    
    
    if strcmp(D{4},'ARNDCQEGHILKMFPSTWYV,')
      count=count+1;
      F = textscan(fid, '%f32',210,'TreatAsEmpty','NA');
      k=0;
      for i=1:20
        for j=1:i
            k=k+1;
            aaindex2(i,j,count)=F{1}(k);
        end
      end
      aaind2_name(count)=C{1};
      aaind2_desc(count)=C{2};

    else 
        D{4}
        disp('shitty format');
        zz=textscan(fid,'%[^////]');
    end
%     elseif strcmp(D{4},'-ARNDCQEGHILKMFPSTWYV,') || strcmp(D{4},'ACDEFGHIKLMNPQRSTVWYJ-,')
%       F = textscan(fid, '%f32',441,'TreatAsEmpty','NA');
%       k=0;
%       for i=1:21
%         for j=1:21
%             k=k+1;
%             aaindex2(i,j,count)=F{1}(k);
%         end
%       end
% %      disp('nono');
%     else 
%         D{4}
%    end

    % assign values
    

    textscan(fid,'////');
end
fclose(fid);

