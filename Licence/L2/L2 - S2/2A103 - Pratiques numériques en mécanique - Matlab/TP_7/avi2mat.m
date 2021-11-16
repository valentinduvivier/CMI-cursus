function avi2mat(filename)

vidObj = VideoReader([filename,'.avi']);

k = 1;

while hasFrame(vidObj)
    Mbuffer = readFrame(vidObj);
    M(:,:,k) = Mbuffer(:,:,1); % on ne récupere qu'une "couleur" R
                               % les autres G et B sont identiques.
    k = k+1;
end

save([filename,'.mat'],'M')

end
