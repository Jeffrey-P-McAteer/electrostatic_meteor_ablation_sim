
function timesteps,filename,sizepertime,nsubdomains,basepath=basepath,double=double

wordsize=4
if(n_elements(double) ne 0) then wordsize=8

if n_elements(basepath) eq 0 then basepath="domain*/"

filename_search_strg = basepath + filename
filename_list=file_search(filename_search_strg)

fileunit=11
openr, fileunit,  filename_list[0]
fsize_bytes=long64((fstat(fileunit)).size)
close,fileunit
mintime = long64(fsize_bytes/wordsize/sizePerTime)
for idomain=1,nsubdomains-1 do begin
    openr, fileunit,  filename_list[idomain]
    fsize_bytes=long64((fstat(fileunit)).size)
    close,fileunit
    thistime = long64(fsize_bytes/wordsize/sizePerTime)
    if mintime gt thistime then mintime = thistime

endfor
return,mintime
end
