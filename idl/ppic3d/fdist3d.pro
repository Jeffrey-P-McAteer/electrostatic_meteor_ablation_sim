; Read in the complete distribution function
@ppic3d.i
;@eppic.i
;if (n_elements(fd0) le 1 and file_exist('domain000/fdist0.bin')) then begin
if file_exist('domain000/fdist0.bin') then begin
    if n_elements(fny0) eq 0 then fny0=1
    if n_elements(fnvy0) eq 0 then fnvy0=1
    if n_elements(fnz0) eq 0 then fnz0=1
    if n_elements(fnvz0) eq 0 then fnvz0=1
    fd0size=long(fnx0)*fnvx0*long(fny0)*fnvy0*long(fnz0)*fnvz0



    fd0=readarray('domain000/fdist0.bin',fd0size, $
	/binary,skip=iskip)
    fd0size=size(fd0)
    if fd0size(0) gt 1 then nfd0=fd0size(2) else nfd0=1
    fd0=reform(fd0,fnvz0,fnvy0,fnvx0,fnz0,fny0,fnx0,nfd0)
endif

if file_exist('domain000/fdist1.bin') then begin
    if n_elements(fny1) eq 0 then fny1=1
    if n_elements(fnvy1) eq 0 then fnvy1=1
    if n_elements(fnz1) eq 0 then fnz1=1
    if n_elements(fnvz1) eq 0 then fnvz1=1
    fd1size=long(fnx1)*fnvx1*long(fny1)*fnvy1*long(fnz1)*fnvz1



    fd1=readarray('domain000/fdist1.bin',fd1size, $
	/binary,skip=iskip)
    fd1size=size(fd1)
    if fd1size(0) gt 1 then nfd1=fd1size(2) else nfd1=1
    fd1=reform(fd1,fnvz1,fnvy1,fnvx1,fnz1,fny1,fnx1,nfd1)
endif


if file_exist('domain000/fdist2.bin') then begin
    if n_elements(fny2) eq 0 then fny2=1
    if n_elements(fnvy2) eq 0 then fnvy2=1
    if n_elements(fnz2) eq 0 then fnz2=1
    if n_elements(fnvz2) eq 0 then fnvz2=1
    fd2size=long(fnx2)*fnvx2*long(fny2)*fnvy2*long(fnz2)*fnvz2



    fd2=readarray('domain000/fdist2.bin',fd2size, $
	/binary,skip=iskip)
    fd2size=size(fd2)
    if fd2size(0) gt 1 then nfd2=fd2size(2) else nfd2=1
    fd2=reform(fd2,fnvz2,fnvy2,fnvx2,fnz2,fny2,fnx2,nfd2)
endif


if file_exist('domain000/fdist3.bin') then begin
    if n_elements(fny3) eq 0 then fny3=1
    if n_elements(fnvy3) eq 0 then fnvy3=1
    if n_elements(fnz3) eq 0 then fnz3=1
    if n_elements(fnvz3) eq 0 then fnvz3=1
    fd3size=long(fnx3)*fnvx3*long(fny3)*fnvy3*long(fnz3)*fnvz3



    fd3=readarray('domain000/fdist3.bin',fd3size, $
	/binary,skip=iskip)
    fd3size=size(fd3)
    if fd3size(0) gt 1 then nfd3=fd3size(2) else nfd3=1
    fd3=reform(fd3,fnvz3,fnvy3,fnvx3,fnz3,fny3,fnx3,nfd3)
endif

; Let's print all the images at the first time step side by side
;s=size(fd0)
;;!p.multi=[0,s(3),s(4)]
;for iy=0,s(3)-1 do begin
;    for ix=0,s(4)-1 do begin
;        fn=transpose(reform(fd0(*,*,iy,ix,0)))
;        image_plot,fn,/zlog
;    endfor
;endfor
;;!p.multi=[0,1,1]

end
