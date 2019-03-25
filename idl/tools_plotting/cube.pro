pro cube, x,y,z

plots,[0,0],[0,y],[0,0],/t3d ; plot y axis

plots,[x,x],[0,y],[0,0],/t3d ; plot y axis
plots,[0,0],[0,y],[z,z],/t3d ; plot y axis
plots,[x,x],[0,y],[z,z],/t3d ; plot y axis

;plot z parallels
plots,[0,0],[0,0],[0,z],/t3d ; 

plots,[x,x],[y,y],[0,z],/t3d ; plot y axis
plots,[x,x],[0,0],[0,z],/t3d ; plot y axis
plots,[0,0],[y,y],[0,z],/t3d ; plot y axis

; plot  x parallels
plots,[0,x],[0,0],[0,0],/t3d ; plot x axis
plots,[0,x],[y,y],[0,0],/t3d ; plot x axis
plots,[0,x],[y,y],[z,z],/t3d ; plot x axis
plots,[0,x],[0,0],[z,z],/t3d ; 

return
end

pro coordgrid,dire,a,b,x,y,z

; dire =1 for x plane
; dire =2 for y plane
; dire =3 for z plane
; a,b how many lines in each direction
; x,y,z size of cube



case dire of 

1: begin 
for i=0,a-1 do begin
    pa=i*y/(a-1)
    plots,[0.5*x,0.5*x],[pa,pa],[0,z],/t3d
endfor

for j=0,b-1 do begin
    pb=j*z/(b-1)
    plots,[0.5*x,0.5*x],[0,y],[pb,pb],/t3d
endfor

end

2: begin 
for i=0,a-1 do begin
    pa=i*x/(a-1)
    plots,[pa,pa],[0.5*y,0.5*y],[0,z],/t3d
endfor

for j=0,b-1 do begin
    pb=j*z/(b-1)
    plots,[0,x],[0.5*y,0.5*y],[pb,pb],/t3d
endfor

end


3: begin 
for i=0,a-1 do begin
    pa=i*x/(a-1)
    plots,[pa,pa],[0,y],[0.5*z,0.5*z],/t3d
endfor

for j=0,b-1 do begin
    pb=j*y/(b-1)
    plots,[0,x],[pb,pb],[0.5*z,0.5*z],/t3d
endfor

end

endcase

end

    
