
.r particles3d

phase_angle=atan(vy0,vx0)
freq=(phase_angle-shift(phase_angle,[0,-1]))/(tvp(1)-tvp(0))
psize=size(phase_angle)

;Correct for the 2pi shifts.
nt=psize[2]
for ip=0, psize[1] -1 do $
  for it=0,nt-2 do $
  if (freq(ip,it)/omega lt 0.0 ) then phase_angle(ip,it+1:nt-1) += 2*!PI

omega=qd0*bz/md0
theoretical_phase_angle=omega*tvp

shift=phase_angle
for i =0, (size(phase_angle))[1] -1 do shift[i,*] += omega*tvp

plot,tvp,freq(0,*),yrange=omega+[-1,1]*0.1
