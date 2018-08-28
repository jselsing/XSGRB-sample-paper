!P.font = 0
!P.MULTI=[0,1,2]
!P.Charsize=1.7
setpl,25,20,test,0,0,/encap
;set_plot,'ps'
device,/color
device,file='magdt_hd.ps'

zstart=0.
zend=8.

loadct,11

tab = rdtab('magdt_highdisp.data')
medtime=median(tab(1,*))
ii = where(tab(1,*) GE medtime)
jj = where(tab(1,*) LT medtime)
print,mean(tab(3,ii)),sigma(tab(1,ii))
print,mean(tab(3,jj)),sigma(tab(1,jj))
print
za = tab(3,*)

plot_oi,tab(1,*),tab(2,*),xra=[0.05,110],psym=3,yrange=[26,16],/ysty,/xsty,$
	xtitle='Hours after the burst',ytitle='[R,i,z]-mag',/nodata,$
	position=[0.10,0.10,0.98,0.80]
t = findgen(20)/19.*2.*!Pi
x = cos(t)
y = sin(t)
usersym,x,y,/fill

for n=0,n_elements(tab(0,*))-1 do begin
	z = tab(3,n)
	colnum = 124+(z-zstart)/(zend-zstart)*124
	oplot,[tab(1,n)],[tab(2,n)],psym=8,symsize=2,color=colnum
endfor

;Now the bursts for which the redshift was determined from the host lines
tab = rdtab('magdt_highdisp_dark.data')
zd = tab(3,*)

t = findgen(5)/4.*2.*!Pi-!Pi/2.
x = cos(t)
y = sin(t)
usersym,x,y,/fill
ii = where(tab(4,*) EQ 1.)
for n=0,n_elements(tab(0,ii))-1 do begin
	z = tab(3,ii(n))
	colnum = 124+(z-zstart)/(zend-zstart)*124
	oplot,[tab(1,ii(n))],[tab(2,ii(n))],psym=8,symsize=3,color=colnum
endfor

t = findgen(4)/3.*2.*!Pi-!Pi/2.
x = cos(t)
y = sin(t)
usersym,x,y,/fill
ii = where(tab(4,*) EQ 0.)
for n=0,n_elements(tab(0,ii))-1 do begin
	z = tab(3,ii(n))
	colnum = 124+(z-zstart)/(zend-zstart)*124
	oplot,[tab(1,ii(n))],[tab(2,ii(n))],psym=8,symsize=3,color=colnum
endfor

;Plot redshift distribution
plot,[zstart,zend],[0,2],/nodata,xtitle='Redshift',$
        position=[0.10,0.88,0.98,0.98],ytickname=[' ',' ',' ',' ',' '],$
        /xst,charsize=1.3;,xtickname=[' ','1','2','3','4','5','6',' ']

for nz=0,49 do begin
      z = zstart+(zend-zstart)*float(nz)/float(50.)
      colnum = 124+(z-zstart)/(zend-zstart)*124
      polyfill,[z,z,z+(zend-zstart)/50.,z+(zend-zstart)/50.],[0.02,1.99,1.99,0.02],color=colnum
endfor

devcl

print,'Mean redshift',mean([[za],[zd]]),median([[za],[zd]])

end
