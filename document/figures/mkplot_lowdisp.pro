!P.font = 0
!P.MULTI=[0,1,2]
!P.Charsize=1.7
setpl,25,20,test,0,0,/encap
;set_plot,'ps'
device,/color
device,file='magdt.ps'

zstart=0.
zend=7.

loadct,11

tab = rdtab('magdt_lowdisp.data')
medtime=median(tab(2,*))
ii = where(tab(2,*) GE medtime)
jj = where(tab(2,*) LT medtime)
print,mean(tab(1,ii)),sigma(tab(1,ii))
print,mean(tab(1,jj)),sigma(tab(1,jj))
print
print,'2005',median(tab(1,0:7))
print,'2006',median(tab(1,8:31))
print,'2007',median(tab(1,32:46))
print,'2008',median(tab(1,47:*))
print
print,'2005',median(tab(2,0:7))
print,'2006',median(tab(2,8:31))
print,'2007',median(tab(2,32:46))
print,'2008',median(tab(2,47:*))
print
print,'2005',mean(tab(2,0:7))
print,'2006',mean(tab(2,8:31))
print,'2007',mean(tab(2,32:46))
print,'2008',mean(tab(2,47:*))
plot_oi,tab(2,*),tab(3,*),xra=[0.05,100],psym=3,yrange=[25,13],/ysty,/xsty,$
	xtitle='Hours after the burst',ytitle='[R,i,z]-mag',/nodata,$
	position=[0.10,0.10,0.98,0.80]
t = findgen(20)/19.*2.*!Pi
x = cos(t)
y = sin(t)
usersym,x,y,/fill

for n=0,n_elements(tab(0,*))-1 do begin
	z = tab(1,n)
	colnum = 124+(z-zstart)/(zend-zstart)*124
	oplot,[tab(2,n)],[tab(3,n)],psym=8,symsize=2,color=colnum
	;xyouts,tab(2,n)*0.90,tab(3,n)*1.01,strmid(tab(1,n),3,7),charsize=0.9,$
		color=colnum
endfor

;oplot,[0.01,1000.],[20.,20.],linestyle=2
;Overplot powerlaw
t = findgen(100000)/1000.
f = t^(-1.)
mag = -2.5*alog10(f)
ii = where(t EQ 1.)
mags = mag-mag(ii(0))+18.5
;oplot,t,mags,linestyle=2

;Plot GRBs with redshift limits
oplot,[2.15],[21.3],psym=8,symsize=2 ; GRB070129
oplot,[43.0],[22.9],psym=8,symsize=2 ; GRB060708
oplot,[9.5],[22.9],psym=8,symsize=2 ; GRB06080
;oplot,[21.1],[18.8],psym=8,symsize=2
oplot,[10.9],[23.0],psym=8,symsize=2; GRB080523
oplot,[5.7],[20.7],psym=8,symsize=2; GRB050801
;oplot,[13.0],[25.1],psym=8,symsize=2
;arrow,13.0,25.1,13.0,25.8,/data
;oplot,[18.8],[24.0],psym=8,symsize=2
;arrow,18.8,24.0,18.8,24.7,/data

;Overplot limits
;GRB061004
;oplot,[0.542*24.],[25.1],psym=8,symsize=1
;arrow,0.542*24.,25.1,0.542*24.,25.6,/data
;;GRB070129
;oplot,[2.15],[21.3],psym=8,symsize=1
;xyouts,2.15,21.0,'z<3'


;Plot redshift distribution
plot,[zstart,zend],[0,2],/nodata,xtitle='Redshift',$
        position=[0.10,0.88,0.98,0.98],ytickname=[' ',' ',' ',' ',' '],$
        /xst,charsize=1.3,xtickname=[' ','1','2','3','4','5','6',' ']

for nz=0,49 do begin
      z = zstart+(zend-zstart)*float(nz)/float(50.)
      colnum = 124+(z-zstart)/(zend-zstart)*124
      polyfill,[z,z,z+(zend-zstart)/50.,z+(zend-zstart)/50.],[0.02,1.99,1.99,0.02],color=colnum
endfor

col = getcolor(/load)
;limits from spectra
oplot,[2.15],[25.3],psym=8,symsize=2,color=colnum-15

devcl

end
