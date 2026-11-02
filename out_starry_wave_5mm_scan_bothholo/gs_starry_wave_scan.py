# Python port of gs_starry_night_great_wave_hologram.m + gs_great_wave_scanned_lines.m (case 3)
# 5x5 mm, 1000x1000 sites @5um, 5um dots/lines, carrier 10/3 sites, 405 nm, 200 GS iters, seed 1.
# GS aperture is yellow + black. Pictures load from the repo root; outputs are written in this folder.
import os
from pathlib import Path
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
import numpy as np, csv, json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from PIL import Image
from scipy.ndimage import uniform_filter, map_coordinates
R = str(HERE.parent) + os.sep
S=1000; P=5.0; D=5.0; SPS=20; WL=0.405; N=2*S; CP=10/3; FR=0.38; IT=200
rng=np.random.default_rng(1)
fft2c=lambda a: np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(a)))
def resample(img,side,sq,odd=True):
    sc=max(img.shape)/side; k=max(2*round((sc-1)/2)+1,1)
    if k>1: img=uniform_filter(img,k,mode='constant')
    off=(np.arange(side)-(side-1)/2)*sc
    yq,xq=np.meshgrid((img.shape[0]-1)/2+off,(img.shape[1]-1)/2+off,indexing='ij')
    return map_coordinates(img,[yq,xq],order=1,cval=0)
# near field classes
pal=np.array([[1,55,175],[0,0,0],[253,239,3]])/255
rgb=np.asarray(Image.open(R+'starry_night_v2.png').convert('RGB'))/255
cls=np.argmin(((rgb[...,None,:]-pal)**2).sum(-1),-1)[::-1]
fr=np.stack([resample((cls==c).astype(float),S,None) for c in range(3)],-1)
sc_=np.argmax(fr,-1)  # 0 sky 1 black 2 yellow
idx=slice(N//2-S//2,N//2+S//2)
written=np.zeros((N,N),bool); written[idx,idx]=True
cl=np.full((N,N),-1); cl[idx,idx]=sc_
yellow=(cl==2)|(cl==1); notblack=written  # both yellow+black carry hologram
# site transforms
po=((np.arange(SPS)+.5)/SPS)*P-P/2; px,py=np.meshgrid(po,po)
f=(np.arange(N)-N/2)/(N*P); pp=np.exp(-2j*np.pi*np.outer(f,po))/SPS
T=lambda m: pp@m.astype(float)@pp.T
disk=np.hypot(px,py)<=D/2; band=np.abs(py)<=D/2
tiles=[disk,disk|(band&(px>=0)),band,disk|(band&(px<=0))]
C=T(np.ones((SPS,SPS))); Dt=T(disk); TT=[T(t) for t in tiles]
# target
g=np.asarray(Image.open(R+'great_wave.png').convert('RGB'))/255; g=g.max(2); g=(g-g.min())/(g.max()-g.min())
side=round(FR*N); sq=slice(N//2-side//2,N//2-side//2+side)
tgtI=np.zeros((N,N)); tgtI[sq,sq]=resample(g[::-1],side,None)
placed=np.zeros((N,N),bool); placed[sq,sq]=True
cs=round(N/CP); target=np.roll(tgtI,cs,0); box=np.roll(placed,cs,0)
den=np.roll(np.abs(Dt),-cs,0); gsT=np.where(tgtI>0,np.sqrt(tgtI)/np.maximum(den,1e-12),0.0)
def gs(t,ap,ph0):
    t=np.fft.ifftshift(t); ap=np.fft.ifftshift(ap); fc=t*np.exp(1j*np.fft.ifftshift(ph0))
    for _ in range(IT):
        h=np.fft.ifft2(fc); h=ap*h/np.maximum(abs(h),1e-300); F=np.fft.fft2(h); fc=t*F/np.maximum(abs(F),1e-300)
    return np.fft.fftshift(np.angle(h))
ph=gs(gsT,yellow.astype(float),2*np.pi*rng.random((N,N)))
Y=np.arange(N)[:,None]-N/2
dots=yellow&(np.cos(2*np.pi*cs/N*Y+ph)<=0)
# runs along x -> tile masks + segment list
d=dots; L=np.zeros_like(d); L[:,1:]=d[:,:-1]; Rn=np.zeros_like(d); Rn[:,:-1]=d[:,1:]
m=[d&~L&~Rn, d&~L&Rn, d&L&Rn, d&L&~Rn]
F_dot=C*fft2c(notblack.astype(float))-Dt*fft2c(d.astype(float))
F_line=C*fft2c(notblack.astype(float))-sum(t*fft2c(mm.astype(float)) for t,mm in zip(TT,m))
def metr(F):
    I=abs(F)**2; eff=I[box].sum()/N**2/S**2
    a=I[box]; b=target[box]; c=np.corrcoef(a,b)[0,1]
    Is=uniform_filter(I,5); cs5=np.corrcoef(Is[box],b)[0,1]; return eff,c,cs5,I
res={}
for k,F in [('dots',F_dot),('lines',F_line)]:
    e,c,c5,I=metr(F); res[k]=dict(efficiency=e,corr=c,corr_5x5=c5)
    np.save(f'farfield_{k}_intensity_crop.npy',I[box.any(1)][:,box.any(0)].astype(np.float32))
    res[k]['I']=I
# site map: rows y (row 0 = y min = bottom), cols x
site=dots[idx,idx].astype(np.uint8); cls_s=sc_
np.save('site_map_1000x1000.npy',site); np.savetxt('site_map_1000x1000.csv',site,fmt='%d',delimiter=',')
np.save('nearfield_class_1000x1000.npy',cls_s.astype(np.uint8))
Image.fromarray(((1-site[::-1])*255).astype(np.uint8)).save('site_map_1000x1000_black_is_dot.png')
# scan segments: per row, line from x_start to x_end (dot centres, mm, origin 0,0 at site (0,0))
segs=[]
for r in range(S):
    row=site[r]; dd=np.diff(np.r_[0,row.astype(int),0]); st=np.where(dd==1)[0]; en=np.where(dd==-1)[0]-1
    for a,b in zip(st,en): segs.append((r,a,b,r*P/1000,a*P/1000,b*P/1000,b-a+1))
with open('scan_segments.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['row','col_start','col_end','y_mm','x_start_mm','x_end_mm','n_sites']); w.writerows(segs)
lens=np.array([s[6] for s in segs])
th=np.degrees(np.arcsin(WL*f)); sinc=WL/(CP*P)
# figures
blue=matplotlib.colors.LinearSegmentedColormap.from_list('b',['k','#3060ff','w'])
br=box.any(1); bc=box.any(0)
fig,ax=plt.subplots(1,3,figsize=(18,6.5))
for a,k in zip(ax[:2],['dots','lines']):
    I=res[k]['I']; a.imshow(I[br][:,bc],origin='lower',cmap=blue,vmax=1.3*np.percentile(I[box],99.5),
      extent=[th[bc][0],th[bc][-1],th[br][0],th[br][-1]])
    a.set_title(f"{k}: eff {100*res[k]['efficiency']:.2f}%, corr {res[k]['corr']:.3f} (5x5 {res[k]['corr_5x5']:.3f})"); a.set_xlabel('θx (°)'); a.set_ylabel('θy (°)')
ax[2].imshow(target[br][:,bc],origin='lower',cmap=blue); ax[2].set_title('target')
fig.suptitle('Starry Night 5x5 mm (1000x1000 @5 µm), Great Wave +1 order, 405 nm, scan lines along x'); fig.savefig('farfield_sim_405nm.png',dpi=110); plt.close()
I=res['lines']['I']; plt.figure(figsize=(8,8)); plt.imshow(I,origin='lower',cmap=blue,vmax=1.3*np.percentile(I[box],99.5),extent=[th[0],th[-1],th[0],th[-1]]); plt.title('full far field (lines): +1 picture, zero order, -1 twin'); plt.savefig('farfield_full_lines.png',dpi=100); plt.close()
dc=np.array([[1,55,178],[0,0,0],[255,240,0]])/255
nf=np.where(cls_s[...,None]==1,dc[1],dc[0]); nf[(site==1)&(cls_s==2)]=dc[2]; nf[(site==1)&(cls_s==1)]=np.array([1,0.2,0.2])
fig,ax=plt.subplots(1,2,figsize=(14,7))
ax[0].imshow(nf,origin='lower',extent=[0,5,0,5]); ax[0].set_title('near field (sites): sky none; yellow=preset14 sites, red=preset1 sites (black region)'); ax[0].set_xlabel('x (mm)')
z=slice(480,520); ax[1].imshow(nf[z,z],origin='lower',interpolation='nearest',extent=[2.4,2.6,2.4,2.6]); ax[1].set_title('40x40-site close-up')
fig.savefig('nearfield.png',dpi=110); plt.close()
out=dict(sites=f'{S}x{S}',size_mm=5,pitch_um=P,dot_um=D,line_width_um=D,wavelength_nm=405,carrier_period_sites=CP,carrier_period_um=CP*P,
 carrier_angle_deg=float(np.degrees(np.arcsin(sinc))),picture_halfspan_deg=float(np.degrees(np.arcsin(WL*FR/(2*P)))),
 gs_iters=IT,seed=1,fft_N=N,n_dot_sites=int(site.sum()),n_aperture_sites=int(yellow.sum()),n_lines=len(segs),
 median_line_sites=float(np.median(lens)),single_site_lines=int((lens==1).sum()),max_line_sites=int(lens.max()),
 total_on_length_mm=float(lens.sum()*P/1000),
 dots={k:v for k,v in res['dots'].items() if k!='I'},lines={k:v for k,v in res['lines'].items() if k!='I'})
json.dump(out,open('summary.json','w'),indent=1); print(json.dumps(out,indent=1))
