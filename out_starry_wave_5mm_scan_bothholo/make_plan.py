# Current both-holo stage file: yellow preset 14 at z=0, black preset 1 at z=+0.5 mm,
# power 40, speed 0.03 mm/s, pause 0. Sky is not written.
# Writes writing_plan_v2_starry_wave_5mm_scan_bothholo_zblack0p5.csv.
# writing_plan_v2_starry_wave_5mm_scan_bothholo.csv is the earlier black focus at z=+1.5 mm; leave it.
import os
from pathlib import Path
os.chdir(Path(__file__).resolve().parent)
import numpy as np, csv
D=np.load('site_map_1000x1000.npy'); C=np.load('nearfield_class_1000x1000.npy')  # C: 0 sky,1 black,2 yellow
P=0.005; V=0.03; PW=40; SRC='gs_starry_wave_5mm_bothholo'
cols="schema_version,operation,group_id,segment_index,laser_state,x_mm,y_mm,z_mm,x2_mm,y2_mm,z2_mm,speed_mm_s,power,dwell_s,pause_s,carbide_preset,source_recipe".split(',')
r=lambda v: round(v,6)
def runs(line):
    d=np.diff(np.r_[0,line.astype(int),0]); return list(zip(np.where(d==1)[0],np.where(d==-1)[0]))
rows=[];g=0;stats={}
for name,mask,z,pre in [('yellow',(D==1)&(C==2),0.0,14),('black',(D==1)&(C==1),0.5,1)]:
    n_on=0;L=[]
    for i in range(mask.shape[0]):
        rr=runs(mask[i])
        if not rr: continue
        g+=1;y=r(i*P);s=0;cur=0
        for a,b in rr:
            if a>cur: s+=1;rows.append([2,'path',g,s,'off',r(cur*P),y,z,r(a*P),y,z,V,PW,'',0,pre,SRC])
            s+=1;rows.append([2,'path',g,s,'on',r(a*P),y,z,r(b*P),y,z,V,PW,'',0,pre,SRC]);L.append(b-a);cur=b
    L=np.array(L);stats[name]=dict(groups=None,on_segments=len(L),on_len_mm=L.sum()*P,min_um=L.min()*5,median_um=np.median(L)*5,max_um=L.max()*5)
with open('writing_plan_v2_starry_wave_5mm_scan_bothholo_zblack0p5.csv','w',newline='') as f:
    w=csv.writer(f,lineterminator='\n');w.writerow(cols);w.writerows(rows)
A=np.array([[x[5],x[6],x[7],x[8]] for x in rows],float)
t=sum(abs(x[8]-x[5]) for x in rows)/V
# self-check
pre=[x[15] for x in rows]; sw=sum(pre[i]!=pre[i-1] for i in range(1,len(pre)))
ok=dict(schema_v2=all(x[0]==2 for x in rows),path=all(x[1]=='path' for x in rows),
 continuity=all(rows[i][5]==rows[i-1][8] for i in range(1,len(rows)) if rows[i][2]==rows[i-1][2]),
 group_ids_monotonic=all(rows[i][2]>=rows[i-1][2] for i in range(1,len(rows))),
 preset_switches=sw, x_bbox=(A[:,0].min(),max(A[:,0].max(),A[:,3].max())),y_bbox=(A[:,1].min(),A[:,1].max()),z_vals=sorted(set(A[:,2])),
 groups=g,segments=len(rows),est_time_h_at_0p03=t/3600)
print(stats);print(ok)

open('plan_check.txt','w').write(str(stats)+'\n'+str(ok)+'\n')
