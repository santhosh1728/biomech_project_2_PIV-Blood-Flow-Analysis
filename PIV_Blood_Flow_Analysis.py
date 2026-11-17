
# PIV analysis for two grayscale BMP images (PIV_1.bmp then PIV_2.bmp)
# Requires: numpy, pandas, scipy, opencv-python, matplotlib, pillow
# Install: pip install numpy pandas scipy opencv-python matplotlib pillow
# Run in the folder containing PIV_1.bmp and PIV_2.bmp.

import os
try:
    import cv2
except ImportError:
    cv2 = None
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
from scipy.signal import fftconvolve
from scipy.ndimage import median_filter, gaussian_filter, distance_transform_edt

DT = 0.300
PIXEL_SIZE = 1.4e-6
RHO_BLOOD = 1060.0
MU_BLOOD = 0.0035
WINDOW = 64
OVERLAP = 32
Y_START = 160  # excludes the bright upper wall/boundary

def preprocess(im):
    if cv2 is not None:
        bg = cv2.GaussianBlur(im, (0, 0), sigmaX=12)
    else:
        bg = gaussian_filter(im, sigma=12)
    hp = im - bg
    lo, hi = np.percentile(hp, [1, 99.8])
    return np.clip((hp-lo)/(hi-lo+1e-12), 0, 1).astype(np.float32)

def subpixel(c, py, px):
    result = []
    for axis, idx in [(0, py), (1, px)]:
        p = [py, px]; pm = p.copy(); pp = p.copy()
        if idx <= 0 or idx >= c.shape[axis]-1:
            result.append(0.0); continue
        pm[axis] -= 1; pp[axis] += 1
        cm, c0, cp = float(c[tuple(pm)]), float(c[tuple(p)]), float(c[tuple(pp)])
        if min(cm,c0,cp) <= 0:
            delta = 0.0
        else:
            den = np.log(cm)-2*np.log(c0)+np.log(cp)
            delta = 0.5*(np.log(cm)-np.log(cp))/den if abs(den)>1e-12 else 0.0
        result.append(float(np.clip(delta,-1,1)))
    return result

def piv(a,b):
    step = WINDOW-OVERLAP
    ys = np.arange(Y_START, a.shape[0]-WINDOW+1, step)
    xs = np.arange(0, a.shape[1]-WINDOW+1, step)
    X,Y = np.meshgrid(xs+WINDOW/2, ys+WINDOW/2)
    du = np.full(X.shape,np.nan); dv = np.full(X.shape,np.nan)
    snr = np.full(X.shape,np.nan)
    hann = np.outer(np.hanning(WINDOW),np.hanning(WINDOW)).astype(np.float32)
    center = WINDOW-1; lim = WINDOW//2
    for j,y in enumerate(ys):
        for i,x in enumerate(xs):
            A=(a[y:y+WINDOW,x:x+WINDOW]-a[y:y+WINDOW,x:x+WINDOW].mean())*hann
            B=(b[y:y+WINDOW,x:x+WINDOW]-b[y:y+WINDOW,x:x+WINDOW].mean())*hann
            if A.std()<1e-6 or B.std()<1e-6: continue
            C=fftconvolve(B,A[::-1,::-1],mode="full")
            C[:center-lim,:]=0; C[center+lim+1:,:]=0
            C[:,:center-lim]=0; C[:,center+lim+1:]=0
            py,px=np.unravel_index(np.argmax(C),C.shape)
            dy_sub,dx_sub=subpixel(C,py,px)
            C2=C.copy()
            C2[max(0,py-2):py+3,max(0,px-2):px+3]=0
            du[j,i]=(px-center)+dx_sub
            dv[j,i]=(py-center)+dy_sub
            snr[j,i]=C[py,px]/(C2.max()+1e-12)
    return X,Y,du,dv,snr

def fill_invalid(arr, mask):
    out=arr.copy()
    out[~mask]=np.nan
    inds=distance_transform_edt(~mask,return_distances=False,return_indices=True)
    out[~mask]=out[tuple(ind[~mask] for ind in inds)]
    return out

file1 = "PIV_1.bmp" if os.path.exists("PIV_1.bmp") else "image1.bmp"
file2 = "PIV_2.bmp" if os.path.exists("PIV_2.bmp") else "Image2.bmp"
im1 = np.asarray(Image.open(file1).convert("L"), dtype=np.float32)
im2 = np.asarray(Image.open(file2).convert("L"), dtype=np.float32)
if im1.shape != im2.shape: raise ValueError("Image dimensions differ.")
a,b=preprocess(im1),preprocess(im2)
H,W=im1.shape
Xpx,Ypx,du,dv,snr=piv(a,b)

valid=np.isfinite(du)&np.isfinite(dv)&(snr>=1.3)
mu_u=median_filter(np.nan_to_num(du,nan=0),size=3)
mu_v=median_filter(np.nan_to_num(dv,nan=0),size=3)
valid &= (np.abs(du-mu_u)<5)&(np.abs(dv-mu_v)<5)
du=gaussian_filter(fill_invalid(du,valid),sigma=0.8)
dv=gaussian_filter(fill_invalid(dv,valid),sigma=0.8)

dx=(WINDOW-OVERLAP)*PIXEL_SIZE
dy=dx
X=Xpx*PIXEL_SIZE
Y=(H-Ypx)*PIXEL_SIZE
U=du*PIXEL_SIZE/DT
V=-dv*PIXEL_SIZE/DT
Wcomp=np.zeros_like(U)
speed=np.sqrt(U**2+V**2+Wcomp**2)

dUdx=np.gradient(U,dx,axis=1)
dUdy=-np.gradient(U,dy,axis=0)
dVdx=np.gradient(V,dx,axis=1)
dVdy=-np.gradient(V,dy,axis=0)

Exx=dUdx; Eyy=dVdy; Exy=0.5*(dUdy+dVdx)
Ezz=np.zeros_like(U)
strain_mag=np.sqrt(Exx**2+Eyy**2+2*Exy**2)

Txx=2*MU_BLOOD*Exx
Tyy=2*MU_BLOOD*Eyy
Txy=2*MU_BLOOD*Exy
Tzz=np.zeros_like(U)
stress_mag=np.sqrt(Txx**2+Tyy**2+2*Txy**2)

out="PIV_Blood_Flow_Output"
os.makedirs(out,exist_ok=True)
data=pd.DataFrame({
 "x_mm":(X*1e3).ravel(),"y_mm":(Y*1e3).ravel(),
 "SNR":snr.ravel(),"valid_vector":valid.ravel(),
 "U_m_s":U.ravel(),"V_m_s":V.ravel(),"W_m_s_assumed":Wcomp.ravel(),
 "speed_m_s":speed.ravel(),"Exx_s-1":Exx.ravel(),"Eyy_s-1":Eyy.ravel(),
 "Exy_s-1":Exy.ravel(),"Ezz_s-1_assumed":Ezz.ravel(),
 "strain_rate_magnitude_s-1":strain_mag.ravel(),
 "tau_xx_Pa":Txx.ravel(),"tau_yy_Pa":Tyy.ravel(),"tau_xy_Pa":Txy.ravel(),
 "tau_zz_Pa_assumed":Tzz.ravel(),"viscous_stress_magnitude_Pa":stress_mag.ravel()
})
data.to_csv(os.path.join(out,"PIV_Blood_Flow_Results.csv"),index=False)

def plot(data,title,label,name,symmetric=False):
    fig,ax=plt.subplots(figsize=(8,6))
    if symmetric:
        lim=max(np.nanpercentile(np.abs(data),98),1e-12)
        im=ax.pcolormesh(X*1e3,Y*1e3,data,shading="auto",cmap="coolwarm",vmin=-lim,vmax=lim)
    else:
        im=ax.pcolormesh(X*1e3,Y*1e3,data,shading="auto",cmap="turbo")
    fig.colorbar(im,ax=ax,label=label)
    ax.set(xlabel="x (mm)",ylabel="y (mm, positive upward)",title=title)
    ax.set_aspect("equal",adjustable="box"); fig.tight_layout()
    fig.savefig(os.path.join(out,name),dpi=220); plt.close(fig)

plot(U,"Horizontal velocity U","U (m/s)","U_velocity.png",True)
plot(V,"Vertical velocity V","V (m/s)","V_velocity.png",True)
plot(Wcomp,"Out-of-plane velocity W (assumed zero)","W (m/s)","W_assumed_zero.png")
plot(speed,"Velocity magnitude","|V| (m/s)","velocity_magnitude.png")
plot(Exx,"Normal strain rate Exx","Exx (s-1)","Exx.png",True)
plot(Eyy,"Normal strain rate Eyy","Eyy (s-1)","Eyy.png",True)
plot(Exy,"Shear strain rate Exy","Exy (s-1)","Exy.png",True)
plot(strain_mag,"Strain-rate magnitude","s-1","strain_magnitude.png")
plot(Txx,"Viscous stress tau_xx","Pa","tau_xx.png",True)
plot(Tyy,"Viscous stress tau_yy","Pa","tau_yy.png",True)
plot(Txy,"Viscous shear stress tau_xy","Pa","tau_xy.png",True)
plot(stress_mag,"Viscous stress magnitude","Pa","stress_magnitude.png")

fig,ax=plt.subplots(figsize=(9,6))
q=ax.quiver(X*1e3,Y*1e3,U,V,speed,cmap="turbo")
fig.colorbar(q,ax=ax,label="Speed (m/s)")
ax.set(xlabel="x (mm)",ylabel="y (mm)",title="PIV velocity vectors")
ax.set_aspect("equal",adjustable="box"); fig.tight_layout()
fig.savefig(os.path.join(out,"velocity_vectors.png"),dpi=220); plt.close(fig)

print("Valid vectors:",valid.sum(),"/",valid.size)
print("Median speed (mm/s):",np.median(speed)*1e3)
print("Maximum speed (mm/s):",np.max(speed)*1e3)
print("Saved results in:",out)
