"""Final-figure parameters: SNR vs rebin; sinusoid upper limit vs period."""

import json

import matplotlib

matplotlib.use("Agg")
import numpy as np

import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--npz', required=True)
parser.add_argument('--out', required=True)
_args = parser.parse_args()
NPZ = _args.npz
OUT = _args.out
from matplotlib import pyplot as plt


d = np.load(NPZ)
t = d["time_ps"]

def band(i, src="photon"):
    return d[f"band_traces_{src}"][i], d[f"band_sem_{src}"][i]

# --- SNR vs rebin for L/R photon traces ---
bs = [band(i) for i in range(4)]
m = (t > 20) & (t < 800)
print("rebin -> median SNR (L, R) over 20-800 ps:")
for k in [2, 4, 8, 16]:
    snrs = []
    for idx in ([0, 2], [1, 3]):
        tr = sum(bs[i][0][m] for i in idx)
        se = np.sqrt(sum(bs[i][1][m] ** 2 for i in idx))
        n = len(tr) // k * k
        tr = tr[:n].reshape(-1, k).mean(1)
        se = se[:n].reshape(-1, k)
        se = np.sqrt((se**2).sum(1) / k**2 * (1 / k)) if False else np.sqrt((se**2).sum(1)) / k
        snrs.append(np.median(tr / se))
    print(f"  {k} ({k*0.984:.1f} ps/bin): L {snrs[0]:.1f}, R {snrs[1]:.1f}")

# --- sinusoid amplitude upper limit from combined z ---
# combine arms: z per arm, weighted average, propagate SEM
def z_arm(i, j):
    tr_i, se_i = bs[i]
    tr_j, se_j = bs[j]
    z = (tr_i - tr_j) / (tr_i + tr_j)
    sz = 2 * np.sqrt((tr_j * se_i) ** 2 + (tr_i * se_j) ** 2) / (tr_i + tr_j) ** 2
    return z, sz

z1, sz1 = z_arm(0, 1)
z2, sz2 = z_arm(2, 3)
w = 1 / sz1[m] ** 2 + 1 / sz2[m] ** 2
z = (z1[m] / sz1[m] ** 2 + z2[m] / sz2[m] ** 2) / w
sz = np.sqrt(1 / w)
tt = t[m]
# remove slow trend
c = np.polyfit(tt, z, 2, w=1 / sz)
zr = z - np.polyval(c, tt)

periods = np.concatenate([np.arange(3, 20, 0.5), np.arange(20, 200, 2.0)])
ul95 = np.zeros_like(periods)
amp = np.zeros_like(periods)
for k, P in enumerate(periods):
    X = np.column_stack([np.cos(2 * np.pi * tt / P), np.sin(2 * np.pi * tt / P)])
    W = 1 / sz**2
    XtW = X.T * W
    cov = np.linalg.inv(XtW @ X)
    beta = cov @ (XtW @ zr)
    amp[k] = np.hypot(*beta)
    ul95[k] = amp[k] + 2 * np.sqrt(max(cov[0, 0], cov[1, 1]))

fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
axes[0].semilogy(periods, amp, ".", ms=3, label="fitted |A|")
axes[0].semilogy(periods, ul95, "-", lw=1.5, label="|A| + 2σ (≈95% UL)")
axes[0].set_xlabel("period [ps]")
axes[0].set_ylabel("z oscillation amplitude")
axes[0].set_title("z(t) sinusoid amplitude vs period (combined arms)")
axes[0].legend(fontsize=8)
axes[1].plot(tt, zr, ".", ms=2, color="magenta")
axes[1].fill_between(tt, -sz, sz, alpha=0.3, color="gray")
axes[1].set_xlabel("t [ps]")
axes[1].set_title("combined detrended z with ±1 SEM")
fig.savefig(OUT, dpi=150)
print(f"z SEM per 1-px bin: median {np.median(sz):.4f}")
print(f"95% UL on z amplitude: median over periods {np.median(ul95):.4f}, worst {ul95.max():.4f}")
print("wrote", OUT)
