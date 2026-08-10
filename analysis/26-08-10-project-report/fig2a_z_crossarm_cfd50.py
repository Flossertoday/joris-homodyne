"""Cross-arm coherence check: do the two HBT arms see the same z(t) wiggle?

z_arm1 = (B1-B2)/(B1+B2), z_arm2 = (B3-B4)/(B3+B4). A real Josephson
oscillation is common to both arms; independent photon noise is not.
"""

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


def rebin(x, y, s, k):
    n = len(y) // k * k
    yr = y[:n].reshape(-1, k).sum(1)
    sr = np.sqrt((s[:n].reshape(-1, k) ** 2).sum(1))
    return x[:n].reshape(-1, k).mean(1), yr / k, sr / k


fig, axes = plt.subplots(3, 2, figsize=(13, 11), constrained_layout=True)
for col, k in enumerate([2, 4]):
    b1, s1 = band(0)
    b2, s2 = band(1)
    b3, s3 = band(2)
    b4, s4 = band(3)
    m = (t > 0) & (t < 840)
    tt, B1, S1 = rebin(t[m], b1[m], s1[m], k)
    _, B2, S2 = rebin(t[m], b2[m], s2[m], k)
    _, B3, S3 = rebin(t[m], b3[m], s3[m], k)
    _, B4, S4 = rebin(t[m], b4[m], s4[m], k)

    z1 = (B1 - B2) / (B1 + B2)
    z2 = (B3 - B4) / (B3 + B4)
    sz1 = 2 * np.sqrt((B2 * S1) ** 2 + (B1 * S2) ** 2) / (B1 + B2) ** 2
    sz2 = 2 * np.sqrt((B4 * S3) ** 2 + (B3 * S4) ** 2) / (B3 + B4) ** 2

    # light detrend: remove best-fit quadratic over the window
    def detrend(y):
        c = np.polyfit(tt, y, 2, w=1.0 / np.maximum(sz1, 1e-9))
        return y - np.polyval(c, tt)

    d1, d2 = detrend(z1), detrend(z2)
    # weighted mean z and its noise floor
    w = 1 / sz1**2 + 1 / sz2**2
    zm = (d1 / sz1**2 + d2 / sz2**2) / w
    szm = np.sqrt(1 / w)

    ax = axes[0, col]
    ax.errorbar(tt, d1, yerr=sz1, fmt=".", ms=3, lw=0.5, alpha=0.6, label="arm1 z")
    ax.errorbar(tt, d2, yerr=sz2, fmt=".", ms=3, lw=0.5, alpha=0.6, label="arm2 z")
    ax.set_title(f"z per arm, rebin {k} ({k*0.984:.1f} ps/bin)")
    ax.legend(fontsize=8)

    # cross-arm Pearson on detrended (unweighted) + null from noise sims
    r = np.corrcoef(d1, d2)[0, 1]
    rng = np.random.default_rng(0)
    null = []
    for _ in range(200):
        n1 = rng.normal(0, sz1)
        n2 = rng.normal(0, sz2)
        null.append(np.corrcoef(n1, n2)[0, 1])
    null = np.array(null)
    p = (np.abs(null) >= abs(r)).mean()
    axes[1, col].plot(tt, zm, ".", ms=3, color="magenta")
    axes[1, col].fill_between(tt, -szm, szm, alpha=0.3, color="gray", label="±1 SEM")
    axes[1, col].set_title(f"combined z (detrended); cross-arm r={r:.3f}, null p={p:.3f}")
    axes[1, col].legend(fontsize=8)

    # FFT of combined z vs noise floor
    freq = np.fft.rfftfreq(len(zm), d=k * 0.984)
    spec = np.abs(np.fft.rfft(zm * np.hanning(len(zm))))
    noise_spec = np.zeros_like(spec)
    for _ in range(50):
        noise_spec += np.abs(np.fft.rfft(rng.normal(0, szm) * np.hanning(len(zm))))
    noise_spec /= 50
    axf = axes[2, col]
    per = 1 / freq[1:]
    axf.semilogy(per, spec[1:], ".", ms=3, label="|FFT(z)|")
    axf.semilogy(per, noise_spec[1:], "-", lw=1, label="noise floor")
    axf.set_xlim(2, 200)
    axf.set_xlabel("period [ps]")
    axf.legend(fontsize=8)
    # peak above floor
    ratio = spec[1:] / np.maximum(noise_spec[1:], 1e-12)
    i = int(np.argmax(ratio[(per > 3) & (per < 150)][3:-3]))
    print(f"rebin {k}: cross-arm r={r:.3f} (null p={p:.3f}); "
          f"FFT peak/noise max ratio {ratio.max():.2f}")

fig.savefig(OUT, dpi=150)
print("wrote", OUT)
