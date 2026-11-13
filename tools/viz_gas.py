# -*- coding: utf-8 -*-
"""
川谷さんの修論（図 3.3.22）と同じ体裁の「二酸化炭素吸収量分布」の図
=====================================================================
左側 : 飛んでいる液滴を液滴径で色分けした散布図
右側 : 気相の二酸化炭素吸収量 C_total [mol]（各セルで、初期濃度に対する
       変化率 β × 初期濃度 × セルの体積。flow.dat の cvolume）
         円筒（Single・対向）: C_total = ∭ β C_g,ini r dr dθ dz（川谷さんの式 (34)）
         平面（並列・xy）    : C_total = β C_g,ini Δx Δy × 奥行き 1 m
- 円筒の計算は、左に液滴（r を負側に折り返す）、右に C_total（川谷さんと同じ）
- 並列の計算は、ノズル間隔の真ん中で分けて、左半分に液滴、右半分に C_total
  （左右対称なので、それぞれ片方のノズルの噴霧に当たる）
- 電極の中（気相の計算をしていないセル）は灰色
- 色の区切りは川谷さんと同じ比（上端の 1, 2, 4, 6, 8 %、10, 20, 40, 60, 80, 100 %）。
  上端は全コマの最大値から自動（--ctop で指定可）

使い方（計算フォルダで）:
  python3 viz_gas.py                    # 86〜100 mm
  python3 viz_gas.py --zmin 0           # 下まで全部
  python3 viz_gas.py --gif              # GIF も作る（flow.dat のコマ数だけ）
出力: viz_gas/
"""
import os, re, sys, argparse, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

EXPFIX = re.compile(r"(\d\.\d+)([+-]\d{2,3})(?![\d.Ee])")

def parse_block(lines, nc):
    txt = "".join(lines)
    a = np.array(txt.split(), dtype=object)
    try:
        v = np.array(txt.split(), dtype=float)
    except ValueError:
        v = np.array(EXPFIX.sub(r"\1E\2", txt).split(), dtype=float)
    return v.reshape(-1, nc)

def read_flow(path, stride=1):
    """flow.dat の全コマ（x, y, cvolume, cgg）を読む。"""
    frames = []
    with open(path, errors="replace") as fh:
        k = -1
        for line in fh:
            if "zone" not in line.lower():
                continue
            m = re.search(r"n=\s*(\d+)", line)
            mi = re.search(r"i=\s*(\d+)", line)
            mj = re.search(r"j=\s*(\d+)", line)
            ni, nj = int(mi.group(1)), int(mj.group(1))
            block = list(itertools.islice(fh, ni * nj))
            k += 1
            if k % stride:
                continue
            v = parse_block(block, 5)
            frames.append(dict(istp=int(m.group(1)), ni=ni, nj=nj,
                               x=v[:, 0].reshape(nj, ni), y=v[:, 1].reshape(nj, ni),
                               cvol=v[:, 3].reshape(nj, ni), cgg=v[:, 4].reshape(nj, ni)))
            print(f"  flow.dat: step {frames[-1]['istp']}", file=sys.stderr)
    return frames

def read_spray(path, want):
    """spray.dat から、ステップ番号が want に入るコマだけ読む。"""
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, errors="replace") as fh:
        for line in fh:
            if "zone" not in line.lower():
                continue
            m = re.search(r"n=\s*(\d+)", line); mi = re.search(r"i=\s*(\d+)", line)
            istp, n = int(m.group(1)), int(mi.group(1))
            block = list(itertools.islice(fh, n))
            if istp in want:
                out[istp] = parse_block(block, 5)
            if len(out) == len(want):
                break
    return out

CT_COLORS = ["#1e50ff", "#14a0ff", "#00d2ff", "#00f0c8", "#00e68c",
             "#00d23c", "#8cf000", "#f0f000", "#ffb400", "#ff6400"]
CT_RATIO = np.array([0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0])

def nice_top(v):
    """川谷さんの図と同じように、上端を 1×10^n にそろえる。"""
    return 10.0 ** np.ceil(np.log10(v) - 1e-9)

def nozzles(xlen):
    """go.out.log からノズル先端の高さと半幅を読む（なければ描かない）。"""
    try:
        t = open("go.out.log", errors="replace").read(20000)
        m = re.search(r"nozzle tip y=\s*([0-9.Ee+-]+)\s+half-width x=\s*([0-9.Ee+-]+)", t)
        return float(m.group(1)) * 1e3, float(m.group(2)) * 1e3
    except Exception:
        return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flow", default="flow.dat")
    ap.add_argument("--spray", default="spray.dat")
    ap.add_argument("--outdir", default="viz_gas")
    ap.add_argument("--zmin", type=float, default=86.0)
    ap.add_argument("--zmax", type=float, default=100.0)
    ap.add_argument("--rmax", type=float, default=8.0, help="1本・対向のときの r の範囲 [mm]")
    ap.add_argument("--ctop", type=float, default=None, help="C_total の色の上端 [mol]")
    ap.add_argument("--stride", type=int, default=1, help="flow.dat のコマを間引く")
    ap.add_argument("--gif", action="store_true")
    ap.add_argument("--layout", default="absorb", choices=["absorb", "auto", "half", "mirror"],
                    help="absorb（既定）: 吸収の図だけを1枚（並列は全幅、1本・対向は r = 0〜rmax）。"
                         "half/mirror: 左に液滴・右に C_total（川谷さんの図 3.3.22 の体裁）")
    ap.add_argument("--full", action="store_true", help="下まで全部の図も viz_gas_full に出す")
    a = ap.parse_args()
    here = os.path.abspath(".")
    single = ("1_single" in here) or ("3_opposed" in here)
    if a.layout == "auto":
        a.layout = "mirror" if single else "half"
    if "3_opposed" in here and a.zmin == 86.0 and a.zmax == 100.0:
        a.zmin, a.zmax = 30.0, 70.0          # 対向は上下のノズルの間
    # flow.dat がなければ D: に移したもの（D:\kawatani_data に同じフォルダ構成）を探す
    if not os.path.exists(a.flow):
        rel = here.split("/kawatani/", 1)[-1]
        alt = os.path.join("/mnt/d/kawatani_data", rel, "flow.dat")
        if os.path.exists(alt):
            print(f"  flow.dat は D: のものを使う: {alt}", file=sys.stderr)
            a.flow = alt
        else:
            sys.exit(f"flow.dat が見つかりません（{alt} にもない）")
    # コマを自動で間引く（1000 ステップごとに書いた計算は 11 コマおき = 2.2 ms ごと）
    if a.stride == 1:
        with open(a.flow, errors="replace") as fh:
            hdr = [next(fh) for _ in range(4)]
        mi = re.search(r"i=\s*(\d+)", hdr[2]); mj = re.search(r"j=\s*(\d+)", hdr[2])
        nfr = os.path.getsize(a.flow) / (int(mi.group(1)) * int(mj.group(1)) * len(hdr[3]))
        if nfr > 20:
            a.stride = 11
    os.makedirs(a.outdir, exist_ok=True)

    fl = read_flow(a.flow, a.stride)
    if not fl:
        sys.exit("flow.dat にコマがありません")
    want = set(f["istp"] for f in fl)
    sp = read_spray(a.spray, want)
    # 円筒か平面か: 円筒は diffusion1.f に 2*pi*x がある
    cyl = False
    if os.path.exists("diffusion1.f"):
        cyl = bool(re.search(r"^\s+vcell\(i,j\) *= *.*2 *\* *pi", open("diffusion1.f", errors="replace").read(), re.M))
    x = fl[0]["x"] * 1e3; y = fl[0]["y"] * 1e3
    xlen = float(x.max() + (x[0, 1] - x[0, 0]) / 2)
    elec = fl[0]["cgg"] >= 0.44   # 電極の中（cgg = 0.445 に固定）
    zsel = (y[:, 0] >= a.zmin - 0.3) & (y[:, 0] <= a.zmax + 0.3)
    top = a.ctop
    if top is None:
        mx = max(float(np.where(elec, 0, f["cvol"])[zsel].max()) for f in fl)
        top = nice_top(mx) if mx > 0 else 1.0
    lv = top * CT_RATIO
    cmap = ListedColormap(CT_COLORS).with_extremes(under="#0000e6", over="#ff0000")
    norm = BoundaryNorm(lv, cmap.N)
    # 液滴径の色
    alld = np.concatenate([s[s[:, 4] == 0, 2] for s in sp.values()]) if sp else np.array([1e-8])
    if np.median(alld) > 1e-7:
        dlv = np.arange(2.0, 8.5, 0.5) * 1e-7; dlab = "d [m]"; dsc = 1.0
    else:
        dlv = np.arange(4.0, 30.0, 2.0); dlab = "d [nm]"; dsc = 1e9
    dnorm = BoundaryNorm(dlv, plt.cm.jet.N, extend="both")

    if a.layout == "absorb":
        xl_d, xl_g = None, ((0.0, a.rmax) if single else (0.0, xlen))
    elif a.layout == "mirror":
        xl_d, xl_g = (-a.rmax, 0.0), (0.0, a.rmax)
    else:
        xl_d, xl_g = (0.0, xlen / 2), (xlen / 2, xlen)
    mirror = a.layout == "mirror"
    zsets = [(a.zmin, a.zmax, a.outdir)]
    if a.full:
        zsets.append((0.0, 100.0, a.outdir + "_full"))
    for zmin_, zmax_, odir in zsets:
        os.makedirs(odir, exist_ok=True)
        render_all(a, fl, sp, cyl, x, y, xlen, elec, cmap, norm, lv, dnorm, dlab, dsc,
                   xl_d, xl_g, zmin_, zmax_, odir)
    print(f"C_total の色の上端 = {top:.1E} mol（{'円筒' if cyl else '平面'}、{a.layout}）")


def render_all(a, fl, sp, cyl, x, y, xlen, elec, cmap, norm, lv, dnorm, dlab, dsc,
               xl_d, xl_g, zmin_, zmax_, odir):
    pngs = []
    if a.layout == "absorb":
        return render_absorb(a, fl, x, y, xlen, elec, cmap, norm, lv, xl_g, zmin_, zmax_, odir, cyl)
    for f in fl:
        fig = plt.figure(figsize=(7.4 if a.layout == "mirror" else 9.0, 5.2))
        ax = fig.add_axes([0.17, 0.11, 0.62, 0.80])
        cax_d = fig.add_axes([0.04, 0.11, 0.022, 0.72])
        cax_g = fig.add_axes([0.88, 0.11, 0.022, 0.72])
        # 右: C_total
        C = np.ma.masked_where(elec, f["cvol"])
        Cp = np.ma.masked_where(elec, np.clip(f["cvol"], lv[0] * 0.5, None))
        m = x >= xl_g[0] - 1e-9
        cols = np.where(m[0])[0]
        ax.pcolormesh(x[:, cols], y[:, cols], Cp[:, cols], cmap=cmap, norm=norm, shading="nearest", zorder=1)
        E = np.ma.masked_where(~elec, np.ones_like(f["cvol"]))
        ax.pcolormesh(x, y, E, cmap=ListedColormap(["0.8"]), shading="nearest", zorder=2)
        if a.layout != "mirror":
            ax.pcolormesh(x[:, :cols[0]], y[:, :cols[0]],
                          np.ma.masked_where(~elec[:, :cols[0]], np.ones_like(f["cvol"][:, :cols[0]])),
                          cmap=ListedColormap(["0.8"]), shading="nearest", zorder=2)
        else:
            ax.pcolormesh(-x, y, E, cmap=ListedColormap(["0.8"]), shading="nearest", zorder=2)
        # 左: 液滴
        s = sp.get(f["istp"])
        if s is not None:
            fly = s[:, 4] == 0
            xp = s[fly, 0] * 1e3; yp = s[fly, 1] * 1e3; dp = s[fly, 2] * dsc
            if a.layout == "mirror":
                xp = -xp; sel = xp >= xl_d[0]
            else:
                sel = xp <= xl_d[1]
            ax.scatter(xp[sel], yp[sel], s=2, c=dp[sel], cmap=plt.cm.jet, norm=dnorm,
                       edgecolors="none", zorder=3)
        nz = nozzles(xlen)
        if nz:
            from matplotlib.patches import Rectangle
            for x0 in ((-nz[1],) if a.layout == "mirror" else (0.0, xlen - nz[1])):
                ax.add_patch(Rectangle((x0, nz[0]), nz[1], 100.0 - nz[0], facecolor="0.8",
                                       edgecolor="0.5", lw=0.6, zorder=5))
            if a.layout == "mirror":
                ax.add_patch(Rectangle((0.0, nz[0]), nz[1], 100.0 - nz[0], facecolor="0.8",
                                       edgecolor="0.5", lw=0.6, zorder=5))
        ax.axvline(xl_g[0], color="k", lw=0.6, zorder=4)
        ax.set_xlim(xl_d[0], xl_g[1]); ax.set_ylim(zmin_, zmax_)
        ax.set_xlabel("r [mm]" if cyl else "x [mm]"); ax.set_ylabel("z [mm]" if cyl else "y [mm]")
        t_ms = f["istp"] * 0.2e-3
        ax.set_title(f"t = {t_ms:.1f} ms  (step {f['istp']})")
        fig.colorbar(plt.cm.ScalarMappable(norm=dnorm, cmap=plt.cm.jet), cax=cax_d)
        cax_d.set_title(dlab, fontsize=9, loc="left"); cax_d.yaxis.set_ticks_position("left")
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax_g, extend="max")
        cb.set_ticks(lv); cb.set_ticklabels([f"{v:.0E}".replace("E-0","E-") for v in lv]); cb.ax.tick_params(labelsize=7)
        cax_g.set_title("$C_{total}$ [mol]", fontsize=9, loc="left")
        fn = os.path.join(odir, f"ctotal_{f['istp']:07d}.png")
        fig.savefig(fn, dpi=130); plt.close(fig); pngs.append(fn)
        print("  ->", fn, file=sys.stderr)
    if a.gif and len(pngs) > 1:
        from PIL import Image
        ims = [Image.open(p).convert("P", palette=Image.ADAPTIVE) for p in pngs]
        ims[0].save(os.path.join(odir, "ctotal_animation.gif"), save_all=True,
                    append_images=ims[1:], duration=250, loop=0)

def render_absorb(a, fl, x, y, xlen, elec, cmap, norm, lv, xl_g, zmin_, zmax_, odir, cyl):
    """吸収の図だけを1枚（今までの visualize.py の吸収の図と同じ配置）。"""
    from matplotlib.patches import Rectangle
    pngs = []
    w = 9.5 if (xl_g[1] - xl_g[0]) > 60 else 7.0
    for f in fl:
        fig = plt.figure(figsize=(w, 4.6))
        ax = fig.add_axes([0.08, 0.12, 0.72, 0.74])
        cax = fig.add_axes([0.84, 0.12, 0.025, 0.70])
        Cp = np.ma.masked_where(elec, np.clip(f["cvol"], lv[0] * 0.5, None))
        ax.pcolormesh(x, y, Cp, cmap=cmap, norm=norm, shading="nearest", zorder=1)
        E = np.ma.masked_where(~elec, np.ones_like(f["cvol"]))
        ax.pcolormesh(x, y, E, cmap=ListedColormap(["0.8"]), shading="nearest", zorder=2)
        nz = nozzles(xlen)
        if nz:
            for x0 in ((0.0,) if xl_g[1] < xlen - 1e-6 else (0.0, xlen - nz[1])):
                ax.add_patch(Rectangle((x0, nz[0]), nz[1], 100.0 - nz[0], facecolor="0.8",
                                       edgecolor="0.5", lw=0.6, zorder=5))
        ax.set_xlim(*xl_g); ax.set_ylim(zmin_, zmax_)
        ax.set_xlabel("r [mm]" if cyl else "x [mm]"); ax.set_ylabel("z [mm]" if cyl else "y [mm]")
        ax.set_title(f"Carbon dioxide absorption (gas phase)   t = {f['istp'] * 0.2e-3:.1f} ms  (step {f['istp']})", fontsize=10)
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, extend="max")
        cb.set_ticks(lv); cb.set_ticklabels([f"{v:.0E}".replace("E-0", "E-") for v in lv]); cb.ax.tick_params(labelsize=7)
        cax.set_title("$C_{total}$ [mol]", fontsize=9, loc="left")
        fn = os.path.join(odir, f"ctotal_{f['istp']:07d}.png")
        fig.savefig(fn, dpi=130); plt.close(fig); pngs.append(fn)
        print("  ->", fn, file=sys.stderr)
    if a.gif and len(pngs) > 1:
        from PIL import Image
        ims = [Image.open(p).convert("P", palette=Image.ADAPTIVE) for p in pngs]
        ims[0].save(os.path.join(odir, "ctotal_animation.gif"), save_all=True,
                    append_images=ims[1:], duration=250, loop=0)


if __name__ == "__main__":
    main()
