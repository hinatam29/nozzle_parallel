# -*- coding: utf-8 -*-
"""
静電噴霧シミュレーション結果の可視化プログラム
===================================================

2 本ノズル並列配置の計算結果を可視化します。

読み込むファイル
  flow.dat  : 格子場データ（Tecplot point 形式）
              変数 = x, y, phi(電位), cvolume(吸収体積), cgg(吸収濃度)
  spray.dat : 噴霧粒子データ（Tecplot point 形式）
              変数 = xp, yp, dp(粒径), clt(粒子ごとの吸収量), pout(状態)

見せるもの
  ・両ノズルからの噴霧粒子（散布図。大きさ=粒径 dp、色=その液滴が吸収した
    二酸化炭素量 clt）… 「液滴が二酸化炭素を吸収している様子」
  ・気相の二酸化炭素吸収分布（cgg のヒートマップ。対向電極の帯は自動でマスク）

生成物
  1) 静止画スナップショット (PNG) : 代表的な時刻を複数枚
  2) アニメーション (MP4 もしくは GIF) : 全時間の変化

使い方（例）
  python visualize.py
  python visualize.py --mode anim --fps 8
  python visualize.py --field cvolume
  python visualize.py --cmax 0.02        # 吸収場の色スケール上限を手動指定

必要ライブラリ : numpy, matplotlib
"""

import argparse
import os
import re
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation


# ----------------------------------------------------------------------
# 日本語フォント（無ければ英語ラベルに自動切替）
# ----------------------------------------------------------------------
def setup_font():
    from matplotlib import font_manager

    candidates = [
        "Yu Gothic", "Meiryo", "MS Gothic", "MS PGothic",
        "Hiragino Sans", "Hiragino Kaku Gothic Pro",
        "Noto Sans CJK JP", "Noto Sans JP",
        "IPAexGothic", "IPAGothic", "TakaoGothic", "VL Gothic",
    ]
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in candidates:
        if name in available:
            matplotlib.rcParams["font.family"] = name
            matplotlib.rcParams["axes.unicode_minus"] = False
            return True
    return False


JP = setup_font()

L = {
    "spray": "噴霧粒子" if JP else "Spray droplets",
    "r": "半径方向 r [m]" if JP else "radial  r [m]",
    "z": "軸方向 z [m]" if JP else "axial  z [m]",
    "cbar_cgg": "気相 吸収濃度 cgg" if JP else "gas absorption cgg",
    "cbar_cvol": "吸収体積 cvolume" if JP else "absorbed volume cvolume",
    "cbar_clt": "液滴が吸収した二酸化炭素量 clt"
                if JP else "droplet absorbed clt",
    "step": "ステップ" if JP else "step",
    "title": "静電噴霧と二酸化炭素吸収"
             if JP else "Electrospray & carbon dioxide absorption",
    "n1": "ノズル1" if JP else "nozzle 1",
    "n2": "ノズル2" if JP else "nozzle 2",
}

# 対向電極の帯では cgg が 0.445 に固定されている（人工値）。この値以上は
# 実際の吸収ではないのでマスクして色スケールから外す。
CGG_ELECTRODE = 0.445


# ----------------------------------------------------------------------
# Tecplot point 形式パーサ
# ----------------------------------------------------------------------
_num_re = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?")
# Fortran は指数が3桁になると 'E' を省略する : 0.35212361-147 -> 0.35212361E-147
_exp_fix = re.compile(r"(\d\.\d+)([+-]\d{2,3})(?![\d.])")


def fix_fortran_exp(s):
    return _exp_fix.sub(r"\1E\2", s)


def _to_float(tok):
    try:
        return float(tok.replace("D", "E").replace("d", "e"))
    except ValueError:
        return np.nan


def parse_variables(line):
    return [n.strip() for n in re.findall(r'"([^"]*)"', line)]


def parse_zone_header(line):
    istp = None
    m = re.search(r"n=\s*(\d+)", line)
    if m:
        istp = int(m.group(1))
    ni = nj = 1
    m = re.search(r"\bi\s*=\s*(\d+)", line, re.I)
    if m:
        ni = int(m.group(1))
    m = re.search(r"\bj\s*=\s*(\d+)", line, re.I)
    if m:
        nj = int(m.group(1))
    return istp, ni, nj


def read_tecplot(path, stride=1, max_frames=None, verbose=True):
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    varnames = None
    frames = []
    zone_index = -1

    with open(path, "r", errors="replace") as fh:
        pending = None
        buf = []
        keep = False

        def flush():
            nonlocal pending, buf
            if pending is None:
                return
            need = pending["npoints"] * (pending["ncols"] or 0)
            if keep and pending["ncols"] and len(buf) >= need:
                arr = np.array(buf[:need], dtype=np.float32)
                arr = arr.reshape(pending["npoints"], pending["ncols"])
                frames.append(dict(istp=pending["istp"], ni=pending["ni"],
                                   nj=pending["nj"], data=arr))
            pending = None
            buf = []

        for line in fh:
            s = line.strip()
            if not s:
                continue
            low = s.lower()
            if low.startswith("title"):
                continue
            if low.startswith("variable"):
                varnames = parse_variables(s)
                continue
            if low.startswith("zone"):
                flush()
                zone_index += 1
                istp, ni, nj = parse_zone_header(s)
                ncols = len(varnames) if varnames else None
                keep = (zone_index % stride == 0) and (
                    max_frames is None or len(frames) < max_frames)
                pending = dict(istp=istp, ni=ni, nj=nj,
                               npoints=ni * nj, ncols=ncols)
                buf = []
                if verbose and zone_index % 10 == 0:
                    print(f"  zone {zone_index} (istp={istp}) ...",
                          file=sys.stderr)
                continue
            if pending is None or not keep:
                continue
            # 空白区切りでトークン化（NaN/Infinity も _to_float で処理）。
            # 3桁指数の E 省略はここで補正する。
            toks = fix_fortran_exp(s).split()
            if pending["ncols"] is None:
                pending["ncols"] = len(toks)
            for t in toks:
                buf.append(_to_float(t))
        flush()

    if varnames is None:
        raise ValueError(f"{path}: variables 行が見つかりません")
    if verbose:
        print(f"  -> {len(frames)} フレーム読み込み ({path})", file=sys.stderr)
    return varnames, frames


# ----------------------------------------------------------------------
# データ取り出し
# ----------------------------------------------------------------------
def field_grid(varnames, frame, field):
    idx = {n: i for i, n in enumerate(varnames)}
    ni, nj = frame["ni"], frame["nj"]
    d = frame["data"]
    xname = "x" if "x" in idx else varnames[0]
    yname = "y" if "y" in idx else varnames[1]
    X = d[:, idx[xname]].reshape(nj, ni)
    Y = d[:, idx[yname]].reshape(nj, ni)
    if field not in idx:
        field = "cvolume" if "cvolume" in idx else varnames[-1]
    C = d[:, idx[field]].reshape(nj, ni).astype(float)
    return X, Y, C, field


def masked_field(C, field):
    """電極帯（cgg≒0.445 の固定値）と非有限値をマスクする。"""
    C = np.array(C, dtype=float)
    bad = ~np.isfinite(C)
    if field == "cgg":
        bad |= C >= CGG_ELECTRODE * 0.999
    return np.ma.array(C, mask=bad)


def spray_points(varnames, frame, active_only=True):
    idx = {n: i for i, n in enumerate(varnames)}
    d = frame["data"]
    xp = d[:, idx.get("xp", 0)].astype(float)
    yp = d[:, idx.get("yp", 1)].astype(float)
    dp = d[:, idx["dp"]].astype(float) if "dp" in idx else np.zeros(len(xp))
    clt = d[:, idx["clt"]].astype(float) if "clt" in idx else np.zeros(len(xp))
    pout = d[:, idx["pout"]].astype(float) if "pout" in idx else np.zeros(len(xp))
    good = np.isfinite(xp) & np.isfinite(yp)
    if active_only:
        good &= (pout == 0)
    return xp[good], yp[good], dp[good], clt[good]


# ----------------------------------------------------------------------
# スケール（全フレーム共通）
# ----------------------------------------------------------------------
def field_vmax(varnames, frames, field, override=None):
    if override is not None:
        return override
    vmax = 0.0
    for fr in frames:
        _, _, C, used = field_grid(varnames, fr, field)
        M = masked_field(C, used).compressed()
        M = M[M > 0]
        if M.size:
            vmax = max(vmax, float(np.percentile(M, 99.0)))
    return vmax if vmax > 0 else 1.0


def clt_vmax(varnames, frames):
    mx = 0.0
    if not frames:
        return 1.0
    idx = {n: i for i, n in enumerate(varnames)}
    if "clt" not in idx:
        return 1.0
    for fr in frames:
        c = fr["data"][:, idx["clt"]].astype(float)
        c = c[np.isfinite(c)]
        if c.size:
            mx = max(mx, float(c.max()))
    return mx if mx > 0 else 1.0


def extent_of(varnames, frame):
    X, Y, _, _ = field_grid(varnames, frame, varnames[-1])
    return [float(X.min()), float(X.max()), float(Y.min()), float(Y.max())]


def marker_sizes(dp):
    dp = np.asarray(dp, dtype=float)
    if dp.size == 0 or dp.max() == dp.min():
        return np.full(dp.shape, 22.0)
    s = (dp - dp.min()) / (dp.max() - dp.min())
    return 10.0 + 70.0 * s


def draw_nozzles(ax, extent):
    xmin, xmax, ymax = extent[0], extent[1], extent[3]
    for xpos, name, ha in ((xmin, L["n1"], "left"), (xmax, L["n2"], "right")):
        ax.scatter([xpos], [ymax], marker="v", s=130, c="red",
                   edgecolors="k", linewidths=0.6, zorder=6, clip_on=False)
        ax.annotate(name, xy=(xpos, ymax), xytext=(0, 8),
                    textcoords="offset points", ha=ha, va="bottom",
                    color="red", fontsize=9, zorder=6, clip_on=False)


# ----------------------------------------------------------------------
# 1 フレーム描画
# ----------------------------------------------------------------------
def draw_frame(ax, cax_f, cax_p, vnf, fframe, vns, sframe,
               field, fvmax, cltmax, extent):
    ax.clear()
    X, Y, C, used = field_grid(vnf, fframe, field)
    Cm = masked_field(C, used)

    cmap_f = plt.cm.YlGnBu.copy()
    cmap_f.set_bad("0.85")  # 電極帯などは薄いグレー
    pcm = ax.pcolormesh(X, Y, Cm, cmap=cmap_f, shading="auto",
                        vmin=0.0, vmax=fvmax)

    sc = None
    if sframe is not None:
        xp, yp, dp, clt = spray_points(vns, sframe)
        if len(xp) > 0:
            sc = ax.scatter(xp, yp, s=marker_sizes(dp), c=clt,
                            cmap="plasma", vmin=0.0, vmax=cltmax,
                            edgecolors="k", linewidths=0.25, alpha=0.9,
                            zorder=5)

    draw_nozzles(ax, extent)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel(L["r"])
    ax.set_ylabel(L["z"])
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f'{L["title"]}   ({L["step"]} = {fframe.get("istp")})')

    cax_f.cla()
    cbl = L["cbar_cgg"] if used == "cgg" else L["cbar_cvol"]
    plt.colorbar(pcm, cax=cax_f, label=cbl)

    cax_p.cla()
    if sc is not None:
        plt.colorbar(sc, cax=cax_p, orientation="horizontal",
                     label=L["cbar_clt"])
    else:
        cax_p.axis("off")
    return pcm, sc


def make_figure():
    fig = plt.figure(figsize=(6.6, 7.2))
    ax = fig.add_axes([0.12, 0.20, 0.70, 0.72])
    cax_f = fig.add_axes([0.85, 0.20, 0.03, 0.72])   # 場（縦）
    cax_p = fig.add_axes([0.20, 0.07, 0.52, 0.022])  # 粒子 clt（横）
    return fig, ax, cax_f, cax_p


# ----------------------------------------------------------------------
def make_snapshots(vnf, ff, vns, sf, field, outdir, nsnap, fvmax, cltmax, ext):
    if not ff:
        print("場データのフレームがありません。", file=sys.stderr)
        return
    picks = np.unique(np.linspace(0, len(ff) - 1, min(nsnap, len(ff))).astype(int))
    for k in picks:
        fig, ax, cax_f, cax_p = make_figure()
        sframe = sf[k] if sf and k < len(sf) else None
        draw_frame(ax, cax_f, cax_p, vnf, ff[k], vns, sframe,
                   field, fvmax, cltmax, ext)
        p = os.path.join(outdir, f"snapshot_{k:03d}_istp{ff[k].get('istp')}.png")
        fig.savefig(p, dpi=140)
        plt.close(fig)
        print(f"  スナップショット保存: {p}", file=sys.stderr)


def make_animation(vnf, ff, vns, sf, field, outdir, fps, fvmax, cltmax, ext):
    if not ff:
        print("場データのフレームがありません。", file=sys.stderr)
        return
    fig, ax, cax_f, cax_p = make_figure()

    def update(k):
        sframe = sf[k] if sf and k < len(sf) else None
        draw_frame(ax, cax_f, cax_p, vnf, ff[k], vns, sframe,
                   field, fvmax, cltmax, ext)
        return []

    anim = animation.FuncAnimation(fig, update, frames=len(ff), blit=False)
    out_mp4 = os.path.join(outdir, "animation.mp4")
    out_gif = os.path.join(outdir, "animation.gif")
    try:
        anim.save(out_mp4, writer=animation.FFMpegWriter(fps=fps, bitrate=2400),
                  dpi=130)
        print(f"  アニメーション保存: {out_mp4}", file=sys.stderr)
    except Exception as e:
        print(f"  MP4 出力不可 ({e}) -> GIF を試します。", file=sys.stderr)
        try:
            anim.save(out_gif, writer=animation.PillowWriter(fps=fps), dpi=100)
            print(f"  アニメーション保存: {out_gif}", file=sys.stderr)
        except Exception as e2:
            print(f"  GIF も失敗: {e2}", file=sys.stderr)
    plt.close(fig)


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="静電噴霧と二酸化炭素吸収の可視化")
    ap.add_argument("--flow", default="flow.dat")
    ap.add_argument("--spray", default="spray.dat")
    ap.add_argument("--outdir", default="viz_out")
    ap.add_argument("--field", default="cgg", choices=["cgg", "cvolume"])
    ap.add_argument("--mode", default="both", choices=["both", "anim", "snap"])
    ap.add_argument("--nsnap", type=int, default=6)
    ap.add_argument("--fps", type=int, default=8)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-frames", type=int, default=None)
    ap.add_argument("--cmax", type=float, default=None,
                    help="吸収場の色スケール上限を手動指定（既定は自動）")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("場データを読み込み中 ...", file=sys.stderr)
    vnf, ff = read_tecplot(args.flow, stride=args.stride,
                           max_frames=args.max_frames)
    vns, sf = [], []
    if os.path.exists(args.spray):
        print("噴霧データを読み込み中 ...", file=sys.stderr)
        vns, sf = read_tecplot(args.spray, stride=args.stride,
                               max_frames=args.max_frames)
    else:
        print(f"注意: {args.spray} が無いので粒子散布は省略します。",
              file=sys.stderr)

    fvmax = field_vmax(vnf, ff, args.field, override=args.cmax)
    cmx = clt_vmax(vns, sf) if sf else 1.0
    ext = extent_of(vnf, ff[0])
    print(f"  吸収場スケール上限={fvmax:.3e}, clt上限={cmx:.3e}", file=sys.stderr)

    if args.mode in ("snap", "both"):
        make_snapshots(vnf, ff, vns, sf, args.field, args.outdir,
                       args.nsnap, fvmax, cmx, ext)
    if args.mode in ("anim", "both"):
        make_animation(vnf, ff, vns, sf, args.field, args.outdir,
                       args.fps, fvmax, cmx, ext)

    print("完了。出力先: " + os.path.abspath(args.outdir), file=sys.stderr)


if __name__ == "__main__":
    main()
