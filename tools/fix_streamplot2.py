"""visualize.py の電気力線（streamplot）を、すべての種類の図で確実に描けるようにする（2026-11-16）。
格子の座標の丸め誤差で「等間隔でない」とエラーになり黙って飛ばされていたので、
streamplot に渡す座標を等間隔に置き直す関数 _eq() を通す。1本・対向・並列・電場の時間変化（efieldt）すべて。
使い方: python3 tools/fix_streamplot2.py"""
import glob, os, re
helper = '''

def _eq(v):
    """streamplot 用に座標を等間隔に置き直す（丸め誤差で失敗しないように）。2026-11-16"""
    v = np.asarray(v)
    return np.linspace(v[0], v[-1], len(v))
'''
pat = re.compile(r"ax\.streamplot\((?!_eq\()([^,]+?) \* 1e3, (?!_eq\()([^,]+?) \* 1e3,")
n = 0
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
for f in glob.glob(os.path.join(root, '**', 'visualize.py'), recursive=True):
    s = open(f).read()
    s2 = pat.sub(lambda m: f"ax.streamplot(_eq({m.group(1)}) * 1e3, _eq({m.group(2)}) * 1e3,", s)
    if s2 != s:
        if "def _eq(" not in s2:
            i = s2.index("\nimport numpy as np")
            j = s2.index("\n", i + 1)
            s2 = s2[:j] + helper + s2[j:]
        open(f, 'w').write(s2); n += 1
print("直した visualize.py:", n)
