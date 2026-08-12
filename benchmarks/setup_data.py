"""Download every dataset the benchmark needs. Idempotent — safe to re-run.

Four sources:

1. ``RANSAC-Tutorial-Data-ValOnly.tar`` (2.2 GB) — only the first 600 MB is
   requested; that range covers all five HDF5 files of the first scene
   (``st_peters_square``). Images are not fetched: the F metric is pose error
   from ground-truth calibration.
2. ``homography.tar.gz`` (798 MB) — EVD + HPatchesSeq, both splits, with
   images (the H metric integrates over the jointly visible image area, so it
   needs image sizes).
3. ``hpatches-sequences-release.zip`` (1.28 GB) — read over HTTP range
   requests, transferring only the 580 ``H_1_N`` files (~200 KB). The tutorial
   ships ``Hgt.h5`` for ``val`` only; this recovers it for ``test``.
4. ``EVD.zip`` (28 MB) — same, for EVD's test split.

The recovery in (3)/(4) is self-checked: the same lookup is applied to the
``val`` split, where the tutorial *does* ship ground truth, and must reproduce
it exactly before any test-split ground truth is written.
"""
import argparse
import io
import shutil
import tarfile
import urllib.request
import zipfile
from pathlib import Path

import h5py
import numpy as np

DATA = Path(__file__).resolve().parent / "data"

VALONLY_URL = ("https://cmp.felk.cvut.cz/~mishkdmy/CVPR-RANSAC-Tutorial-2020/"
               "RANSAC-Tutorial-Data-ValOnly.tar")
HOMOGRAPHY_URL = ("http://cmp.felk.cvut.cz/~mishkdmy/CVPR-RANSAC-Tutorial-2020/"
                  "homography.tar.gz")
HPATCHES_URL = ("https://huggingface.co/datasets/vbalnt/hpatches/resolve/main/"
                "hpatches-sequences-release.zip")
EVD_URL = "https://cmp.felk.cvut.cz/wbs/datasets/EVD.zip"

#: Bytes of the ValOnly tar to request. st_peters_square's HDF5 files end at
#: ~364 MB; the margin covers a re-ordered rebuild of the archive.
VALONLY_RANGE = 600_000_000

F_SCENE = "st_peters_square"
F_FILES = ("matches.h5", "match_conf.h5", "K1_K2.h5", "R.h5", "T.h5")


def _log(msg):
    print(msg, flush=True)


def _download(url, dest, headers=None):
    """Stream a URL to ``dest`` via a temporary file, so a killed download
    never leaves a truncated file that a later run would treat as complete."""
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers=headers or {})
    with urllib.request.urlopen(req) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f, length=1 << 20)
    tmp.rename(dest)


class _RemoteFile(io.RawIOBase):
    """Seekable read-only file over HTTP range requests.

    Lets ``zipfile`` read a remote archive's central directory and then pull
    single members, instead of downloading the whole thing. Reads are served
    from aligned cached blocks: ``zipfile`` issues two tiny reads per member
    (local header, then payload), and the members we want are ~100-byte text
    files clustered per directory, so block caching turns hundreds of
    round-trips into a few dozen.
    """

    BLOCK = 1 << 19  # 512 KB

    def __init__(self, url):
        self.url = url
        self.pos = 0
        self._blocks = {}
        with urllib.request.urlopen(
                urllib.request.Request(url, method="HEAD")) as r:
            self.size = int(r.headers["Content-Length"])
            if "bytes" not in r.headers.get("Accept-Ranges", ""):
                raise RuntimeError(f"{url} does not support range requests")

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=io.SEEK_SET):
        base = {io.SEEK_SET: 0, io.SEEK_CUR: self.pos, io.SEEK_END: self.size}
        self.pos = base[whence] + off
        return self.pos

    def _fetch(self, start, end):
        req = urllib.request.Request(
            self.url, headers={"Range": f"bytes={start}-{end}"})
        with urllib.request.urlopen(req) as r:
            return r.read()

    def _block(self, idx):
        if idx not in self._blocks:
            start = idx * self.BLOCK
            self._blocks[idx] = self._fetch(
                start, min(start + self.BLOCK, self.size) - 1)
        return self._blocks[idx]

    def read(self, n=-1):
        if n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        if n <= 0:
            return b""
        # Bypass the cache for reads larger than a block (the central
        # directory) — caching those would just waste memory.
        if n > self.BLOCK:
            data = self._fetch(self.pos, self.pos + n - 1)
            self.pos += len(data)
            return data
        out = bytearray()
        while len(out) < n:
            idx, off = divmod(self.pos + len(out), self.BLOCK)
            chunk = self._block(idx)[off:]
            if not chunk:
                break
            out += chunk[:n - len(out)]
        self.pos += len(out)
        return bytes(out)


# --------------------------------------------------------------------------
# F data
# --------------------------------------------------------------------------

def setup_f(force=False):
    dest = DATA / "f_data" / F_SCENE
    if not force and all((dest / n).exists() for n in F_FILES):
        _log(f"F data already present: {dest}")
        return
    dest.mkdir(parents=True, exist_ok=True)
    head = DATA / "valonly_head.tar"
    if force or not head.exists():
        _log(f"fetching first {VALONLY_RANGE // 10**6} MB of {VALONLY_URL}")
        _download(VALONLY_URL, head,
                  headers={"Range": f"bytes=0-{VALONLY_RANGE - 1}"})
    # The truncated tar ends mid-member; tarfile raises at the tail, which is
    # expected — everything we need is decoded before that point.
    prefix = f"RANSAC-Tutorial-Data-ValOnly/val/{F_SCENE}/"
    got = set()
    with tarfile.open(head, "r|") as tf:
        try:
            for member in tf:
                name = member.name
                if name.startswith(prefix) and name[len(prefix):] in F_FILES:
                    src = tf.extractfile(member)
                    with open(dest / name[len(prefix):], "wb") as f:
                        shutil.copyfileobj(src, f)
                    got.add(name[len(prefix):])
                    if got == set(F_FILES):
                        break
        except tarfile.TarError:
            pass
    missing = set(F_FILES) - got
    if missing:
        raise RuntimeError(
            f"missing {sorted(missing)} in the first {VALONLY_RANGE} bytes of "
            f"the ValOnly tar — raise VALONLY_RANGE")
    head.unlink()
    _log(f"F data ready: {dest}")


# --------------------------------------------------------------------------
# H data
# --------------------------------------------------------------------------

def setup_h(force=False):
    root = DATA / "homography"
    if not force and (root / "HPatchesSeq" / "val" / "Hgt.h5").exists():
        _log(f"H data already present: {root}")
    else:
        arch = DATA / "homography.tar.gz"
        if force or not arch.exists():
            _log(f"fetching {HOMOGRAPHY_URL} (798 MB)")
            _download(HOMOGRAPHY_URL, arch)
        _log("extracting")
        DATA.mkdir(parents=True, exist_ok=True)
        with tarfile.open(arch, "r:gz") as tf:
            tf.extractall(DATA)
        arch.unlink()
        _log(f"H data ready: {root}")
    recover_test_gt(force=force)


def _hpatches_upstream(sequences):
    """{sequence: {n: 3x3 H}} for the requested sequences, via range reads."""
    zf = zipfile.ZipFile(_RemoteFile(HPATCHES_URL))
    out = {}
    for name in zf.namelist():
        parts = name.split("/")
        if len(parts) == 3 and parts[2].startswith("H_1_") \
                and parts[1] in sequences:
            H = np.loadtxt(io.StringIO(zf.read(name).decode()))
            out.setdefault(parts[1], {})[parts[2][-1]] = H
    missing = sequences - set(out)
    if missing:
        raise RuntimeError(f"sequences absent from {HPATCHES_URL}: "
                           f"{sorted(missing)}")
    return out


def _evd_upstream():
    """{name: 3x3 H} from EVD.zip."""
    path = DATA / "EVD.zip"
    if not path.exists():
        _log(f"fetching {EVD_URL} (28 MB)")
        DATA.mkdir(parents=True, exist_ok=True)
        _download(EVD_URL, path)
    zf = zipfile.ZipFile(path)
    return {n.split("/")[-1][:-4]: np.loadtxt(io.StringIO(zf.read(n).decode()))
            for n in zf.namelist()
            if n.startswith("EVD/h/") and n.endswith(".txt")}


def _lookup(dataset, upstream, key):
    """Ground-truth homography for a tutorial pair key.

    HPatches keys are ``{sequence}_1_{n}`` -> ``{sequence}/H_1_{n}``;
    EVD keys are the image-pair stem -> ``h/{stem}.txt``.
    """
    if dataset == "HPatchesSeq":
        return upstream[key[:-4]][key[-1]]
    return upstream[key]


def _pair_keys(dataset, split):
    with h5py.File(DATA / "homography" / dataset / split / "matches.h5",
                   "r") as f:
        return sorted(f.keys())


def _verify_val(dataset, upstream):
    """The tutorial ships Hgt.h5 for val; the recovery must reproduce it."""
    with h5py.File(DATA / "homography" / dataset / "val" / "Hgt.h5", "r") as f:
        keys = sorted(f.keys())
        gt = {k: np.asarray(f[k][()], float) for k in keys}
    worst = 0.0
    for k in keys:
        a, b = gt[k], np.asarray(_lookup(dataset, upstream, k), float)
        worst = max(worst, float(np.abs(a / a[2, 2] - b / b[2, 2]).max()))
    if worst > 1e-9:
        raise RuntimeError(
            f"{dataset}: recovered val ground truth disagrees with the "
            f"shipped Hgt.h5 (max |diff| = {worst:g}) — the upstream mapping "
            f"is wrong, refusing to write test ground truth")
    _log(f"  {dataset}: val self-check passed over {len(keys)} pairs "
         f"(max |diff| = {worst:g})")


def recover_test_gt(force=False):
    """Write ``test/Hgt.h5`` for both H datasets from the source datasets."""
    targets = [(ds, DATA / "homography" / ds / "test" / "Hgt.h5")
               for ds in ("HPatchesSeq", "EVD")]
    if not force and all(p.exists() for _, p in targets):
        _log("test-split ground truth already present")
        return
    _log("recovering test-split ground truth from the source datasets")
    for dataset, out in targets:
        keys = _pair_keys(dataset, "test")
        if dataset == "HPatchesSeq":
            # val keys too: the self-check below needs their homographies.
            seqs = {k[:-4] for k in keys + _pair_keys(dataset, "val")}
            upstream = _hpatches_upstream(seqs)
        else:
            upstream = _evd_upstream()
        _verify_val(dataset, upstream)
        with h5py.File(out, "w") as f:
            for k in keys:
                f.create_dataset(
                    k, data=np.asarray(_lookup(dataset, upstream, k), float))
        _log(f"  {dataset}: wrote {len(keys)} test homographies -> {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", choices=["f", "h"], help="fetch one side only")
    ap.add_argument("--force", action="store_true", help="re-download")
    args = ap.parse_args()
    DATA.mkdir(parents=True, exist_ok=True)
    if args.only != "h":
        setup_f(force=args.force)
    if args.only != "f":
        setup_h(force=args.force)


if __name__ == "__main__":
    main()
