"""Check the native RNG against libc and the legacy Windows stream."""
import ctypes
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import platform
import threading

import pytest


@pytest.fixture(scope="module", params=[False, True], ids=["native", "windows"])
def rng(request, tmp_path_factory):
    # Build only the RNG so these tests also exercise the Windows branch on Unix.
    ccompiler = pytest.importorskip("setuptools._distutils.ccompiler")
    sysconfig = pytest.importorskip("setuptools._distutils.sysconfig")
    directory = tmp_path_factory.mktemp("rng")
    compiler = ccompiler.new_compiler()
    sysconfig.customize_compiler(compiler)
    source = Path(__file__).resolve().parents[1] / "src/pydegensac/degensac/bsd_random.c"
    objects = compiler.compile(
        [str(source)], output_dir=str(directory),
        macros=[("WIN32", "1")] if request.param else [],
    )
    library = str(directory / ("rng" + compiler.shared_lib_extension))
    compiler.link_shared_object(objects, library, export_symbols=[
        "degensac_srand", "degensac_rand", "degensac_srandom", "degensac_random",
    ])
    lib = ctypes.CDLL(library)
    lib.degensac_srand.argtypes = [ctypes.c_uint]
    lib.degensac_srandom.argtypes = [ctypes.c_uint]
    lib.degensac_random.restype = ctypes.c_long
    return lib, request.param or platform.system() == "Windows"


SEEDS = [0, 1, 42, 2024, 2147483646, 2147483647, 2147483648, 4294967295]


@pytest.mark.parametrize("seed", SEEDS)
def test_rand_matches_libc(rng, seed):
    lib, windows = rng
    if windows and platform.system() != "Windows":
        pytest.skip("Windows libc unavailable on this host")
    if not windows and platform.system() != "Darwin" and platform.libc_ver()[0] != "glibc":
        pytest.skip("libc compatibility is defined for glibc, Darwin and Windows")
    libc = ctypes.CDLL("ucrtbase" if windows else None)
    libc.srand(ctypes.c_uint(seed))
    expected = [libc.rand() for _ in range(1000)]
    lib.degensac_srand(seed)
    assert [lib.degensac_rand() for _ in expected] == expected


@pytest.mark.parametrize("seed", SEEDS)
def test_windows_sampling_and_degeneracy_share_stream(rng, seed):
    lib, windows = rng
    if not windows:
        pytest.skip("Windows stream only")
    # The CRT recurrence is also a portable oracle when running on Unix.
    state = seed
    expected = []
    for _ in range(1000):
        state = (state * 214013 + 2531011) & 0xffffffff
        expected.append((state >> 16) & 0x7fff)
    lib.degensac_srand(seed)
    lib.degensac_srandom(seed)
    actual = [lib.degensac_random() if i % 2 else lib.degensac_rand()
              for i in range(len(expected))]
    assert actual == expected


def test_rng_state_is_thread_local(rng):
    lib, _ = rng
    def draws(seed):
        lib.degensac_srand(seed)
        lib.degensac_srandom(seed)
        return [lib.degensac_random() if i % 2 else lib.degensac_rand()
                for i in range(1000)]
    expected = [draws(seed) for seed in SEEDS]
    barrier = threading.Barrier(len(SEEDS))
    def concurrent(seed):
        lib.degensac_srand(seed)
        lib.degensac_srandom(seed)
        barrier.wait()
        return [lib.degensac_random() if i % 2 else lib.degensac_rand()
                for i in range(1000)]
    with ThreadPoolExecutor(max_workers=len(SEEDS)) as executor:
        assert list(executor.map(concurrent, SEEDS)) == expected
