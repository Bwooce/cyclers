from pathlib import Path

from cyclerfinder.verify.spice_kernels import ensure_jup365_kernel


def test_jup365_kernel_path_exists() -> None:
    p = ensure_jup365_kernel()
    assert Path(p).is_file()
    assert p.endswith("jup365.bsp")
